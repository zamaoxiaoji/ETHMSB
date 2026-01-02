# TFHEpp 代码导读（本次对话总结）

- TRGSW/RGSW：gadget decomposition、gadget matrix（对角线嵌入）、external product
- CMUX、BlindRotate（盲旋转）
- KeySwitch（密钥切换）
---

## 1. 参数与类型：对象长什么样

核心类型别名集中在 `thirdparty/TFHEpp/include/params.hpp`：

- `TLWE<P>`：`std::array<T, k*n + 1>`（`a` 向量 + 最后一项 `b`），见 `thirdparty/TFHEpp/include/params.hpp:35`
- `Polynomial<P>`：`std::array<T, n>`（环元素系数表示），见 `thirdparty/TFHEpp/include/params.hpp:38`
- `TRLWE<P>`：`std::array<Polynomial<P>, k+1>`（`a_0..a_{k-1}, b`），见 `thirdparty/TFHEpp/include/params.hpp:53`
- `TRGSW<P>`：`std::array<TRLWE<P>, (k+1)*l>`（把“GSW 矩阵的行”展平为 `(k+1)*l` 行），见 `thirdparty/TFHEpp/include/params.hpp:64`
- `TRGSWFFT<P>`：TRGSW 的 FFT 域表示，见 `thirdparty/TFHEpp/include/params.hpp:66`
- `BootstrappingKeyFFT<P>`：`std::array<TRGSWFFT<targetP>, domainP::n>`（每个 `i` 一把 TRGSWFFT），见 `thirdparty/TFHEpp/include/params.hpp:73-75`

参数 `P`（如 `lvl1param`）提供 `n/nbit/k/l/Bgbit/Bg/T/α/μ` 等（例如 `thirdparty/TFHEpp/include/params/128bit.hpp`）。

---

## 2. TRGSW/RGSW 的核心：gadget decomposition + gadget matrix

### 2.1 gadget 向量 `h[i]`（权重）

`hgen<P>()` 在 `thirdparty/TFHEpp/include/trgsw.hpp:34`：

- `h[i] = 1 << (w - (i+1)*Bgbit)`，其中 `w = digits(P::T)`
- 作用：把 torus 元素近似写成 `x ≈ Σ d_i * h[i]`，其中 `d_i` 是“居中 digit”（balanced digit）

### 2.2 `DecompositionPolynomial`：把多项式每个系数拆成 digit

实现：`thirdparty/TFHEpp/src/trgsw.cpp:21-39`

关键常量：

- `offset = offsetgen<P>()`：用于把 digit 平移到居中表示，见 `thirdparty/TFHEpp/src/trgsw.cpp:9-18`
- `roundoffset`：用于四舍五入（减少纯截断偏差），见 `thirdparty/TFHEpp/src/trgsw.cpp:25-27`
- `mask = (1<<Bgbit)-1`、`halfBg = 1<<(Bgbit-1)`，见 `thirdparty/TFHEpp/src/trgsw.cpp:28-30`

核心一行（对每个系数 `poly[i]`）：

- `decpoly[i] = ((((poly[i] + offset + roundoffset) >> (w - (digit+1)*Bgbit)) & mask) - halfBg)`

含义：取出第 `digit` 段（每段 `Bgbit` 位）的高位片段，并转换为范围 `[-Bg/2, Bg/2)` 的居中 digit。

FFT/NTT 版本只是把分解结果变换到乘法友好域：

- `DecompositionPolynomialFFT`：`thirdparty/TFHEpp/src/trgsw.cpp:48-54`
- `DecompositionPolynomialNTT`：`thirdparty/TFHEpp/src/trgsw.cpp:63-69`

### 2.3 TRGSW 的矩阵形状与索引 `row = i + k*l`

TRGSW 的“行”天然是二维索引 `(kblk, i)`：

- `kblk ∈ [0..k]`：对应 TRLWE 向量的第几个分量（`a_0..a_{k-1}, b` 共 `k+1` 列）
- `i ∈ [0..l-1]`：gadget decomposition 的 digit（每列块展开成 `l` 行）

把二维 `(kblk, i)` 扁平化为一维行号的最自然方式就是：

- `row = kblk * l + i`

代码里写成 `i + k * P::l`，用于：

- 对角线嵌入：`thirdparty/TFHEpp/src/trgsw.cpp:214-219`
- external product 取行：`thirdparty/TFHEpp/src/trgsw.cpp:93-99`

### 2.4 `trgswSymEncrypt`：先 Enc(0)，再加对角 gadget

实现：`thirdparty/TFHEpp/src/trgsw.cpp:207-223`

两步：

1) 每一行先设为 `trlweSymEncryptZero`（加密 0 的 TRLWE），见 `thirdparty/TFHEpp/src/trgsw.cpp:213`

2) 再把 `p*h[i]` 写到“块对角线”位置：

- 代码：`trgsw[i + k * P::l][k][j] += p[j]*h[i]`（`thirdparty/TFHEpp/src/trgsw.cpp:214-219`）
- 含义：第 `(kblk=k, digit=i)` 这行，只往第 `k` 列（对角）加入 `p*h[i]`。

当 `k=1,l=3`（TFHEpp 常见配置）时，TRGSW 可以看成 `6×2` 的多项式矩阵（忽略 Enc(0) 的噪声行，只看 gadget 结构）：

- 行 0..2（`kblk=0`）：只在列 0 放 `p*h[0..2]`
- 行 3..5（`kblk=1`）：只在列 1 放 `p*h[0..2]`

## 把 TRGSW 看成一个“矩阵”，每个单元格都是一个多项式（环元素）：
$\mathrm{TRGSW}(p)\in R^{(k+1)l\times (k+1)}$







---

## 3. External Product：`trgswfftExternalProduct` 在算什么

实现：`thirdparty/TFHEpp/src/trgsw.cpp:78-103`

输入：

- `trlwe`：一个 `TRLWE<P>`（`k+1` 个多项式）
- `trgswfft`：一个 `TRGSWFFT<P>`（`(k+1)*l` 行 × `k+1` 列）

输出：

- `res`：一个新的 `TRLWE<P>`

语义（省略噪声项）：

1) 对每个输入分量多项式 `trlwe[kblk]` 做 gadget decomposition，得到 `D_i(trlwe[kblk])`（`i=0..l-1`）。

2) 用这些 digit 多项式作为系数，对 TRGSW 的各行做线性组合：

`res = Σ_{kblk=0..k} Σ_{i=0..l-1} D_i(trlwe[kblk]) ⊗ trgsw[row=i+kblk*l]`

代码结构对应这个求和：

- 初始化：`trlwe[0]` 的 digit 0 用 `MulInFD`（`thirdparty/TFHEpp/src/trgsw.cpp:83-86`）
- 累加：`trlwe[0]` 的 digit 1..l-1 用 `FMAInFD`（`thirdparty/TFHEpp/src/trgsw.cpp:87-92`）
- 再累加：`trlwe[1..k]` 的每个 digit，按 `row=i+k*l` 取行继续 `FMAInFD`（`thirdparty/TFHEpp/src/trgsw.cpp:93-101`）
- 回到系数域：`TwistFFT`（`thirdparty/TFHEpp/src/trgsw.cpp:102`）

为什么它实现了“明文乘法”？

- `trgswSymEncrypt(p)` 的对角线 gadget 让 external product 的输出满足近似：`res ≈ p * trlwe`
- 从而解密相位也满足：`phase(res) ≈ p * phase(trlwe)`

这正是 CMUX / BlindRotate 的数学基础。

---

## 4. CMUX：external product 在哪里用到

### 4.1 `CMUXFFT(res, cs, c1, c0)`

实现：`thirdparty/TFHEpp/src/detwfa.cpp:4-13`

流程：

1) `res = c1 - c0`
2) `res = ExternalProduct(res, cs)`（`thirdparty/TFHEpp/src/detwfa.cpp:10`）
3) `res += c0`

当 `cs`（TRGSWFFT）明文是 0/1（或 ±1）时，这就实现 `cs?c1:c0`。

### 4.2 BlindRotate 用的“无拷贝 CMUX”等价式

`CMUXFFTwithPolynomialMulByXaiMinusOne`：`thirdparty/TFHEpp/src/detwfa.cpp:20-30`

它把 CMUX 写成：

- `temp = (X^a - 1) * acc`（`thirdparty/TFHEpp/src/detwfa.cpp:25-26`）
- `temp = ExternalProduct(temp, cs)`（`thirdparty/TFHEpp/src/detwfa.cpp:27`）
- `acc += temp`（`thirdparty/TFHEpp/src/detwfa.cpp:28-29`）

因此：

- `trgswfftExternalProduct<P>(temp, temp, cs);`（`thirdparty/TFHEpp/src/detwfa.cpp:27`）

只对应 **一次** `cs` 的 external product（一次 CMUX 的核心乘法步骤），并不是“对 bkfft 的所有位都做 CMUX”。

对所有位的 CMUX，是 BlindRotate 的 `for (i...)` 循环逐次触发出来的（见下一节）。

---

## 5. BlindRotate：`b̄`、`ā` 到底是什么？哪里体现了 `b - Σ a_i s_i`

实现：`thirdparty/TFHEpp/include/gatebootstrapping.hpp:15-42`

BlindRotate 的目标：把 TLWE 的相位（phase）信息编码为 accumulator（TRLWE）里 testvector 的“旋转位置”。

记：

- `N = P::targetP::n = 2^{P::targetP::nbit}`
- 取 `2N`（而不是 `N`），因为 negacyclic 环满足 `X^N = -1`，指数自然按 `2N` 周期工作；`PolynomialMulByXai` 也按 `a < N`/`a ≥ N` 分支实现符号翻转（`thirdparty/TFHEpp/include/utils.hpp:111-124`）。

### 5.1 `b̄`（由 TLWE 的 `b` 得到的初始旋转步数）

代码：`thirdparty/TFHEpp/include/gatebootstrapping.hpp:22-26`

- `b` 是 TLWE 的最后一项：`tlwe[P::domainP::k * P::domainP::n]`
- `b̄ = 2N - (quantize(b) << bitwidth)`

含义：

- `quantize(b)`：从 torus 的 `b` 取高位映射到 `[0, 2N)` 的整数
- `2N - x`：模 `2N` 意义下的“取负”，即 `b̄ ≡ -quantize(b) (mod 2N)`

随后执行：

- `acc = X^{b̄} * testvector`（`thirdparty/TFHEpp/include/gatebootstrapping.hpp:28`）

### 5.2 `ā`（由每个 TLWE 系数 `a_i` 得到的条件旋转步数）

代码：`thirdparty/TFHEpp/include/gatebootstrapping.hpp:29-41`

- `tlwe[i]` 是 TLWE 的第 `i` 个 `a_i`
- `roundoffset` 用于四舍五入（减少截断偏差），见 `thirdparty/TFHEpp/include/gatebootstrapping.hpp:30-32`
- `ā = (quantize(a_i) << bitwidth)`，是对齐到 `2^{bitwidth}` 倍数的旋转步数（见 `thirdparty/TFHEpp/include/gatebootstrapping.hpp:33-37`）
- `if (ā == 0) continue`：量化为 0 就跳过

然后对每个 `i` 做一次条件旋转：

- `CMUXFFTwithPolynomialMulByXaiMinusOne(acc, bkfft[i], ā)`（`thirdparty/TFHEpp/include/gatebootstrapping.hpp:40-41`）

其中 `bkfft[i]` 明文携带 `s_i`（见下节 5.3），因此效果近似：

- 若 `s_i = 0`：acc 不变
- 若 `s_i = 1`：acc 变为 `X^{ā} * acc`

### 5.3 `bkfft[i]` 里是什么：它让 `s_i` 在密文里参与运算

bootstrapping key 生成在 `thirdparty/TFHEpp/src/cloudkey.cpp:23-30`：

- `plainpoly[0] = sk_domain[i]`
- `bkfft[i] = trgswfftSymEncrypt(plainpoly, ..., sk_target)`

所以 `bkfft[i]` 是一个 TRGSWFFT，加密了（domain 密钥的）第 `i` 个 secret 值 `s_i`。

### 5.4 “哪里计算了 `b - Σ a_i s_i`？”

代码里 **不会显式计算** `b - Σ a_i s_i`（那需要明文密钥 `s_i`）。  
BlindRotate 通过“效果等价”实现：

- 初始旋转 `X^{-quantize(b)}`（由 `b̄` 实现）
- 逐项条件旋转 `X^{quantize(a_i) * s_i}`（由 `bkfft[i]` 的 CMUX 决定是否生效）

因此总旋转指数等价于 `-(quantize(b - Σ a_i s_i))`（差一个符号约定），把 TLWE 相位绑定到 accumulator 的旋转位置上。

`bitwidth` 的作用：当 `num_out>1` 时，让 `b̄/ā` 对齐到 `2^{bitwidth}`，从而一次 BlindRotate 后可以用多个 `SampleExtractIndex(..., i)` 取出多个输出槽位（见 `GateBootstrappingManyLUT`：`thirdparty/TFHEpp/include/gatebootstrapping.hpp:121-131`）。

---

## 6. Bootstrapping：BlindRotate + SampleExtract（以及包裹的 keyswitch）

### 6.1 基本门自举：`GateBootstrappingTLWE2TLWEFFT`

实现：`thirdparty/TFHEpp/src/gatebootstrapping.cpp:15-23`

1) `BlindRotate(acc, tlwe, bkfft, testvector)`  
2) `SampleExtractIndex(res, acc, 0)`（从 TRLWE 抽样提取为 TLWE）

### 6.2 `GateBootstrapping` 为何前后夹 keyswitch

在 `thirdparty/TFHEpp/include/gatebootstrapping.hpp:141-159`，TFHEpp 固定了一条常用通道 `lvl0 ↔ lvl1`：

- 输入 `TLWE<lvl0param>`：先 bootstrap 到 `lvl1`，再 keyswitch 回 `lvl0`
- 输入 `TLWE<lvl1param>`：先 keyswitch 到 `lvl0`，再 bootstrap 回 `lvl1`

---

## 7. KeySwitch：`IdentityKeySwitch` 的分解与查表减法

实现：`thirdparty/TFHEpp/src/keyswitch.cpp:14-70`

核心步骤：

1) 把输入 `b`（TLWE 最后一项）搬到输出（必要时做位宽缩放），见 `thirdparty/TFHEpp/src/keyswitch.cpp:27-37`
2) 对每个 `a_i` 做 base `2^{basebit}` 的 `t` 位分解（`prec_offset` 做近似四舍五入），见 `thirdparty/TFHEpp/src/keyswitch.cpp:18-49`
3) 对每个 digit（若非 0），减掉预计算表项 `ksk[i][j][aij-1]`（整条 TLWE 向量逐项相减），见 `thirdparty/TFHEpp/src/keyswitch.cpp:50-62`

KeySwitchingKey 的生成在 `thirdparty/TFHEpp/src/cloudkey.cpp:55-70`（`ikskgen`）：

- `ksk[i][j][u] = Enc_target( s_i * (u+1) * 2^{w-(j+1)*basebit} )`

因此 keyswitch 的本质：用 target key 下的密文把 domain key 的影响“抵消掉”，使相位在 target key 下保持一致。

---




