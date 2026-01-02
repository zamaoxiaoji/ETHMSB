# IKS 数学过程解析

*导出时间: 2025/12/26 11:34:31*

---

0. 先把代码里的对象换成数学记号
0.1 TLWE<P> 是一个长度为
$N + 1,\quad N = P::k \cdot P::n$
的数组（最后一项是 $b$），也就是：


tlwe[0..N-1] 对应向量 $\mathbf a \in \mathbb Z_{2^w}^{N}$


tlwe[N] 对应标量 $b \in \mathbb Z_{2^w}$


用下标区分 domain/target：

domain 参数：P::domainP
$N_d = k_d n_d,\quad w_d = \texttt{digits(domainT)}$

target 参数：P::targetP
$N_t = k_t n_t,\quad w_t = \texttt{digits(targetT)}$

其中位宽 $w=\texttt{numeric\_limits<T>::digits}$，对应模数
$Q_d = 2^{w_d},\qquad Q_t = 2^{w_t}.$
0.2 Torus 表示
在 TFHE 里用整数环 $\mathbb Z_{2^w}$ 表示 torus（实数模 1）：

一个 torus 实数 $x \in \mathbb R/\mathbb Z$ 用整数
$\lfloor x\cdot 2^w\rceil \bmod 2^w$
表示。

代码里 1ULL << (w_t - (j+1)*basebit) 这种“巨大移位”，其实对应的是 torus 上的分数：
$\frac{2^{w_t - (j+1)b}}{2^{w_t}} = 2^{-(j+1)b} = \frac{1}{B^{j+1}},
\quad B=2^b.$

1. TLWE 加密
对某个密钥 $\mathbf s$（通常是二值向量），TLWE 样本：
$\mathbf c = (\mathbf a, b)\in\mathbb Z_{2^w}^{N}\times \mathbb Z_{2^w}$
满足（对称加密的典型形式）：
$b = \langle \mathbf a,\mathbf s\rangle + \mu + e \pmod{2^w}$
其中：


$\mu$：明文（也是 torus 元素）


$e$：噪声


解密看的是相位 phase：
$\phi_{\mathbf s}(\mathbf c)=b-\langle \mathbf a,\mathbf s\rangle
= \mu + e \pmod{2^w}.$

2. KeySwitch（IKS）的目标
一个 domain 密钥 $\mathbf s_d$ 下的 TLWE：
$\mathbf c_d=(\mathbf a, b)\ \text{under }\mathbf s_d.$
得到一个 target 密钥 $\mathbf s_t$ 下的 TLWE：
$\mathbf c_t=(\mathbf a', b')\ \text{under }\mathbf s_t,$
并且保持“明文不变”（身份切换）：
$\phi_{\mathbf s_t}(\mathbf c_t)\approx \mu \quad(\text{噪声变大但可控}).$

3. KeySwitchingKey代码：
   using KeySwitchingKey = std::array<
    std::array std::array<TLWE<targetP, (1<<basebit)-1>, t>,
    domainN>;

展开就是三重索引：


$i \in [0, N_d-1]$（domain 的每个 secret key 分量）


$j \in [0, t-1]$（分解层数）


$k \in [0, B-2]$，其中 $B=2^{b}$，basebit=b
实际 digit 值是 $d=k+1 \in [1,B-1]$

所以存的是：
$\mathrm{KSK}[i][j][d] \in \text{TLWE}_{tgt}.$

4. ikskgen：KSK 里每个条目到底加密了什么
   ksk[i][j][k] = Enc_(
    s_d[i] * (k+1) * 2^{w_t - (j+1)*b},
    α, s_t);

4.1 数学表达
记 $d=k+1$，且 $g_j = 2^{w_t-(j+1)b}$。则：
$\boxed{
\mathrm{KSK}[i][j][d] = \mathrm{Enc}_{\mathbf s_t}\Big( s_{d,i}\cdot d\cdot g_j \Big)
}$
其中 $s_{d,i}$ 是 domain 密钥 $\mathbf s_d$ 的第 $i$ 个分量。
4.2 gadget（分解基底）
因为：
$\frac{d\cdot g_j}{2^{w_t}}
= d\cdot 2^{-(j+1)b}
= \frac{d}{B^{j+1}}.$
所以 KSK 预先准备了：


对每个 $i$：加密了 $s_{d,i}$ 的各种“分数倍”


分数的基底是 $\left(\frac{1}{B},\frac{1}{B^2},\dots,\frac{1}{B^t}\right)$

5. IdentityKeySwitch：真正的 keyswitch 如何计算
代码顺序拆成三段：
(A) 处理常数项 $b$ → (B) 分解每个 $a_i$ → (C) 查 KSK 并相减

5.A 先把输入的 $b$ 映射到 target 位宽
代码：
if (w_d == w_t) b' = b
else if (w_d > w_t) b' = (b + 2^{(w_d-w_t-1)}) >> (w_d-w_t)
else b' = b << (w_t - w_d)


若 $w_t>w_d$：左移等价于乘 $2^{w_t-w_d}$


若 $w_t<w_d$：右移并加半个单位实现“最近整数舍入”

可以理解为：
$b' \approx \left\lfloor b\cdot 2^{w_t-w_d}\right\rceil \pmod{2^{w_t}}.$

注意：只转换了 $b$，没转换 $\mathbf a$。

5.B 分解（decomposition）：把每个 $a_i$ 拆成 $t$ 个 base-$B$ digit
代码：
prec_offset = 2^{w_d - (1 + b*t)}
aibar = a_i + prec_offset
aij = (aibar >> (w_d - (j+1)*b)) & (B-1)

5.B.1 prec_offset 是什么
$\texttt{prec\_offset} = 2^{w_d - (1+b t)}.$
把它除以 $2^{w_d}$ 看 torus 意义：
$\frac{2^{w_d-(1+bt)}}{2^{w_d}}
=2^{-(1+bt)}
=\frac{1}{2B^t}.$
它的作用：
让“截断分解”变成“近似四舍五入分解”，把分解误差控制在 $\le \frac{1}{2B^t}$ 量级（对应整数误差 $\le 2^{w_d-(bt+1)}$）。
5.B.2 digit 提取的数学形式
令
$\bar a_i = a_i + 2^{w_d-(bt+1)} \pmod{2^{w_d}}.$
代码里的 digit：
$\boxed{
a_{i,j}
=
\left(
\left\lfloor
\frac{\bar a_i}{2^{w_d-(j+1)b}}
\right\rfloor
\right)\bmod B
\in\{0,1,\dots,B-1\}.
}$
这就是：


>> (w_d - (j+1)b)：把目标那一组 b 位移到最低位


& (B-1)：mask 出 b 位

5.C 查表并相减：输出密文是如何组出来的
代码：
res = 0
res.b = b'
for i:
  for j:
    if aij != 0:
       res -= ksk[i][j][aij-1]

写成数学形式（忽略 aij=0 的跳过优化）：
$\boxed{
\mathbf c_t
=
( \mathbf 0, b')
-
\sum_{i=0}^{N_d-1}\sum_{j=0}^{t-1}
\mathrm{KSK}[i][j][a_{i,j}].
}$
注意你 KSK 的第三维只存 $1..B-1$，所以代码用 aij-1 做索引。

6. 相减完成 key switching
  6.1 写出每个 KSK 条目的 TLWE 形式
  记选中的 KSK 样本（在 target 参数下）为：
  $\mathrm{KSK}[i][j][a_{i,j}]
  =
  (\mathbf A_{i,j},\, B_{i,j})$
  它是对明文
  $m_{i,j} = s_{d,i}\cdot a_{i,j}\cdot g_j,\qquad g_j=2^{w_t-(j+1)b}$
  在 target key $\mathbf s_t$ 下的 TLWE 加密，所以：
  $B_{i,j}
  =
  \langle \mathbf A_{i,j}, \mathbf s_t\rangle
  +
  m_{i,j}
  +
  e_{i,j}
  \pmod{2^{w_t}}.$
  6.2 输出密文 $\mathbf c_t=(\mathbf a', b'')$
  由相减得到：
  $\mathbf a' = -\sum_{i,j}\mathbf A_{i,j},
  \qquad
  b'' = b' - \sum_{i,j} B_{i,j}.$
  6.3 看输出在 target key 下的相位
  $\begin{aligned}
  \phi_{\mathbf s_t}(\mathbf c_t)
  &= b'' - \langle \mathbf a', \mathbf s_t\rangle\\
  &=
  \left(b' - \sum_{i,j}B_{i,j}\right)
  -
  \left\langle -\sum_{i,j}\mathbf A_{i,j},\, \mathbf s_t\right\rangle\\
  &=
  b' - \sum_{i,j}\Big(\langle \mathbf A_{i,j}, \mathbf s_t\rangle + m_{i,j}+e_{i,j}\Big)
  +
  \sum_{i,j}\langle \mathbf A_{i,j}, \mathbf s_t\rangle\\
  &=
  \boxed{
  b' - \sum_{i,j} m_{i,j} - \sum_{i,j} e_{i,j}.
  }
  \end{aligned}$

  $\langle \mathbf A_{i,j}, \mathbf s_t\rangle$ 完全抵消了
  用加密的方式把旧密钥信息搬到新密钥下，并在相位里消掉新密钥
  6.4 $\sum m_{i,j}$ 与旧的 $\langle \mathbf a,\mathbf s_d\rangle$ 的关系
  $\sum_{i,j} m_{i,j}
  =
  \sum_i s_{d,i}\sum_j a_{i,j}\, 2^{w_t-(j+1)b}
  \approx
  \sum_i s_{d,i}\cdot \left\lfloor \frac{a_i}{2^{w_d}}2^{w_t}\right\rceil
  \approx
  \left\lfloor \frac{\langle \mathbf a,\mathbf s_d\rangle}{2^{w_d}}2^{w_t}\right\rceil.$
  而 $b'$ 也是 $b$ 在位宽上的等价映射，因此：$\phi_{\mathbf s_t}(\mathbf c_t)
  \approx
  \left\lfloor \frac{b-\langle \mathbf a,\mathbf s_d\rangle}{2^{w_d}}2^{w_t}\right\rceil\sum e_{i,j}$
- 但输入密文相位$b-\langle \mathbf a,\mathbf s_d\rangle = \mu + e$
结论：输出确实是同一明文 $\mu$ 在新密钥下的 TLWE，只是噪声增大。

7. rec_offset 为啥要加？不加会怎样？
   不加的话 digit 提取就是纯截断，会导致分解误差偏向一侧；加上 $\frac{1}{2B^t}$ 的 offset 本质上是做“最近舍入”，把误差限制在：
   $\le \frac{1}{2B^t}$量级，保证 keyswitch 后噪声不会炸得太厉害。
