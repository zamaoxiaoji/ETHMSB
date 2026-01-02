# KeySwitching数学讲解

    TRGSW<P> trgsw;
    for (TRLWE<P> &trlwe : trgsw) trlwe = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::l; i++) {
        for (int k = 0; k < P::k + 1; k++) {
            for (int j = 0; j < P::n; j++) {
                trgsw[i + k * P::l][k][j] +=
                    static_cast<typename P::T>(p[j]) * h[i];
            }
        }
    }
    return trgsw;


1) 先把对象对齐：TRLWE vs TRGSW（RGSW）
在 TFHE/环 LWE 里，典型的（k 维）TRLWE 密文可以写成：
$\mathbf{c}=(a_0,\dots,a_{k-1},\,b)\in R^{k+1}$
其中 $R$ 是多项式环（比如 $R=\mathbb{Z}[X]/(X^n+1)$ 或 torus 版本），密钥是 $s=(s_0,\dots,s_{k-1})$。
对称加密一个明文多项式 $m\in R$ 的标准形式是：
$b=\sum_{t=0}^{k-1} a_t\cdot s_t \;+\; m \;+\; e$
因此解密（忽略噪声）就是：
$m \approx b-\sum_{t=0}^{k-1} a_t s_t$
加密 0（代码里的 trlweSymEncryptZero）就是取 $m=0$：
$b=\sum a_t s_t + e$

TRGSW / RGSW 可以把它理解成“一大块矩阵密文”
GSW 系列（包括 RGSW）不是把明文塞进一个向量里，而是塞进一个矩阵/块矩阵里，核心结构是：
$\mathbf{C} \;=\; \text{Enc}(0)\;+\; m\cdot \mathbf{G}$


$\text{Enc}(0)$：很多行的“加密 0”的 RLWE/TRLWE 样本堆起来


$\mathbf{G}$：gadget matrix（分解矩阵），用来做 gadget decomposition / external product


$m\cdot \mathbf{G}$：把明文按某种固定“模板”填进矩阵（通常是“对角块”）

函数正是在做这一句：先生成 Enc(0)，再把 $p\cdot G$ 加到指定位置。

拆成三步：

Step A：生成 gadget 向量 h
constexpr std::array<typename P::T, P::l> h = hgen<P>();

这对应理论里的 gadget decomposition 里那组权重：
$h=(h_0,h_1,\dots,h_{l-1})$
一般它会是某个基 $B_g$ 的幂次缩放（TFHE 里常见形式类似 $h_i \propto B_g^{l-1-i}$，在 torus 表示下还会带 $2^{-w}$ 的缩放），作用是：


让任意环元素 $x$ 可以被“分解”为小系数 digits：
$x \approx \sum_{i=0}^{l-1} x^{(i)}\cdot h_i$


这样外积（external product）时噪声不会炸得太快（因为 digits 小）

Step B：把整个 TRGSW 先初始化成“很多条 TRLWE(0)”
for (TRLWE<P> &trlwe : trgsw) trlwe = trlweSymEncryptZero<P>(α, key);

这一步在数学上就是构造一个矩阵（很多行），每一行都是一个 TRLWE 加密 0：
$\forall r,\quad \mathbf{c}_r = \text{Enc}_{\text{TRLWE}}(0)$
每行长得像：
$\mathbf{c}_r=(a_{r,0},\dots,a_{r,k-1},b_r),
\quad b_r=\sum_t a_{r,t}s_t + e_r$
注意：这一步还没有塞入明文 $p$。这就是 GSW 的 “Enc(0)”。

Step C：把明文 $p$ 按 gadget 模板加到特定位置
核心这句：
trgsw[i + k * P::l][k][j] += p[j] * h[i];

i：gadget level（第几层分解，0..l-1）


k：列索引（0..k，注意这里是 P::k+1，包含最后那个 b 分量）


j：多项式系数下标（0..n-1）


trgsw[...] 里每个元素其实是一条 TRLWE（一个长度 k+1 的多项式向量），所以：


trgsw[row][col]：这一行 TRLWE 的第 col 个多项式分量


[j]：该多项式的第 j 个系数


所以这句就是在做：

在第 row = i + k*l 行的 TRLWE 里，把“第 k 个分量”加上多项式 $p\cdot h_i$

TRGSW 在 TFHE 的常见存法可以看成一个矩阵：


行数：$(k+1)\cdot l$


列数：$k+1$


每个元素是一个环元素（多项式）

并且按块排列：每个“列 k”有一组 l 行（对应 i=0..l-1）。
索引 row = i + k*l 恰好就是：第 k 个块的第 i 行。
所以明文被加在：
$\text{row block}=k \quad \text{and}\quad \text{column}=k$
这就是一个标准的 block diagonal（块对角） 结构。


假设 $k=1$（所以列数 k+1=2），$l=3$（三层 gadget），那么矩阵是 $6\times 2$。
行按 row=i + k*l 排：


k=0 块：row 0,1,2


k=1 块：row 3,4,5


每一行是一个 TRLWE 向量 (a0, b)（因为 k=1 只有一个 a 分量和 b）
那么“填入明文”位置就是：


row 0..2 的 col 0（a0）加 $p\cdot h_i$


row 3..5 的 col 1（b）加 $p\cdot h_i$

画成表大概这样（每格是一个多项式）：
$\begin{array}{c|cc}
 & \text{col 0} & \text{col 1} \\
\hline
\text{row }0\ (i=0,k=0) & a_{0,0}+p h_0 & b_0 \\
\text{row }1\ (i=1,k=0) & a_{1,0}+p h_1 & b_1 \\
\text{row }2\ (i=2,k=0) & a_{2,0}+p h_2 & b_2 \\
\hline
\text{row }3\ (i=0,k=1) & a_{3,0} & b_3+p h_0 \\
\text{row }4\ (i=1,k=1) & a_{4,0} & b_4+p h_1 \\
\text{row }5\ (i=2,k=1) & a_{5,0} & b_5+p h_2 \\
\end{array}$

让加密的 $p$ 能作为一个“乘法器”，去乘另一个 TRLWE 密文，而噪声增长可控。

这就需要一种结构，使得：
$\text{ExternalProduct}(\text{TRGSW}(p),\ \text{TRLWE}(m)) \approx \text{TRLWE}(p\cdot m)$
4.1 先引入 gadget matrix $G$
用 gadget 向量 $h$，构造块对角 gadget 矩阵：
$\mathbf{G} = \mathrm{diag}(h,h,\dots,h) \quad ((k+1)\text{ 次})$
也就是说：每个列 $k$ 有 l 行，对应 $h_0,\dots,h_{l-1}$。
那么对任意 TRLWE 向量 $\mathbf{c}\in R^{k+1}$，都有一种分解（digits）：
$\mathrm{Decomp}(\mathbf{c}) \in R^{(k+1)l}
\quad \text{s.t.}\quad
\mathrm{Decomp}(\mathbf{c})\cdot \mathbf{G} \approx \mathbf{c}$
这就是 TFHE 里 gadget decomposition 的本质。

4.2 目标：让 TRGSW(p) 满足 “Enc(0)+pG”
对任意 $\mathbf{c}$：
$\mathrm{Decomp}(\mathbf{c})\cdot \mathbf{C}
\approx \mathrm{Enc}(0) + p\cdot \mathbf{c}$
如果 $\mathbf{C}$ 取成：
$\mathbf{C}=\mathrm{Enc}(0)+p\mathbf{G}$
那么：
$\mathrm{Decomp}(\mathbf{c})\cdot \mathbf{C}
=
\underbrace{\mathrm{Decomp}(\mathbf{c})\cdot \mathrm{Enc}(0)}_{\text{仍是 Enc(0)，线性组合}}
+
p\cdot \underbrace{\mathrm{Decomp}(\mathbf{c})\cdot \mathbf{G}}_{\approx \mathbf{c}}
\approx \mathrm{Enc}(0) + p\cdot \mathbf{c}$
而 $p\cdot \mathbf{c}$ 对 TRLWE 来说就是把每个分量多项式都乘上 $p$，解密时：
$(p\cdot b)-\sum (p\cdot a_t)s_t
=
p\cdot (b-\sum a_t s_t)
\approx p\cdot m$
于是外积就实现了“密文×明文”的效果，但这里的 $p$ 本身又是加密的（在 $\mathbf{C}$ 里），这就把乘法提升成了可同态的一步。
因为 $\mathbf{G}$ 是块对角矩阵：每个列 $k$ 的 gadget 行只在那一列有 $h_i$，其它列是 0。
把 $p\mathbf{G}$ 加进去，就意味着：


第 $k$ 个列对应的那组 l 行，在第 $k$ 列出现 $p h_i$

其它位置都是 0
代码用 trgsw[i + k*l][k] += p*h[i] 逐项把这块对角线构造出来。
它在生成一个 TRGSW/RGSW 密文，它不是“单条 RLWE 密文”，而是一个“GSW 矩阵密文”，等价于：
$\boxed{
\mathbf{C}=\mathrm{Enc}(0)\;+\;p\mathbf{G}
}$
其中 $\mathrm{Enc}(0)$ 由很多个 trlweSymEncryptZero(α,key) 组成，$p\mathbf{G}$ 通过对角位置加 p[j]*h[i] 构造。




取 $P::k=1$ ⇒ 一条 TRLWE 有 $k+1=2$ 个分量：col 0 和 col 1


col 0 = $a_0(X)$


col 1 = $b(X)$


取 $P::l=3$ ⇒ gadget 层数 $l=3$：$h_0,h_1,h_2$


TRGSW 的“矩阵”尺寸就是：行数 $(k+1)l = 2\cdot 3 = 6$，列数 $k+1=2$


行号正是：
$\text{row} = i + \kappa\cdot l
\quad
(i=0..l-1,\;\kappa=0..k)$

0) 


每一行是一条 TRLWE 向量：$(a,\;b)$


每一列是 TRLWE 的一个分量（不是行列式那种线代列）



1) 初始化阶段：先把整个 TRGSW 变成「Enc(0) 的堆叠」
   代码对应：
   for (TRLWE<P> &trlwe : trgsw)
    trlwe = trlweSymEncryptZero<P>(α, key);

也就是：6 行里每一行都是一条 TRLWE(0)。
对于 $k=1$ 的 TRLWE(0)，
$b_r(X)=a_{r,0}(X)\cdot s(X)\;+\;e_r(X)$
于是初始化后的“Enc(0) 矩阵” $E$ 可以画成：
$E=\text{Enc}(0)=
\begin{array}{c|cc}
 & \text{col 0} & \text{col 1} \\
\hline
\text{row }0 & a_{0,0} & b_0 \\
\text{row }1 & a_{1,0} & b_1 \\
\text{row }2 & a_{2,0} & b_2 \\
\hline
\text{row }3 & a_{3,0} & b_3 \\
\text{row }4 & a_{4,0} & b_4 \\
\text{row }5 & a_{5,0} & b_5 \\
\end{array}$
并且这些 $b_r$ 都满足：
$b_r=a_{r,0}\cdot s + e_r
\quad(\text{这里每个 }a_{r,0},b_r,e_r\text{ 都是多项式})$

先用噪声参数 $\alpha$ 和密钥 key 生成一堆“加密 0”的 TRLWE 行，堆成一个大矩阵。


2) gadget 向量 $h$
  代码对应：
  constexpr std::array<typename P::T, P::l> h = hgen<P>();
  $h=\big(h_0,\;h_1,\;h_2\big)$
  它的意义是“分解基底权重”：用来让外积（external product）里做 gadget decomposition：
  $x \approx \sum_{i=0}^{l-1} x^{(i)}\cdot h_i$

3) G 矩阵
在这个实现/存储方式下，$G$ 就是一个 $(k+1)l\times (k+1)$ 的“块对角 gadget 矩阵”。
对我们这个例子（2 列，6 行），它是：
$G=
\begin{array}{c|cc}
 & \text{col 0} & \text{col 1} \\
\hline
\text{row }0\ (i=0,\kappa=0) & h_0 & 0 \\
\text{row }1\ (i=1,\kappa=0) & h_1 & 0 \\
\text{row }2\ (i=2,\kappa=0) & h_2 & 0 \\
\hline
\text{row }3\ (i=0,\kappa=1) & 0 & h_0 \\
\text{row }4\ (i=1,\kappa=1) & 0 & h_1 \\
\text{row }5\ (i=2,\kappa=1) & 0 & h_2 \\
\end{array}$
看重点：


前 3 行（$\kappa=0$ 块）只在 col 0 放 $h_i$


后 3 行（$\kappa=1$ 块）只在 col 1 放 $h_i$


这就是标准的 block diagonal（块对角） 结构


用一个一句话概括：
$G_{\text{row}=i+\kappa l,\;\text{col}=t}=
\begin{cases}
h_i,& t=\kappa\\
0,& t\neq \kappa
\end{cases}$

4) 填充阶段：把明文 $p$ 按 $pG$ 的样子填进去
trgsw[i + k * P::l][k][j] += p[j] * h[i];

把它提升到“多项式层面”就是：
$\text{在 }(\text{row}=i+\kappa l,\;\text{col}=\kappa)\text{ 这个格子里，加上 }p(X)\cdot h_i$
所以明文嵌入矩阵 $pG$ 画出来就是：
$pG=
\begin{array}{c|cc}
 & \text{col 0} & \text{col 1} \\
\hline
\text{row }0\ (i=0,\kappa=0) & p h_0 & 0 \\
\text{row }1\ (i=1,\kappa=0) & p h_1 & 0 \\
\text{row }2\ (i=2,\kappa=0) & p h_2 & 0 \\
\hline
\text{row }3\ (i=0,\kappa=1) & 0 & p h_0 \\
\text{row }4\ (i=1,\kappa=1) & 0 & p h_1 \\
\text{row }5\ (i=2,\kappa=1) & 0 & p h_2 \\
\end{array}$

5) 最终“加密”结果：$\boxed{C=E+pG}$
  这就是 RGSW / TRGSW 的标准结构：
  $\boxed{
  C=\text{Enc}(0)+pG
  }$
  把两张表直接相加，就得到最终密文矩阵 $C$：
  $C=
  \begin{array}{c|cc}
   & \text{col 0} & \text{col 1} \\
  \hline
  \text{row }0\ (i=0,\kappa=0) & a_{0,0}+p h_0 & b_0 \\
  \text{row }1\ (i=1,\kappa=0) & a_{1,0}+p h_1 & b_1 \\
  \text{row }2\ (i=2,\kappa=0) & a_{2,0}+p h_2 & b_2 \\
  \hline
  \text{row }3\ (i=0,\kappa=1) & a_{3,0} & b_3+p h_0 \\
  \text{row }4\ (i=1,\kappa=1) & a_{4,0} & b_4+p h_1 \\
  \text{row }5\ (i=2,\kappa=1) & a_{5,0} & b_5+p h_2 \\
  \end{array}$
  这张就是你那张“最终效果图”（我补全了它在理论里的解释：它其实就是 $E+pG$）。

---


维度 $K=P::k=1$


所以一条 TRLWE 是 $(a,b)$ 两列（col 0 是 $a$，col 1 是 $b$）


gadget 层数 $l=3$，gadget 向量 $h=(h_0,h_1,h_2)$


TRGSW 行数 $(K+1)l = 2\cdot 3=6$



1) 先把“输入 TRLWE 密文”画出来：$c=\mathrm{TRLWE}(m)$
设有一条 TRLWE 密文 $c$（两列）：
$c=
\begin{array}{c|cc}
 & \text{col 0} & \text{col 1}\\
\hline
\text{(single row)} & a & b
\end{array}$
它满足（对称加密的典型关系）：
$b = a\cdot s + m + e$
所以解密时：
$b-a\cdot s \approx m$

2) gadget decomposition：把 $a,b$ 都分解成“digits × $h_i$”（这一步决定了外积的“权重”）
external product 的第一步不是直接拿 $a,b$ 去乘 TRGSW，而是要先把它们做 gadget 分解。
2.1 把每个分量（多项式）写成 $h_i$ 的线性组合
$a \approx a^{(0)}h_0 + a^{(1)}h_1 + a^{(2)}h_2$
$b \approx b^{(0)}h_0 + b^{(1)}h_1 + b^{(2)}h_2$
这里：


$a^{(i)}, b^{(i)}$ 是“digits”（通常系数都比较小，来自某个基 $B_g$ 的取整/舍入分解）


“$\approx$”是因为分解会有舍入误差（在 torus/定点表示里很常见）



2.2 把这些 digits 按 TRGSW 的行顺序排成一个长度 6 的向量（这就是 Decomp）
你的 TRGSW 行顺序是 row = i + κ*l，所以自然的分解向量就是：
$\mathrm{Decomp}(c)=
\begin{array}{c|c}
\text{row idx} & d_{\text{row}} \\
\hline
0\ (i=0,\kappa=0) & a^{(0)} \\
1\ (i=1,\kappa=0) & a^{(1)} \\
2\ (i=2,\kappa=0) & a^{(2)} \\
\hline
3\ (i=0,\kappa=1) & b^{(0)} \\
4\ (i=1,\kappa=1) & b^{(1)} \\
5\ (i=2,\kappa=1) & b^{(2)} \\
\end{array}$
external product 代码里：

for row: out += d[row] * trgsw[row];

对应的就是这张表的“row idx”。

3. $G$（6×2）还是这个块对角
   $G=
   \begin{array}{c|cc}
    & \text{col 0} & \text{col 1} \\
   \hline
   \text{row }0 & h_0 & 0 \\
   \text{row }1 & h_1 & 0 \\
   \text{row }2 & h_2 & 0 \\
   \hline
   \text{row }3 & 0 & h_0 \\
   \text{row }4 & 0 & h_1 \\
   \text{row }5 & 0 & h_2 \\
   \end{array}$
   3.2 把每行乘上对应 digit，再把 6 行“加起来”，你就得到 $(a,b)$
   $\begin{array}{c|cc}
    & \text{col 0} & \text{col 1} \\
   \hline
   a^{(0)}\times & a^{(0)}h_0 & 0 \\
   a^{(1)}\times & a^{(1)}h_1 & 0 \\
   a^{(2)}\times & a^{(2)}h_2 & 0 \\
   \hline
   b^{(0)}\times & 0 & b^{(0)}h_0 \\
   b^{(1)}\times & 0 & b^{(1)}h_1 \\
   b^{(2)}\times & 0 & b^{(2)}h_2 \\
   \hline
   \text{sum} & \sum_i a^{(i)}h_i & \sum_i b^{(i)}h_i
   \end{array}$
   因此：
   $\big(\sum_i a^{(i)}h_i,\ \sum_i b^{(i)}h_i\big)\approx (a,b)=c$
   external product：$\boxed{\mathrm{ExtProd}(\mathrm{TRGSW}(p),\mathrm{TRLWE}(m))\approx \mathrm{TRLWE}(p\cdot m)}$
   现在进入你最关心的“外积到底怎么做”。
   4.1 先把 TRGSW(p)（你前面那张最终表）搬过来
   $C=\mathrm{TRGSW}(p)=
   \begin{array}{c|cc}
    & \text{col 0} & \text{col 1} \\
   \hline
   \text{row }0 & a_{0,0}+p h_0 & b_0 \\
   \text{row }1 & a_{1,0}+p h_1 & b_1 \\
   \text{row }2 & a_{2,0}+p h_2 & b_2 \\
   \hline
   \text{row }3 & a_{3,0} & b_3+p h_0 \\
   \text{row }4 & a_{4,0} & b_4+p h_1 \\
   \text{row }5 & a_{5,0} & b_5+p h_2 \\
   \end{array}$
   你可以把它记成：
   $C = E + pG$
   其中 $E=\mathrm{Enc}(0)$ 是那堆 TRLWE(0) 的行堆叠。

4.2 external product 的“图形化算法”：用 Decomp 的 digits 当权重，把 TRGSW 的行加权求和
这一步就是：每一行乘一个 digit，然后把 6 行加起来，结果是一条新的 TRLWE（还是两列）。
把它画成一个“乘权重并求和”的表最直观：
$\mathrm{ExtProd}(C,c)=
\begin{array}{c|cc}
 & \text{col 0} & \text{col 1} \\
\hline
a^{(0)}\times & a^{(0)}(a_{0,0}+p h_0) & a^{(0)}b_0 \\
a^{(1)}\times & a^{(1)}(a_{1,0}+p h_1) & a^{(1)}b_1 \\
a^{(2)}\times & a^{(2)}(a_{2,0}+p h_2) & a^{(2)}b_2 \\
\hline
b^{(0)}\times & b^{(0)}a_{3,0} & b^{(0)}(b_3+p h_0) \\
b^{(1)}\times & b^{(1)}a_{4,0} & b^{(1)}(b_4+p h_1) \\
b^{(2)}\times & b^{(2)}a_{5,0} & b^{(2)}(b_5+p h_2) \\
\hline
\text{sum} & A_{\text{out}} & B_{\text{out}}
\end{array}$
因此输出是一条 TRLWE：
$\mathrm{out}=
\begin{array}{c|cc}
 & \text{col 0} & \text{col 1}\\
\hline
\text{(single row)} & A_{\text{out}} & B_{\text{out}}
\end{array}$

4.3 把 $A_{\text{out}},B_{\text{out}}$ 分成“噪声块” + “$p\cdot(a,b)$”
把上面 sum 展开，先写 $A_{\text{out}}$：
$A_{\text{out}}
=
\underbrace{\sum_{i} a^{(i)}a_{i,0} + \sum_{i} b^{(i)}a_{i+3,0}}_{\text{这堆来自 }E=\mathrm{Enc}(0)}
\;+\;
p\cdot\underbrace{\sum_i a^{(i)}h_i}_{\approx a}$
同理 $B_{\text{out}}$：
$B_{\text{out}}
=
\underbrace{\sum_i a^{(i)}b_i + \sum_i b^{(i)}b_{i+3}}_{\text{这堆来自 }E=\mathrm{Enc}(0)}
\;+\;
p\cdot\underbrace{\sum_i b^{(i)}h_i}_{\approx b}$
于是得到：
$(A_{\text{out}},B_{\text{out}})
\approx
(\underbrace{\star,\star}_{\text{仍然是 Enc(0) 的线性组合}})
\;+\;
p\cdot(a,b)$
也就是一句非常关键的等式（GSW 的灵魂）：
$\boxed{
\mathrm{Decomp}(c)\cdot (E+pG)
=
\underbrace{\mathrm{Decomp}(c)\cdot E}_{\text{Enc(0) 仍是 Enc(0)}}
+
p\cdot \underbrace{\mathrm{Decomp}(c)\cdot G}_{\approx c}
\approx
\mathrm{Enc}(0) + p\cdot c
}$

4.4 为什么这就等于“加密了 $p\cdot m$”？
因为如果原来：
$c=(a,b),\quad b=a s+m+e$
那么 $p\cdot c=(pa,pb)$ 解密时：
$pb-(pa)s = p(bs-as) \approx p(m+e) = pm + p e$
所以 external product 的输出就是一条新的 TRLWE，明文近似是：
$\boxed{pm}$
噪声来源主要两块：


$\mathrm{Decomp}(c)\cdot E$：把很多条 Enc(0) 做线性组合得到的新噪声


分解误差：$\mathrm{Decomp}(c)\cdot G \approx c$ 不是严格等号，会多一点误差项，再乘上 $p$


这也是为什么 digits 要“够小”、$l$ 要“够深”：噪声可控。

---

设：


numeric_limits<...>::digits = w = 32


l = 3


bg = 6（所以每个 digit 是 6 bit）


于是每个 digit 的基数 B = 2^bg = 64


3 个 digit 一共占 l*bg = 18 bit


所以要丢掉的低位 bit 数：
$t = w - l\cdot bg = 32 - 3\cdot 6 = 14$


roundoffset = 2^{t-1} = 2^{13} = 8192

1) 对任意 32-bit 值 x，把它写成：
$x = Q\cdot 2^{14} + R,\quad 0 \le R < 2^{14}$


Q：就是右移 14 位后留下的“高 18 bit”（将来会被拆成 3 个 base-64 digit）


R：就是被你丢掉的低 14 bit（残差）

**截断（不四舍五入）**就是：
$Q_{\text{trunc}} = x >> 14 = Q$
四舍五入就是（先加半单位再右移）：
$Q_{\text{round}} = (x + 2^{13}) >> 14$
所以它只做一件事：根据残差 $R$ 是否 ≥ 8192，决定要不要给 $Q$ 加 1：


若 $R < 8192$：不加 1（round down）


若 $R \ge 8192$：加 1（round up）


而“加 1”发生在 拆 digit 之前，所以这个 +1 的进位会在随后拆 digit 时体现出来，并且会通过 carry 传到更高 digit。

2) 用一个能看出“carry 贯穿多个 digit”的例子
我们挑一个 $Q$ 的 3 个 digit（base 64）是：


$d_2=10$, $d_1=63$, $d_0=63$


写成 6-bit 一组就是：


$10 = 001010$


$63 = 111111$


$63 = 111111$


所以 Q 的 18 bit 是：
$Q = \underbrace{001010}_{d_2}\ \underbrace{111111}_{d_1}\ \underbrace{111111}_{d_0}$
这在十进制里是：
$Q = 10\cdot 64^2 + 63\cdot 64 + 63 = 45055$
然后构造 32-bit 的 x：
$x = Q\cdot 2^{14} + R$
也就是（概念上）：
$x = [\text{高18bit = }Q]\ [\text{低14bit = }R]$
情况 A：残差不够半单位（不会进位）
取 $R=7000$，因为 $7000 < 8192$。


截断：$Q_{\text{trunc}} = Q = 45055$


四舍五入：$(x+8192)>>14$ 仍然不会让高位多 1，所以
$Q_{\text{round}} = 45055$


拆成 digit（base64）还是：
$(10,\ 63,\ 63)$
情况 B：残差达到半单位（触发 round up，并且 carry 连锁）
取 $R=9000$，因为 $9000 \ge 8192$。
这时：
$Q_{\text{round}} = Q + 1 = 45056$
关键来了：看这个 +1 加在 18-bit 的最低位上，会发生什么：
$001010\ 111111\ 111111 \;+\; 1
=
001011\ 000000\ 000000$
你看到没？


最低 digit：111111 + 1 -> 000000，产生 carry


中间 digit：111111 + carry -> 000000，再产生 carry


最高 digit：001010 + carry -> 001011


所以拆成 base64 digit 变成：
$(11,\ 0,\ 0)$

所以整体效果就是：对“截断后的整体值”做四舍五入，而不是每个 digit 各自独立截断。



