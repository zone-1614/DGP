# 作业4
* 作业4是参数化算法, Tutte's embedding method
* 自己实现了一个LSCM方法
> 实现LSCM的感想: 这个算法不论是理解还是实现对我来说都挺难的. 算法上的难度在于把问题建模成复数的形式, 最后又回到实数的形式. 并且这篇论文是2002年的SIGGRAPH, 他对柯西黎曼方程的引入和我学的不一样. 
> 再感慨一下, 这篇LSCM的论文包括的东西真多, 太牛了 (不过感觉符号挺乱的?)

## Tutte
这个是从图论来的定理, 我的数学知识不够, 没法证明. 其结论应用到参数化中就是: 把边界点映射到凸多边形上, 就一定能得到满足以下三点的参数化结果
* flip-free
* bijective
* local injective

## LSCM
LSCM基于*共形映射*(conformal mapping, 也叫保角映射), 共形映射要求函数是*全纯*(holomorphic)的, 也就是需要满足*Cauchy-Riemann*方程:
$$
\frac{\partial u}{\partial x} = \frac{\partial v}{\partial y}, \frac{\partial u}{\partial y} = -\frac{\partial v}{\partial x}
$$

## 基本理论
下面是探讨用矩阵的形式来表示共形映射. 首先抛开这些名词, 保角也就是相似(就是中学学的相似三角形), 所以这个变换相当于对三角形做scale和rotation, 所以矩阵表示应该如下:
$$
s\begin{bmatrix}
    \cos \theta && -\sin \theta \\
    \sin \theta && \cos \theta
\end{bmatrix} = 
\begin{bmatrix}
    a && -b \\
    b && a
\end{bmatrix}
$$ 
### 问题建模
三角形是二维的图形, 这里是在每个三角形面片上建立局部坐标系

接下来就是老套路了, 网格上没办法所有顶点都满足这些约束, 所以只能尽可能的去满足这些约束, 所以问题转化为一个优化问题. 再通过某些技巧转化为线性优化的问题, 从而建模成一个稀疏矩阵, 这样就可以求解了.

$$
E_{LSCM} = \sum_t A_t ((\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y})^2 + (\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x})^2) 
$$

### 具体计算
在网格上计算上述的几个微分量:
$$
\begin{bmatrix}
    \frac{\partial u}{\partial x} \\
    \frac{\partial u}{\partial y}
\end{bmatrix}
= \frac{1}{2A} \begin{bmatrix}
    y_j - y_k && y_k - y_i && y_i - y_j \\
    x_k - x_j && x_i - x_k && x_j - x_i 
\end{bmatrix}
\begin{bmatrix}
    u_i \\ u_j \\ u_k
\end{bmatrix} \\
\begin{bmatrix}
    \frac{\partial v}{\partial x} \\
    \frac{\partial v}{\partial y}
\end{bmatrix}
= \frac{1}{2A} \begin{bmatrix}
    y_j - y_k && y_k - y_i && y_i - y_j \\
    x_k - x_j && x_i - x_k && x_j - x_i 
\end{bmatrix}
\begin{bmatrix}
    v_i \\ v_j \\ v_k
\end{bmatrix}
$$

再把上述形式改为复数形式:
$$
\frac{\partial u}{\partial x} + i\frac{\partial u}{\partial y}
= \frac{i}{2A} \begin{bmatrix}
    W_i && W_j && W_k    
\end{bmatrix}
\begin{bmatrix}
    u_i \\ u_j \\ u_k
\end{bmatrix} \\
W_i = (x_k-x_j) + i(y_k-y_j), W_j = (x_i-x_k) + i(y_i-y_k), W_k = (x_j - x_i) + i(y_j-y_i)
$$
因为三角形面片上的导数是常数, 并且是两个点相减的形式, 所以这一步之后局部坐标系的影响就被消除了
每个三角形的能量为
$$
E(t) = \frac{1}{2A_t}|(W_{j1,t}, W_{j2,t}, W_{j3, t}) \begin{bmatrix}U_{j1}\\ U_{j2} \\ U_{j3}\end{bmatrix}|^2
$$
用T表示所有三角形集合, 则整个网格的能量为
$$
E_{LSCM}=\sum_{t\in T} E(t)
$$
因为$E(t)$是用复数表示的, 所以可以做以下变换
$$
E(t) = \frac{1}{2A_t} U^TW^TWU =  U^T C U, \\
\text{C可以表示为: }C = \frac{1}{2A_t}W^TW \\
\text{再把系数塞进去: }C = M^TM \\
\text{其中, } M = (m_{ij}) \\
m_{ij} = \left\{\begin{aligned} \frac{W_{j, t_i}}{\sqrt{2A_{t_i}}} ,\ j属于三角形t_i \\ 0, \ j不属于三角形t_i \end{aligned} \right.
$$
这样就得到了矩阵M, 矩阵M有 **nf**(number of faces) 行, **nv**(number of vertices) 列, 为了避免平凡解, 一般会固定两个点. 在矩阵M中, 把这两个固定的点所对应的列移动到最后两列, 也就是
$$
M = (M_f, M_p)
$$
其中$M_p$就是固定点对应的那两列. (下标f表示free, 下标p表示pinned)

所以可以把能量表示为
$$
E_{LSCM} = ||M_fU_f + M_pU_p||^2
$$
其中$U_f, U_p$分别是free和pinned的顶点坐标.
最后, 把复矩阵变为实矩阵可得:
$$
E_{LSCM} = ||Ax-b||^2
$$
其中$x$就是待求解的坐标, 三个量的表达式为:
$$
A = \begin{bmatrix}
    M_f^r & -M_f^i \\
    M_f^i & M_f^r
\end{bmatrix}, 
x = \begin{bmatrix}
u_1 \\
u_2 \\
\vdots \\
u_n \\
v_1 \\
v_2 \\
\vdots \\
v_n
\end{bmatrix},
b = -\begin{bmatrix}
    M_p^r & -M_p^i \\
    M_p^i & -M_p^r
\end{bmatrix}\begin{bmatrix}
    U_p^r \\
    U_p^i
\end{bmatrix}
$$
其中, 上标r表示real, i表示imag