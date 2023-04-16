# 作业4
* 作业4是参数化算法, Tutte's embedding method
* 自己实现了一个LSCM方法
> 实现LSCM的感想: 这个算法不论是理解还是实现对我来说都挺难的. 算法上的难度在于把问题建模成复数的形式, 最后又回到实数的形式. 并且这篇论文是2002年的SIGGRAPH, 他对柯西黎曼方程的引入和我学的不一样. 
> 让我比较惊讶的是他这篇论文不仅包含了参数化的算法, 还包含了segmentation, packing的内容

## LSCM
LSCM基于*共形映射*(conformal mapping, 也叫保角映射), 共形映射要求函数是*全纯*(holomorphic)的, 也就是需要满足*Cauchy-Riemann*方程:
$$
\frac{\partial u}{\partial x} = \frac{\partial v}{\partial y}, \frac{\partial u}{\partial y} = -\frac{\partial v}{\partial x}
$$

## 基本理论
下面是探讨用矩阵的形式来表示共形映射. 首先抛开这些名词, 保角也就是相似(就是中学学的相似三角形), 所以这个变换相当于对三角形做scale和rotation, 所以矩阵表示应该如下:
$$
s\begin{pmatrix}
    \cos \theta && -\sin \theta \\
    \sin \theta && \cos \theta
\end{pmatrix} = 
\begin{pmatrix}
    a && -b \\
    b && a
\end{pmatrix}
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
\begin{pmatrix}
    \frac{\partial u}{\partial x} \\
    \frac{\partial u}{\partial y}
\end{pmatrix}
= \frac{1}{2A} \begin{pmatrix}
    y_j - y_k && y_k - y_i && y_i - y_j \\
    x_k - x_j && x_i - x_k && x_j - x_i 
\end{pmatrix}
\begin{pmatrix}
    u_i \\ u_j \\ u_k
\end{pmatrix} \\
\begin{pmatrix}
    \frac{\partial v}{\partial x} \\
    \frac{\partial v}{\partial y}
\end{pmatrix}
= \frac{1}{2A} \begin{pmatrix}
    y_j - y_k && y_k - y_i && y_i - y_j \\
    x_k - x_j && x_i - x_k && x_j - x_i 
\end{pmatrix}
\begin{pmatrix}
    v_i \\ v_j \\ v_k
\end{pmatrix}
$$

再把上述形式改为复数形式:
$$
\frac{\partial u}{\partial x} + i\frac{\partial u}{\partial y}
= \frac{i}{2A} \begin{pmatrix}
    W_i && W_j && W_k    
\end{pmatrix}
\begin{pmatrix}
    u_i \\ u_j \\ u_k
\end{pmatrix} \\
W_i = (x_k-x_j) + i(y_k-y_j), W_j = (x_i-x_k) + i(y_i-y_k), W_k = (x_j - x_i) + i(y_j-y_i)
$$
因为三角形面片上的导数是常数, 并且是两个点相减的形式, 所以这一步之后局部坐标系的影响就被消除了
每个三角形的能量为
$$

$$