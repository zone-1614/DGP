# 作业5
作业5也是参数化, 代码一起放在[hw4.cpp](../hw4/hw4.cpp)里面了

论文: [A Local/Global Approach to Mesh Parameterization](./Gortler_LocalGlobal.pdf)

# 原理
> 这篇论文介绍的方法其实是在已经参数化完成的uv面上进行优化, 所以首先需要用别的方法进行参数化(我用了作业4的Tutte的结果, 所以代码写在[hw4.cpp](../hw4/hw4.cpp)里了).

## 符号约定
* 将mesh的三角形进行编号: $t = 1, 2, \cdots, T$
* 对于三角形t, 其面积$A_t$
* 原始坐标$\bold{x}_t=\{x_t^0, x_t^1, x_t^2\}$, 变化后的坐标$\bold{u}_t=\{u_t^0, u_t^1, u_t^2\}$. 这里的$x_t^i, u_t^i$都是2维列向量
* 对于三角形t, 当前参数化的Jacobian: $J_t(\bold{u})$
* 对每个三角形t, 定义辅助矩阵$L_t\in SO(2)$
* Frobenius范数: $||\cdot||_F$

## 定义能量
要把参数化的结果往$SO(2)$靠拢, 所以定义能量如下
$$
E(\bold{u}, \bold{L}) = \sum_{t=1}^{T}A_t||J_t(\bold{u}_t)-L_t||_F^2
$$
其中$\bold{u} = \{ \bold{u}_t | t = 1, 2,\cdots, T \}$是所有坐标的集合, $\bold{L}=\{L_t | t = 1, 2,\cdots, T\}$

上面这个能量可以用三角形的坐标显示的表达如下(请自动把上标i模2):
$$
E(\bold{u}, \bold{L}) = \frac{1}{2} \sum_{t=1}^{T} \sum_{i=0}^{2} \cot(\theta_t^i)||(u_t^i - u_t^{i+1}) - L_t (x_t^i-x_t^{i+1})||^2
$$
其中$\theta_t^i$是边$(x_t^i-x_t^{i+1})$的对角. 注意是$(x_t^i-x_t^{i+1})$的对角, 不是$(u_t^i-u_t^{i+1})$的对角, x在整个算法过程中是不变的!

## 问题求解
定义完能量之后, 问题就转化为如下
$$
\argmin_{\bold{u}, \bold{L}}E(\bold{u}, \bold{L}) \\
s.t. \ \ L_t \in SO(2)
$$
虽然同时求解出$\bold{u}, \bold{L}$, 不过我们只关心$\bold{u}$. $\bold{L}$只是求解过程的辅助变量

### 找L
$$
d(J, L) = ||J - L||_F^2 = tr ((J-L)^T(J-L))
$$
其中J是Jacobian, $L \in SO(2)$

对J用SVD, 有
$$
J = U\Sigma V^T \\ 
\Sigma = \begin{pmatrix}
    \sigma_1 & 0 \\
    0 & \sigma_2
\end{pmatrix}
$$
为了使得$L\in SO(2)$, 而不仅仅是$L\in O(2)$, 这里的SVD是signed SVD, 也就是
$$
\det (UV^T) = 1, \sigma_1 > 0 \\
\sigma_2 > 0 \ or \  \sigma_2 < 0
$$
然后, 令$L=UV^T$即可. 所以, 只要知道$\bold{u}$, 就可以求出J, 对J做SVD, 就可以得到L, 也就是说只需要$\bold{u}$这一个变量即可.
> 为什么这样赋值给L就行呢? 这里的思想是: 对于J这样一个可能性质比较差的矩阵, 我们从$SO(2)$中找一个比较接近他的, 而signed SVD分解后, $U, V\in SO(2)$, 而$\Sigma$代表一个scale的作用, 所以$UV^T$就是最接近J的那个L

在论文的Appendix中, 证明了能量可以表示为如下 
$$
\sum_{t=1}^{T} A_t [(\sigma_{1, t} - 1)^2 + (\sigma_{2, t} - 1)^2]
$$

## 算法
论文中对求解上诉问题使用了 Local/Global Algorithm

### Local phase
(固定J_t 求 L_t)
对每个三角形t, 利用$\bold{x}, \bold{u}$求出$J_t$, 对$J_t$做SVD $J_t=U\Sigma V^T$, 得$L_t=UV^T$

### Global phase
(固定$\bold{L}=\{L_t | t = 1, 2,\cdots, T\}$, 求 $\bold{u}$)
把能量用半边数据结构表示
$$
E(\bold{u}, \bold{L}) = \frac{1}{2} \sum_{(i, j) \in he} \cot(\theta_{ij}) || (u_i-u_j) - L_{t}(x_i - x_j) ||^2
$$
令$\nabla E=0$有
$$
\sum_{j\in N(i)} (\cot(\theta_{ij}) + \cot(\theta_{ji}))(u_i - u_j) = \sum_{j\in N(i)}(\cot(\theta_{ij}) L_t + \cot(\theta_{ji}) L_t) (x_i - x_j)
$$
