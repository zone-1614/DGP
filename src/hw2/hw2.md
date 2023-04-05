# 作业2
在离散曲面上计算微分量, 从而把微分几何的那套方法搬到三角网格上

计算并可视化平均曲率, 高斯曲率

$$
H_i = 0.5 * ||\Delta x|| \\
K_i = \frac{1}{A_i}(2\pi - \sum_{j\in \Omega(i)} \theta_j)
$$

涉及到的内容: 离散拉普拉斯算子(cot 权), local averaging region(mixed Voronoi 方法)