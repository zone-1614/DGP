# Digital Geometry Processing

## 资料
* [网课](https://www.bilibili.com/video/BV1B54y1B7Uc)
* [课程主页](https://ustc-gcl-f.github.io/course/2020_Spring_DGP/index.html), 里面有PPT
* [我用的作业框架](https://github.com/pmp-library/pmp-library), [官方的作业框架](https://ustc-gcl-f.github.io/code/index.html#sec_surface_framework)

## 作业截图
### 作业1
作业1是在mesh上实现单源最短路径算法--dijkstra算法和最小生成树算法--prim算法, 没有做可视化, 所以没有图.

### 作业2
作业2是可视化平均曲率和高斯曲率

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw2_curvature_pig1.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw2_curvature_pig2.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw2_curvature_pig3.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw2_curvature_pig4.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw2_curvature_multicube1.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw2_curvature_multicube2.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw2_curvature_bunny1.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw2_curvature_bunny2.png)


### 作业3
作业3是曲面去噪, 本来设计了一个误差, 根据误差去调整迭代次数, 但是不同曲面这个误差不太一样, 所以把迭代次数和高斯函数都作为可以调的参数.

下面两张图是基本的去噪效果, 感觉效果还是不错的

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw3_denoising_1.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw3_denoising_3.png)


下面三张图展示的是不同参数对结果的影响. 去噪前:


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw3_denoising_4.png)


去噪后:

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw3_denoising_5.png)


修改参数后的效果:

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw3_denoising_6.png)


最后两张图展示的是对较大噪声的去噪效果, 去噪前:

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw3_denoising_7.png)


去噪后:

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw3_denoising_8.png)


### 作业4
Tutte's embedding method

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw4_00.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw4_01.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw4_10.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw4_11.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw4_20.png)


![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw4_21.png)

自己实现了LSCM, 下面分别是balls, bunny, Nefertiti_face的参数化结果

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw4_LSCM_balls.png)

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw4_LSCM_bunny.png)

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw4_LSCM_face.png)

## 作业5
作业5是ARAP参数化, 因为原理还没完全搞明白, 所以写了好久. 下面分别是balls, bunny, Nefertiti_face的参数化结果

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw5_ARAP_0.png)

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw5_ARAP_1.png)

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw5_ARAP_2.png)

好久没写DGP的作业, 这段时间主要去复习考研和提升一些工程能力了. 回过头来看这个ARAP应该是还有bug, 效果不太对, 等以后再回来改改

## 作业6
作业6是ARAP deformation

自己用blender做的模型, 变形前:

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw6_ARAP_deformation6.png)

变形后:

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw6_ARAP_deformation7.png)


spot模型, 变形前: 

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw6_ARAP_deformation0.png)

变形后(迈步):

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw6_ARAP_deformation4.png)

抬头: 

![](https://raw.githubusercontent.com/zone-1614/pic/main/img/hw6_ARAP_deformation5.png)