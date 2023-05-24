# WaveEequation
这里尝试使用python搭建一个计算一维波动方程的差分方法发求解器

现在尝试搭建一个计算线性双曲型偏微分方程的差分求解器

方程：

$u_{tt}-a^{2}\Delta u=f(x,t)$

这里考虑一维的情况：

$u_{tt}-a^2 u_{xx}=f(x,t) x\in (0,1)$

$u(x,0)=\varphi (x),\frac{\partial u}{\partial t}(x,0)=\psi(x),x\in [0,1]$

$u(0,t)=\alpha (t),u(1,t)=\beta(t)$

为此我们构造其差分格式，记

$\Omega=|\{(x,t)|x\in [0,1],t\in [0,T]\}|$

Step1 使用空间步长$\Delta x=\frac{1}{N}$，时间步长为$\Delta t=\frac{T}{M}$的网格线划分$\Omega$得到网格

$\mathbb{J}=\{(x_i,t_j)|(x_i,t_j)\in\bar{\Omega}\}$

$\mathbb{J}$上$u,f,\varphi,\psi,\alpha,\beta$的取值分别记为$u_{ij},f_{ij},\varphi_{ij},\psi_{ij},\alpha_{ij},\beta_{ij}$

Step2 采用中心差分公式：

$\frac{\partial^2 u}{\partial x^2}=\frac{u_{i-1,j}-2u_{i,j}+u_{i+1,j}}{\Delta x^2}$

$\frac{\partial^2 u}{\partial t^2}=\frac{u_{i,j-1}-2u_{i,j}+u_{i,j+1}}{\Delta t^2}$

由此得到波动方程的差分格式：

$u_{i,j+1}=\frac{a^2 \Delta t^2}{\Delta x^2}(u_{i+1,j}-2u_{ij}+u_{i-1,j})-u_{i,j-1}+2u_{ij}+\Delta t^2 f_{ij}$

$u_{i,0}=\varphi_{i},u_{i,1}=\varphi_{i}+\Delta t\psi_{i}$

$u_{0,j}=\alpha_{j},u_{N,j}=\beta_{j}$
