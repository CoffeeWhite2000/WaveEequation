'''
现在尝试搭建一个计算线性双曲型偏微分方程的差分求解器
方程：
$u_{tt}-a^{2}\Detla u=f(x,t)$
这里考虑一维的情况：
$u_{tt}-a^2 u_{xx}=f(x,t) x\in (0,1)$
$u(x,0)=\varphi (x),\frac{\partial u}{\partial t}(x,0)=\psi(x),x\in [0,1]$
$u(0,t)=\alpha (t),u(1,t)=\beta(t)$
为此我们构造其差分格式，记$\Omega=|\{(x,t)|0<x<1,0<t<T\}|$
Step1 使用空间步长$\Delta x=\frac{1}{N}$，时间步长为$\Delta t=\frac{T}{M}$的网格线划分$\Omega$得到网格$\mathbb{J}=\{(x_i,t_j)|(x_i,t_j)\in\bar{\Omega}\}$
    $\mathbb{J}$上$u,f,\varphi,\psi,\alpha,\beta$的取值分别记为$u_{ij},f_{ij},\varphi_{ij},\psi_{ij},\alpha_{ij},\beta_{ij}$
Step2 采用中心差分公式：
    $\frac{\partial^2 u}{\partial x^2}|_{(x_i,t_j)}=\frac{u_{i-1,j}-2u_{i,j}+u_{i+1,j}}{\Delta x^2}$
    $\frac{\partial^2 u}{\partial t^2}|_{(x_i,t_j)}=\frac{u_{i,j-1}-2u_{i,j}+u_{i,j+1}}{\Delta t^2}$
由此得到波动方程的差分格式：
    $u_{i,j+1}=\frac{a^2 \Delta t^2}{\Delta x^2}(u_{i+1,j}-2u_{ij}+u_{i-1,j})-u_{i,j-1}+2u_{ij}+\Delta t^2 f_{ij}$
    $u_{i,0}=\varphi_{i},u_{i,1}=\varphi_{i}+\Delta t\psi_{i}$
    $u_{0,j}=\alpha_{j},u_{N,j}=\beta_{j}$
'''





import numpy as np
import numpy.matlib
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def var(x):
    #return cos(pi*x)**2*cos(pi*y)**2
    return 0
def psi(x):
    return x
def f(x,t):
    #return sin(pi*x)*sin(pi*y)
    return (t**2-x**2)*sin(x*t)
def alp(t):
    return 0
def beta(t):
    return sin(t)

def get_f_vector(j,N,h,s):
    vector=np.matlib.zeros((N-1,1), dtype=float)
    for i in range(0,N-1):
        vector[i]=(s**2)*f((i+1)*h,j*s)
    return vector


l=1 #x轴长度
N=20 #x轴等分
T=1 #时间秒
M=100 #时间等分
a=1 #方程中的系数,波动系数
s=T/M
h=l/N
r=(a**2)*(s**2)/(h**2) # r必须小于0.5 此差分方程才稳定
print(r)
x=np.linspace(0,l,N+1, endpoint=True,dtype=float)#
t=np.linspace(0,T,M+1, endpoint=True,dtype=float)
u=np.matlib.zeros((N+1,M+1), dtype=float)
A=np.matlib.identity(N-1, dtype=float)
A=2*(1-r)*A
for i in range(0,N-2):
    A[i,i+1]=r
    A[i+1,i]=r

for i in range(0,N+1):
    u[i,0]=var(i*h)
    u[i,1]=var(i*h)+s*psi(i*h)
for j in range(0,M+1):
    u[0,j]=alp(j*s)
    u[N,j]=beta(j*s)
for k in range(2,M+1):
    u[1:N,k]=A*u[1:N,k-1]-u[1:N,k-2]+get_f_vector(k,N,h,s)
    u[1,k]=u[1,k]+r*u[0,k-1]
    u[N-1,k]=u[N-1,k]+r*u[N,k-1]


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
X, Y = np.meshgrid(x, t)
ax.plot_surface(X, Y,u.T , rstride=1, cstride=1, cmap='hot')
plt.show() 

