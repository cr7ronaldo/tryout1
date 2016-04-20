# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 20:02:01 2016

@author: DP BABAR
"""
"""Solve the 2D biharmonic equation DW4/Dx4+2DW4/Dx2Dy2+DW4/Dy4=f(x)
       subject to the boundary conditions 
       u(-1)=u'(-1)=u(+1)=u'(+1)=0
       Inputs:
       n - number of grid points
       f - forcing function
       
       Outputs:
       x - grid points
       u - solution at grid points """
import numpy as np
from scipy.linalg import solveh_banded
import matplotlib.pyplot as plt

def plate4(n,ffun,a):

    
    x = np.linspace(0,1,n+1)
    h = 1.0/n
    
    # Allocate bands of the symmetric matrix
    stencil = np.array((1,-4,6))
    B = np.outer(stencil,np.ones(n))
    B[1,-1] = -2;    B[2,0] = 7;    B[2,-1] = 1;    B[2,-2] = 5    
    
    f = ffun(x)
    f *= h**4;    f[-1] *= 0.5;    f[-1] -= a*h**3
    
    u = np.zeros(n+1)    
    u[1:] = solveh_banded(B,f[1:],lower=False)
    
    return x,u

if __name__ == '__main__':
    
    n = 100    
    
    a = -(np.pi/2)**3
    
    def ffun(y):
        return  -(np.pi/2)**4*np.cos(np.pi*y/2)

    x,u = plate4(n,ffun,a)        

    # Exact solution    
    uex = 1-np.cos(np.pi*x/2)
    plt.plot(x,u)
    plt.plot(x,uex)
    
    plt.show()