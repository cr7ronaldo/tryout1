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
       u - solution at grid points
       To test the code out, we can choose an f(x) and where we know the exact solution. 
       In particular if f(x)=−(π2)4cos(π2x)f(x)=−(π2)4cos⁡(π2x) and a=−(π2)3a=−(π2)3, 
       then the exact solution is u(x)=1−cos(π2x)u(x)=1−cos⁡(π2x).
       """
       
#finite difference scheme
from numpy import matrix
import math
def f(x):
    return math.sinh(x)
w0=0    #y(0)=0
w4=3.62686         #y(2)=3.62686
h1=0.5
W=f(0.5)
#arranging and calculating the values
A=matrix([[-9, 4, 0],[4, -9, 4],[0, 4, -9]])
C=matrix([[0],[0],[-14.50744]])
X=A.I*C
print "computed value with h=%f of y(0.5) is %f\n" %(h1,X[0][0])
print "error in the result with actual value %f\n" %(abs(W-X[0][0]))
h2=1.0
w0=0    #y(0)=0
w2=3.62686   #y(2)=3.62686
w1=(w0+w2)/3
W=(4*X[1][0]-w1)/3
print "with better approximation error is reduced to %f" %(abs(W-f(1.0)))       

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