import numpy as np
import sympy as sp
#constants
m=1
omega=1
h=1
beta=5
L=10

def sinBasis(x,n):
    return np.sqrt(2/L)*np.sin(2*n*np.pi*(x-L/2)/L)

n=np.arange(1,50)
m=np.arange(1,50)
HMatrix=np.