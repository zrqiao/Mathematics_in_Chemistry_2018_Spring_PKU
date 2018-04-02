import numpy as np
import math
raw_hessian = open('hess.dat')
mass = list[1.0008, 15.9997, 1.0008]
hess = list(map(lambda line: list(map(float,
    [num for num in line.split()])), raw_hessian.readlines()))
for i, j in range(9):
    hess[i][j] = hess[i][j]/(math.sqrt(mass(math.floor(i/3)))*math.sqrt(mass(math.floor(j/3))))
freq, mode=np.linalg.eig(hess)


formatNum="{:.8f}".format(Coor).append("   ")