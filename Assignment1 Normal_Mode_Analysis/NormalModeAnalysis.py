import numpy as np
import math
raw_hessian = open('materials/hess.dat')
mass = [1.0008,1.0008,1.0008, 15.9997, 15.9997,15.9997,1.0008,1.0008, 1.0008]
hess = list(map(lambda line: list(map(float,
    [num for num in line.split()])), raw_hessian.readlines()))
for i in range(9):
    for j in range(9):
        hess[i][j] = hess[i][j]/(math.sqrt(mass[i])*math.sqrt(mass[j]))
freq2, mode=np.linalg.eig(hess)
#print(freq2,mode)
freq=[0 for i in range(9)]
Coor0=[0 for i in range(9)]
Coor=[[0 for i in range(9)]for t in range (1000)]
iniCoor = open('materials/config_opt.dat')
L=1
for i in range (3):
    temp=iniCoor.readline().split()
    for j in range (3):
        Coor0[3*i+j]=float(temp[j+1])
FreqW=open('Frequencies.dat','w')
for i in range(9):
    if freq2[i]>0:
        freq[i]=math.sqrt(freq2[i])
        freqr=freq[i]/(2*math.pi)/(137*5.2917721*10**(-9)*math.sqrt(1836))
        FreqW.write(str(freqr)+'\n')
        Coor=[[Coor0[j]+L*mode[j][i]*math.sin(freq[i]*0.05*t)/(math.sqrt(mass[j])) for j in range(9)] for t in range (500)]
        traj=open('materials/traj'+str(i)+'.mdcrd','w')
        for t in range(500):
            traj.write('   '.join(map(lambda x: "{:.8f}".format(x),Coor[t]))+'\n')
        traj.close()
FreqW.close()