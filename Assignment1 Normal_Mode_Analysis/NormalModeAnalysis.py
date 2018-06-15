import numpy as np
import math
hess = np.loadtxt('materials/hess.dat')
mass = np.array([1.0008,1.0008,1.0008, 15.9997, 15.9997,15.9997,1.0008,1.0008, 1.0008])

for i in range(9):
    for j in range(9):
        hess[i][j] = hess[i][j]/(np.sqrt(mass[i])*np.sqrt(mass[j]))
freq2, mode=np.linalg.eig(hess)
#print(freq2,mode)
freq=np.zeros(9)
Coor0=np.zeros((9,500))
Coor=np.zeros((500,9))
iniCoor = open('materials/config_opt.dat')
L=1
for i in range (3):
    temp=iniCoor.readline().split()
    for j in range (3):
        Coor0[3*i+j]=float(temp[j+1])
Coor0=Coor0.transpose()
FreqW=open('Frequencies.dat','w')
for i in range(9):
    if freq2[i]>0:
        T=np.linspace(0,50,num=500).repeat(9)
        T=np.reshape(T,(500,9))
        freq[i]=np.sqrt(freq2[i])
        freqr=freq[i]/(2*np.pi)/(137*5.2917721*10**(-9)*np.sqrt(1836))
        FreqW.write(str(freqr)+'\n')
        Coor=Coor0+L*mode[:][i]*np.sin(freq[i],T)*np.sqrt(1/mass)
        traj=open('materials/traj'+str(i)+'.mdcrd','w')
        np.savetxt(traj,Coor,fmt='%.8f')
        traj.close()
FreqW.close()