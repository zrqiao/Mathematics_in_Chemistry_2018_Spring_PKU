import numpy as np

C=299792458
Planck_Constant=6.62607004e-34
Boltzmann_Constant=1.38064852e-23

hess = np.loadtxt('materials/hess.dat')
mass = np.array([1.0008,1.0008,1.0008, 15.9997, 15.9997,15.9997,1.0008,1.0008, 1.0008])*1836

for i in range(9):
    for j in range(9):
        hess[i][j] = hess[i][j]/(np.sqrt(mass[i])*np.sqrt(mass[j]))
freq2, mode=np.linalg.eig(hess)#得到质量约化坐标下的mode
#print(freq2,mode)
freq=np.zeros(9)
freqr=np.zeros(9)
Coor0=np.zeros((9,500))
Coor=np.zeros((500,9))
part=np.zeros(3)

iniCoor = open('materials/config_opt.dat')

out=open('log.txt','w')

for i in range (3):
    mode[:,i]=mode[:,i]/np.sqrt(np.dot(mode[:,i],mode[:,i]))#Unitary Transformation
    temp=iniCoor.readline().split()
    for j in range (3):
        Coor0[3*i+j]=float(temp[j+1])
Coor0=Coor0.transpose()
Coor_S=Coor0[0]
out.write('Calculating Modes...')
ni=[]
for i in range(9):
    if freq2[i]>=0:
        ni.append(i)
        T = np.linspace(0.0, 50.0, num=500)
        T = np.reshape(np.repeat(T, 9), (500, 9))
        freq[i] = np.sqrt(freq2[i])
        freqr[i] = freq[i]/(2*np.pi)/(137.036*5.291772108*10**(-9))
        Coor = Coor0+np.sin(freq[i]*T)*(mode[:,i]*(1/np.sqrt(mass)))
        traj = open('materials/traj'+str(i)+'.mdcrd','w')
        np.savetxt(traj,Coor,fmt='%.8f')
        traj.close()

out.write('Finished\n\n')
out.write('Calculating partition functions...\n')

#deltaX2=(i+1/2)/freq[i] #原子单位制下质量约化归一化的dx,hb=1，需要乘以mode
mode_cart=mode.transpose()/np.sqrt(mass)
D0_OH1=np.linalg.norm(Coor_S[3:6]-Coor_S[:3])
D0_OH2=np.linalg.norm(Coor_S[3:6]-Coor_S[6:])
Theta=np.arccos(np.dot(Coor_S[:3]-Coor_S[3:6], Coor_S[6:]-Coor_S[3:6])/(D0_OH1*D0_OH2))
dD0_OH1_m1 = np.linalg.norm(mode_cart[0][:3] - mode_cart[0][3:6])
dD0_OH2_m1 = np.linalg.norm(mode_cart[0][6:9] - mode_cart[0][3:6])
dD0_OH1_m2 = np.linalg.norm(mode_cart[1][:3] - mode_cart[1][3:6])
dD0_OH2_m2 = np.linalg.norm(mode_cart[1][6:9] - mode_cart[1][3:6])
dTheta0 = (np.linalg.norm(mode_cart[2][:3] - mode_cart[2][3:6]) + np.linalg.norm(mode_cart[2][6:] - mode_cart[2][3:6]) )/(D0_OH1+D0_OH2)
out.write(
    'Average O_H1 length: ' + np.str(D0_OH1)+'\nAverage O_H2 length: ' + np.str(D0_OH2) + '\nAverage bond angle: ' + np.str(Theta/np.pi*180)+'\n')

StateList=np.arange(1/2,10000+1/2,1)


T=300
EkO=np.zeros(3)
EkH=np.zeros(3)
Boltzmann_List = np.zeros((3,10000))
Fluct=np.zeros(3)#归一化的谐振子涨落幅
out.write('Temperature: 300K\n')
for i in range(3):
    n=ni[i]
    Boltzmann_List[i] = np.exp(-(StateList - 1 / 2) * freq[i] * 3.1577464e5 / (T))
    part[i]=np.sum(Boltzmann_List[i])
    Fluct[i]= (np.dot(Boltzmann_List[i],StateList))/freq[i]/part[i]
    EkO[i]=np.dot(Boltzmann_List[i],StateList)/2 *np.linalg.norm(mode[i][3:6])**2/np.linalg.norm(mode[i])**2
    EkH[i] = (np.dot(Boltzmann_List[i], StateList) / 2 - EkO[i])/2
    out.write('Mode ' + np.str(i + 1) + '    Frequency: ' + np.str(
        freqr[i]) + 'cm-1    Partition function per mode: ' + np.str(part[n]) + '\n')
dD1 = np.sqrt(Fluct[0] * dD0_OH1_m1 ** 2 + Fluct[1] * dD0_OH1_m2 ** 2)
dD2 = np.sqrt(Fluct[0] * dD0_OH2_m1 ** 2 + Fluct[1] * dD0_OH2_m2 ** 2)
dTheta = np.sqrt(Fluct[2]) * dTheta0

out.write('Partition function: ' + np.str(np.prod(part)) + '\n'
          + 'O-H1 bond length fluctuation: '+np.str(dD1)+'\n'
          + 'O-H2 bond length fluctuation: ' + np.str(dD2) + '\n'
          + 'H-O-H bond angle fluctuation: ' + np.str(dTheta) + '\n')

out.write('Kinetic Energy:\n' + 'O '+np.str(np.sum(EkO))+'\n'
          + 'H ' + np.str(np.sum(EkH)) + '\n')




T=0#Ground state only
out.write('Temperature: 0K\n')
for i in range(3):
    n=ni[i]
    part[i]=1
    Fluct[i] = StateList[0] / freq[i] / part[i]
    EkO[i] = StateList[0] / 2 * np.linalg.norm(mode[i][3:6]) ** 2 / np.linalg.norm(mode[i]) ** 2
    EkH[i] = (StateList[0] / 2 - EkO[i]) / 2

    out.write('Mode '+np.str(i+1)+'    Frequency: '+np.str(freqr[i])+'cm-1    Partition function per mode: '+np.str(part[n])+'\n')
dD1=np.sqrt((1/2)/freq[0]*dD0_OH1_m1**2+(1/2)/freq[1]*dD0_OH1_m2**2)
dD2=np.sqrt((1/2)/freq[0]*dD0_OH2_m1**2+(1/2)/freq[1]*dD0_OH2_m2**2)
dTheta=np.sqrt((1/2)/freq[2])*dTheta0
out.write('Partition function: ' + np.str(np.prod(part)) + '\n'
          + 'O-H1 bond length fluctuation: '+np.str(dD1)+'\n'
          + 'O-H2 bond length fluctuation: ' + np.str(dD2) + '\n'
          + 'H-O-H bond angle fluctuation: ' + np.str(dTheta) + '\n')
out.write('Kinetic Energy:\n' + 'O '+np.str(np.sum(EkO))+'\n'
          + 'H ' + np.str(np.sum(EkH)) + '\n')
out.write('Wigner Distribution')


out.close()