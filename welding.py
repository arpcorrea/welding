# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 16:26:41 2022

@author: Admin
"""
import numpy as np
import matplotlib.pyplot as plt

a0      =   10.
w0      =   50.
b0      =   50.
Dic0    =   1. / (1. + w0/b0)
C1      =   5. * (1. + (w0/b0)) * (a0/b0)**2
P       =   800000
u       =   72730.9
tR      =   1.
Tg      =   280
t_total =   30.
dt      =   0.1
nstep   =   int(t_total/dt)

time    =   np.zeros((nstep+1))
Dic     =   np.ones((nstep+1))*Dic0
Dic2     =   np.ones((nstep+1))*Dic0
Dh      =   np.zeros((nstep+1))
Db      =   np.zeros((nstep+1))
T      =   300
D_int      =   np.ones((nstep+1))*Dic0
DDic      =   np.zeros((nstep+1))
for i in range(1,nstep+1):
    time[i] = i*dt

n=0
b_i=b0
Dic_sum=0
for i in range (1,nstep+1): 
    if n ==0:
        if b_i<b0+w0 and T>Tg:
            Dic_sum=Dic_sum+((P/u)*dt)
            b_i=((b0**5)+5*(b0+w0)*((a0*b0)**2)*Dic_sum)**0.2       
        Dic[i]=b_i/(b0+w0)
        DDic[i]=Dic[i]-Dic[i-1]
        if Dic[i] > 1.:
            # import pdb; pdb.set_trace()
            DDic[i]=DDic[i] - (Dic[i]-1)
            Dic[i]=Dic[i-1]+DDic[i]
            n=1
    else:
        Dic[i]=1.
        
    Dh[i]   = ((1/tR)*time[i])**(1./4.)
    if Dh[i] > 1.:
        Dh[i] = 1.
    aux=0
    
    for j in range (1,i):
        aux = aux+((Dic[j]-Dic[j-1]))*Dh[i-j]
    Db[i]   = Dic[0]*Dh[i] + aux
    if Db[i] >1.:
        Db[i] = 1.
        
plt.plot(time, Dh, label="Dh")
plt.plot(time, Dic, label="Dic")
plt.plot(time, Db, label="Db")
plt.xlabel('time, s')
plt.ylabel('Degree of, -')
plt.legend()

def calculate_Dic(total_time=30, dt=0.1, a0=10, b0=50, w0=50, P=800000):
    n=0
    nstep=int(t_total/dt)
    b_i=b0
    Dic_sum=0
    DDic      =   np.zeros((nstep+1))
    for i in range (1,nstep+1): 
        if n ==0:
            if b_i<b0+w0 and T>Tg:
                Dic_sum=Dic_sum+((P/u)*dt)
                b_i=((b0**5)+5*(b0+w0)*((a0*b0)**2)*Dic_sum)**0.2       
            Dic[i]=b_i/(b0+w0)
            DDic[i]=Dic[i]-Dic[i-1]
            if Dic[i] > 1.:
                # import pdb; pdb.set_trace()
                DDic[i]=DDic[i] - (Dic[i]-1)
                Dic[i]=Dic[i-1]+DDic[i]
                n=1
        else:
            Dic[i]=1.
            
        Dh[i]   = ((1/tR)*time[i])**(1./4.)
        if Dh[i] > 1.:
            Dh[i] = 1.
    plt.plot(time, Dic, label="Dic")
    return time, Dh

def calculate_Dh0(total_time=30, dt=0.1, tR=10.5):
    nstep=int(total_time/dt)
    Dh0=[0]
    time=[0]
    for i in range (1,nstep+1): 
        time.append(i*dt)
        Dh0.append(((1/tR)*time[i])**(1./4.))
        if Dh0[i] > 1.:
            Dh0[i] = 1.      
    plt.plot(time, Dh0, label="Dh0")
    return time,Dh0
    
    