#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 13:49:13 2019

@author: u301023
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.close('all')


#### Settings ####
x=np.linspace(-80000,80000,1601)
z=np.linspace(0,10000,1601)
l_x = 10000          #m
l_z = 300           #m
g = 9.81            #m/s2
U = 10              #m/s
theta_0 = 290       #K
theta_dz = 0.005    #K/m
z_s=l_z * l_x**2/(l_x**2 + x**2)
N=np.sqrt((g/theta_0)*theta_dz)
delta=np.zeros((len(x),len(z)))
theta=np.zeros((len(x),len(z)))



for i in range(len(z)):
    delta[i,:]= z_s[:] * np.cos((N/U)*(z[i]-z_s[:])) + ((-x*z_s[:])/l_x) * np.sin((N/U)*(z[i]-z_s[:]))
    theta[i,:]= theta_0 * (1 + (N**2*(z[i]-delta[i,:]))/g)

lambda_b = 2*np.pi * (U/N)
print(lambda_b)


#### Plot ####
plt.figure(1)
plt.plot(x,z_s)
plt.ylim(0,10000)
plt.contour(x,z,delta,colors = ['black'])
plt.contourf(x,z,delta,cmap = cm.jet)
plt.colorbar()

plt.figure(2)
plt.plot(x,z_s)
plt.ylim(0,10000)
plt.contour(x,z,theta,colors = ['black'])
plt.contourf(x,z,theta,cmap = cm.jet)
plt.colorbar()
