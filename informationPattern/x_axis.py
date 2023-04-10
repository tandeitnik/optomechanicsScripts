# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 13:46:11 2023

@author: tandeitnik
"""

import numpy as np

#This script evaluates and plots the information radiation pattern for the x direction.
#It is based on the paper Optimal position detection of a dipolar scatterer in a focused field by Felix Tebbenjohanns, Martin Frimmer, and Lukas Novotny

N = 300 #resolution

theta = np.linspace(0,np.pi,N)
phi = np.linspace(0,2*np.pi,N)

r = np.zeros([N,N])
x = np.zeros([N,N])
y = np.zeros([N,N])
z = np.zeros([N,N])

for i in range(N):
    
    for j in range(N):
        
        r[i,j] = (15/(8*np.pi))*np.sin(theta[i])**2*np.cos(phi[j])**2*(1-np.sin(theta[i])**2*np.cos(phi[j])**2)
        x[i,j] = r[i,j]*np.sin(theta[i])*np.cos(phi[j])
        y[i,j] = r[i,j]*np.sin(theta[i])*np.sin(phi[j])
        z[i,j] = r[i,j]*np.cos(theta[i])
        
import matplotlib.pyplot as plt
from matplotlib import cm

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)


ax.set_xlabel('$X$', fontsize=20)
ax.set_ylabel('$Y$', fontsize=20)
ax.set_zlabel(r'$Z$', fontsize=20)


plt.show()
