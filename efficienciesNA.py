import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from matplotlib import cm

#This script evaluates the detection efficiency for foward and baward detection of a levitated dieletic in a focused EM field.
#It is based on the paper Optimal position detection of a dipolar scatterer in a focused field by Felix Tebbenjohanns, Martin Frimmer, and Lukas Novotny

NA_tl = 0.78 #trapping lens NA

NA_cl_array = np.linspace(0.01,1,100)  
NA_tl_array = np.linspace(0.01,1,100)

n_fw = np.zeros([3,len(NA_cl_array)])
n_bw = np.zeros([3,len(NA_tl_array)])

#forward
for i in range(len(NA_cl_array)):
    
    if NA_cl_array[i] > NA_tl:
        
        n_fw[0,i] = n_fw[0,i-1]
        n_fw[1,i] = n_fw[1,i-1]
        n_fw[2,i] = n_fw[2,i-1]
        
    else:
    
        #evaluating the geometric parameter A
        theta_tl = np.arcsin(NA_tl)
        C = 2*(8/15 - np.cos(theta_tl)**(3/2)/3 - np.cos(theta_tl)**(5/2)/5)
        D = 2*(12/35 - np.cos(theta_tl)**(5/2)/5 - np.cos(theta_tl)**(7/2)/7)
        A = D/C
        
        #evaluating B_fw to calculate the foward detection efficienty
        theta_cl = np.linspace(0,np.arcsin(NA_cl_array[i]),1000)
        s = np.sin(theta_cl)
        c = np.cos(theta_cl)
        Bx_integrand = s*np.sqrt(c)*s*(1+2*c)/3 
        By_integrand = s*np.sqrt(c)*s*(2+c)/3
        Bz_integrand = s*np.sqrt(c)*np.pi*(c-A)*(1+c)/4
        
        Bx_fw = integrate.trapezoid(Bx_integrand,theta_cl)
        By_fw = integrate.trapezoid(By_integrand,theta_cl)
        Bz_fw = integrate.trapezoid(Bz_integrand,theta_cl)
        
        n_fw[0,i] = (30*Bx_fw**2)/(np.pi**2*NA_cl_array[i]**2)
        n_fw[1,i] = (15*By_fw**2)/(np.pi**2*NA_cl_array[i]**2)
        n_fw[2,i] = (1/(1+2.5*A**2))*((15*Bz_fw**2)/(np.pi**2*NA_cl_array[i]**2))


#backward
for i in range(len(NA_tl_array)):
    
    #evaluating the geometric parameter A
    theta_tl = np.arcsin(NA_tl_array[i])
    C = 2*(8/15 - np.cos(theta_tl)**(3/2)/3 - np.cos(theta_tl)**(5/2)/5)
    D = 2*(12/35 - np.cos(theta_tl)**(5/2)/5 - np.cos(theta_tl)**(7/2)/7)
    A = D/C
    
    #evaluating B_bw to calculate the backward detection efficienty
    integrationLimit = np.linspace(0,theta_tl,1000)
    s = np.sin(integrationLimit)
    c = np.cos(integrationLimit)
    Bx_integrand = s*np.sqrt(c)*s*(1+2*c)/3 
    By_integrand = s*np.sqrt(c)*s*(2+c)/3
    Bz_integrand = (np.pi/2)*s*np.sqrt(c)*(1+c)*c
    
    Bx_bw = integrate.trapezoid(Bx_integrand,integrationLimit)
    By_bw = integrate.trapezoid(By_integrand,integrationLimit)
    Bz_bw = integrate.trapezoid(Bz_integrand,integrationLimit)
    
    n_bw[0,i] = (30*Bx_bw**2)/(np.pi**2*NA_tl_array[i]**2)
    n_bw[1,i] = (15*By_bw**2)/(np.pi**2*NA_tl_array[i]**2)
    n_bw[2,i] = (1/(1+2.5*A**2))*((15*Bz_bw**2)/(np.pi**2*NA_tl_array[i]**2))

#determining max value of foward z efficiency and for which NA_cl value it is achieved
max_n_fw_z = np.max(n_fw[2,:])

for i in range(len(NA_cl_array)):
    
    if n_fw[2,i] >= max_n_fw_z:
        
        NA_cl_max = NA_cl_array[i]
        break
    
evenly_spaced_interval = np.linspace(0, 1, 10)
colors = [cm.viridis(x) for x in evenly_spaced_interval]

fig, ax = plt.subplots(1,2, figsize=(18,11), sharex=False)
    
plt.rcParams.update({'font.size': 20})
plt.rcParams["axes.linewidth"] = 1

ax[0].plot(NA_cl_array,n_fw[0,:], label = r"$\eta_x^{fw}$", alpha = 1,lw = 4)
ax[0].plot(NA_cl_array,n_fw[1,:], label = r"$\eta_y^{fw}$", alpha = 1,lw = 4)
ax[0].plot(NA_cl_array,n_fw[2,:]*100, label = r"$\eta_z^{fw}\cdot 100$", alpha = 1,lw = 4)
ax[0].plot(np.linspace(0,NA_cl_max,100),np.ones(100)*max_n_fw_z*100, "--", color = 'k', lw = 2)
ax[0].plot(np.ones(100)*NA_cl_max,np.linspace(0,max_n_fw_z,100)*100, "--", color = 'k', lw = 2)
ax[0].set(xlabel=r'$NA_{cl}$')
ax[0].set(ylabel= r'detection efficiency $\eta$')
ax[0].grid(alpha = 0.4)
ax[0].legend(loc = 'upper left')
ax[0].set(title = "")
ax[0].set_xlim(0,1)
ax[0].set_ylim(0,0.6)
ax[0].text(0.05, 0.35, r"for $NA_{tl} = $"+"${:2.4g}$".format(NA_tl)+\
" the optimal $NA_{cl}$ is "+"${:2.4g}$".format(NA_cl_max)+"\nwith "+\
"$\eta_z^{fw} = $"+"${:2.2g}$".format(max_n_fw_z*100)+"%", fontsize=20)
#ax.ticklabel_format(axis="x", style="sci", scilimits=(0,0))

ax[1].plot(NA_tl_array,n_bw[0,:], label = r"$\eta_x^{bw}$", alpha = 1,lw = 4)
ax[1].plot(NA_tl_array,n_bw[1,:], label = r"$\eta_y^{bw}$", alpha = 1,lw = 4)
ax[1].plot(NA_tl_array,n_bw[2,:], label = r"$\eta_z^{bw}$", alpha = 1,lw = 4)
ax[1].set(xlabel=r'$NA_{tl}$')
#ax[1].set(ylabel= r'detection efficiency $\eta$')
ax[1].grid(alpha = 0.4)
ax[1].legend(loc = 'upper left')
ax[1].set(title = "")
ax[1].set_xlim(0,1)
ax[1].set_ylim(0,1)

fig.tight_layout()
fig.subplots_adjust(hspace=0.1)
