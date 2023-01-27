#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 09:35:30 2019

@author: valeria
"""

import numpy as np

import matplotlib.pyplot as plt
from scipy.integrate import odeint, ode


####################################################
#============= Time array ===================
####################################################

# --- Initial time, hours
t0 = 0.
# --- End time, hours
tf = 480.
# --- Time step
dt = 2.4
# --- Time vector, hours
t = np.arange(t0, tf, dt)
# --- Relaxation time
relax_t = 5

######################################################
# ============ Model-related parameters =============
######################################################

# --------------- Parameters -----------------

# --- Cell wall extensibility, MPa-1 . h-1
phi_w= 0.019

# --- Plasma membrane hydraulic permeability, MPa-1 . h-1
L_r_m = 569.

# --- Yield turgor pressure, MPa
y = 0.1

# --- Plasmodesmal hydraulic permeability, mm3. MPa-1. h-1
mu = 8e-3

# --- Source of solutes, MPa-1. h-1
alpha_s = 0.5


# --- Turgor pressure of seed cell, MPa
P_seed = 0.11

# --- Osmotic pressure of seed cell, MPa
pi_seed = 1.1

# --- Water potential of seed cell, MPa
Psi_ext= P_seed - pi_seed

# --------------- Initial conditions -----------------

# --- Volume at t0, mm^3
v_0 = 1.88e-4
# --- Osmotic pressure of the fiber at t0, MPa
pi_fiber_0 =  pi_seed
# --- Vector of initial conditions
z0 = [v_0, pi_fiber_0]

# ============== Dynamical parameters ===================

# --------------- Function for plasmodesmal permeability -----------------
def pd(t):
    return mu
#    if 240 <= t <= 275:
#        return mu/2
#    else:
#        return mu

# --------------- Function for source of solutes -----------------   
def alpha(t):
#    return alpha_s
    if 240 < t < 275:
        return alpha_s
    else:
        return 0.
    
    
#def alpha(t):
#    return alpha_s
#    if 0 <= t < 370:
#        return alpha_s * 2
#    else:
#        return alpha_s
#    

# --------------- Function for extensibility of the cell wall -----------------

def phi(t):
    return phi_w
#    return phi_w*(1-t/tf)


# --------------- Function for yield turgor pressure -----------------

def Y(t):
    return y
#    return y*(1+t/tf)


# --------------- Function aquaporins -----------------   
def L_r(t):
    return L_r_m
#    if 170 < t < 290:
#        return L_r_m * 100
#    else:
#        return L_r_m 
    

# ============== Condition for turgor pressure ===================   

def P(v, pi_fiber, t):
    h = pd(t)/L_r(t)/v
    if (Psi_ext + pi_fiber - Y(t) * (1+ h) + P_seed * h) > 0:
        P_fiber = (Psi_ext + pi_fiber + P_seed * h + phi(t) * Y(t) /L_r(t)) / (1+h+phi(t)/L_r(t))
    else:
        P_fiber = (Psi_ext + pi_fiber + P_seed * h)/ (1+h)
    return P_fiber


# ============== Differential equation system ===================

def cotton_elongation(z, t, *args):
    v = z[0]
    pi_fiber = z[1]
    
    # --- Equation for volume
    dv =  L_r(t) * v * (Psi_ext + pi_fiber - P(v,pi_fiber, t)) - pd(t) *(P(v,pi_fiber, t) - P_seed)
    
    # --- Equation for osmotic pressure
    
    dpi = alpha(t) - pi_fiber * 1/v * dv - pd(t)/v * (pi_fiber * np.maximum((P(v,pi_fiber, t) - P_seed),0) - pi_seed * np.maximum((P_seed - P(v,pi_fiber, t)),0))
    
    #--- Save the values of volume and osmostic pressure into array
    return [dv, dpi]


# ============== Numerical solution ===================

output = odeint(cotton_elongation, z0, t)   

# --- Numerical solution of volume is: 
model_volume = output[:,0]

# --- Numerical solution of osmotic pressure is: 
model_pi = output[:,1]

# --- Turgor pressure computation
model_turgor = np.array([P(model_volume[i],model_pi[i], t[i]) for i in range(len(t))])



# ------- Subplot for variables ----------


fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, sharex=True)
ax1.plot(t,model_volume, c= 'orange', linewidth= 3.0, label='V')
ax1.grid(True)
ax1.legend(loc='lower right', bbox_to_anchor=(1., 0.7), bbox_transform=plt.gcf().transFigure)

ax2.plot(t,model_pi, c= 'blue', linewidth= 3.0, label='PI')
ax2.grid(True)
ax2.legend(loc='lower right', bbox_to_anchor=(1., 0.45), bbox_transform=plt.gcf().transFigure)

ax3.plot(t,model_turgor, c= 'green', linewidth= 3.0, label='P')
ax3.grid(True)
#plt.legend(bbox_to_anchor=(1.1, 0.9))
ax3.legend(loc='lower right', bbox_to_anchor=(1., 0.2), bbox_transform=plt.gcf().transFigure)

#plt.show()
plt.savefig("0Toalpha_dynamic.png")

b = np.vectorize(alpha)


plt.figure()
plt.plot(t, b(t), linewidth= 3.0)
#plt.title('Cell wall estensibility')
plt.title('Source of solutes')
#plt.title('Yield turgor pressure')
#plt.title('Aquaporins *100')
#plt.title('Plasmodesmal permeability (half)')

plt.xlabel('Time (hours)')
#plt.ylabel()
plt.ylim(0., 0.55)
#plt.show()
plt.savefig('0Toalpha_profile.png', bbox_inches = 'tight')