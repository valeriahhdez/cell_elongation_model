#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 16:11:42 2019

@author: valeria
"""

import numpy as np

import matplotlib.pyplot as plt
import time
from scipy.integrate import odeint, ode
import csv


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
PHI= np.array([0.09, 2.])
#PHI= np.array([2e-4, 4e-3, 0.019, 0.09, 2.])

# --- Plasma membrane hydraulic permeability, MPa-1 . h-1
L_R = np.array([0.018, 18. ])
#L_R = np.array([0.018, 18., 569., 18e3, 18e6])

# --- Yield turgor pressure, MPa
Y = np.array([0.018, 0.06])
#Y = np.array([0.018, 0.06, 0.1, 0.2, 0.7])

# --- Plasmodesmal hydraulic permeability, mm3. MPa-1. h-1
MU = np.array([8e-3, 0.23])
#MU = np.array([0, 4e-7, 3e-4, 8e-3, 0.23, 176])

# --- Source of solutes, MPa-1. h-1
ALPHA = np.array([0.1, , 0.3])
#ALPHA = np.array([0.1, 0.3, 0.5, 0.7, 1.6])

# --- Turgor pressure of seed cell, MPa
P_SEED = np.array([0.11, 0.18])
#P_SEED = np.array([0.023, 0.07, 0.11, 0.18, 0.54])

# --- Osmotic pressure of seed cell, MPa
PI_SEED = np.array([0.8, 1.])
#PI_SEED = np.array([0.8, 1., 1.1, 1.3, 1.7])

#=====================================================

# --------------- Initial conditions -----------------

# --- Volume at t0, mm^3
v_0 = 1.88e-4
# --- Osmotic pressure of the fiber at t0, MPa
PI_FIBER_0 =  PI_SEED
# --- Vector of initial conditions
z0 = [v_0, PI_FIBER_0]


######################################################
# =================== Output lists ==================
######################################################

# ------------ Lists with numerical outputs ------------

# -- R1 = dvT1/dv
R1 = []
# -- R2 = dvT2/dv
R2 = []
# -- R3 = Vf/V0
R3 = []
# -- Maximum volume
V_max = []
# -- Minimum volume
V_min = []
# --  % (+) V
pos_v = []
# -- % (-) V
neg_v = []
# -- % (constant) V
const_v = []
# -- R4 = dpiT2/alpha
R4 = []
# -- R5 = dpiT3/alpha
R5 = []
# -- R6 PImax/PImin
R6 = []
# -- Maximum osmotic pressure
PI_max = []
# -- Minimum osmotic pressure
PI_min = []
# --  % (+) PI
pos_pi = []
# --  % (-) PI
neg_pi = []
# -- % (constant) PI
const_pi = []
# -- R7 Pmax/Pmin
R7 = []
# -- Maximum turgor
P_max = []
# -- minimum turgor
P_min = []
# -- Final turgor
F_p = []
# -- Initial turgor
I_p = []
# -- % (+) P
pos_p = []
# -- % (-) P
neg_p =[]
# -- % (constant) P
const_p = []


# ------------ Lists with values of parameters ------------
Phi = []
Lr = []
Y_ = []
mu_ = []
alpha_ = []
#Psi_e = []
P_s = []
Pi_s = []
Pi_f = []



######################################################
# =============== Model equations =================
######################################################

start_time = time.time()

for phi in PHI:
    for L_r in L_R:
        for y in Y:
            for mu in MU:
                for alpha in ALPHA:
                    for P_seed in P_SEED:
                        for pi_seed in PI_SEED:
                            Psi_ext = P_seed - pi_seed


                            # ============== Condition for turgor pressure ===================
                            def P(v, pi_fiber,t):
                                h = mu / L_r /v
                                if (Psi_ext + pi_fiber - y * (1+ h) + P_seed * h) > 0:
                                #if (Psi_ext + pi_fiber + P_seed * h) > (1 + h) * Y :
                                    P_fiber = (Psi_ext + pi_fiber + P_seed * h + phi * y / L_r) /(1 + h + phi/L_r)
                                # otherwise, use this equation
                                else:
                                    P_fiber = (Psi_ext + pi_fiber + P_seed * h)/(1+h)
                                return P_fiber


                            # ============== Differential equation system ===================
                            def cotton_elongation(z, t, *args):

                                    v = z[0]
                                    pi_fiber = z[1]

                                    # --- Equation for volume
                                    dv =  L_r * v * (Psi_ext + pi_fiber - P(v,pi_fiber, t)) - mu *(P(v,pi_fiber, t) - P_seed)

                                    # --- Equation for osmotic pressure

                                    dpi = alpha - pi_fiber * 1/v * dv - mu/v * (pi_fiber * np.maximum((P(v,pi_fiber, t) - P_seed),0) - pi_seed * np.maximum((P_seed - P(v,pi_fiber, t)),0))

                                    #--- Saves the values volume and ormostic pressure into array
                                    return [dv, dpi]

                            # --- Loop over the initial conditions for PI_fiber
                            for pi in PI_SEED:
                                cond_ini = [v_0, pi]
                                # --- Calculate numerical solution
                                output = odeint(cotton_elongation, cond_ini, t)
                                Pi_f.append(pi)

                                # --- Numerical solution of volume is saved into variable:
                                model_volume = output[:,0]
                                # --- Get the maximum value from numerical solution
                                v_mx = np.max(model_volume)
                                # --- Append to list
                                V_max.append(v_mx)
                                # --- Get the minimum value from numerical solution
                                v_mn = np.min(model_volume)
                                # --- Append to list
                                V_min.append(v_mn)


                                # --- Numerical solution of osmotic pressure is saved into variable:
                                model_pi = output[:,1]
                                # --- Get the maximum value from numerical solution
                                pi_mx = np.max(model_pi[relax_t:])
                                # --- Append to list
                                PI_max.append(pi_mx)
                                # --- Get the minimum value from numerical solution
                                pi_mn = np.min(model_pi[relax_t:])
                                # --- Append to list
                                PI_min.append(pi_mn)
                                r6 = pi_mx/pi_mn
                                R6.append(r6)


                                # --- Turgor pressure computation
                                model_turgor = np.array([P(model_volume[i],model_pi[i], t[i]) for i in range(len(t))])
                                # --- Get the maximum value from numerical solution
                                p_mx = np.max(model_turgor[relax_t:])
                                # --- Append to list
                                P_max.append(p_mx)
                                # --- Get the minimum value from numerical solution
                                p_mn = np.min(model_turgor[relax_t:])
                                # --- Append to list
                                P_min.append(p_mn)
                                r7 = p_mx/p_mn
                                R7.append(r7)
                                # --- Initial and final values of P
                                Init_p = model_turgor[relax_t]
                                I_p.append(Init_p)
                                Final_p = model_turgor[-1]
                                F_p.append(Final_p)



                                # ==================== Functions to compute individual terms from equations of V and PI_fiber  ================

                                # ------------- Terms from volume equation ------------------

                                #--- First term
                                def f_dv(v, pi, P, t):
                                    a = L_r * v * (Psi_ext + pi - P)
                                    return a

                                # --- Second term
                                def sec_dv(P, t):
                                    b = - mu *(P - P_seed)
                                    return b

                                # ----------- Terms from osmotic pressure equation ------------
                                # --- Second term
                                def sec_pi(v, pi, P, t):
                                    g = L_r * (Psi_ext + pi - P) - mu/v *(P - P_seed)

                                    a = -pi * g
                                    return a

                                # --- Third term
                                def third_pi(v, pi, P, t):
                                    b= - mu/v * (pi * np.maximum((P-P_seed),0)- pi_seed*np.maximum((P_seed-P),0))
                                    return b


                                # ======================= Save solutions into variables ===================
#
                                # ------------- Terms from volume equation ------------------
                                # --- First term
                                dv_1 = np.array([f_dv(model_volume[i],model_pi[i], model_turgor[i], t[i]) for i in range(len(t))])

                                # --- Second term term
                                dv_2 = np.array([sec_dv(model_turgor[i], t[i]) for i in range(len(t))])

                                # --- Total dv/dt
                                dv = dv_1 + dv_2


                                # ------------- Terms from osmotic pressure equation ------------------
                                # --- Second term
                                dpi_2 =  np.array([sec_pi(model_volume[i], model_pi[i], model_turgor[i], t[i]) for i in range(len(t))])

                                # --- Third term
                                dpi_3 =  np.array([third_pi(model_volume[i], model_pi[i], model_turgor[i], t[i]) for i in range(len(t))])

                                # --- Total dPI/dt
                                dpi = alpha + dpi_2 + dpi_3

                               # ======================= Computation of numerical outputs ===================


                               # --------- Function to compute the contributions of terms 1 and 2 of dV/dt---------

                                def vol_contribution(x, y):
                                    t_y = []
                                    for i in range(len(x)):
                                        # --- Do not include the first point
                                        if i == 0:
                                            continue

                                        if x[i] > 0:
                                            ty_y = y[i]/x[i]
                                            t_y.append(ty_y)

                                    return t_y

                                # --- Ratio 1: contribution of first term dvT1/dv
                                r1 = vol_contribution(dv, dv_1)
                                r1_1 = np.asarray(r1)
                                # -- Calculate the mean
                                r1_11 = np.mean(r1_1)
                                # -- Save values into list R1
                                R1.append(r1_11)

                                # --- Ratio 2: contribution of second term dvT2/dv
                                r2 = vol_contribution(dv, dv_2)
                                r2_2 = np.asarray(r2)
                                # -- Calculate the mean
                                r2_22 = np.mean(r2_2)
                                # -- Save values into list R2
                                R2.append(r2_22)


                                # --------- Function to compute the ratio of Vf/V0 ---------
                                def Vf_V0(x):
                                    F_0 = x[-1] /x[0]

                                    return F_0

                                r3 = Vf_V0(model_volume)
                                # -- Save values into list R5
                                R3.append(r3)




                                # --------- Function to compute the contributions of terms 2 and 3 of dPI/dt---------
                                def pi_contribution(x, y):
                                    t_y = []
                                    # --- Skip the first point
                                    for i in range(1, len(x)):
                                        t_yy = x[i]/y
                                    t_y.append(t_yy)
                                    return t_y

                                # --- R4: contribution of second term, dPIT2/alpha
                                r4 = pi_contribution(dpi_2, alpha)
                                r4_4 = np.asarray(r4)
                                # -- Calculate the mean
                                r4_44 = np.mean(r4_4)
                                # -- Save values into list R4
                                R4.append(r4_44)

                                # --- R5: contribution of third term, dPIT3/alpha
                                r5 = pi_contribution(dpi_3, alpha)
                                r5_5 = np.asarray(r5)
                                # -- Calculate the mean
                                r5_55 = np.mean(r5_5)
                                # -- Save values into list R5
                                R5.append(r5_55)

                                # --------- Function to estimate the number of positive and negative points in V, PI and P ------
                                thresh = 0.01

                                def F_points(x):
                                    pos = 0.
                                    neg = 0.
                                    cons = 0.
                                    relax = relax_t + 1

                                    for i in range(relax_t, len(x)-1):
                                        a = ( x[i+1]-x[i] ) / x[i]
                                        if a < -thresh:
                                            neg += 1.
                                        elif a > thresh:
                                            pos += 1.
                                        else:
                                            cons +=1.
                                    pos = pos/ (len(x)-relax) * 100.
                                    neg = neg/ (len(x)-relax) * 100.
                                    cons = cons/ (len(x)-relax) * 100.

                                    return np.array([pos, neg, cons])



                                points_pi = F_points(model_pi)
                                # -- % of positive points in PI
                                posiv_Pi = points_pi[0]
                                # -- Save values into list
                                pos_pi.append(posiv_Pi)

                                # -- % of negative points in PI
                                Nega_pi = points_pi[1]
                                # -- Save values into list
                                neg_pi.append(Nega_pi)

                                # -- % of nconstant points in PI
                                Const_pi = points_pi[2]
                                # -- Save values into list
                                const_pi.append(Const_pi)


                                points_P = F_points(model_turgor)
                                # -- Calculate the % of positive points in P
                                PosP = points_P[0]
                                # -- Save values into list R8
                                pos_p.append(PosP)

                                # -- Calculate the % of negative points in P
                                NegP = points_P[1]
                                # -- Save values into list R9
                                neg_p.append(NegP)

                                # -- % of constant points in P
                                Const_p = points_P[2]
                                # -- Save values into list
                                const_p.append(Const_p)


                                points_v = F_points(model_volume)
                                # -- Calculate the % of positive points in P
                                PosV = points_v[0]
                                # -- Save values into list R8
                                pos_v.append(PosV)

                                # -- Calculate the % of negative points in P
                                NegV = points_v[1]
                                # -- Save values into list R9
                                neg_v.append(NegV)

                                # -- % of constant points in P
                                Const_V = points_v[2]
                                # -- Save values into list
                                const_v.append(Const_V)




                                # ----- Add values from current iteration to lists
                                Phi.append(phi)
                                Lr.append(L_r)
                                Y_.append(y)
                                mu_.append(mu)
                                alpha_.append(alpha)
                                #Psi_e.append(Psi_ext)
                                P_s.append(P_seed)
                                Pi_s.append(pi_seed)



# --- Header for columns
header = [ "Phi", "Lr", "Y", "Mu", "Alpha", "P_s", "Pi_s", "Pi_f", "R1", "R2", "R3", "V_max", "V_min", "% (+) V", "% (-) V", "% (c) V", "R4",  "R5", "PImax/PImin", "PI maximum", "PI minimum", "% (+) dpi", "%(-) dpi", "%(c) dpi", "Pmx/Pmin","P maximum", "P minimum", "Initial P", "Final P", "%(+) P", "%(-) P", "%(c) P"]


## ---------------- Expor data to csv file -------------------
with open('28.08.19_xx.csv', 'w') as out_file:
    writer = csv.DictWriter(out_file, fieldnames = header)
    writer.writeheader()
    for i in range(len(mu_)):
        out_string=""
        out_string += str(Phi[i])
        out_string += "," + str(Lr[i])
        out_string += "," + str(Y_[i])
        out_string += "," + str(mu_[i])
        out_string += "," + str(alpha_[i])
        out_string += "," + str(P_s[i])
        out_string += "," + str(Pi_s[i])
        out_string += "," + str(Pi_f[i])
        out_string += "," + str(R1[i])
        out_string += "," + str(R2[i])
        out_string += "," + str(R3[i])
        out_string += "," + str(V_max[i])
        out_string += "," + str(V_min[i])
        out_string += "," + str(pos_v[i])
        out_string += "," + str(neg_v[i])
        out_string += "," + str(const_v[i])
        out_string += "," + str(R4[i])
        out_string += "," + str(R5[i])
        out_string += "," + str(R6[i])
        out_string += "," + str(PI_max[i])
        out_string += "," + str(PI_min[i])
        out_string += "," + str(pos_pi[i])
        out_string += "," + str(neg_pi[i])
        out_string += "," + str(const_pi[i])
        out_string += "," + str(R7[i])
        out_string += "," + str(P_max[i])
        out_string += "," + str(P_min[i])
        out_string += "," + str(I_p[i])
        out_string += "," + str(F_p[i])
        out_string += "," + str(pos_p[i])
        out_string += "," + str(neg_p[i])
        out_string += "," + str(const_p[i])
        out_string += "\n"
        out_file.write(out_string)







elapsed_time = time.time() - start_time
