# -*- coding: UTF-8 -*-

import getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd
import seaborn as sns
from scipy.integrate import odeint, ode
from pandas import DataFrame, read_csv

# Simpler model of one cell:
# one cell that exchanges water with its extracellular medium

#parameters of the cell
area = 2615 # in micrometers^3
permeability = 0.46 # in micrometers s⁻1
solutes_concentration = 2.25e-16 # in mol micrometers^-3
yield_strain = 0.1 #dimensionless
yield_turgor = 0.28 # MPa
wall_extensibility = 6e-3 # MPa s⁻1
elastic_modulus = 6 # MPa
gas_constant = 8.314e-6 # gas constant in J mol-1 K-1
temperature =  293 # in K
water_potential = 0.8 # MPa
dt = 1e-2

# initial conditions
volume_0 = 8432 # microm^3
turgor_0 = 0.5 # MPa
y0 = np.array([volume_0, turgor_0]) # vector of initial values for volume and turgor
time = np.arange(0, 200, dt) # time vector


# Script parameters
# name of the run
#output_basis = 'onemodel_test'
#STORE = True

# variables of the system:
def osmotic(temperature, gas_constant, volume, solutes_concentration):
    osmotic_pressure = gas_constant * temperature * solutes_concentration
    return osmotic_pressure

print osmotic
# Differential equation system to be solved:
def f(y, t):
    volume = y[0]
    turgor = y[1]
    f0 = area * permeability * (water_potential + (osmotic(temperature, gas_constant, volume, solutes_concentration) \
            - turgor))
    f1 = (yield_turgor/yield_strain * volume)* (area * permeability * water_potential * \
            (osmotic(temperature, gas_constant, volume, solutes_concentration) \
            - turgor)) - wall_extensibility * elastic_modulus * (max(0, turgor -  yield_turgor))
    d = np.array([f0, f1])
    return d



if __name__ == '__main__':
    # Resolution of the differential system
    solflux = odeint(f, y0, time)
    # f0 and f1 are vectors that contain the solution values for each time step; transpose these vectors to plot them
    # f0 = y.T[0]
    # f1 = y.T[1]
    time.shape =(1, time.size)
    #vec = np.concatenate((time, y.T)).T
    vec = pd.DataFrame(np.concatenate((time, solflux.T)).T, columns=['time','volume','turgor'])
vec.to_csv('datos.csv')
