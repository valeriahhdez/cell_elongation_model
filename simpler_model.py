# -*- coding: UTF-8 -*-

import getopt
import numpy as np
import sys
import pandas as pd

from scipy.integrate import odeint, ode
from pandas import DataFrame, read_csv

# Simpler model of one cell:
# one cell that exchanges water with its extracellular medium

#parameters of the cell
area = 1.79e-8 # in m^3
permeability = 1.7e-8 # in m s⁻1
solutes_quantity = 225 # in mol m^3
yield_strain = 0.1 #dimensionless
yield_turgor = 0.28 # MPa
wall_extensibility = 6e-3 # MPa s⁻1
elastic_modulus = 6 # MPa
gas_constant = 8.32 # gas constant in J mol-1 K-1
temperature =  293 # in K
water_potential = - 0.8 # MPa
dt = 1e-2

# initial conditions
volume_0 = 8.63e-14 # m^3
turgor_0 = 0.5 # MPa
i_0 = [volume_0, turgor_0] # vector of initial values for volume and turgor
time = np.arange(0, 200, dt) # time vector

# Script parameters
# name of the run
#output_basis = 'onemodel_test'
#STORE = True

# variables of the system:
def osmotic(temperature, gas_constant, volume, solutes_quantity):
    osmotic_pressure = (gas_constant * temperature * solutes_quantity)/volume
    return osmotic_pressure

# Differential equation system to be solved:
def f(y, t):
    volume = y[0]
    turgor = y[1]
    f0 = [area * permeability * water_potential * (osmotic(temperature, gas_constant, volume, solutes_quantity) \
            - turgor)]
    f1 = [(yield_turgor/yield_strain * volume)* (area * permeability * water_potential * \
            (osmotic(temperature, gas_constant, volume, solutes_quantity) \
            - turgor)) - wall_extensibility * elastic_modulus * (max(0, turgor -  yield_turgor))]
    return np.array([f0, f1])

# or try data = np.array([dvdt, dTdt])
# and then return data


# Resolution of the differential system
solflux = odeint(f, i_0, time)

# f0 and f1 are vectors that contain the solution values for each time step; transpose these vectors to plot them
#f0 = y.T[0]
#f1 = y.T[1]
