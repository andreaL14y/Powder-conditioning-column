import numpy as np
import matplotlib.pyplot as plt
import os, sys
import matplotlib

import numpy as np
import math
from define_functions import*
from input_parameters import*
from copy import deepcopy

matplotlib.rc('font', size=18)
matplotlib.rc('font', family='Arial')

#definition of numerical parameters
NN = 101 #number of grid points
L = float(0.2) #size of grid
dx = L/(NN-1) #grid spacing
dt = dx / superficial_velocity #time step
nsteps = 1000 #number of time steps

r1 = dt/dx**2 #assuming heat diffusion coefficient == 1
r2 = dt/(2*dx)

############### INITIAL CONDITIONS ############################################################
moisture_particle = np.zeros(NN) + moisture_particle_initial
temperature_particle = np.zeros(NN) + temp_initial
moisture_gas = np.zeros(NN) + moisture_gas_initial_in
moisture_gas_2 = np.zeros(NN) + moisture_gas_initial_in                   
temperature_gas = np.zeros(NN) + temp_initial

############### INITIAL PARAMETERS THAT CHANGE DIRECTLY WITH THE TEMPERATURE ###################
#moisture_density = np.zeros(NN)+moisture_density                                    # TODO: for simplifications this will be considered as constant
pressure_saturated = np.zeros(NN)+pressure_saturated_initial

############### INITIAL PARAMETERS THAT CHANGE INDIRECTLY WITH THE TEMPERATURE #################
molar_concentration_moisture = np.zeros(NN)+molar_concentration_moisture_initial #dep. on moisture_density
relative_humidity = np.zeros(NN)+relative_humidity_gas_initial #dep. on pressure saturated
partial_pressure_moisture = np.zeros(NN)+partial_pressure_moisture_initial #dep. on pressure saturated and RH_gas

################ INITIAL PARAMETERS THAT CHANGE WITH MOISTURE & TEMPERATURE (PARTICLE/GAS)#####

constant = np.zeros(NN)+constant_initial #dep. on pressure saturated, k_GP
k_GP = np.zeros(NN)+k_GP_initial #dep. on molar_concentration_moisture

# MOISTURE PARTICLE
b_mp = np.zeros(NN)
r_mp = np.zeros(NN)

#initialize matrices A, B and b array
A_tp = np.zeros((NN,NN))
B_tp = np.zeros((NN,NN))
b = np.zeros((NN))
#define matrices A, B and b array
for i in range(NN):
    b_mp[i] = -constant[i]
    r_mp[i] = constant[i]*relative_humidity[i]

    if i==0:
        A_tp[i,:] = [1 if j==0 else 0 for j in range(NN)]
        B_tp[i,:] = [1 if j==0 else 0 for j in range(NN)]
        b[i] = 0. #boundary condition at i=1
    elif i==N-1:
        A_tp[i,:] = [1-dt*b_mp[i]*0.5 if j==i else 0 for j in range(NN)]
        B_tp[i,:] = [1+dt*b_mp[i]*0.5 if j==i else 0 for j in range(NN)]
        b[i] = dt*r_mp[i] #boundary condition at i=NN
    else:
        A_tp[i,:] = [1-dt*b_mp[i]*0.5 if j==i else 0 for j in range(NN)]
        B_tp[i,:] = [1+dt*b_mp[i]*0.5 if j==i else 0 for j in range(NN)]
        b[i] = 0.5*dt*(r_mp[i]+r_mp[i-1])
    
    pressure_saturated[i]=compute_p_saturated(A, B, temperature_gas[i], C)
        
    relative_humidity[i]=compute_relative_humidity_from_Y(molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[i], pressure_saturated[i])
       
    molar_concentration_moisture[i]=compute_molar_concentration(relative_humidity[i], pressure_saturated[i], R_gas_constant, temperature_gas[i])
        
    partial_pressure_moisture[i]=compute_partial_pressure_moisture(molar_concentration_moisture[i], R_gas_constant, temperature_gas[i])  

    k_GP[i] = compute_mass_transfer_coefficient(
            moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, 
            flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i])[3]
        
    constant[i]=k_GP[i] * specific_surface_area * pressure_saturated[i] / pressure_ambient
    


#initialize grid
x = np.linspace(0,1,NN)
#initial condition
u = np.zeros(NN)+moisture_particle_initial
#evaluate right hand side at t=0
bb = B_tp.dot(u) + b

for j in range(nsteps):
    #find solution inside domain
    u = np.linalg.solve(A_tp,bb)
    #update right hand side
    bb = B_tp.dot(u) + b
print(u)
