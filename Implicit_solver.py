import numpy as np
import matplotlib.pyplot as plt
import os, sys
import matplotlib

import numpy as np
import math
from define_functions import*
from define_parameters_OLD import*
from copy import deepcopy

matplotlib.rc('font', size=18)
matplotlib.rc('font', family='Arial')

#definition of numerical parameters
NN = 101 #number of grid points
L = float(0.2) #size of grid
dx = L/(NN-1) #grid spacing
dt = dx / superficial_velocity #time step
nsteps = 20 #number of time steps

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

#initialize grid
x = np.linspace(0,1,NN)
#initial condition
X_P = np.zeros(NN)+moisture_particle_initial
#evaluate right hand side at t=0
bb = B_tp.dot(X_P) + b

for t in range(nsteps):
    position_gas = t * gas_velocity
    reached_points = position_gas/dx

    for j in range(int(min(NN, reached_points))):
        b_mp[j] = -constant[j]
        r_mp[j] = constant[j]*relative_humidity[j]

        A_mp = np.diagflat([1.]+[1-dt*b_mp[j]*0.5 for j in range(NN-1)])
        B_mp = np.diagflat([1.]+[1+dt*b_mp[j]*0.5 for j in range(NN-1)])
        b = np.array([0.]+[0.5*dt*(r_mp[j]+r_mp[j-1]) for j in range(NN-2)]+[dt*r_mp[j]])

        pressure_saturated[j]=compute_p_saturated(A, B, temperature_gas[j], C)
        
        relative_humidity[j]=compute_relative_humidity_from_Y(molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[j], pressure_saturated[j])
       
        molar_concentration_moisture[j]=compute_molar_concentration(relative_humidity[j], pressure_saturated[j], R_gas_constant, temperature_gas[j])
        
        partial_pressure_moisture[j]=compute_partial_pressure_moisture(molar_concentration_moisture[j], R_gas_constant, temperature_gas[j])  

        k_GP[j] = compute_mass_transfer_coefficient(
            moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, 
            flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[j])[3]
        
        constant[j]=k_GP[j] * specific_surface_area * pressure_saturated[j] / pressure_ambient
        
        #find solution inside domain
        X_P [:j] = (np.linalg.solve(A_mp,bb))[:j]
        #update right hand side
        print('b: ', b)
        bb[:j] = (B_mp.dot(X_P) + b)[:j]
#print(X_P)


    #define matrices A, B and b array
    # for i in range(NN):
    #     b_mp[i] = -constant[i]
    #     r_mp[i] = constant[i]*relative_humidity[i]

    #     A_mp = np.diagflat([1.]+[1-dt*b_mp[i]*0.5 for i in range(NN-1)])
    #     B_mp = np.diagflat([1.]+[1+dt*b_mp[i]*0.5 for i in range(NN-1)])
    #     b = np.array([0.]+[0.5*dt*(r_mp[i]+r_mp[i-1]) for i in range(NN-2)]+[dt*r_mp[i]])

    #     pressure_saturated[i]=compute_p_saturated(A, B, temperature_gas[i], C)
        
    #     relative_humidity[i]=compute_relative_humidity_from_Y(molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[i], pressure_saturated[i])
       
    #     molar_concentration_moisture[i]=compute_molar_concentration(relative_humidity[i], pressure_saturated[i], R_gas_constant, temperature_gas[i])
        
    #     partial_pressure_moisture[i]=compute_partial_pressure_moisture(molar_concentration_moisture[i], R_gas_constant, temperature_gas[i])  

    #     k_GP[i] = compute_mass_transfer_coefficient(
    #         moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, 
    #         flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i])[3]
        
    #     constant[i]=k_GP[i] * specific_surface_area * pressure_saturated[i] / pressure_ambient