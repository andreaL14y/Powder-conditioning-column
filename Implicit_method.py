import numpy as np
import math
from define_functions import*
from test_main import*

#from control_volume_computation_andrea import *
time_steps = 100
number_measure_points = 100
space_step = bed_length/number_measure_points
time_step = space_step / superficial_velocity

############### INITIAL CONDITIONS ############################################################
moisture_particle = np.zeros(number_measure_points) + moisture_particle_initial
temperature_particle = np.zeros(number_measure_points) + temp_initial
moisture_gas = np.zeros(number_measure_points) + moisture_gas_initial_in
moisture_gas_2 = np.zeros(number_measure_points) + moisture_gas_initial_in                   
temperature_gas = np.zeros(number_measure_points) + temp_initial

############### INITIAL PARAMETERS THAT CHANGE DIRECTLY WITH THE TEMPERATURE ###################
#moisture_density = np.zeros(number_measure_points)+moisture_density                                    # TODO: for simplifications this will be considered as constant
pressure_saturated = np.zeros(number_measure_points)+pressure_saturated_initial

############### INITIAL PARAMETERS THAT CHANGE INDIRECTLY WITH THE TEMPERATURE #################
molar_concentration_moisture = np.zeros(number_measure_points)+molar_concentration_moisture_initial #dep. on moisture_density
relative_humidity = np.zeros(number_measure_points)+relative_humidity_gas_initial #dep. on pressure saturated
partial_pressure_moisture = np.zeros(number_measure_points)+partial_pressure_moisture_initial #dep. on pressure saturated and RH_gas

################ INITIAL PARAMETERS THAT CHANGE WITH MOISTURE & TEMPERATURE (PARTICLE/GAS)#####

constant = np.zeros(number_measure_points)+constant_initial #dep. on pressure saturated, k_GP
k_GP = np.zeros(number_measure_points)+k_GP_initial #dep. on molar_concentration_moisture


################ SYSTEM OF LIN EQ FOR du/dt=a*Laplacian*u-b*u+r ################################
################ Temperature Particle                           ################################
M=np.zeros((number_measure_points, number_measure_points))
r=np.zeros(number_measure_points)
for t in range(time_steps):
    #at_most_bed_length = min(number_measure_points, t+1)
    temperature_particle_current = temperature_particle
    a =  conductivity_particle/(particle_density*particle_heat_capacity)
    b = -heat_transfer_coefficient_initial*specific_surface_area*1/particle_heat_capacity
    for i in range(number_measure_points):
        r[i] = constant[i]*(relative_humidity[i]-compute_equilibrium_moisture(alpha_parameter, moisture_particle[i],N)) \
            *heat_of_vaporization/particle_heat_capacity+heat_transfer_coefficient_initial*specific_surface_area \
                *temperature_gas[i]/particle_heat_capacity
        if i==0:
            M[i,i]=1/time_step+2*a/(space_step**2)-b
            M[i,i+1]=-a/(space_step**2)
        elif i==(number_measure_points-1):
            M[i,i-1]=-a/(space_step**2)
            M[i,i]=1/time_step+2*a/(space_step**2)-b
        else: 
            M[i,i-1]=-a/(space_step**2)
            M[i,i]=1/time_step+2*a/(space_step**2)-b
            M[i,i+1]=-a/(space_step**2)

    ############ UPDATE PARAMETERS 
        pressure_saturated[i]=compute_p_saturated(A, B, temperature_gas[i], C)
        
        relative_humidity[i]=compute_relative_humidity_from_Y(molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[i], pressure_saturated[i])
        
        molar_concentration_moisture[i]=compute_molar_concentration(relative_humidity[i], pressure_saturated[i], R_gas_constant, temperature_gas[i])
        
        partial_pressure_moisture[i]=compute_partial_pressure_moisture(molar_concentration_moisture[i], R_gas_constant, temperature_gas[i])  

        k_GP[i] = compute_mass_transfer_coefficient(
            moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, 
            flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i])[3]
        
        constant[i]=k_GP[i] * specific_surface_area * pressure_saturated[i] / pressure_ambient
    #print(M)
    temperature_particle = np.linalg.solve(M, temperature_particle_current)+time_step*r
    # UPDATE PARAMETERS
    print(temperature_particle[0])
