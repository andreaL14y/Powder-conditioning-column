import numpy as np
import math
from define_functions import*
from test_main import*

time_step=0.1
number_measure_points=100
end_time=10

############### INITIAL CONDITIONS ############################################################
moisture_particle = np.zeros(number_measure_points) + moisture_particle_initial
temperature_particle = np.zeros(number_measure_points) + temp_initial
moisture_gas = np.zeros(int(end_time/time_step)) + relative_humidity_gas_initial         #TODO: Check initially RH_0=Y_G,0
temperature_gas = np.zeros(int(end_time/time_step)) + temp_initial

############### INITIAL PARAMETERS THAT CHANGE DIRECTLY WITH THE TEMPERATURE ################### #TODO: temperature of the gas? or some mean?
moisture_density = moisture_density                                                 #TODO: for simplifications this will be considered as constant
pressure_saturated = pressure_saturated_initial

############### INITIAL PARAMETERS THAT CHANGE INDIRECTLY WITH THE TEMPERATURE #################
molar_concentration_moisture = molar_concentration_moisture_initial #dep. on moisture_density
relative_humidity = relative_humidity_gas_initial #dep. on pressure saturated
partial_pressure_moisture = partial_pressure_moisture_initial #dep. on pressure saturated and RH_gas

################ INITIAL PARAMETERS THAT CHANGE WITH MOISTURE & TEMPERATURE (PARTICLE/GAS)#####

constant = constant_initial #dep. on pressure saturated, k_GP                       #TODO: CHECK THIS, this formula is simplified
k_GP = k_GP_initial #dep. on molar_concentration_moisture
heat_transfer_coefficient = heat_transfer_coefficient_initial #dep. on molar_concentration_moisture

############### UPDATE OF PARTICLE PARAMETERS THAT THE GAS HAS YET REACHED #####################
time=0
gas_velocity=gas_velocity = compute_velocity(volumetric_flow_rate_liters_per_minute, bed_length, column_diameter, porosity_powder) #TODO: dont know why this import didnt work
while time<=end_time: 
    position_gas = time*gas_velocity
    #print('pos_gas: ', position_gas)
    for i in range(number_measure_points): 
        position_particle = i*(bed_length/number_measure_points)
        if position_particle <= position_gas: 
            # Update moisture of the particles
            moisture_particle[i] += compute_moisture_particle(moisture_particle[i], alpha_parameter, N, relative_humidity, time_step, constant)
            # Update temperature of the particles
            temperature_particle[i] += compute_temperature_particle(temperature_particle[i], constant, time_step, conductivity_particle, laplacian, particle_density, alpha_parameter, moisture_particle[i], relative_humidity, N,
                                                                    heat_of_vaporization, heat_transfer_coefficient, specific_surface_area, temperature_gas[i], particle_heat_capacity, 1)
            # TODO: Update moisture of the gas inex until where it is changing? only first entry maybe
            # compute_moisture_gas(moisture_particle, moisture_gas, alpha, N, relative_humidity, dt, constant, velocity, diffusivity,
                         #gradient_moisture, laplacian, density_gas, density_particle, porosity, x)
            # TODO: Update temperature of the gas
    time += time_step
    if math.floor(time) % 2 == 0:
        print('time: ', time)
        print('temp particle: ', temperature_particle)
        print('moisture particle: ', moisture_particle)
