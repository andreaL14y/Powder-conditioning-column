import numpy as np
import math
from define_functions import*
from test_main import*

time_step=0.01
number_measure_points=int(bed_length/(time_step*superficial_velocity))+1 #number of measure points
end_time=20

############### INITIAL CONDITIONS ############################################################
moisture_particle = np.zeros(number_measure_points) + moisture_particle_initial
temperature_particle = np.zeros(number_measure_points) + temp_initial
moisture_gas = np.zeros(number_measure_points) + relative_humidity_gas_initial                      #TODO: Check initially RH_0=Y_G,0
temperature_gas = np.zeros(number_measure_points) + temp_initial

############### INITIAL PARAMETERS THAT CHANGE DIRECTLY WITH THE TEMPERATURE ###################    #TODO: temperature of the gas? or some mean?
moisture_density = np.zeros(number_measure_points)+moisture_density                                 #TODO: for simplifications this will be considered as constant
pressure_saturated = np.zeros(number_measure_points)+pressure_saturated_initial

############### INITIAL PARAMETERS THAT CHANGE INDIRECTLY WITH THE TEMPERATURE #################
molar_concentration_moisture = np.zeros(number_measure_points)+molar_concentration_moisture_initial #dep. on moisture_density
relative_humidity = np.zeros(number_measure_points)+relative_humidity_gas_initial #dep. on pressure saturated
partial_pressure_moisture = np.zeros(number_measure_points)+partial_pressure_moisture_initial #dep. on pressure saturated and RH_gas

################ INITIAL PARAMETERS THAT CHANGE WITH MOISTURE & TEMPERATURE (PARTICLE/GAS)#####

constant = np.zeros(number_measure_points)+constant_initial #dep. on pressure saturated, k_GP         #TODO: CHECK THIS, this formula is simplified
k_GP = np.zeros(number_measure_points)+k_GP_initial #dep. on molar_concentration_moisture
heat_transfer_coefficient = np.zeros(number_measure_points)+heat_transfer_coefficient_initial #dep. on molar_concentration_moisture

############### HELPER VARIABLES ###############################################################
superficial_mass_velocity=np.zeros(number_measure_points)
particle_surface_area=np.zeros(number_measure_points)
reynolds_number=np.zeros(number_measure_points)


############### UPDATE OF PARTICLE PARAMETERS THAT THE GAS HAS YET REACHED #####################

    # position_gas[k]=k*time_step*superficial_velocity = approximately measure point
for l in range(int(end_time/time_step)):
    # print('l:', l)
    if l*time_step*superficial_velocity >= bed_length:
        print('in the if loop where l is: ', l)
        for i in range(time_step):
            # Update moisture of the particles
            moisture_particle[i] = compute_moisture_particle(moisture_particle[i], alpha_parameter, N, relative_humidity[i], time_step, constant[i])
            # Update temperature of the particles
            temperature_particle[i] = compute_temperature_particle(temperature_particle[i], constant[i], time_step, conductivity_particle, laplacian, particle_density, alpha_parameter, moisture_particle[i], relative_humidity[i], N,
                                                                heat_of_vaporization, heat_transfer_coefficient[i], specific_surface_area, temperature_gas[i], particle_heat_capacity, 1)
            # Update moisture of the gas
            moisture_gas[i] = compute_moisture_gas(moisture_particle[i], moisture_gas[i], alpha_parameter, N, relative_humidity[i], time_step, constant[i], superficial_velocity, moisture_diffusivity,
                                                gradient_moisture, laplacian, gas_density, particle_density, porosity_powder, 1)
            # Update temperature of the gas
            temperature_gas[i] = compute_temperature_gas(temperature_particle[i], constant[i], time_step, conductivity_gas, laplacian, gas_density, alpha_parameter, moisture_gas[i], 
                                    N, moisture_vapor_heat_capacity, relative_humidity[i], heat_transfer_coefficient[i], specific_surface_area, temperature_gas[i], gas_heat_capacity, superficial_velocity, temp_gradient, porosity_powder, particle_density, 1)
            # Update all other parameters
            #moisture_density[i]=moisture_density[i] #for now since we dont know how the vapor density develops with temp
            pressure_saturated[i]=compute_p_saturated(A, B, temperature_gas[i], C)
            molar_concentration_moisture[i]=moisture_density[i] / molar_mass_moisture
            partial_pressure_moisture[i]=compute_partial_pressure_moisture(molar_concentration_moisture[i], R_gas_constant, temperature_gas[i])  
            relative_humidity[i]=partial_pressure_moisture[i]/pressure_saturated[i]
        
            superficial_mass_velocity[i], particle_surface_area[i], reynolds_number[i], k_GP[i] = compute_mass_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i])
            # constant[j]=k_GP[j] * specific_surface_area * pressure_saturated[j] / pressure_ambient        #TODO: doesn't change for now
            #heat_transfer_coefficient[i]=compute_heat_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i], gas_heat_capacity) #NOT changing unless K_GP changes
        print('i: ', i)
        print('temp particle: ', temperature_particle)
        print('moisture particle: ', moisture_particle)
    elif l==int(end_time/time_step): 
        print('temp particle: ', temperature_particle)
        print('moisture particle: ', moisture_particle)
    else: 
        for i in range(l):
            # Update moisture of the particles
            moisture_particle[i] = compute_moisture_particle(moisture_particle[i], alpha_parameter, N, relative_humidity[i], time_step, constant[i])
            # Update temperature of the particles
            temperature_particle[i] = compute_temperature_particle(temperature_particle[i], constant[i], time_step, conductivity_particle, laplacian, particle_density, alpha_parameter, moisture_particle[i], relative_humidity[i], N,
                                                                heat_of_vaporization, heat_transfer_coefficient[i], specific_surface_area, temperature_gas[i], particle_heat_capacity, 1)
            # Update moisture of the gas
            moisture_gas[i] = compute_moisture_gas(moisture_particle[i], moisture_gas[i], alpha_parameter, N, relative_humidity[i], time_step, constant[i], superficial_velocity, moisture_diffusivity,
                                                gradient_moisture, laplacian, gas_density, particle_density, porosity_powder, 1)
            # Update temperature of the gas
            temperature_gas[i] = compute_temperature_gas(temperature_particle[i], constant[i], time_step, conductivity_gas, laplacian, gas_density, alpha_parameter, moisture_gas[i], 
                                    N, moisture_vapor_heat_capacity, relative_humidity[i], heat_transfer_coefficient[i], specific_surface_area, temperature_gas[i], gas_heat_capacity, superficial_velocity, temp_gradient, porosity_powder, particle_density, 1)
            # Update all other parameters
            #moisture_density[i]=moisture_density[i] #for now since we dont know how the vapor density develops with temp
            pressure_saturated[i]=compute_p_saturated(A, B, temperature_gas[i], C)
            molar_concentration_moisture[i]=moisture_density[i] / molar_mass_moisture
            partial_pressure_moisture[i]=compute_partial_pressure_moisture(molar_concentration_moisture[i], R_gas_constant, temperature_gas[i])  
            relative_humidity[i]=partial_pressure_moisture[i]/pressure_saturated[i]
        
            superficial_mass_velocity[i], particle_surface_area[i], reynolds_number[i], k_GP[i] = compute_mass_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i])
            # constant[j]=k_GP[j] * specific_surface_area * pressure_saturated[j] / pressure_ambient #TODO: doesn't change for now
            #heat_transfer_coefficient[i]=compute_heat_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i], gas_heat_capacity) #NOT changing unless K_GP changes

            # if i % 500 == 1:
            #     print('i: ', i)
            #     print('temp particle: ', temperature_particle)
            #     print('moisture particle: ', moisture_particle)