import numpy as np
import math
from define_functions import*
from test_main import*

# time_step=0.1
# number_measure_points=int(bed_length/(time_step*superficial_velocity))+1 #number of measure points
# end_time=1
# space_step=bed_length/number_measure_points
time_steps = 6
number_measure_points = 1000
space_step = bed_length/number_measure_points           # 0.2 mm w 1000 divisions
time_step = space_step / superficial_velocity
print('time step: ', time_step)
print('Length per step: ', space_step)
print('Gas velocity: ', superficial_velocity)       # 2 mm per s

############### INITIAL CONDITIONS ############################################################
moisture_particle = np.zeros(number_measure_points) + moisture_particle_initial
temperature_particle = np.zeros(number_measure_points) + temp_initial
moisture_gas = np.zeros(number_measure_points) + relative_humidity_gas_initial                      #TODO: Check initially RH_0=Y_G,0
moisture_gas_2 = np.zeros(number_measure_points) + relative_humidity_gas_initial                      #TODO: Check initially RH_0=Y_G,0
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

############### UPDATE OF PARTICLE PARAMETERS THAT THE GAS HAS YET REACHED #####################
moisture_difference_y_absorption=np.zeros(number_measure_points)
moisture_gradient_gas=np.zeros(number_measure_points)

    # position_gas[k]=k*time_step*superficial_velocity = approximately measure point
for t in range(time_steps):
    at_most_bed_length = min(number_measure_points, t+1)
    moisture_gas_old = moisture_gas
    for i in range(at_most_bed_length):

        ###### SAVE OLD PARAMETERS
        m_p_old = moisture_particle[i]
        t_p_old = temperature_particle[i]
        m_G_old = moisture_gas[i]
        t_G_old = temperature_gas[i]

        ###### COMPUTE GRADIENTS AND LAPLACIAN
        # print('t: ', t)
        # print('i: ', i)
        # print('moisture_gas[i]: ', moisture_gas[i])
        # print('moisture_gas[i+1]: ', moisture_gas[i+1])
        moisture_gradient_gas[i] = compute_gradient(moisture_gas, i, space_step)
        print('moisture_gradient_gas i: ', moisture_gradient_gas[i])
        #print('moisture_gradient_gas: ', moisture_gradient_gas)
        
        ###### UPDATE MOISTURE PARTICLE
        moisture_particle[i], moisture_gas[i], moisture_difference_y_absorption[i]= compute_moisture_particle(
            m_p_old, alpha_parameter, N, relative_humidity[i], time_step, constant[i],  m_G_old, particle_density, porosity_powder, gas_density)
        
        ###### UPDATE TEMP PARTICLE
        temperature_particle[i] = compute_temperature_particle(
            t_p_old, constant[i], time_step, conductivity_particle, laplacian, particle_density, alpha_parameter, 
            m_p_old, relative_humidity[i], N, heat_of_vaporization, heat_transfer_coefficient_initial, 
            specific_surface_area, t_G_old, particle_heat_capacity, 1)
        
        ###### UPDATE MOISTURE GAS
        moisture_gas_2[i] = (particle_density*porosity_powder*moisture_particle[i])/(gas_density*(1-porosity_powder))
            # compute_moisture_gas(
            # m_p_old, m_G_old, alpha_parameter, N, relative_humidity[i], time_step, constant[i], superficial_velocity, 
            # moisture_diffusivity, gradient_moisture, laplacian, gas_density, particle_density, porosity_powder, 1)
        #print('moisture Gas computed: ', moisture_gas[i])
        
        ###### UPDATE TEMP GAS
        temperature_gas[i] = compute_temperature_gas(
            t_p_old, constant[i], time_step, conductivity_gas, laplacian, gas_density, alpha_parameter, m_G_old, N, 
            moisture_vapor_heat_capacity, relative_humidity[i], heat_transfer_coefficient_initial, 
            specific_surface_area, t_G_old, gas_heat_capacity, superficial_velocity, temp_gradient, porosity_powder, 
            particle_density, 1)
        
        ###### UPDATE PARAMETERS 
        #moisture_density[i]=moisture_density[i] #for now since we dont know how the vapor density develops with temp
        pressure_saturated[i]=compute_p_saturated(A, B, temperature_gas[i], C)

        # molar_concentration_moisture[i]=moisture_density[i] / molar_mass_moisture

        partial_pressure_moisture[i]=compute_partial_pressure_moisture(
            molar_concentration_moisture[i], R_gas_constant, temperature_gas[i])  

        relative_humidity[i]=partial_pressure_moisture[i]/pressure_saturated[i]

        k_GP[i] = compute_mass_transfer_coefficient(
            moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, 
            flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i])[3]
        
        # constant[j]=k_GP[j] * specific_surface_area * pressure_saturated[j] / pressure_ambient #TODO: doesn't change for now

        if i % 500 == 1:
        #     print('i: ', i)
        #     print('temp particle: ', temperature_particle)
            print('moisture particle: ', moisture_particle[0])
            print('moisture gas: ', moisture_gas[0])

        # print('t: ', t)
        # print('moisture-gas gradient', (moisture_gas[i]-m_G_old)/space_step)

# print('Change in temp particles:\n', (temperature_particle - temp_initial)[0:5])
# print('Change in temp gas:\n', (temperature_gas - temp_initial)[0:5])

print('Change in moisture particles:\n', (moisture_particle - moisture_particle_initial)[0:5])
print('Change in moisture gas:\n', (moisture_gas - relative_humidity_gas_initial)[0:5])



################# OLD COMPUTE MOISTURE PARTICLE
# def compute_moisture_particle(moisture_particle, alpha, N, relative_humidity, dt, constant):
#     change_moisture = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle, N))
#     moisture_particle_current = moisture_particle + change_moisture * dt
#     return moisture_particle_current