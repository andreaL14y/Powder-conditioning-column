import numpy as np
from define_functions import*
from test_main import*
number_of_time_steps = 6

number_of_divisions = 10000
d_length = bed_length/number_of_divisions
dt = d_length / superficial_velocity
print('time step: ', dt)
specific_surface_area = 400

################################### SET INITIAL CONDITIONS #############################################################
relative_humidity = np.zeros(number_of_divisions) + relative_humidity_gas_initial
pressure_saturated = np.zeros(number_of_divisions) + pressure_saturated_initial
partial_pressure_moisture = np.zeros(number_of_divisions) + partial_pressure_moisture_initial
molar_concentration_moisture = np.zeros(number_of_divisions) + molar_concentration_moisture_initial

mass_transfer_coefficient = np.zeros(number_of_divisions) + k_GP_initial
heat_transfer_coefficient = np.zeros(number_of_divisions) + heat_transfer_coefficient_initial
print('h_GP: ', heat_transfer_coefficient_initial)
print('k_GP: ', mass_transfer_coefficient)
constant = np.zeros(number_of_divisions) + constant_initial

moisture_particle = np.zeros((1, number_of_divisions)) + molar_concentration_moisture_initial

temp_particle = np.zeros((1, number_of_divisions)) + temp_initial
moisture_gas = np.zeros((1, number_of_divisions)) + relative_humidity_gas_initial
temp_gas = np.zeros((1, number_of_divisions)) + temp_initial

for t in range(number_of_time_steps):
    number_x = min(number_of_divisions, t+1)
    for y in range(1):
        for x in range(number_x):
            # print('x is: ', x)
            temp_particle_current = temp_particle[y, x]
            temp_gas_current = temp_gas[y, x]
            moisture_particle_current = moisture_particle[y, x]
            moisture_gas_current = moisture_gas[y, x]

            moisture_particle[y, x] = compute_moisture_particle(
                moisture_particle[y, x], alpha_parameter, N, relative_humidity[x], dt, constant[x])

            temp_particle[y, x] = compute_temperature_particle(
                temp_particle_current, constant[x], dt, conductivity_particle, laplacian, particle_density, alpha_parameter,
                moisture_particle_current, relative_humidity[x], N, heat_of_vaporization, heat_transfer_coefficient[x],
                specific_surface_area, temp_gas_current, particle_heat_capacity, x)

            temp_gas[y, x] = compute_temperature_gas(
                temp_particle_current, constant[x], dt, conductivity_gas, laplacian, gas_density, alpha_parameter,
                moisture_particle_current, N, moisture_vapor_heat_capacity, relative_humidity[x], heat_transfer_coefficient[x],
                specific_surface_area, temp_gas_current, gas_heat_capacity, superficial_velocity, temp_gradient, porosity_powder,
                particle_density, x)

            moisture_gas[y, x] = compute_moisture_gas(
                moisture_particle_current, moisture_gas_current, alpha_parameter, N, relative_humidity[x], dt, constant[x], superficial_velocity,
                moisture_diffusivity, gradient_moisture, laplacian_moisture, gas_density, particle_density, porosity_powder, x)

            ############################ UPDATE PARAMETERS #############################################################
            pressure_saturated[x] = compute_p_saturated(A, B, temp_gas[y, x], C)
            partial_pressure_moisture[x] = compute_partial_pressure_moisture(molar_concentration_moisture[x],
                                                                             R_gas_constant, temp_gas[y, x])

            relative_humidity[x] = compute_relative_humidity(partial_pressure_moisture[x], pressure_saturated[x])

            mass_transfer_coefficient[x] = compute_mass_transfer_coefficient(
                moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density,
                flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[x])[3]

            constant[x] = mass_transfer_coefficient[x] * specific_surface_area * pressure_saturated[x]/pressure_ambient  # just some simplification

            heat_transfer_coefficient[x] = compute_heat_transfer_coefficient(
                moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density,
                flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[x],
                gas_heat_capacity, conductivity_gas)




print('Change in temp particles:\n', (temp_particle - temp_initial)[0, 0:5])
print('Change in temp gas:\n', (temp_gas - temp_initial)[0, 0:5])

print('Change in moisture particles:\n', (moisture_particle - moisture_particle_initial)[0, 0:5])
print('Change in moisture gas:\n', (moisture_gas - relative_humidity_gas_initial)[0, 0:5])
# print('Moisture gas:', (moisture_gas)[0, 0:5])
# TODO: moisture in gas keeps going down, becoming negative. This is since RH is not changing, since gas temp and molar
#  concentration are not changing. How does c change?!

# print(relative_humidity_gas_initial - moisture_gas)
# print(temp_particle - temp_initial)
# print(temp_gas[0, 0:5])
# print(relative_humidity)