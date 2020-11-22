import numpy as np
from define_functions import*
from test_main import*
number_of_time_steps = 10

number_of_divisions = 100
d_length = bed_length/number_of_divisions
dt = d_length / gas_velocity
print(d_length, dt)

################################### SET INITIAL CONDITIONS #############################################################
relative_humidity = np.zeros(number_of_divisions) + relative_humidity_gas_initial
pressure_saturated = np.zeros(number_of_divisions) + pressure_saturated_initial
partial_pressure_moisture = np.zeros(number_of_divisions) + partial_pressure_moisture_initial
molar_concentration_moisture = np.zeros(number_of_divisions) + molar_concentration_moisture_initial

mass_transfer_coefficient = np.zeros(number_of_divisions) + k_GP_initial
heat_transfer_coefficient = np.zeros(number_of_divisions) + heat_transfer_coefficient_initial
constant = np.zeros(number_of_divisions) + constant_initial

moisture_particle = np.zeros((1, number_of_divisions)) + initial_moisture_particle

temp_particle = np.zeros((1, number_of_divisions)) + temp_initial
moisture_gas = np.zeros((1, number_of_divisions)) + relative_humidity_gas_initial
temp_gas = np.zeros((1, number_of_divisions)) + temp_initial

for t in range(number_of_time_steps):
    number_x = min(number_of_divisions, t)
    for y in range(1):
        for x in range(1, number_x):

            moisture_particle[y, x] = compute_moisture_particle(
                moisture_particle[y, x], alpha_parameter, N, relative_humidity[x], dt, constant[x])

            temp_particle[y, x] = compute_temperature_particle(
                temp_particle[y, x], constant[x], dt, conductivity_particle, laplacian, particle_density, alpha_parameter,
                moisture_particle[y, x], relative_humidity[x], N, heat_of_vaporization, heat_transfer_coefficient[x],
                specific_surface_area, temp_gas[y, x], particle_heat_capacity)

            # temp_gas[y, x] = compute_temperature_gas(
            #     temp_particle[y, x], constant[x], dt, conductivity_gas, laplacian, gas_density, alpha_parameter,
            #     moisture_particle[y, x], N, moisture_vapor_heat_capacity, relative_humidity[x], heat_transfer_coefficient[x],
            #     specific_surface_area, temp_gas[y, x], gas_heat_capacity, gas_velocity, temp_gradient, porosity_powder,
            #     particle_density)

            moisture_gas[y, x] = compute_moisture_gas(
                moisture_particle[y, x], alpha_parameter, N, relative_humidity[x], dt, constant[x], gas_velocity,
                moisture_diffusivity, gradient_moisture, laplacian_moisture, gas_density, particle_density, porosity_powder)

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
                gas_heat_capacity)

print(moisture_particle - moisture_particle_initial)
# print(relative_humidity_gas_initial - moisture_gas)
# print(temp_particle - temp_initial)
# print(temp_gas)
print(relative_humidity)