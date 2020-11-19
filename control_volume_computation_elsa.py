import numpy as np
from define_functions import*
from test_main import*
number_of_time_steps = 1000

number_of_divisions = 100
d_length = bed_length/number_of_divisions
dt = d_length / gas_velocity
# print(d_length, dt)
# print(velocity)
print(k_GP)
moisture_particle = np.zeros((1, number_of_divisions))
temp_particle = np.zeros((1, number_of_divisions)) + initial_temp

moisture_gas = np.zeros((1, number_of_divisions)) + relative_humidity
temp_gas = np.zeros((1, number_of_divisions)) + initial_temp

for t in range(100):
    number_x = min(100, t)
    # print(number_x)
    for y in range(1):
        # print('y is', y)
        for x in range(1, number_x):
            # print('x is', x)
            moisture_particle[y, x] = compute_moisture_particle(moisture_particle[y, x], alpha_parameter, N,
                                                                relative_humidity, dt, constant)
            temp_particle[y, x] = compute_temperature_particle(temp_particle[y, x], constant, dt, conductivity_particle,
                                                               laplacian, particle_density, alpha_parameter, moisture_particle[y, x],
                                                               relative_humidity, N, heat_of_vaporization, heat_transfer_coefficient,
                                                               specific_surface_area, temp_gas[y, x], particle_heat_capacity)
            temp_gas[y, x] = compute_temperature_gas(temp_particle[y, x], constant, dt, conductivity_gas, laplacian, gas_density,
                                                     alpha_parameter, moisture_particle[y, x], moisture_vapor_heat_capacity,
                                                     relative_humidity, N, heat_transfer_coefficient, specific_surface_area,
                                                     temp_gas[y, x], gas_heat_capacity, gas_velocity, temp_gradient, porosity_powder, particle_density)
            # moisture_gas[y, x] = compute_moisture_gas()
# print(moisture_particle)
print(temp_particle)
print(temp_gas)


