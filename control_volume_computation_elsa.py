import numpy as np
import matplotlib.pyplot as plt
from define_functions import*
from test_main import*

number_of_time_steps = 1000

number_of_divisions = 100000
d_length = bed_length/number_of_divisions           # 0.2 mm w 1000 divisions
dt = d_length / superficial_velocity
print('time step: ', dt)
print('Length per step: ', d_length)
print('Gas velocity: ', gas_velocity)       # 5 mm per s
specific_surface_area = 400

################################### SET INITIAL CONDITIONS #############################################################
relative_humidity = np.zeros(number_of_divisions) + relative_humidity_gas_initial
pressure_saturated = np.zeros(number_of_divisions) + pressure_saturated_initial
partial_pressure_moisture = np.zeros(number_of_divisions) + partial_pressure_moisture_initial
molar_concentration_moisture = np.zeros(number_of_divisions) + molar_concentration_moisture_initial

mass_transfer_coefficient = np.zeros(number_of_divisions) + k_GP_initial
heat_transfer_coefficient = heat_transfer_coefficient_initial
constant = np.zeros(number_of_divisions) + constant_initial

moisture_particle = np.zeros((1, number_of_divisions)) + moisture_particle_initial
moisture_gas = np.zeros((1, number_of_divisions)) + moisture_gas_initial_bed
moisture_gas[0, 0] = moisture_gas_initial_in

temp_particle = np.zeros((1, number_of_divisions)) + temp_initial
temp_gas = np.zeros((1, number_of_divisions)) + temp_initial

gradient_moisture = gradient_moisture_initial
laplacian_moisture = laplacian_initial

print('Pressure sat initial: ', pressure_saturated_initial)
print('Moisture particle initial at x = 0: ', moisture_particle[0,0])
print('Moisture gas initial at x = 0: ', moisture_gas[0,0], '\n\n')

plot_value = np.zeros(number_of_time_steps)

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
                moisture_particle_current, alpha_parameter, N, relative_humidity[x], dt, constant[x], x)
            moisture_gas[y, x + 1] = compute_moisture_gas(
                moisture_particle_current, moisture_gas_current, alpha_parameter, N, relative_humidity[x], dt,
                constant[x], gas_velocity,
                moisture_diffusivity, gradient_moisture, laplacian_moisture, gas_density, particle_density,
                porosity_powder, x)

            # if x == 0:
            #     print('Moisture particle at x = 0: ', moisture_particle[y, x])
            #     print('Moisture gas at x = 1: ', moisture_gas[y, x + 1])

            # temp_particle[y, x] = compute_temperature_particle(
            #     temp_particle_current, constant[x], dt, conductivity_particle, laplacian, particle_density, alpha_parameter,
            #     moisture_particle_current, relative_humidity[x], N, heat_of_vaporization, heat_transfer_coefficient,
            #     specific_surface_area, temp_gas_current, particle_heat_capacity, x)
            #
            # temp_gas[y, x] = compute_temperature_gas(
            #     temp_particle_current, constant[x], dt, conductivity_gas, laplacian, gas_density, alpha_parameter,
            #     moisture_particle_current, N, moisture_vapor_heat_capacity, relative_humidity[x], heat_transfer_coefficient,
            #     specific_surface_area, temp_gas_current, gas_heat_capacity, gas_velocity, temp_gradient, porosity_powder,
            #     particle_density, x)

            ############################ UPDATE PARAMETERS #############################################################
            # Use temperature to find saturated pressure
            pressure_saturated[x] = compute_p_saturated(A, B, temp_gas[y, x], C)

            # Use moisture gas and p_sat to find RH
            relative_humidity[x] = compute_relative_humidity_from_Y(
                molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[y, x], pressure_saturated[x])

            # Use RH, p_sat and T_gas to find c
            molar_concentration_moisture[x] = compute_molar_concentration(relative_humidity[x], pressure_saturated[x],
                                                                          R_gas_constant, temp_gas[y, x])

            # Use c and T_gas to find partial pressure
            partial_pressure_moisture[x] = compute_partial_pressure_moisture(molar_concentration_moisture[x],
                                                                             R_gas_constant, temp_gas[y, x])

            if x == 0:
                # print('RH computed: ', relative_humidity[x])
                # print(' \n')
                plot_value[t] = moisture_particle[y, x]

            mass_transfer_coefficient[x] = compute_mass_transfer_coefficient(
                moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density,
                flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[x])[3]

            constant[x] = mass_transfer_coefficient[x] * specific_surface_area * pressure_saturated[x]/pressure_ambient  # just some simplification

# print('Change in temp particles:\n', (temp_particle - temp_initial)[0, 0:5])
# print('Change in temp gas:\n', (temp_gas - temp_initial)[0, 0:5])

print('Change in moisture particles:\n', (moisture_particle - moisture_particle_initial)[0, 0:5])
print('Change in moisture gas:\n', (moisture_gas - moisture_gas_initial_in)[0, 0:5])
print('Moisture gas:', (moisture_gas)[0, 0:5])

# print(relative_humidity_gas_initial - moisture_gas)
# print(temp_particle - temp_initial)
# print(temp_gas[0, 0:5])
# print(relative_humidity)

plt.plot(np.arange(number_of_time_steps), plot_value)
plt.show()