import numpy as np
import matplotlib.pyplot as plt
from define_functions import*
from test_main import*

space_steps = 10
space_step = bed_length / space_steps


time_step = space_step ** 2 / 2.1            # max /2, here /2.1

time_steps_per_space_step = int((space_step/gas_velocity)/time_step)+1
time_steps = time_steps_per_space_step + 55

print('time step: ', time_step)
print('space step: ', space_step)
print('time steps per space step: ', time_steps_per_space_step)
print('Length per step: ', gas_velocity*time_step)
print('Gas velocity: ', gas_velocity)       # 5 mm per s
specific_surface_area = 400

################################### SET INITIAL CONDITIONS #############################################################
relative_humidity = np.zeros((time_steps, space_steps)) + relative_humidity_bed_initial #+0.5
relative_humidity[:, 0] = relative_humidity_gas_initial                    # whole first section of tube is filled with more humid gas
print(relative_humidity[0:5, :])

moisture_particle = np.zeros((time_steps, space_steps)) + moisture_particle_initial
moisture_gas = np.zeros((time_steps, space_steps)) + moisture_gas_initial_bed

temp_particle = np.zeros((time_steps, space_steps)) + temp_initial
temp_gas = np.zeros((time_steps, space_steps)) + temp_initial


pressure_saturated = np.zeros(space_steps) + pressure_saturated_initial
partial_pressure_moisture = np.zeros(space_steps) + partial_pressure_moisture_initial
molar_concentration_moisture = np.zeros(space_steps) + molar_concentration_moisture_initial

mass_transfer_coefficient = np.zeros(space_steps) + k_GP_initial
heat_transfer_coefficient = heat_transfer_coefficient_initial
constant = np.zeros(space_steps) + constant_initial

gradient_moisture = gradient_moisture_initial
gradient_temp = gradient_moisture_initial
laplacian_moisture = laplacian_initial
laplacian_temp = laplacian_initial

# print('Pressure sat initial: ', pressure_saturated_initial)
# print('Moisture particle initial at x = 0: ', moisture_particle[0,0])
# print('Moisture gas initial at x = 1: ', moisture_gas[0, 1])
# print('Moisture gas initial at x = 0: ', moisture_gas[0, 0], '\n\n')

plot_value = np.zeros(time_steps)

for t in range(time_steps):
    gas_distance = gas_velocity * t * time_step
    gas_position = int(gas_distance / space_step)+1
    number_x = min(space_steps, gas_position)
    # print('number_x', number_x)
    # print('Gas position: ', gas_distance)

    for x in range(1, number_x):
        # print('\nx is: ', x)

        temp_particle_current = temp_particle[t-1, x]
        temp_gas_current = temp_gas[t-1, x]
        moisture_particle_current = moisture_particle[t-1, x]
        moisture_gas_current = moisture_gas[t-1, x]
        relative_humidity_current = relative_humidity[t-1, x-1]
        # print('RH is: ', relative_humidity_current)

        laplacian_moisture = compute_laplacian(moisture_gas[t, :], x, space_step)
        gradient_moisture = compute_gradient(moisture_gas[t, :], x, space_step)

        moisture_particle[t, x] = compute_moisture_particle(
            moisture_particle_current, alpha_parameter, N, relative_humidity_current, time_step, constant[x], x,
            verbose=False)

        moisture_gas[t, x] = compute_moisture_gas(
            moisture_particle_current, moisture_gas_current, alpha_parameter, N, relative_humidity_current, time_step,
            constant[x], gas_velocity,
            moisture_diffusivity, gradient_moisture, laplacian_moisture, gas_density, particle_density,
            porosity_powder, x, verbose=True)

        # if x == 0:
        #     print('Moisture particle at x = 0: ', moisture_particle[y, x])
        #     print('Moisture gas at x = 1: ', moisture_gas[y, x + 1])

        temp_particle[t, x] = compute_temperature_particle(
            temp_particle_current, constant[x], time_step, conductivity_particle, laplacian_temp, particle_density,
            alpha_parameter,
            moisture_particle_current, relative_humidity_current, N, heat_of_vaporization, heat_transfer_coefficient,
            specific_surface_area, temp_gas_current, particle_heat_capacity, x, verbose=False)

        temp_gas[t, x] = compute_temperature_gas(
            temp_particle_current, constant[x], time_step, conductivity_gas, laplacian_temp, gas_density,
            alpha_parameter, moisture_particle_current, N, moisture_vapor_heat_capacity, relative_humidity_current,
            heat_transfer_coefficient, specific_surface_area, temp_gas_current, gas_heat_capacity, gas_velocity,
            gradient_temp, porosity_powder, particle_density, x, verbose=True)

        ############################ UPDATE PARAMETERS #############################################################
        # Use temperature to find saturated pressure
        pressure_saturated[x] = compute_p_saturated(A, B, temp_gas[t, x], C)
        # print(pressure_saturated[x])

        # Use moisture gas and p_sat to find RH
        relative_humidity[t, x] = compute_relative_humidity_from_Y(
            molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[t, x], pressure_saturated[x])

        # Use RH, p_sat and T_gas to find c
        molar_concentration_moisture[x] = compute_molar_concentration(relative_humidity[t, x], pressure_saturated[x],
                                                                      R_gas_constant, temp_gas[t, x])

        # Use c and T_gas to find partial pressure
        partial_pressure_moisture[x] = compute_partial_pressure_moisture(molar_concentration_moisture[x],
                                                                         R_gas_constant, temp_gas[t, x])

        if x == 1:
            print('RH computed: ', relative_humidity[t, x+1])
            # print(' \n')
            plot_value[t] = moisture_particle[t, x]

        mass_transfer_coefficient[x] = compute_mass_transfer_coefficient(
            moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density,
            flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[x])[3]

        constant[x] = mass_transfer_coefficient[x] * specific_surface_area * pressure_saturated[x] / pressure_ambient

print('\nRESULTS:')
# print('Change in temp particles:\n', (temp_particle - temp_initial)[0, 0:5])
# print('Change in temp gas:\n', (temp_gas - temp_initial)[0, 0:5])

t = time_steps_per_space_step
# print('RH\n %.2f' % relative_humidity[t-1:(t+50), 0:3])
print('RH\n', np.around(relative_humidity[t-1:(t+25), 0:3],3))
# print('Change in moisture particles:\n', (moisture_particle - moisture_particle_initial)[t:t+10, 1])
# print('Change in moisture gas:\n', (moisture_gas - moisture_gas_initial_in)[t:t+10, 1])
# print('Moisture particle:\n', (moisture_particle)[t:t+10, 1])
# print('Moisture gas:\n', (moisture_gas)[t:t+10, 1])

# print(relative_humidity_gas_initial - moisture_gas)
# print(temp_particle - temp_initial)
# print('Temp gas is: ', temp_gas)
# print(relative_humidity)

# plt.plot(np.arange(number_of_time_steps), plot_value)
# plt.show()