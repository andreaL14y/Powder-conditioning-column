import numpy as np
import matplotlib.pyplot as plt
from define_functions import*
from test_main import*

space_steps = 10
space_step = bed_length / space_steps

time_step = space_step ** 2 / 2                                               # max /2, here /2.1

time_steps_per_space_step = int((space_step/gas_velocity)/time_step)+1
time_steps = 15000

print('time step: ', time_step)
print('space step: ', space_step)
print('time steps per space step: ', time_steps_per_space_step)
print('Length per step: ', gas_velocity*time_step)
print('Gas velocity: ', gas_velocity)
print('Total time computed: ', time_steps * time_step, 'seconds.')
specific_surface_area = 400

################################### SET INITIAL CONDITIONS #############################################################
relative_humidity = np.zeros((time_steps, space_steps)) + relative_humidity_bed_initial
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

moisture_full_humidity = compute_initial_moisture_particle(alpha_parameter, N, 0.9)
print('Particle moisture content at 0.9 RH: ', moisture_full_humidity)

################################### SET BOUNDARY CONDITIONS ############################################################
relative_humidity[:, 0] = relative_humidity_gas_initial                 # x = 0 is before tube starts, humid gas
# relative_humidity[0, 1] = relative_humidity_gas_initial                 # at the start, RH is higher in x = 1
relative_humidity[:, space_steps-1] = 0                                 # 0 humidity at end of tube
moisture_gas[:, 0] = moisture_gas_initial_in
# moisture_gas[0, 1] = moisture_gas_initial_in                            # at the start, Y is higher in x = 1
moisture_gas[:, space_steps-1] = 0

# print('Pressure sat initial: ', pressure_saturated_initial)
print('Moisture particle initial at x = 0: ', moisture_particle[0,0])
print('Moisture particle initial at x = 1: ', moisture_particle[0,1])
print('Moisture gas initial at x = 0: ', moisture_gas[0, 0])
print('Moisture gas initial at x = 1: ', moisture_gas[0, 1], '\n\n')

plot_particle = np.zeros((time_steps, space_steps))
plot_gas = np.zeros((time_steps, space_steps))

for t in range(time_steps-1):
    gas_distance = gas_velocity * t * time_step
    gas_position = int(gas_distance / space_step)+1
    number_x = min(space_steps, gas_position)
    # print('number_x', number_x)
    # print('Gas position: ', gas_distance)

    # for x in range(1, number_x+1):
    for x in range(1, space_steps):
        temp_particle_current = temp_particle[t, x]
        temp_gas_current = temp_gas[t, x]
        moisture_particle_current = moisture_particle[t, x]
        moisture_gas_current = moisture_gas[t, x]
        relative_humidity_current = relative_humidity[t, x]
        laplacian_moisture = compute_laplacian(moisture_gas[t, :], x, space_step)
        gradient_moisture = compute_gradient(moisture_gas[t, :], x, space_step)

        # if x ==1:
        # print('\nx is: ', x)
        #     print('Laplacian: ', laplacian_moisture)
        #     print('Gradient: ', gradient_moisture)

        moisture_particle[t+1, x] = compute_moisture_particle(
            moisture_particle_current, alpha_parameter, N, relative_humidity_current, time_step, constant[x], x,
            verbose=False)

        moisture_gas[t+1, x] = compute_moisture_gas(
            moisture_particle_current, moisture_gas_current, alpha_parameter, N, relative_humidity_current, time_step,
            constant[x], gas_velocity,
            moisture_diffusivity, gradient_moisture, laplacian_moisture, gas_density, particle_density,
            porosity_powder, x, verbose=False)

        # if x == 0:
        #     print('Moisture particle at x = 0: ', moisture_particle[y, x])
        #     print('Moisture gas at x = 1: ', moisture_gas[y, x + 1])

        # temp_particle[t, x] = compute_temperature_particle(
        #     temp_particle_current, constant[x], time_step, conductivity_particle, laplacian_temp, particle_density,
        #     alpha_parameter,
        #     moisture_particle_current, relative_humidity_current, N, heat_of_vaporization, heat_transfer_coefficient,
        #     specific_surface_area, temp_gas_current, particle_heat_capacity, x, verbose=False)
        #
        # temp_gas[t, x] = compute_temperature_gas(
        #     temp_particle_current, constant[x], time_step, conductivity_gas, laplacian_temp, gas_density,
        #     alpha_parameter, moisture_particle_current, N, moisture_vapor_heat_capacity, relative_humidity_current,
        #     heat_transfer_coefficient, specific_surface_area, temp_gas_current, gas_heat_capacity, gas_velocity,
        #     gradient_temp, porosity_powder, particle_density, x, verbose=True)

        ############################ UPDATE PARAMETERS #############################################################
        # Use temperature to find saturated pressure
        pressure_saturated[x] = compute_p_saturated(A, B, temp_gas[t+1, x], C)
        # print(pressure_saturated[x])

        # Use moisture gas and p_sat to find RH
        relative_humidity[t+1, x] = compute_relative_humidity_from_Y(
            molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[t+1, x], pressure_saturated[x])

        # Use RH, p_sat and T_gas to find c
        molar_concentration_moisture[x] = compute_molar_concentration(relative_humidity[t+1, x], pressure_saturated[x],
                                                                      R_gas_constant, temp_gas[t+1, x])

        # Use c and T_gas to find partial pressure
        partial_pressure_moisture[x] = compute_partial_pressure_moisture(molar_concentration_moisture[x],
                                                                         R_gas_constant, temp_gas[t+1, x])

        plot_particle[t, x] = moisture_particle[t, x]
        plot_gas[t, x] = moisture_gas[t, x]

        mass_transfer_coefficient[x] = compute_mass_transfer_coefficient(
            moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density,
            flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[x])[3]

        constant[x] = mass_transfer_coefficient[x] * specific_surface_area * pressure_saturated[x] / pressure_ambient

print('\nRESULTS:')
# print('Change in temp particles:\n', (temp_particle - temp_initial)[0, 0:5])
# print('Change in temp gas:\n', (temp_gas - temp_initial)[0, 0:5])

# t = time_steps_per_space_step
# print('RH\n %.2f' % relative_humidity[t-1:(t+50), 0:3])
# print('RH\n', np.around(relative_humidity[t-1:(t+25), 0:3],3))
# print('Change in moisture particles:\n', (moisture_particle - moisture_particle_initial)[t:t+10, 1])
# print('Change in moisture gas:\n', (moisture_gas - moisture_gas_initial_in)[t:t+10, 1])
print('Moisture particle:\n', (moisture_particle)[0:t, 1])
# print('Moisture gas x=0:\n', (moisture_gas)[t+950:t+1000, 0])
print('Moisture gas x=1:\n', (moisture_gas)[0:t, 1])

# print(relative_humidity_gas_initial - moisture_gas)
# print(temp_particle - temp_initial)
# print('Temp gas is: ', temp_gas)
print('RH: \n', relative_humidity[0:t, 1])
x = len(plot_gas)
plt.plot(np.arange(x)[:t-1], plot_particle[:t-1, 1])
plt.plot(np.arange(x)[:t-1], plot_particle[:t-1, 2])
# plt.plot(np.arange(x)[:t-1], plot_gas[:t-1])
plt.show()