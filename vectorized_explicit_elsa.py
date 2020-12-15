import numpy as np
import matplotlib.pyplot as plt
from define_functions import*
from vectorized_define_functions import gas_velocity, moisture_particle_initial, moisture_gas_initial_bed, \
    pressure_saturated_initial, partial_pressure_moisture_initial, molar_concentration_moisture_initial, k_GP_initial, \
    heat_transfer_coefficient, constant_initial, moisture_gas_initial_in, flow_rate, superficial_velocity
from vectorized_define_functions_OLD import*
from input_parameters import*

space_steps = 5
space_step = bed_length / space_steps

time_step = space_step ** 2 / 2                                               # max /2, here /2.1

time_steps_per_space_step = int((space_step/gas_velocity)/time_step)+1
time_steps = 20

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
heat_transfer_coefficient = heat_transfer_coefficient
constant = np.zeros(space_steps) + constant_initial

moisture_full_humidity = compute_initial_moisture_particle(alpha_parameter, N, 0.9)
print('Particle moisture content at 0.9 RH: ', moisture_full_humidity)

################################### SET BOUNDARY CONDITIONS ############################################################
relative_humidity[:, 0] = relative_humidity_gas_initial                 # x = 0 is before tube starts, humid gas
# relative_humidity[0, 1] = relative_humidity_gas_initial                 # at the start, RH is higher in x = 1
relative_humidity[:, space_steps - 1] = 0  # 0 humidity at end of tube
moisture_gas[:, 0] = moisture_gas_initial_in
# moisture_gas[0, 1] = moisture_gas_initial_in                            # at the start, Y is higher in x = 1
moisture_gas[:, space_steps - 1] = 0

# print('Pressure sat initial: ', pressure_saturated_initial)
print('Moisture particle initial at x = 0: ', moisture_particle[0, 0])
print('Moisture particle initial at x = 1: ', moisture_particle[0, 1])
print('Moisture gas initial at x = 0: ', moisture_gas[0, 0])
print('Moisture gas initial at x = 1: ', moisture_gas[0, 1], '\n\n')

plot_particle = np.zeros((time_steps, space_steps))
plot_gas = np.zeros((time_steps, space_steps))
plot_particle_temp = np.zeros((time_steps, space_steps))
plot_gas_temp = np.zeros((time_steps, space_steps))

plot_RH = np.zeros((time_steps, space_steps))

print('Total time computed: ', time_steps * time_step, 'seconds.')
for t in range(time_steps - 1):
    if t*time_step % 20 == 0:
        print('Time elapsed:', t*time_step)
    temp_particle_current = temp_particle[t, 1:]
    temp_gas_current = temp_gas[t, 1:]
    moisture_particle_current = moisture_particle[t, 1:]
    moisture_gas_current = moisture_gas[t, 1:]
    relative_humidity_current = relative_humidity[t, 1:]

    laplacian_moisture = compute_laplacian_moisture_vector(moisture_gas[t, :], space_step, moisture_gas_initial_in)
    gradient_moisture = compute_gradient_moisture_vector(moisture_gas[t, :], space_step, moisture_gas_initial_in)

    moisture_particle[t + 1, 1:] = compute_moisture_particle_vector(
        moisture_particle_current, alpha_parameter, N, relative_humidity_current, time_step, constant[1:], verbose=False)

    moisture_gas[t + 1, 1:] = compute_moisture_gas_vector(
        moisture_particle_current, moisture_gas_current, alpha_parameter, N, relative_humidity_current, time_step,
        constant[1:], gas_velocity, moisture_diffusivity, gradient_moisture[1:], laplacian_moisture[1:], gas_density,
        particle_density, porosity_powder, verbose=False)

    laplacian_temp = compute_laplacian_temp_vector(temp_gas[t, :], space_step)
    gradient_temp = compute_gradient_temp_vector(temp_gas[t, :], space_step)

    temp_particle[t + 1, 1:] = compute_temperature_particle_vector(
        temp_particle_current, constant[1:], time_step, conductivity_particle, laplacian_temp[1:], particle_density,
        alpha_parameter,
        moisture_particle_current, relative_humidity_current, N, heat_of_vaporization, heat_transfer_coefficient,
        specific_surface_area, temp_gas_current, particle_heat_capacity, verbose=False)

    temp_gas[t + 1, 1:] = compute_temperature_gas_vector(
        temp_particle_current, constant[1:], time_step, conductivity_gas, laplacian_temp[1:], gas_density,
        alpha_parameter, moisture_particle_current, N, moisture_vapor_heat_capacity, relative_humidity_current,
        heat_transfer_coefficient, specific_surface_area, temp_gas_current, gas_heat_capacity, gas_velocity,
        gradient_temp[1:], porosity_powder, particle_density, verbose=True)

    ############################ UPDATE PARAMETERS #############################################################
    # Use temperature to find saturated pressure
    pressure_saturated[1:] = compute_p_saturated_vector(A, B, temp_gas[t + 1, 1:], C)
    # print(pressure_saturated[x])

    # Use moisture gas and p_sat to find RH
    relative_humidity[t + 1, 1:] = compute_relative_humidity_from_Y_vector(
        molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[t + 1, 1:], pressure_saturated[1:])

    # Use RH, p_sat and T_gas to find c
    molar_concentration_moisture[1:] = compute_molar_concentration_vector(
        relative_humidity[t + 1, 1:], pressure_saturated[1:], R_gas_constant, temp_gas[t + 1, 1:])

    # Use c and T_gas to find partial pressure
    partial_pressure_moisture[1:] = compute_partial_pressure_moisture_vector(
        molar_concentration_moisture[1:], R_gas_constant, temp_gas[t + 1, 1:])

    plot_particle[t, :] = moisture_particle[t, :]
    plot_gas[t, :] = moisture_gas[t, :]

    plot_particle_temp[t, :] = temp_particle[t, :]
    plot_gas_temp[t, :] = temp_gas[t, :]

    plot_RH[t, :] = relative_humidity[t, :]

    mass_transfer_coefficient[1:] = compute_mass_transfer_coefficient_vector(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density,
        flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[1:])[3]

    constant = mass_transfer_coefficient * specific_surface_area * pressure_saturated / pressure_ambient

print('\nRESULTS:')
# print('Change in temp particles:\n', (temp_particle - temp_initial)[0, 0:5])
print('Change in temp gas:\n', (temp_gas - temp_initial)[0:50, 1])

# print('Change in moisture particles:\n', (moisture_particle - moisture_particle_initial)[t:t+10, 1])
# print('Change in moisture gas:\n', (moisture_gas - moisture_gas_initial_in)[t:t+10, 1])
print('Moisture particle:\n', (moisture_particle)[0:t, 1])
print('Moisture gas x=1:\n', (moisture_gas)[0:t, 1])

# print(relative_humidity_gas_initial - moisture_gas)
# print(temp_particle - temp_initial)
# print('Temp gas is: ', temp_gas)
print('RH: \n', relative_humidity[0:t, 1])

x = len(plot_gas)
fig, ax = plt.subplots(1, 4, figsize=(20, 13))
fig.suptitle(f'Time: {time_steps * time_step} s', fontsize=16)
# ax[0].plot(np.arange(x)[:t-1], plot_particle[:t-1, 2], label = 'Particle moisture x2')
ax[0].plot(np.arange(x)[:t-1], plot_particle[:t-1, 1], label = 'Particle moisture x1')

ax[1].plot(np.arange(x)[:t-1], plot_gas[:t-1, 1], label = 'Gas moisture x1')
# ax[1].plot(np.arange(x)[:t-1], plot_gas[:t-1, 2], label = 'Gas moisture x2')

ax[2].plot(np.arange(x)[:t-1], plot_particle_temp[:t-1, 1], label = 'Particle temp x1')
# ax[2].plot(np.arange(x)[:t-1], plot_particle_temp[:t-1, 2], label = 'Particle temp x2')
ax[2].set_ylim(0, 500)
# ax1.set_ylim([0, 5])

ax[3].plot(np.arange(x)[:t-1], plot_gas_temp[:t-1, 1], label = 'Gas temp x1')
# ax[3].plot(np.arange(x)[:t-1], plot_gas_temp[:t-1, 1], label = 'Gas temp x1')
ax[3].set_ylim(0, 500)

ax[0].legend(loc="upper left")
ax[1].legend(loc="upper left")
ax[2].legend(loc="upper left")
ax[3].legend(loc="upper left")

# plt.title(f'Time: {time_steps * time_step} s')
plt.xlabel('Time t')

# plt.plot(np.arange(x)[:t-1], plot_gas[:t-1])
plt.savefig('explicit_method_issues.pdf')
plt.show()