print('\n############################################ PROGRAM STARTED ############################################')
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use('TkAgg')

from old_files.tam_sim import *
import time

################################## CHOOSE DISCRETIZATION ###############################################################
max_time = 4 * 60 * 60                                # hours to seconds
resolution = 5000                                       # Number of outputs, TODO: more for nn, k estimation?
n_features = 4                                          # X_C, X_A, Y, T_G, T_P, AM

######################################### SETUP ########################################################################
discrete_time       = np.linspace(0, max_time, resolution)
tabs = 40

moisture_gas_initial            = moisture_gas_initial_bed
moisture_particle_initial_cryst = moisture_cryst_particle_initial
print('Old moisture cryst initial'.ljust(tabs), '{:.5f}'.format(moisture_particle_initial_cryst))
# moisture_particle_initial_cryst = compute_GAB_equilibrium_moisture_cryst(relative_humidity_bed_initial)
# print('New moisture cryst initial'.ljust(tabs), '{:.5f}'.format(moisture_particle_initial_cryst))
# moisture_particle_initial_cryst = 0

moisture_particle_initial_am    = moisture_am_particle_initial
# moisture_particle_initial_am    = 0

# moisture_cryst_particle_saturated = compute_GAB_equilibrium_moisture_cryst(relative_humidity_gas_inlet)

temp_gas                        = temp_initial
temp_particle_initial           = temp_initial
amorphous_material_initial      = amorphous_material_initial
amorphous_material_initial      = 1
relative_humidity               = relative_humidity_gas_inlet

########################################## COMPUTE #####################################################################
initial_system = moisture_particle_initial_cryst, moisture_particle_initial_am, temp_particle_initial, amorphous_material_initial
                # change_m_particle_cryst, change_m_particle_am,                change_temp_particle,   change_amorphous_material
computed_system = np.zeros([resolution, n_features])


print('Initial moisture cryst powder:'.ljust(tabs), '{:.5f}'.format(moisture_cryst_particle_initial))
print('Initial moisture am powder:'.ljust(tabs), '{:.4f}'.format(moisture_am_particle_initial))
print(f'Saturated moisture cryst powder:'.ljust(tabs), '{:.5f}'.format(moisture_cryst_particle_saturated))
print(f'Saturated moisture am powder:'.ljust(tabs), '{:.4f}'.format(moisture_am_particle_saturated))

print('\n       ***        STARTING COMPUTATION       ***        ')
run_time_start = time.time()

computed_system = odeint(conditioning_column, initial_system, discrete_time, args=(temp_gas, relative_humidity))

elapsed = time.time() - run_time_start
print(f'       ***        COMPUTATION COMPLETE IN {elapsed:.2f} SECONDS       ***        \n')

########################################## SPLIT #######################################################################
moisture_particle_cryst_vector  = computed_system[:, 0]
moisture_particle_am_vector     = computed_system[:, 1]
temp_particle_vector            = computed_system[:, 2]
amorphous_material_vector       = computed_system[:, 3]

total_moisture_vector = moisture_particle_cryst_vector * (1 - amorphous_material_vector) + moisture_particle_am_vector * amorphous_material_vector

# amorphous_material_vector       = normalize_data(amorphous_material_vector)
cryst_start_index = np.where(amorphous_material_vector != 1)[0][0]
# index = cryst_start_index
print(cryst_start_index)

glass_temp_vector = compute_glass_temp_mix( 1-moisture_particle_am_vector, glass_temp_lactose, glass_temp_water_1 )
temp_diff = temp_particle_vector - glass_temp_vector

diff_heat_flow_powder = (temp_particle_vector - temp_initial) * heat_capacity_particle

############################################ PLOT ######################################################################
# Convert to easier-to-read units
seconds = max_time
hours = seconds / 3600
minutes = seconds / 60
discrete_time /= 60

fig, ax = plt.subplots(2, 3, figsize=(20, 10))
fig.suptitle(f'Moisture, temperature & amorphous material over time. '
             f'Total time: {int(minutes)} minutes. RH: {relative_humidity}, T: {temp_initial-kelvin}', fontsize=16)
ax[0, 0].set_title('Moisture Particle Am')
ax[0, 1].set_title('Moisture Particle Cryst')
ax[0, 2].set_title('Moisture Particle Total')
ax[1, 0].set_title('Amount amorphous')
ax[1, 1].set_title('Temp Particle')

# Set styles
moisture_color = 'navy'
temp_color = 'orangered'
initial_line = 'dashed'
initial_color = 'gray'
saturated_color = 'lightcoral'
amorph_color = 'green'

eps = 0.01

ax[0, 0].set_ylabel(f'M_am', rotation=0.5, size='large')
ax[0, 0].plot(discrete_time, moisture_particle_am_vector, c=moisture_color)
# patch = mpatches.Patch(color=moisture_color, label=f'Y {step}')
# ax[0, 0].set_ylim(moisture_particle_initial_am-eps, moisture_am_particle_saturated+eps)
ax[0, 0].hlines(moisture_particle_initial_am, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)
ax[0, 0].text(minutes * 4 / 5, moisture_gas_initial_bed, ('{:.4f}'.format(moisture_particle_initial_am)), ha='left', va='center')
ax[0, 0].hlines(moisture_am_particle_saturated, 0, discrete_time[-1], colors=saturated_color, linestyles=initial_line)
ax[0, 0].vlines(discrete_time[cryst_start_index], moisture_particle_initial_am-eps, moisture_am_particle_saturated+eps, colors=saturated_color, linestyles=initial_line)
# ax[0, 0].text(1, moisture_gas_initial_in, ('{:.4f}'.format(moisture_gas_initial_in)), ha='left', va='center')

# ax[0, 1].set_ylabel(f'M_cr', rotation=0, size='large')
# ax[0, 1].plot(discrete_time, moisture_particle_cryst_vector, c=moisture_color)
# # patch = mpatches.Patch(color=moisture_color, label=f'Y {step}')
# ax[0, 1].set_ylim(moisture_particle_initial_cryst-eps*0.001, moisture_cryst_particle_saturated+eps*0.001)
# ax[0, 1].hlines(moisture_particle_initial_cryst, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)
# ax[0, 1].text(minutes * 4 / 5, moisture_gas_initial_bed, ('{:.4f}'.format(moisture_particle_initial_cryst)), ha='left', va='center')
# ax[0, 1].hlines(moisture_cryst_particle_saturated, 0, discrete_time[-1], colors=saturated_color, linestyles=initial_line)
# ax[0, 1].vlines(discrete_time[cryst_start_index], 0, 1, colors=saturated_color, linestyles=initial_line)

ax[0, 1].set_ylabel(f'Heat flow', rotation=0, size='large')
ax[0, 1].plot(discrete_time, diff_heat_flow_powder, c=temp_color)
# ax[0, 1].set_ylim(moisture_particle_initial_cryst-eps*0.001, moisture_cryst_particle_saturated+eps*0.001)
# ax[0, 1].hlines(moisture_particle_initial_cryst, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)

ax[0, 2].set_ylabel(f'M_tot', rotation=0, size='large')
ax[0, 2].plot(discrete_time, total_moisture_vector, c=moisture_color)
ax[0, 2].set_ylim(moisture_particle_initial_cryst-eps*0.001, moisture_am_particle_saturated+eps*0.001)
ax[0, 2].hlines(moisture_particle_initial_cryst, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)
ax[0, 2].hlines(moisture_cryst_particle_saturated, 0, discrete_time[-1], colors=saturated_color, linestyles=initial_line)
ax[0, 2].hlines(moisture_am_particle_saturated, 0, discrete_time[-1], colors=saturated_color, linestyles=initial_line)
ax[0, 2].vlines(discrete_time[cryst_start_index], 0, 1, colors=saturated_color, linestyles=initial_line)

ax[1, 0].set_ylabel(f'Am g/g', rotation=0, size='large')
ax[1, 0].plot(discrete_time, amorphous_material_vector, c=amorph_color)
ax[1, 0].set_ylim(0-eps, 1+eps)
ax[1, 0].hlines(amorphous_material_initial, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)
ax[1, 0].vlines(discrete_time[cryst_start_index], -1, 1, colors=saturated_color, linestyles=initial_line)

ax[1, 1].set_ylabel(f'T', rotation=0, size='large')
ax[1, 1].plot(discrete_time, temp_particle_vector-kelvin, c=temp_color)
# ax[1, 1].set_ylim(temp_initial-eps - kelvin, np.max(temp_particle_vector) + eps - kelvin)
ax[1, 1].hlines(temp_initial-kelvin, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)
ax[1, 1].vlines(discrete_time[cryst_start_index], 0, 30, colors=saturated_color, linestyles=initial_line)

ax[1, 2].set_ylabel(f'T - Tg', rotation=0, size='large')
ax[1, 2].plot(discrete_time, temp_diff, c=temp_color)
# ax[1, 2].set_ylim(temp_initial-eps - kelvin, np.max(temp_particle_vector) + eps - kelvin)
ax[1, 2].hlines(temp_initial-kelvin, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)
ax[1, 2].vlines(discrete_time[cryst_start_index], 0, 30, colors=saturated_color, linestyles=initial_line)


plt.show()

print('\n############################################ PROGRAM ENDED ############################################')