print('\n############################################ PROGRAM STARTED ############################################')
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from COND_main_functions import *
from COND_plot import *
import time
from scipy.integrate import odeint

################################## CHOOSE DISCRETIZATION ###############################################################
max_time = 1.5 * 60 * 60                                # hours to seconds
n_space_steps = 10                                      # MUST BE EVEN NUMBER
n_height_steps = int(n_space_steps/2)
resolution = 5000                                       # Number of outputs, TODO: more for nn, k estimation?
height_of_interest = 3                                  # Where do we want output
n_features = 6                                          # X_C, X_A, Y, T_G, T_P, AM
n_rotations = int(max_time/rotation_time_interval)
n_rotations = 0

######################################### SETUP ########################################################################
values_per_feature  = n_space_steps * n_height_steps
space_step_size     = bed_length / n_space_steps
discrete_time       = np.linspace(0, max_time, resolution)

moisture_gas_initial_all            = np.zeros((n_height_steps, n_space_steps)) + moisture_gas_initial_bed
moisture_particle_initial_cryst_all = np.zeros((n_height_steps, n_space_steps)) + moisture_cryst_particle_initial
moisture_particle_initial_cryst_all = np.zeros((n_height_steps, n_space_steps))
moisture_particle_initial_am_all    = np.zeros((n_height_steps, n_space_steps)) + moisture_am_particle_initial
moisture_particle_initial_am_all    = np.zeros((n_height_steps, n_space_steps))
temp_gas_initial_all                = np.zeros((n_height_steps, n_space_steps)) + temp_initial
temp_particle_initial_all           = np.zeros((n_height_steps, n_space_steps)) + temp_initial
amorphous_material_initial_all      = np.zeros((n_height_steps, n_space_steps)) + amorphous_material_initial

########################################## COMPUTE #####################################################################
initial_system = np.concatenate(
    [moisture_gas_initial_all.flatten(), moisture_particle_initial_cryst_all.flatten(), moisture_particle_initial_am_all.flatten(),
     temp_gas_initial_all.flatten(), temp_particle_initial_all.flatten(), amorphous_material_initial_all.flatten()])

if n_rotations > 0:
    n_time_outputs_per_rotation = int(resolution / n_rotations)
    discrete_time = discrete_time[0:(n_rotations*n_time_outputs_per_rotation)]
    computed_system = np.zeros([n_rotations, n_time_outputs_per_rotation, n_features * values_per_feature])
else:
    computed_system = np.zeros([resolution, n_features * values_per_feature])

tabs = 40
print('Initial moisture cryst powder:'.ljust(tabs), '{:.3e}'.format(moisture_cryst_particle_initial))
print('Initial moisture am powder:'.ljust(tabs), '{:.3f}'.format(moisture_am_particle_initial))
print(f'Saturated moisture powder:'.ljust(tabs), '{:.3e}'.format(moisture_cryst_particle_saturated))
print(f'Saturated moisture am powder:'.ljust(tabs), '{:.3f}'.format(moisture_am_particle_saturated))

print(f'Initial moisture gas:'.ljust(tabs), '{:.3f}'.format(moisture_gas_initial_bed))
print(f'Saturated moisture gas:'.ljust(tabs), '{:.3f}'.format(moisture_gas_initial_in))

print('\n       ***        STARTING COMPUTATION       ***        ')
run_time_start = time.time()

if n_rotations > 0:
    for rotation in range(n_rotations):
        print('Computing rotation:'.ljust(tabs), rotation + 1, '/', n_rotations)

        start_time = rotation * rotation_time_interval
        discrete_time_division = np.linspace(start_time, start_time + rotation_time_interval, n_time_outputs_per_rotation)
        computed_system[rotation, :] = odeint(conditioning_column, initial_system, discrete_time_division,
                                              args=(space_step_size, n_space_steps, n_height_steps))

        ################### initial_system = average ###################
        initial_system = computed_system[rotation, -1, :]
        for feature in range(n_features):
            avg = np.average(computed_system[rotation, -1, (feature * values_per_feature):((feature + 1) * values_per_feature)])
            initial_system[(feature * values_per_feature):((feature + 1) * values_per_feature)] = avg
else:
    computed_system = odeint(conditioning_column, initial_system, discrete_time,
                             args=(space_step_size, n_space_steps, n_height_steps))#, mxstep=5000)
    for feature in range(n_features):
        avg = np.average(computed_system[-1, (feature * values_per_feature):((feature + 1) * values_per_feature)])

elapsed = time.time() - run_time_start
print(f'       ***        COMPUTATION COMPLETE IN {elapsed:.2f} SECONDS       ***        \n')

########################################## SPLIT #######################################################################
if n_rotations > 0:
    moisture_gas_vector =               computed_system[:, :, 0:values_per_feature]
    moisture_particle_cryst_vector =    computed_system[:, :, values_per_feature:(values_per_feature * 2)]
    moisture_particle_am_vector =       computed_system[:, :, (values_per_feature * 2):(values_per_feature * 3)]
    temp_gas_vector =                   computed_system[:, :, (values_per_feature * 3):(values_per_feature * 4)]
    temp_particle_vector =              computed_system[:, :, (values_per_feature * 4):(values_per_feature * 5)]
    amorphous_material_vector =         computed_system[:, :, (values_per_feature * 5):(values_per_feature * 6)]
else:
    moisture_gas_vector =               computed_system[:, 0:values_per_feature]
    moisture_particle_cryst_vector =    computed_system[:, values_per_feature:(values_per_feature * 2)]
    moisture_particle_am_vector =       computed_system[:, (values_per_feature * 2):(values_per_feature * 3)]
    temp_gas_vector =                   computed_system[:, (values_per_feature * 3):(values_per_feature * 4)]
    temp_particle_vector =              computed_system[:, (values_per_feature * 4):(values_per_feature * 5)]
    amorphous_material_vector =         computed_system[:, (values_per_feature * 5):(values_per_feature * 6)]

moisture_gas_vector =               moisture_gas_vector.reshape(-1, n_height_steps, n_space_steps)
moisture_particle_cryst_vector =    moisture_particle_cryst_vector.reshape(-1, n_height_steps, n_space_steps)
moisture_particle_am_vector =       moisture_particle_am_vector.reshape(-1, n_height_steps, n_space_steps)
temp_gas_vector =                   temp_gas_vector.reshape(-1, n_height_steps, n_space_steps)
temp_particle_vector =              temp_particle_vector.reshape(-1, n_height_steps, n_space_steps)
amorphous_material_vector =         amorphous_material_vector.reshape(-1, n_height_steps, n_space_steps)
glass_temp_vector = compute_glass_temp_mix( 1-moisture_particle_am_vector, glass_temp_lactose, glass_temp_water )
temp_diff = temp_particle_vector - glass_temp_vector

if n_rotations > 0:
    avg_amorphous_material = np.average(computed_system[-1, -1, (5 * values_per_feature):(6 * values_per_feature)])
else:
    avg_amorphous_material = np.average(computed_system[ -1, (5 * values_per_feature):(6 * values_per_feature)])

diff_heat_flow_powder = (temp_particle_vector[1:, :, :] - temp_particle_vector[:-1, :, :]) * particle_heat_capacity

############################################ PLOT ######################################################################
# Convert to easier-to-read units
seconds = max_time
hours = seconds / 3600
discrete_time /= 3600

max_temp_gas = np.max(temp_gas_vector)
max_temp_gas_index = np.where(temp_gas_vector == max_temp_gas)

max_temp_particle = np.max(temp_particle_vector)
max_moisture_gas = np.max(moisture_gas_vector)
max_moisture_particle_cryst = np.max(moisture_particle_cryst_vector)
max_moisture_particle_cryst_index = np.where(moisture_particle_cryst_vector == max_moisture_particle_cryst)

max_moisture_particle_am = np.max(moisture_particle_am_vector)
max_moisture_particle_am_index = np.where(moisture_particle_am_vector == max_moisture_particle_am)


print('\n############################################ RESULTS ############################################')
print('Time computed:'.ljust(tabs), '{:.0f} seconds = {:.2f} hours = {:.2f} days'.format(seconds, hours, hours/24))
print('Relative humidity:'.ljust(tabs), '{:.0f} %'.format(relative_humidity_gas_inlet * 100))
print('Temp:'.ljust(tabs), '{:.2f} degrees Celsius'.format(temp_initial - kelvin))
print('Temp walls:'.ljust(tabs), '{:.2f} degrees Celsius'.format(temp_walls - kelvin))
print('Mass powder:'.ljust(tabs), '{:.2f} g\n'.format(mass_powder))

print('Max humidity gas:'.ljust(tabs), '{:.4f}'.format(max_moisture_gas))
print('Saturated humidity gas at inlet:'.ljust(tabs), '{:.4f}'.format(moisture_gas_initial_in), '\n')

print('Max humidity in cryst particles:'.ljust(tabs), '{:.4f}'.format(max_moisture_particle_cryst))
print('At time, height, length:'.ljust(tabs), max_moisture_particle_cryst_index[0],
      max_moisture_particle_cryst_index[1], max_moisture_particle_cryst_index[2])
print('Saturated humidity cryst particles:'.ljust(tabs), '{:.4f}'.format(moisture_cryst_particle_saturated))

print('Max humidity in am particles:'.ljust(tabs), '{:.4f}'.format(max_moisture_particle_am))
print('At time, height, length:'.ljust(tabs), max_moisture_particle_am_index[0],
      max_moisture_particle_am_index[1], max_moisture_particle_am_index[2])
print('Here, temp in particles is:'.ljust(tabs), temp_particle_vector[max_moisture_particle_am_index] - kelvin)

print('Saturated humidity am particles:'.ljust(tabs), '{:.4f} \n'.format(moisture_am_particle_saturated))

print('Max temperature gas:'.ljust(tabs), '{:.2f} degrees Celcius'.format(max_temp_gas - kelvin))
print('Max temperature particles:'.ljust(tabs), '{:.2f} degrees Celcius\n'.format(max_temp_particle - kelvin))

print('Avg amorphous material is:'.ljust(tabs), '{:.2f} %\n'.format(avg_amorphous_material * 100))

#print('Avg:', init_test.avg_temp_particle_vector)


#plot_average_temperature(average_temp_particles, kelvin, discrete_time)

#plot_average_moisture(average_moisture_particles, discrete_time)

# plot_heat_flow(diff_heat_flow_powder, discrete_time, hours)
#plot_heat_flow_slider(diff_heat_flow_powder, discrete_time, hours)

# plot_sections_over_time(
#     moisture_gas_vector, moisture_particle_cryst_vector, moisture_particle_am_vector, temp_gas_vector,
#     temp_particle_vector, amorphous_material_vector, height_of_interest, n_space_steps, discrete_time,
#     moisture_gas_initial_bed, moisture_gas_initial_in, moisture_cryst_particle_initial,
#     moisture_cryst_particle_saturated, moisture_am_particle_initial, moisture_am_particle_saturated, temp_min, kelvin,
#     hours, max_temp_gas, max_temp_particle)

# plot_heatmap(
#     moisture_gas_vector, moisture_particle_cryst_vector, temp_gas_vector, temp_particle_vector, height_of_interest,
#     n_space_steps, discrete_time, moisture_gas_initial_bed, moisture_gas_initial_in, moisture_am_particle_initial,
#     moisture_am_particle_saturated, temp_min, kelvin, hours, max_temp_gas, max_temp_particle)

slide_heat_map(
    moisture_gas_vector, moisture_particle_am_vector, temp_gas_vector, temp_particle_vector, amorphous_material_vector,
    temp_diff, temp_min, max_temp_particle, max_temp_gas, moisture_am_particle_initial,
    moisture_cryst_particle_saturated, moisture_am_particle_saturated, moisture_gas_initial_bed,
    moisture_gas_initial_in, amorphous_material_initial, hours)

############################################ SAVE DATA #################################################################
answer = 'N'
if answer == 'Y':
    np.save(f'RH-{relative_humidity_gas_inlet}_T-{temp_initial-kelvin}_TTg', temp_diff)
    np.save(f'RH-{relative_humidity_gas_inlet}_T-{temp_initial-kelvin}_AmMaterial', amorphous_material_vector)

print('\n############################################ PROGRAM ENDED ############################################')