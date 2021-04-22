from scipy.integrate import odeint
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from v_ode_functions_2 import *
from plot_ode_functions_2 import *

################################## CHOOSE DISCRETIZATION ###############################################################
max_time = 400000                       # seconds
max_time = 2.5 * 60 * 60                       # seconds
n_space_steps = 10                      # MUST BE EVEN NUMBER
n_height_steps = int(n_space_steps/2)
resolution = 1000
height_of_interest = 3
n_features = 4

######################################### SETUP ########################################################################
values_per_feature = n_space_steps * n_height_steps
space_step_size = bed_length / n_space_steps
discrete_time = np.linspace(0, max_time, resolution)

moisture_gas_initial_all = np.zeros((n_height_steps, n_space_steps)) + moisture_gas_initial_bed
moisture_particle_initial_all = np.zeros((n_height_steps, n_space_steps)) + moisture_particle_initial
temp_gas_initial = np.zeros((n_height_steps, n_space_steps)) + temp_initial
temp_particle_initial = np.zeros((n_height_steps, n_space_steps)) + temp_initial

########################################## COMPUTE #####################################################################
initial_system = np.concatenate(
    [moisture_gas_initial_all.flatten(), moisture_particle_initial_all.flatten(),
     temp_gas_initial.flatten(), temp_particle_initial.flatten()])

n_rotations = int(max_time/rotation_time_interval)
n_time_outputs_per_rotation = int(resolution / n_rotations)
discrete_time = discrete_time[0:(n_rotations*n_time_outputs_per_rotation)]
print("Number of rotations are:", n_rotations)
print("Number of time outputs per rotation is:", n_time_outputs_per_rotation)
print("Shape of time is:", discrete_time.shape)
#print("Time is:", discrete_time)

computed_system = np.zeros([n_rotations, n_time_outputs_per_rotation, 4 * values_per_feature])

for rotation in range(n_rotations):
    print(rotation)
    start_time = rotation * rotation_time_interval
    discrete_time_division = np.linspace(start_time, start_time + rotation_time_interval, n_time_outputs_per_rotation)
    computed_system[rotation, :] = odeint(conditioning_column, initial_system, discrete_time_division, args=(space_step_size, n_space_steps, n_height_steps))

    #initial_system = average
    initial_system = computed_system[rotation, -1, :]
    for feature in range(n_features):
        avg = np.average(computed_system[rotation, -1, (feature * values_per_feature):((feature + 1) * values_per_feature)])
        initial_system[(feature * values_per_feature):((feature + 1) * values_per_feature)] = avg

########################################## SPLIT #######################################################################
moisture_gas_vector = computed_system[:, :, 0:values_per_feature]
moisture_particle_vector = computed_system[:, :, values_per_feature:(values_per_feature * 2)]
temp_gas_vector = computed_system[:, :, (values_per_feature * 2):(values_per_feature * 3)]
temp_particle_vector = computed_system[:, :, (values_per_feature * 3):(values_per_feature * 4)]

moisture_gas_vector = moisture_gas_vector.reshape(-1, n_height_steps, n_space_steps)
print("Shape moisture gas vector:", moisture_gas_vector.shape)
moisture_particle_vector = moisture_particle_vector.reshape(-1, n_height_steps, n_space_steps)
temp_gas_vector = temp_gas_vector.reshape(-1, n_height_steps, n_space_steps)
temp_particle_vector = temp_particle_vector.reshape(-1, n_height_steps, n_space_steps)

############################################ PLOT ######################################################################
# Convert to easier-to-read units
seconds = max_time
hours = seconds / 3600
discrete_time /= 3600

max_temp_gas = np.max(temp_gas_vector)
max_temp_gas_index = np.where(temp_gas_vector == max_temp_gas)

max_temp_particle = np.max(temp_particle_vector)
max_moisture_gas = np.max(moisture_gas_vector)
max_moisture_particle = np.max(moisture_particle_vector)
max_moisture_particle_index = np.where(moisture_particle_vector == max_moisture_particle)

print(f'Time computed is: {seconds} seconds = {int(hours)} hours = {int(hours/24)} days\n')
print('Max humidity in gas is: {:.4f}'.format(max_moisture_gas))
print('Saturated humidity in gas at inlet is: {:.4f}\n'.format(moisture_gas_initial_in))

print('Max humidity in particles is: {:.4f}'.format(max_moisture_particle))
print('This happens at time, height, length: ', max_moisture_particle_index[0], max_moisture_particle_index[1], max_moisture_particle_index[2])
print('Saturated humidity in particles is: {:.4f} \n'.format(moisture_particle_saturated))

print('Max temperature in gas is: {:.4f} degrees Celcius'.format(max_temp_gas - kelvin))
print('Max temperature in particles is: {:.4f} degrees Celcius\n'.format(max_temp_particle - kelvin))

plot_sections_over_time(
    moisture_gas_vector, moisture_particle_vector, temp_gas_vector, temp_particle_vector, height_of_interest,
    n_space_steps, discrete_time, moisture_gas_initial_bed, moisture_gas_initial_in, moisture_particle_initial,
    moisture_particle_saturated, temp_min, kelvin, hours, max_temp_gas, max_temp_particle)

# plot_heatmap(
#     moisture_gas_vector, moisture_particle_vector, temp_gas_vector, temp_particle_vector, height_of_interest,
#     n_space_steps, discrete_time, moisture_gas_initial_bed, moisture_gas_initial_in, moisture_particle_initial,
#     moisture_particle_saturated, temp_min, kelvin, hours, max_temp_gas, max_temp_particle)
#
# slide_heat_map(
#     moisture_gas_vector, moisture_particle_vector, temp_gas_vector, temp_particle_vector, temp_min, max_temp_particle,
#     max_temp_gas, moisture_particle_initial, moisture_particle_saturated, moisture_gas_initial_bed, moisture_gas_initial_in, hours)