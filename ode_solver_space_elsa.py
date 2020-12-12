from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from vectorized_define_functions import *

################################## CHOOSE DISCRETIZATION ###############################################################
max_time = 1000000
n_space_steps = 3
n_height_steps = 2
resolution = 100

######################################### SETUP ########################################################################
values_per_feature = n_space_steps * n_height_steps
space_step_size = bed_length / n_space_steps
discrete_time = np.linspace(0, max_time, resolution)

# moisture_gas_initial_all = np.zeros(n_space_steps) + moisture_gas_initial_bed
moisture_gas_initial_all = np.zeros((n_height_steps, n_space_steps)) + moisture_gas_initial_bed

# moisture_particle_initial_all = np.zeros(n_space_steps) + moisture_particle_initial
moisture_particle_initial_all = np.zeros((n_height_steps, n_space_steps)) + moisture_particle_initial

# temp_gas_initial = np.zeros(n_space_steps) + temp_initial
temp_gas_initial = np.zeros((n_height_steps, n_space_steps)) + temp_initial

# temp_particle_initial = np.zeros(n_space_steps) + temp_initial
temp_particle_initial = np.zeros((n_height_steps, n_space_steps)) + temp_initial

########################################## COMPUTE #####################################################################
initial_system = np.concatenate(
    [moisture_gas_initial_all.flatten(), moisture_particle_initial_all.flatten(),
     temp_gas_initial.flatten(), temp_particle_initial.flatten()])

computed_system = odeint(
    conditioning_column, initial_system, discrete_time, args=(space_step_size, n_space_steps, n_height_steps))
print(np.shape(computed_system))
# shape is (resolution, space_steps * height_steps * 4)     features

moisture_gas_vector = computed_system[:, 0:values_per_feature]
moisture_particle_vector = computed_system[:, values_per_feature:(values_per_feature * 2)]
temp_gas_vector = computed_system[:, (values_per_feature * 2):(values_per_feature * 3)]
temp_particle_vector = computed_system[:, (values_per_feature * 3):(values_per_feature * 4)]

moisture_gas_vector = np.reshape(moisture_gas_vector, (resolution, n_height_steps, n_space_steps))
moisture_particle_vector = np.reshape(moisture_particle_vector, (resolution, n_height_steps, n_space_steps))
temp_gas_vector = np.reshape(temp_gas_vector, (resolution, n_height_steps, n_space_steps))
temp_particle_vector = np.reshape(temp_particle_vector, (resolution, n_height_steps, n_space_steps))

print(np.shape(moisture_gas_vector))
print(np.shape(moisture_particle_vector))
print(np.shape(temp_gas_vector))
print(np.shape(temp_particle_vector))

############################################ PLOT ######################################################################
# Convert to easier-to-read units
seconds = max_time
hours = seconds / 3600
discrete_time /= 3600
# computed_system[:, n_space_steps * 2:n_space_steps * 3] -= kelvin
temp_gas_vector -= kelvin
temp_particle_vector -= kelvin
max_temp_gas = np.max(temp_gas_vector)
max_temp_gas_index = np.where(temp_gas_vector == max_temp_gas)
# print(max_temp_gas_index, system[max_temp_gas_index])

# computed_system[:, n_space_steps * 3:n_space_steps * 4] -= kelvin
# max_temp_particle = np.max(computed_system[:, (n_space_steps * 3):(n_space_steps * 4)])
max_temp_particle = np.max(temp_particle_vector)

print(f'Time computed is: {seconds} seconds = {int(hours)} hours = {int(hours/24)} days.')
print('Max temperature in gas is: {:.4f} degrees Celcius'.format(max_temp_gas))
print('Max temperature in particles is: {:.4f} degrees Celcius'.format(max_temp_particle))

fig, ax = plt.subplots(n_space_steps, 4)
fig.suptitle(f'Moisture & temperature in cylinder sections over time. Total time: {int(hours)} hours.', fontsize=14)
ax[0, 0].set_title('Moisture Gas')
ax[0, 1].set_title('Moisture Particle')
ax[0, 2].set_title('Temp Gas')
ax[0, 3].set_title('Temp Particle')

# Set styles
moisture_color = 'navy'
temp_color = 'orangered'
initial_line = 'dashed'
initial_color = 'gray'
saturated_color = 'lightcoral'

# counter = 0
epsilon = 0.001
height_of_interest = 1
for feature in range(4):
    if feature == 0 or feature == 1:
        current_color = moisture_color
    else:
        # print('yo')
        current_color = temp_color

    for step in range(n_space_steps):
        counter = 0
        if feature == 0:
            ax[step, feature].set_ylabel(f'Section {step+1}                      ', rotation=0, size='large')
        # ax[step, feature].plot(discrete_time, computed_system[:, counter], c=current_color)

        if feature == 0:
            ax[step, feature].plot(discrete_time, moisture_gas_vector[:, height_of_interest, step], c=current_color)

            patch = mpatches.Patch(color=current_color, label=f'Y {step}')
            ax[step, feature].set_ylim(moisture_gas_initial_bed - epsilon, moisture_gas_initial_in + epsilon)
            ax[step, feature].hlines(moisture_gas_initial_bed, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)
            ax[step, feature].text(hours*4/5, moisture_gas_initial_bed+epsilon, ('{:.4f}'.format(moisture_gas_initial_bed)), ha='left', va='center')

            ax[step, feature].hlines(moisture_gas_initial_in, 0, discrete_time[-1], colors=saturated_color, linestyles=initial_line)
            ax[step, feature].text(1, moisture_gas_initial_in - epsilon, ('{:.4f}'.format(moisture_gas_initial_in)), ha='left', va='center')

            counter += 1
        elif feature == 1:
            ax[step, feature].plot(discrete_time, moisture_particle_vector[:, height_of_interest, step], c=current_color)

            patch = mpatches.Patch(color=current_color, label=f'X {step}')
            ax[step, feature].set_ylim(moisture_particle_initial - epsilon, moisture_particle_saturated + epsilon)
            ax[step, feature].hlines(moisture_particle_initial, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)
            ax[step, feature].text(hours * 4 / 5, moisture_particle_initial + 2*epsilon,
                                   ('{:.4f}'.format(moisture_particle_initial)), ha='left', va='center')

            ax[step, feature].hlines(moisture_particle_saturated, 0, discrete_time[-1], colors=saturated_color, linestyles=initial_line)
            ax[step, feature].text(1, moisture_particle_saturated - 2*epsilon, ('{:.4f}'.format(moisture_particle_saturated)),
                                   ha='left', va='center')
            counter += 1
        elif feature == 2:
            ax[step, feature].plot(discrete_time, temp_gas_vector[:, height_of_interest, step], c=current_color)

            patch = mpatches.Patch(color=current_color, label=f'{step}')
            ax[step, feature].set_ylim(temp_initial - (kelvin + epsilon), max_temp_gas + 1)
            ax[step, feature].hlines(temp_initial - kelvin, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)
            ax[step, feature].hlines(max_temp_gas, 0, discrete_time[-1], colors=saturated_color, linestyles=initial_line)
            ax[step, feature].text(hours * 4 / 5, max_temp_gas - 0.7,
                                   ('{:.2f}'.format(max_temp_gas)), ha='left', va='center')
            counter += 1
        else:
            ax[step, feature].plot(discrete_time, temp_particle_vector[:, height_of_interest, step], c=current_color)

            patch = mpatches.Patch(color=current_color, label=f'{step}')
            ax[step, feature].set_ylim(temp_initial - (kelvin + epsilon), max_temp_particle + 1)
            ax[step, feature].hlines(temp_initial - kelvin, 0, discrete_time[-1], colors=initial_color, linestyles=initial_line)

            ax[step, feature].hlines(max_temp_particle, 0, discrete_time[-1], colors=saturated_color, linestyles=initial_line)
            ax[step, feature].text(hours * 4 / 5, max_temp_particle - 0.7,
                                   ('{:.2f}'.format(max_temp_particle)), ha='left', va='center')
            # counter += 1

        # ax[step, feature].legend(handles=[patch], loc="lower center")
        ax[step, feature].grid()
        # ax[step, feature].xaxis.grid()
        # counter += 1
# fig.tight_layout()
# plt.savefig('system_over_time.pdf')
plt.show()
