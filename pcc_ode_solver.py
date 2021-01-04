from scipy.integrate import odeint
from pcc_define_functions import *
from pcc_plot_functions import *
from input_class import*

def computation_of_system(in_param: input_parameter, dis_param: discretization_parameter, init_param: initial_conditions):
    ######################################### SETUP ########################################################################
    values_per_feature = dis_param.n_space_steps * dis_param.n_height_steps
    space_step_size = in_param.bed_length / dis_param.n_space_steps
    discrete_time = np.linspace(0, dis_param.max_time, dis_param.resolution)

    moisture_gas_initial_all = np.zeros((dis_param.n_height_steps, dis_param.n_space_steps)) + init_param.moisture_gas_initial_bed
    moisture_particle_initial_all = np.zeros((dis_param.n_height_steps, dis_param.n_space_steps)) + init_param.moisture_particle_initial
    temp_gas_initial = np.zeros((dis_param.n_height_steps, dis_param.n_space_steps)) + in_param.temp_initial
    temp_particle_initial = np.zeros((dis_param.n_height_steps, dis_param.n_space_steps)) + in_param.temp_initial

    ########################################## COMPUTE #####################################################################
    initial_system = np.concatenate(
        [moisture_gas_initial_all.flatten(), moisture_particle_initial_all.flatten(),
        temp_gas_initial.flatten(), temp_particle_initial.flatten()])

    computed_system = odeint(
        conditioning_column, initial_system, discrete_time, args=(space_step_size, dis_param.n_space_steps, dis_param.n_height_steps, in_param, init_param))

    moisture_gas_vector = computed_system[:, 0:values_per_feature]
    moisture_particle_vector = computed_system[:, values_per_feature:(values_per_feature * 2)]
    temp_gas_vector = computed_system[:, (values_per_feature * 2):(values_per_feature * 3)]
    temp_particle_vector = computed_system[:, (values_per_feature * 3):(values_per_feature * 4)]

    moisture_gas_vector = np.reshape(moisture_gas_vector, (dis_param.resolution, dis_param.n_height_steps, dis_param.n_space_steps))
    moisture_particle_vector = np.reshape(moisture_particle_vector, (dis_param.resolution, dis_param.n_height_steps, dis_param.n_space_steps))
    temp_gas_vector = np.reshape(temp_gas_vector, (dis_param.resolution, dis_param.n_height_steps, dis_param.n_space_steps))
    temp_particle_vector = np.reshape(temp_particle_vector, (dis_param.resolution, dis_param.n_height_steps, dis_param.n_space_steps))

    ############################################ PLOT ######################################################################
    # Convert to easier-to-read units
    seconds = dis_param.max_time
    hours = seconds / 3600
    discrete_time /= 3600

    # temp_gas_vector -= kelvin
    # temp_particle_vector -= kelvin
    max_temp_gas = np.max(temp_gas_vector)
    max_temp_gas_index = np.where(temp_gas_vector == max_temp_gas)

    max_temp_particle = np.max(temp_particle_vector)
    max_moisture_gas = np.max(moisture_gas_vector)
    max_moisture_particle = np.max(moisture_particle_vector)
    max_moisture_particle_index = np.where(moisture_particle_vector == max_moisture_particle)

    print(f'Time computed is: {seconds} seconds = {int(hours)} hours = {int(hours/24)} days\n')
    print('Max humidity in gas is: {:.4f} degrees Celcius'.format(max_moisture_gas))
    print('Saturated humidity in gas at inlet is: {:.4f} degrees Celcius\n'.format(init_param.moisture_gas_initial_in))

    print('Max humidity in particles is: {:.4f} degrees Celcius'.format(max_moisture_particle))
    print('This happens at time, height, length: ', max_moisture_particle_index[0], max_moisture_particle_index[1], max_moisture_particle_index[2])
    print('Saturated humidity in particles is: {:.4f} degrees Celcius\n'.format(init_param.moisture_particle_saturated))

    print('Max temperature in gas is: {:.4f} degrees Celcius'.format(max_temp_gas))
    print('Max temperature in particles is: {:.4f} degrees Celcius\n'.format(max_temp_particle))

    plot_sections_over_time(
        moisture_gas_vector, moisture_particle_vector, temp_gas_vector, temp_particle_vector, dis_param.height_of_interest,
        dis_param.n_space_steps, discrete_time, init_param.moisture_gas_initial_bed, init_param.moisture_gas_initial_in, init_param.moisture_particle_initial,
        init_param.moisture_particle_saturated, init_param.temp_min, init_param.kelvin, hours, max_temp_gas, max_temp_particle)

    plot_heatmap(
        moisture_gas_vector, moisture_particle_vector, temp_gas_vector, temp_particle_vector, dis_param.height_of_interest,
        dis_param.n_space_steps, discrete_time, init_param.moisture_gas_initial_bed, init_param.moisture_gas_initial_in, init_param.moisture_particle_initial,
        init_param.moisture_particle_saturated, init_param.temp_min, init_param.kelvin, hours, max_temp_gas, max_temp_particle)

    slide_heat_map(
        moisture_gas_vector, moisture_particle_vector, temp_gas_vector, temp_particle_vector, init_param.temp_min, max_temp_particle,
        max_temp_gas, init_param.moisture_particle_initial, init_param.moisture_particle_saturated, init_param.moisture_gas_initial_bed, init_param.moisture_gas_initial_in, hours)