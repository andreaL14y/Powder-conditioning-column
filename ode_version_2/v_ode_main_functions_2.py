from input_parameters import*
from v_ode_functions_2 import*
from find_k_from_data import crystallization_parameters, crystallization_speed_curves
from glass_transition_curve import glass_temp_mix, glass_temp_lactose, glass_temp_water_1, glass_temp_water_2
counter = 0

###################################### MAIN EQUATIONS (1-4) ############################################################
def conditioning_column(moisture_matrix, t, space_step, n_space_steps, n_height_steps):
    values_per_feature = n_space_steps * n_height_steps
    global counter
    counter += 1

    ################################### CREATE ALL MATRICES ############################################################
    moisture_gas_vector =               moisture_matrix[0:values_per_feature]
    moisture_particle_cryst_vector =    moisture_matrix[values_per_feature:(values_per_feature * 2)]
    moisture_particle_am_vector =       moisture_matrix[(values_per_feature * 2):(values_per_feature * 3)]
    temp_gas_vector =                   moisture_matrix[(values_per_feature * 3):(values_per_feature * 4)]
    temp_particle_vector =              moisture_matrix[(values_per_feature * 4):(values_per_feature * 5)]
    amorphous_material_vector =         moisture_matrix[(values_per_feature * 5):(values_per_feature * 6)]

    moisture_gas_vector =               np.reshape(moisture_gas_vector, (n_height_steps, n_space_steps))
    moisture_particle_cryst_vector =    np.reshape(moisture_particle_cryst_vector, (n_height_steps, n_space_steps))
    moisture_particle_am_vector =       np.reshape(moisture_particle_am_vector, (n_height_steps, n_space_steps))
    temp_gas_vector =                   np.reshape(temp_gas_vector, (n_height_steps, n_space_steps))
    temp_particle_vector =              np.reshape(temp_particle_vector, (n_height_steps, n_space_steps))
    amorphous_material_vector =         np.reshape(amorphous_material_vector, (n_height_steps, n_space_steps))
    amorphous_material_vector[amorphous_material_vector < 0] = 0

    ############################# CHECK GLASS TRANSITION TEMP ##########################################################
    glass_transition_temp_vector = glass_temp_mix(1 - moisture_particle_am_vector, glass_temp_lactose, glass_temp_water_1)
    temp_glass_temp_diff = temp_particle_vector - glass_transition_temp_vector

    ##################################### UPDATE PARAMETERS ############################################################
    pressure_saturated = compute_p_saturated_vector(temp_gas_vector)
    relative_humidity = compute_relative_humidity_from_Y_vector(moisture_gas_vector, pressure_saturated)

    equilibrium_state_cryst = compute_equilibrium_moisture_vector(moisture_particle_cryst_vector)
    equilibrium_state_am = compute_equilibrium_moisture_am_vector(relative_humidity)

    molar_concentration_moisture = compute_molar_concentration_vector(relative_humidity, pressure_saturated, temp_gas_vector)

    mass_transfer_coefficient = compute_mass_transfer_coefficient_vector(molar_concentration_moisture)[3]
    constant = mass_transfer_coefficient * specific_surface_area * pressure_saturated / pressure_ambient

    laplacian_moisture_gas = compute_laplacian_vector(
        moisture_gas_vector, space_step, moisture_gas_initial_in, boundary_condition_L=moisture_gas_end,
        boundary_condition_wall=0, temperature=False)
    gradient_moisture_gas = compute_gradient_vector(
        moisture_gas_vector, space_step, moisture_gas_initial_in, boundary_condition_L=moisture_gas_end)

    laplacian_temp_gas = compute_laplacian_vector(
        temp_gas_vector, space_step, temp_initial, boundary_condition_L=temp_initial,
        boundary_condition_wall=temp_walls, temperature=True)
    gradient_temp_gas = compute_gradient_vector(
        temp_gas_vector, space_step, temp_initial, boundary_condition_L=temp_initial)

    laplacian_temp_particle = compute_laplacian_vector(
        temp_particle_vector, space_step, temp_initial, boundary_condition_L=temp_initial,
        boundary_condition_wall=temp_walls, temperature=True)

    ##################################### UPDATE MOISTURE ##############################################################
    change_m_diffusion_gas = moisture_diffusivity * gas_density * (1 - porosity_powder) * laplacian_moisture_gas
    change_m_absorption_gas = - constant * particle_density * porosity_powder * (relative_humidity - equilibrium_state_cryst)

    change_m_gas = (change_m_diffusion_gas + change_m_absorption_gas) / \
                   (gas_density * (1 - porosity_powder)) - gas_velocity * gradient_moisture_gas

    change_m_particle = constant * (relative_humidity - equilibrium_state_cryst)
    change_m_particle_cryst = change_m_particle * (1 - amorphous_material_vector)

    change_m_particle_am = constant * (equilibrium_state_am - moisture_particle_am_vector)
    change_m_particle_am = change_m_particle_am * amorphous_material_vector


    ##################################### PRINT STUFF ##################################################################
    if counter % 2500 == 0:
        print("\nCounter:", counter)
        print(f"T, Tg, T - Tg in middle row, first col: {temp_particle_vector[2, 0] - kelvin:.2f}, "
              f"{glass_transition_temp_vector[2, 0] - kelvin:.2f}, {temp_glass_temp_diff[2, 0]:.2f}"
              f"\nMoisture am {moisture_particle_am_vector[2, 0]:.4f}\n"
              f"Moisture eq am {equilibrium_state_am[2, 0]:.4f}\n"
              f"RH {relative_humidity[2, 0]:.2}")

    ##################################### UPDATE TEMP ##################################################################
    conduction_gas = conductivity_gas * (1 - porosity_powder) * laplacian_temp_gas
    heat_of_sorption_gas = particle_density * porosity_powder * constant * (relative_humidity - equilibrium_state_cryst) * \
                       moisture_vapor_heat_capacity * (temp_gas_vector - temp_particle_vector)
    heat_transfer_gas = -heat_transfer_coefficient * particle_density * porosity_powder * specific_surface_area * \
                    (temp_gas_vector - temp_particle_vector)

    change_temp_gas = (conduction_gas + heat_of_sorption_gas + heat_transfer_gas) / \
                      (gas_density * (1 - porosity_powder) * gas_heat_capacity) - gas_velocity * gradient_temp_gas

    conduction_particle = conductivity_particle * laplacian_temp_particle / particle_density
    heat_of_sorption_particle = constant * (relative_humidity - equilibrium_state_cryst) * heat_of_vaporization
    heat_transfer_particle = heat_transfer_coefficient * specific_surface_area * (temp_gas_vector-temp_particle_vector)

    change_temp_particle = (conduction_particle + heat_of_sorption_particle + heat_transfer_particle) / \
                           particle_heat_capacity

    ############################## UPDATE AMORPHOUS MATERIAL ###########################################################
    change_amorphous_material = compute_amorphicity(
        moisture_particle_am_vector, temp_particle_vector, amorphous_material_vector, glass_transition_temp_vector)

    ##################################### CONCATENATE ##################################################################
    moisture_matrix_updated = np.concatenate([change_m_gas.flatten(), change_m_particle_cryst.flatten(),
                                              change_m_particle_am.flatten(), change_temp_gas.flatten(),
                                              change_temp_particle.flatten(), change_amorphous_material.flatten()])
    return moisture_matrix_updated


######################################### RECURRENT ####################################################################
def compute_amorphicity(
        moisture_particle_am_vector, temp_particle_vector, amorphous_material_vector, glass_transition_temp_vector):

    # Only change if above Tg
    k_parameter = compute_k_vector(temp_particle_vector, moisture_particle_am_vector, glass_transition_temp_vector)[:, :, 1]
    k_parameter[glass_transition_temp_vector > temp_particle_vector] = 0
    change_amorphous_material = - k_parameter * amorphous_material_vector

    return change_amorphous_material


def compute_k_vector(temperature_vector, humidity_vector, glass_transition_temp_vector):
    height_steps, length_steps = temperature_vector.shape
    parameters_vector = np.zeros((height_steps, length_steps, 3))
    for h in range(height_steps):
        for l in range(length_steps):
            temp = temperature_vector[h, l]
            rh = humidity_vector[h, l]
            if temp < glass_transition_temp_vector[h, l]:
                starting_point, k_parameter, rest = 0, 0, 0
            elif temp < 40 + kelvin:
                if rh < 55:
                    starting_point, k_parameter, rest = crystallization_parameters[0]
                else:
                    starting_point, k_parameter, rest = crystallization_parameters[1]
            else:
                if rh < 55:
                    starting_point, k_parameter, rest = crystallization_parameters[2]
                else:
                    starting_point, k_parameter, rest = crystallization_parameters[3]
            k_parameter /= 60 * 60 #* 7 * 24                                     # Unit transformed from 1/h to 1/s

            parameters_vector[h, l, :] = starting_point, k_parameter, rest
    return parameters_vector


def compute_gradient_vector(vector, space_step, boundary_condition_0, boundary_condition_L):
    rows, cols = np.shape(vector)
    gradient = np.zeros((rows, cols))

    gradient[:, 0] = (vector[:, 1] - boundary_condition_0) / (2 * space_step)
    gradient[:, cols - 1] = (vector[:, cols - 1] - vector[:, cols - 2]) / (space_step)

    for c in range(1, cols-1):
        gradient[:, c] = (vector[:, c+1] - vector[:, c-1]) / (2*space_step)
    return gradient


def compute_laplacian_vector(vector, space_step, boundary_condition_0, boundary_condition_L, boundary_condition_wall,
                             temperature=True):
    rows, cols = np.shape(vector)
    rows_padded = rows + 2
    cols_padded = cols + 2
    # Create the matrix with zeros all around it.
    padded_vector = np.zeros((rows_padded, cols_padded))
    padded_vector[1:rows+1, 1:cols+1] = vector

    # Fill with boundary conditions at start and end of tube.
    padded_vector[0:rows+2, 0] = boundary_condition_0
    # padded_vector[0:rows+2, cols+1] = boundary_condition_L
    padded_vector[1:rows+1, cols+1] = vector[0:rows, cols-1]

    # Fill with boundary conditions at walls for temperature
    if temperature:
        padded_vector[0, 0:cols_padded] = boundary_condition_wall
        padded_vector[rows_padded-1, 0:cols_padded] = boundary_condition_wall

    laplacian = np.zeros((rows, cols))
    for c in range(0, cols):
        laplacian[:, c] = (padded_vector[1:(rows_padded-1), c] - 2 * padded_vector[1:(rows_padded-1), c+1] +
                           padded_vector[1:(rows_padded-1), c + 2])
    if temperature:
        for r in range(0, rows):
            laplacian[r, :] += (padded_vector[r, 1:(cols_padded-1)] - 2 * padded_vector[r+1, 1:(cols_padded-1)] +
                                padded_vector[r + 2, 1:(cols_padded-1)])
    else:
        for r in range(1, rows-1):
            laplacian[r, :] += (padded_vector[r, 1:(cols_padded-1)] - 2 * padded_vector[r+1, 1:(cols_padded-1)] +
                                padded_vector[r + 2, 1:(cols_padded-1)])
        laplacian[0, :] += ( - padded_vector[1, 1:(cols_padded - 1)] + padded_vector[2, 1:(cols_padded - 1)])
        laplacian[rows-1, :] += (padded_vector[rows_padded-3, 1:(cols_padded - 1)] - padded_vector[rows_padded-2, 1:(cols_padded - 1)])
    laplacian /= (space_step ** 2)
    return laplacian

############################################## TESTING #################################################################
# r, c = 2, 5
# temp_vector = np.zeros((r, c)) + temp_initial
# rh_vector = np.zeros((r, c)) + relative_humidity_bed_initial
# am_vector = np.zeros((r, c)) + 1
#
# # parameters_vector = compute_k_vector(temp_vector, rh_vector)
# # print(parameters_vector[0,0,:])
#
# test = compute_amorphicity(temp_vector, rh_vector, am_vector)
# print(test)