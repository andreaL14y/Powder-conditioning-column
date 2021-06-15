from input_parameters import*
from v_ode_functions_2 import*
# from find_k_from_data import crystallization_parameters, crystallization_speed_curves
from glass_transition_curve import compute_glass_temp_mix, glass_temp_lactose, glass_temp_water_1, glass_temp_water_2
counter = 0
from art_functions import compute_entropy_of_activation_S, compute_enthalpy_of_activation_H, compute_crystal_growth_rate
from moisture_content_equilibrium import compute_air_equilibrium, compute_GAB_equilibrium_moisture_am, compute_GAB_equilibrium_moisture_cryst

###################################### MAIN EQUATIONS (1-4) ############################################################
def conditioning_column(moisture_matrix, t, space_step, n_space_steps, n_height_steps):
    values_per_feature = n_space_steps * n_height_steps
    global counter
    counter += 1

    ################################### CREATE ALL MATRICES ############################################################
    moisture_gas_vector, moisture_particle_cryst_vector, moisture_particle_am_vector, temp_gas_vector, \
    temp_particle_vector, amorphous_material_vector = create_matrices(
        moisture_matrix, values_per_feature, n_height_steps, n_space_steps)


    ##################################### CREATE T_G MATRIX ############################################################
    glass_transition_temp_vector = compute_glass_temp_mix(1 - moisture_particle_am_vector, glass_temp_lactose, glass_temp_water_1)
    temp_glass_temp_diff = temp_particle_vector - glass_transition_temp_vector


    ##################################### UPDATE PARAMETERS ############################################################
    pressure_saturated = compute_p_saturated_vector(temp_gas_vector)
    relative_humidity = compute_relative_humidity_from_Y_vector(moisture_gas_vector, pressure_saturated)

    equilibrium_state_cryst = compute_air_equilibrium(moisture_particle_cryst_vector, alpha_parameter_cryst, N_cryst)
    # equilibrium_state_cryst_new = compute_equilibrium_air_cryst(relative_humidity)          #TODO: CHANGED, eq moisture content powder

    equilibrium_state_am = compute_air_equilibrium(moisture_particle_am_vector, alpha_parameter_am, N_am)
    # equilibrium_state_am_new = compute_GAB_equilibrium_moisture_am(relative_humidity)

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


    ############################## UPDATE AMORPHOUS MATERIAL ###########################################################
    change_amorphous_material = compute_amorphicity(
        moisture_particle_am_vector, temp_particle_vector, amorphous_material_vector, glass_transition_temp_vector)     # will be negative

    ##################################### UPDATE MOISTURE ##############################################################
    total_equilibrium_term = ((1 - amorphous_material_vector) * (relative_humidity - equilibrium_state_cryst) +
                              amorphous_material_vector * (relative_humidity - equilibrium_state_am))

    # total_equilibrium_term_new = ((1 - amorphous_material_vector) * (equilibrium_state_cryst_new - moisture_particle_cryst_vector) +
    #                           amorphous_material_vector * (equilibrium_state_am - moisture_particle_am_vector)) #TODO: CHANGED

    #TODO: stop sorbing after start of crystallization

    # Moisture in am material is only rising, due to sorption. When T_G is reached, this material crystallizes, and is
    # no longer am. During crystallization, moisture is expelled into the gas.

    # amount crystallizing x difference saturated moisture
    change_m_particle_am_crystallization = change_amorphous_material * \
                                           (moisture_particle_am_vector - moisture_particle_cryst_vector)   # negative, losing moisture
    if np.where(change_m_particle_am_crystallization > 0)[0] != []:        #np.where(arr < 6)[0]
        print('I should be negative...')
        print(change_m_particle_am_crystallization)

    change_m_particle_cryst = constant * (relative_humidity - equilibrium_state_cryst)  # sorbing from gas
    # change_m_particle_cryst_new = constant * (equilibrium_state_cryst_new - moisture_particle_cryst_vector)  # TODO: CHANGED
    # if counter & 10000 == 0:
    #     print('\nOriginal', change_m_particle_cryst[0, 0])
    #     print('New', change_m_particle_cryst_new[0, 0])
    #     print('Eq state or', equilibrium_state_cryst[0, 0])
    #     print('Eq state new', equilibrium_state_cryst_new[0, 0])
    #     print('Current m', moisture_particle_cryst_vector[0, 0])

    change_m_particle_am_sorption = constant * (relative_humidity - equilibrium_state_am)  # sorbing from gas
    # change_m_particle_am_sorption = constant * (equilibrium_state_am - moisture_particle_am_vector)  # sorbing from gas
    change_m_particle_am = change_m_particle_am_sorption

    change_m_diffusion_gas = moisture_diffusivity * gas_density * (1 - porosity_powder) * laplacian_moisture_gas
    change_m_absorption_gas = - constant * particle_density * porosity_powder * total_equilibrium_term
    change_m_gas = (change_m_diffusion_gas + change_m_absorption_gas + change_m_particle_am_crystallization) / \
                   (gas_density * (1 - porosity_powder)) - gas_velocity * gradient_moisture_gas     # Added crystallization


    ##################################### UPDATE TEMP ##################################################################
    conduction_gas          = conductivity_gas * (1 - porosity_powder) * laplacian_temp_gas
    heat_of_sorption_gas    = particle_density * porosity_powder * constant * moisture_vapor_heat_capacity * \
                              (temp_gas_vector - temp_particle_vector) * total_equilibrium_term

    heat_transfer_gas       = -heat_transfer_coefficient * particle_density * porosity_powder * \
                              specific_surface_area * (temp_gas_vector - temp_particle_vector)

    change_temp_gas         = (conduction_gas + heat_of_sorption_gas + heat_transfer_gas) / \
                              (gas_density * (1 - porosity_powder) * gas_heat_capacity) - gas_velocity * gradient_temp_gas

    conduction_particle     = conductivity_particle * laplacian_temp_particle / particle_density
    heat_of_sorption_particle = constant * total_equilibrium_term * heat_of_sorption
    heat_transfer_particle  = heat_transfer_coefficient * specific_surface_area * (temp_gas_vector-temp_particle_vector)

    # enthalpy_of_activation = compute_enthalpy_of_activation_H(moisture_particle_am_vector)

    # heat_of_cryst_particle  = heat_of_crystallization * (-change_amorphous_material * mass_powder/(n_space_steps * n_height_steps) )
    heat_of_cryst_particle  = heat_of_crystallization * constant * (-change_amorphous_material)     # kind of copied from heat of sorption

    change_temp_particle    = (conduction_particle + heat_of_sorption_particle + heat_transfer_particle + heat_of_cryst_particle) / \
                           particle_heat_capacity


    ##################################### PRINT STUFF ##################################################################
    tabs = 50
    if counter % 10000 == 0:
        print("Counter:".ljust(tabs), counter)
        print("T - Tg:".ljust(tabs), "{:.2f}".format(temp_glass_temp_diff[2, 0]))
        print(f"Moisture am:".ljust(tabs), "{:.4f}".format(moisture_particle_am_vector[2, 0]))
        print(f"Moisture eq am:".ljust(tabs), "{:.4f}".format(equilibrium_state_am[2, 0]))
        print("RH:".ljust(tabs), "{:.2}".format(relative_humidity[2, 0]))
        print("Amorphous material change: ".ljust(tabs), "{:.3e}".format(change_amorphous_material[2, 0]))
        print("Moisture change:".ljust(tabs), "{:.3e}".format(change_m_particle_am[2, 0]) )
        print("Time computed:".ljust(tabs), "{:.3f}".format(t/(3600)), "h \n")

        # print(heat_of_cryst_particle[2, 0])
        # print(heat_of_sorption_particle[2, 0])
    # if heat_of_cryst_particle[2, 0] > 0.119:
    #     print('\n', heat_of_cryst_particle[2, 0])
    #     print(heat_of_sorption_particle[2, 0])

    ##################################### CONCATENATE ##################################################################
    moisture_matrix_updated = np.concatenate([change_m_gas.flatten(), change_m_particle_cryst.flatten(),
                                              change_m_particle_am.flatten(), change_temp_gas.flatten(),
                                              change_temp_particle.flatten(), change_amorphous_material.flatten()])
    return moisture_matrix_updated


######################################### RECURRENT ####################################################################
def compute_amorphicity(
        moisture_particle_am_vector, temp_particle_vector, amorphous_material_vector, glass_transition_temp_vector):

    # Only change if above Tg
    k_vector = compute_k_vector(temp_particle_vector, moisture_particle_am_vector, glass_transition_temp_vector)
    change_amorphous_material = - k_vector * amorphous_material_vector
    return change_amorphous_material


def compute_k_vector(temperature_vector, moisture_content_vector, glass_transition_temp_vector):
    k_vector = compute_crystal_growth_rate(moisture_content_vector, temperature_vector)
    diff_temp = temperature_vector - glass_transition_temp_vector

    # k_vector *= 1/(1 + np.exp(-3 * diff_temp - 3))          # Sigmoid to ramp up, TODO: change this, nothing <0
    k_vector *= 1/(1 + np.exp(-5 * diff_temp + 3))
    k_vector = np.where(diff_temp < -5, 0, k_vector)
    return k_vector


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


def create_matrices(moisture_matrix, values_per_feature, n_height_steps, n_space_steps):
    # n_features = len(moisture_matrix)/values_per_feature
    # for f in range(n_features):
    #     matrix = moisture_matrix[f * values_per_feature:(1 + f)*values_per_feature]
    #     matrix = np.reshape(matrix, n_height_steps, n_space_steps)

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
    return moisture_gas_vector, moisture_particle_cryst_vector, moisture_particle_am_vector, temp_gas_vector, \
           temp_particle_vector, amorphous_material_vector
    
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