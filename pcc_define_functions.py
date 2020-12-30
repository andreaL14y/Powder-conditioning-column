from input_class import input_parameter, initial_conditions
import numpy as np

###################################### MAIN EQUATIONS (1-4) ############################################################
def conditioning_column(moisture_matrix, t, space_step, n_space_steps, n_height_steps, in_param: input_parameter, init_param: initial_conditions):
    values_per_feature = n_space_steps * n_height_steps

    moisture_gas_vector = moisture_matrix[0:values_per_feature]
    moisture_particle_vector = moisture_matrix[values_per_feature:(values_per_feature * 2)]
    temp_gas_vector = moisture_matrix[(values_per_feature * 2):(values_per_feature * 3)]
    temp_particle_vector = moisture_matrix[(values_per_feature * 3):(values_per_feature * 4)]

    moisture_gas_vector = np.reshape(moisture_gas_vector, (n_height_steps, n_space_steps))
    moisture_particle_vector = np.reshape(moisture_particle_vector, (n_height_steps, n_space_steps))
    temp_gas_vector = np.reshape(temp_gas_vector, (n_height_steps, n_space_steps))
    temp_particle_vector = np.reshape(temp_particle_vector, (n_height_steps, n_space_steps))

    ##################################### UPDATE PARAMETERS ############################################################
    equilibrium_state = compute_equilibrium_moisture_vector(moisture_particle_vector, in_param)
    pressure_saturated = compute_p_saturated_vector(temp_gas_vector, in_param, init_param)
    relative_humidity = compute_relative_humidity_from_Y_vector(moisture_gas_vector, pressure_saturated, in_param)
    molar_concentration_moisture = compute_molar_concentration_vector(
        relative_humidity, pressure_saturated, temp_gas_vector, in_param)

    mass_transfer_coefficient = compute_mass_transfer_coefficient_vector(molar_concentration_moisture, in_param, init_param)[3]
    constant = mass_transfer_coefficient * init_param.specific_surface_area * pressure_saturated / in_param.pressure_ambient

    laplacian_moisture_gas = compute_laplacian_vector(
        moisture_gas_vector, space_step, init_param.moisture_gas_initial_in, boundary_condition_L=init_param.moisture_gas_end,
        boundary_condition_wall=0, temperature=False)
    gradient_moisture_gas = compute_gradient_vector(
        moisture_gas_vector, space_step, init_param.moisture_gas_initial_in, boundary_condition_L=init_param.moisture_gas_end)

    laplacian_temp_gas = compute_laplacian_vector(
        temp_gas_vector, space_step, in_param.temp_initial, boundary_condition_L=in_param.temp_initial,
        boundary_condition_wall=in_param.temp_walls, temperature=True)
    gradient_temp_gas = compute_gradient_vector(
        temp_gas_vector, space_step, in_param.temp_initial, boundary_condition_L=in_param.temp_initial)

    laplacian_temp_particle = compute_laplacian_vector(
        temp_particle_vector, space_step, in_param.temp_initial, boundary_condition_L=in_param.temp_initial,
        boundary_condition_wall=in_param.temp_walls, temperature=True)

    ##################################### UPDATE MOISTURE ##############################################################
    change_m_diffusion_gas = in_param.moisture_diffusivity * in_param.gas_density * (1 - in_param.porosity_powder) * laplacian_moisture_gas
    change_m_absorption_gas = - constant * in_param.particle_density * in_param.porosity_powder * (relative_humidity - equilibrium_state)

    change_m_gas = (change_m_diffusion_gas + change_m_absorption_gas) / \
                   (in_param.gas_density * (1 - in_param.porosity_powder)) - init_param.gas_velocity * gradient_moisture_gas

    change_m_particle = constant * (relative_humidity - equilibrium_state)

    ##################################### UPDATE TEMP ##################################################################
    conduction_gas = in_param.conductivity_gas * (1 - in_param.porosity_powder) * laplacian_temp_gas
    heat_of_sorption_gas = in_param.particle_density * in_param.porosity_powder * constant * (relative_humidity - equilibrium_state) * \
                       in_param.moisture_vapor_heat_capacity * (temp_gas_vector - temp_particle_vector)
    heat_transfer_gas = -init_param.heat_transfer_coefficient * in_param.particle_density * in_param.porosity_powder * init_param.specific_surface_area * \
                    (temp_gas_vector - temp_particle_vector)

    change_temp_gas = (conduction_gas + heat_of_sorption_gas + heat_transfer_gas) / \
                      (in_param.gas_density * (1 - in_param.porosity_powder) * in_param.gas_heat_capacity) - init_param.gas_velocity * gradient_temp_gas

    conduction_particle = in_param.conductivity_particle * laplacian_temp_particle / in_param.particle_density
    heat_of_sorption_particle = constant * (relative_humidity - equilibrium_state) * in_param.heat_of_vaporization
    heat_transfer_particle = init_param.heat_transfer_coefficient * init_param.specific_surface_area * (temp_gas_vector-temp_particle_vector)

    change_temp_particle = (conduction_particle + heat_of_sorption_particle + heat_transfer_particle) / \
                           in_param.particle_heat_capacity
    # change_temp_particle[:,:] = 0
    # change_temp_gas[:,:] = 0
    return np.concatenate([change_m_gas.flatten(), change_m_particle.flatten(), change_temp_gas.flatten(), change_temp_particle.flatten()])

    # return np.concatenate([change_m_gas, change_m_particle, change_temp_gas, change_temp_particle])


######################################### ONE-TIME USE #################################################################
def volumetric_flow_rate_m3_per_second(in_param: input_parameter):
    volumetric_flow_rate = in_param.volumetric_flow_rate / (60 * 10 ** 3)  # from l/min -> cubic meters per s
    return volumetric_flow_rate


def compute_velocity(in_param: input_parameter):
    volumetric_flow_rate = volumetric_flow_rate_m3_per_second(in_param)
    area_column = (in_param.column_diameter / 2) ** 2 * np.pi
    fraction_gas = 1 - in_param.porosity_powder
    velocity = volumetric_flow_rate / (area_column * fraction_gas)                  # only area with gas, not powder
    return velocity


def compute_specific_surface_area(in_param: input_parameter):
    r = in_param.particle_diameter / 2
    specific_surface_area = 3 / (r * in_param.particle_density)
    return specific_surface_area


def compute_moisture_particle_from_RH(relative_humidity, in_param: input_parameter):
    moisture_particle = relative_humidity/in_param.alpha_parameter
    return moisture_particle


def compute_heat_transfer_coefficient(molar_concentration_moisture, in_param: input_parameter, init_param: initial_conditions):
    reynolds_number = compute_mass_transfer_coefficient_vector(molar_concentration_moisture, in_param, init_param)[2]
    j_m = 0.61 * reynolds_number ** -0.41
    h_GP = (j_m * in_param.gas_density * in_param.gas_heat_capacity * init_param.superficial_velocity) / (
            (in_param.gas_heat_capacity * in_param.gas_viscosity / in_param.conductivity_gas) ** (2 / 3))
    return h_GP


######################################### RECURRENT ####################################################################
def compute_equilibrium_moisture_vector(moisture_particle_vector, in_param: input_parameter):
    rows, cols = np.shape(moisture_particle_vector)
    f_of_x = np.zeros((rows, cols))
    for r in range(rows):
        for c in range(cols):
            if moisture_particle_vector[r, c] < 1/in_param.alpha_parameter:
                f_of_x[r, c] = in_param.alpha_parameter * moisture_particle_vector[r, c]
            else:
                f_of_x[r, c] = 1 - np.exp(-in_param.alpha_parameter * moisture_particle_vector[r, c] ** in_param.N_parameter)
    return f_of_x


def compute_p_saturated_vector(temp_kelvin_vector, in_param: input_parameter, init_param: initial_conditions):
    temp_celsius = temp_kelvin_vector - init_param.kelvin
    p_saturated = 10 ** (in_param.antoine_constant_A - in_param.antoine_constant_B / (temp_celsius + in_param.antoine_constant_C))
    p_saturated_pascal_vector = p_saturated * 133.322                                           # Torr to Pascal
    return p_saturated_pascal_vector


def compute_molar_concentration_vector(relative_humidity_vector, pressure_saturated_vector, temp_vector, in_param: input_parameter):
    molar_concentration = relative_humidity_vector * pressure_saturated_vector / (in_param.R_gas_constant * temp_vector)
    return molar_concentration


def compute_partial_pressure_moisture_vector(molar_concentration_vector, temperature_vector, in_param: input_parameter):
    partial_pressure_moisture = molar_concentration_vector * in_param.R_gas_constant * temperature_vector
    return partial_pressure_moisture


def compute_mass_transfer_coefficient_vector(molar_concentration_vector, in_param: input_parameter, init_param: initial_conditions):
    superficial_mass_velocity = (4 * in_param.gas_density * init_param.flow_rate) / (np.pi * in_param.column_diameter ** 2)  # G_0
    particle_surface_area = in_param.porosity_powder * in_param.particle_density * init_param.specific_surface_area          # a

    reynolds_number = superficial_mass_velocity / (particle_surface_area * in_param.gas_viscosity)       # Re
    denominator = in_param.gas_viscosity / (in_param.gas_density * in_param.moisture_diffusivity)
    j_m = 0.61 * reynolds_number ** -0.41

    k_gp_vector = j_m * (molar_concentration_vector * in_param.molar_mass_moisture) * init_param.superficial_velocity / \
                  (denominator ** (2 / 3))
    return superficial_mass_velocity, particle_surface_area, reynolds_number, k_gp_vector


def compute_Y_from_RH(relative_humidity, pressure_saturated, in_param: input_parameter):
    Y_initial = 1 / (1 + in_param.molar_mass_dry_air / in_param.molar_mass_moisture * (
            in_param.pressure_ambient / relative_humidity - pressure_saturated) / pressure_saturated)
    return Y_initial


def compute_relative_humidity_from_Y_vector(Y_current_vector, pressure_saturated_vector, in_param: input_parameter):
    denominator = Y_current_vector * pressure_saturated_vector - in_param.molar_mass_moisture / in_param.molar_mass_dry_air * \
                  pressure_saturated_vector * (Y_current_vector - 1)
    relative_humidity = Y_current_vector * in_param.pressure_ambient / denominator
    return relative_humidity


def compute_gradient_vector(vector, space_step, boundary_condition_0, boundary_condition_L):
    rows, cols = np.shape(vector)
    gradient = np.zeros((rows, cols))

    gradient[:, 0] = (vector[:, 1] - boundary_condition_0) / (2 * space_step)
    gradient[:, cols - 1] = (vector[:, cols - 1] - vector[:, cols - 2]) / (space_step)

    # for c in range(1, cols):
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


    # print(padded_vector)
    laplacian = np.zeros((rows, cols))
    for c in range(0, cols):
        laplacian[:, c] = (padded_vector[1:(rows_padded-1), c] - 2 * padded_vector[1:(rows_padded-1), c+1] +
                           padded_vector[1:(rows_padded-1), c + 2])
    # print('\n', laplacian)
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
    # print(laplacian)
    laplacian /= (space_step ** 2)
    return laplacian