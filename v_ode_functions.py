from input_parameters import*
###################################### MAIN EQUATIONS (1-4) ############################################################
def conditioning_column(moisture_matrix, t, space_step, n_space_steps, n_height_steps):
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
    equilibrium_state = compute_equilibrium_moisture_vector(moisture_particle_vector)
    pressure_saturated = compute_p_saturated_vector(temp_gas_vector)
    relative_humidity = compute_relative_humidity_from_Y_vector(moisture_gas_vector, pressure_saturated)
    molar_concentration_moisture = compute_molar_concentration_vector(
        relative_humidity, pressure_saturated, temp_gas_vector)

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
    change_m_absorption_gas = - constant * particle_density * porosity_powder * (relative_humidity - equilibrium_state)

    change_m_gas = (change_m_diffusion_gas + change_m_absorption_gas) / \
                   (gas_density * (1 - porosity_powder)) - gas_velocity * gradient_moisture_gas
    # if t ==0:
    #     print(change_m_diffusion_gas[0,0])
    #     print(change_m_absorption_gas[0,0])
    change_m_particle = constant * (relative_humidity - equilibrium_state)

    ##################################### UPDATE TEMP ##################################################################
    conduction_gas = conductivity_gas * (1 - porosity_powder) * laplacian_temp_gas
    heat_of_sorption_gas = particle_density * porosity_powder * constant * (relative_humidity - equilibrium_state) * \
                       moisture_vapor_heat_capacity * (temp_gas_vector - temp_particle_vector)
    heat_transfer_gas = -heat_transfer_coefficient * particle_density * porosity_powder * specific_surface_area * \
                    (temp_gas_vector - temp_particle_vector)

    change_temp_gas = (conduction_gas + heat_of_sorption_gas + heat_transfer_gas) / \
                      (gas_density * (1 - porosity_powder) * gas_heat_capacity) - gas_velocity * gradient_temp_gas

    conduction_particle = conductivity_particle * laplacian_temp_particle / particle_density
    heat_of_sorption_particle = constant * (relative_humidity - equilibrium_state) * heat_of_vaporization
    heat_transfer_particle = heat_transfer_coefficient * specific_surface_area * (temp_gas_vector-temp_particle_vector)

    change_temp_particle = (conduction_particle + heat_of_sorption_particle + heat_transfer_particle) / \
                           particle_heat_capacity
    # print(t)
    # if t == 1:
    #     print(conduction_gas[0,0])
    #     print(heat_of_sorption_gas[0,0])
    #     print(heat_transfer_gas[0,0])
    return np.concatenate([change_m_gas.flatten(), change_m_particle.flatten(), change_temp_gas.flatten(), change_temp_particle.flatten()])

######################################### ONE-TIME USE #################################################################
def volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute):
    volumetric_flow_rate = volumetric_flow_rate_liters_per_minute / (60 * 10 ** 3)  # from l/min -> cubic meters per s
    return volumetric_flow_rate


def compute_velocity(volumetric_flow_rate_liters_per_minute):
    volumetric_flow_rate = volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute)
    area_column = (column_diameter / 2) ** 2 * np.pi
    fraction_gas = 1 - porosity_powder
    velocity = volumetric_flow_rate / (area_column * fraction_gas)                  # only area with gas, not powder
    return velocity


def compute_specific_surface_area():
    r = particle_diameter / 2
    specific_surface_area = 3 / (r * particle_density)
    return specific_surface_area


def compute_moisture_particle_from_RH(relative_humidity):
    moisture_particle = relative_humidity/alpha_parameter
    return moisture_particle


def compute_heat_transfer_coefficient(molar_concentration_moisture):
    reynolds_number = compute_mass_transfer_coefficient_vector(molar_concentration_moisture)[2]
    j_m = 0.61 * reynolds_number ** -0.41
    h_GP = (j_m * gas_density * gas_heat_capacity * superficial_velocity) / (
            (gas_heat_capacity * gas_viscosity / conductivity_gas) ** (2 / 3))
    return h_GP


######################################### RECURRENT ####################################################################
def compute_equilibrium_moisture_vector(moisture_particle_vector):
    rows, cols = np.shape(moisture_particle_vector)
    f_of_x = np.zeros((rows, cols))
    for r in range(rows):
        for c in range(cols):
            if moisture_particle_vector[r, c] < 1/alpha_parameter:
                f_of_x[r, c] = alpha_parameter * moisture_particle_vector[r, c]
            else:
                f_of_x[r, c] = 1 - np.exp(-alpha_parameter * moisture_particle_vector[r, c] ** N)
    return f_of_x


def compute_p_saturated_vector(temp_kelvin_vector):
    temp_celsius = temp_kelvin_vector - kelvin
    p_saturated = 10 ** (A - B / (temp_celsius + C))
    p_saturated_pascal_vector = p_saturated * 133.322                                           # Torr to Pascal
    return p_saturated_pascal_vector


def compute_molar_concentration_vector(relative_humidity_vector, pressure_saturated_vector, temp_vector):
    molar_concentration = relative_humidity_vector * pressure_saturated_vector / (R_gas_constant * temp_vector)
    return molar_concentration


def compute_partial_pressure_moisture_vector(molar_concentration_vector, temperature_vector):
    partial_pressure_moisture = molar_concentration_vector * R_gas_constant * temperature_vector
    return partial_pressure_moisture


def compute_mass_transfer_coefficient_vector(molar_concentration_vector):
    superficial_mass_velocity = (4 * gas_density * flow_rate) / (np.pi * column_diameter ** 2)  # G_0
    particle_surface_area = porosity_powder * particle_density * specific_surface_area          # a

    reynolds_number = superficial_mass_velocity / (particle_surface_area * gas_viscosity)       # Re
    denominator = gas_viscosity / (gas_density * moisture_diffusivity)
    j_m = 0.61 * reynolds_number ** -0.41

    k_gp_vector = j_m * (molar_concentration_vector * molar_mass_moisture) * superficial_velocity / \
                  (denominator ** (2 / 3))
    return superficial_mass_velocity, particle_surface_area, reynolds_number, k_gp_vector


def compute_Y_from_RH(relative_humidity, pressure_saturated):
    Y_initial = 1 / (1 + molar_mass_dry_air / molar_mass_moisture * (
            pressure_ambient / relative_humidity - pressure_saturated) / pressure_saturated)
    return Y_initial


def compute_relative_humidity_from_Y_vector(Y_current_vector, pressure_saturated_vector):
    denominator = Y_current_vector * pressure_saturated_vector - molar_mass_moisture / molar_mass_dry_air * \
                  pressure_saturated_vector * (Y_current_vector - 1)
    relative_humidity = Y_current_vector * pressure_ambient / denominator
    return relative_humidity


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

################################## INITIAL CONDITIONS ##################################################################
gas_velocity = compute_velocity(volumetric_flow_rate_liters_per_minute)
specific_surface_area = compute_specific_surface_area()                                 # m2/kg, SSA
flow_rate = volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute)  # m3/s
superficial_velocity = flow_rate/(np.pi*(column_diameter/2)**2)                         # superficial velocity U in m/s

pressure_saturated_initial = compute_p_saturated_vector(temp_initial)
partial_pressure_moisture_initial = pressure_saturated_initial * relative_humidity_gas_initial
molar_concentration_moisture_initial = compute_molar_concentration_vector(
    relative_humidity_gas_initial, pressure_saturated_initial, temp_initial)

moisture_gas_initial_bed = compute_Y_from_RH(relative_humidity_bed_initial, pressure_saturated_initial)
moisture_gas_initial_in = compute_Y_from_RH(relative_humidity_gas_initial, pressure_saturated_initial)
moisture_gas_end = compute_Y_from_RH(relative_humidity_gas_end, pressure_saturated_initial)
moisture_particle_initial = compute_moisture_particle_from_RH(relative_humidity_bed_initial)
moisture_particle_saturated = compute_moisture_particle_from_RH(relative_humidity_gas_initial)

k_GP_initial = compute_mass_transfer_coefficient_vector(molar_concentration_moisture_initial)[3]
constant_initial = k_GP_initial * specific_surface_area * pressure_saturated_initial / pressure_ambient

heat_transfer_coefficient = compute_heat_transfer_coefficient(molar_concentration_moisture_initial)
temp_min = min(temp_initial, temp_walls)