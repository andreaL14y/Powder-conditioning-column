from input_parameters import*

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


def compute_heat_transfer_coefficient(molar_concentration_moisture):
    reynolds_number = compute_mass_transfer_coefficient_vector(molar_concentration_moisture)[2]
    j_m = 0.61 * reynolds_number ** -0.41
    h_GP = (j_m * gas_density * gas_heat_capacity * superficial_velocity) / (
            (gas_heat_capacity * gas_viscosity / conductivity_gas) ** (2 / 3))
    return h_GP


def compute_moisture_particle_from_RH(relative_humidity):
    moisture_particle = relative_humidity/alpha_parameter_cryst
    return moisture_particle


def compute_moisture_particle_from_RH_cryst(relative_humidity):                   # TODO: added this instead of above?!
    moisture_particle = (-np.log(1-relative_humidity)/alpha_parameter_cryst)**(1/N_cryst)
    return moisture_particle


def compute_moisture_particle_from_RH_am(relative_humidity):                   # TODO: added this instead of above?!
    moisture_particle = (-np.log(1-relative_humidity)/alpha_parameter_am)**(1/N_am)
    return moisture_particle

######################################### RECURRENT ####################################################################
def compute_p_saturated_vector(temp_kelvin_vector):
    temp_celsius = temp_kelvin_vector - kelvin
    p_saturated = 10 ** (A - B / (temp_celsius + C))
    p_saturated_pascal_vector = p_saturated * 133.322                                           # Torr to Pascal
    return p_saturated_pascal_vector


def compute_molar_concentration_vector(relative_humidity_vector, pressure_saturated_vector, temp_vector):
    molar_concentration = relative_humidity_vector * pressure_saturated_vector / (r_gas_constant * temp_vector)
    return molar_concentration


def compute_partial_pressure_moisture_vector(molar_concentration_vector, temperature_vector):
    partial_pressure_moisture = molar_concentration_vector * r_gas_constant * temperature_vector
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


################################## INITIAL CONDITIONS ##################################################################
gas_velocity = compute_velocity(volumetric_flow_rate_liters_per_minute)
specific_surface_area = compute_specific_surface_area()                                 # m2/kg, SSA
flow_rate = volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute)  # m3/s
superficial_velocity = flow_rate/(np.pi*(column_diameter/2)**2)                         # superficial velocity U in m/s

pressure_saturated_initial = compute_p_saturated_vector(temp_initial)
partial_pressure_moisture_initial = pressure_saturated_initial * relative_humidity_gas_inlet
molar_concentration_moisture_initial = compute_molar_concentration_vector(
    relative_humidity_gas_inlet, pressure_saturated_initial, temp_initial)

moisture_gas_initial_bed            = compute_Y_from_RH(relative_humidity_bed_initial, pressure_saturated_initial)
moisture_gas_initial_in             = compute_Y_from_RH(relative_humidity_gas_inlet, pressure_saturated_initial)
moisture_gas_end                    = compute_Y_from_RH(relative_humidity_gas_end, pressure_saturated_initial)

# moisture_cryst_particle_initial     = compute_moisture_particle_from_RH(relative_humidity_bed_initial)
# moisture_cryst_particle_saturated   = compute_moisture_particle_from_RH(relative_humidity_gas_inlet)

moisture_cryst_particle_initial     = compute_GAB_equilibrium_moisture_cryst(relative_humidity_bed_initial)
moisture_cryst_particle_saturated   = compute_GAB_equilibrium_moisture_cryst(relative_humidity_gas_inlet)

moisture_am_particle_initial        = compute_GAB_equilibrium_moisture_am(relative_humidity_bed_initial)
moisture_am_particle_saturated      = compute_GAB_equilibrium_moisture_am(relative_humidity_gas_inlet)

k_GP_initial = compute_mass_transfer_coefficient_vector(molar_concentration_moisture_initial)[3]
constant_initial = k_GP_initial * specific_surface_area * pressure_saturated_initial / pressure_ambient

heat_transfer_coefficient = compute_heat_transfer_coefficient(molar_concentration_moisture_initial)
temp_min = min(temp_initial, temp_walls)

# EQUILIBRIUM DATA
relative_humidities_eq = np.linspace(0, 0.9, 10000)
moistures_am_eq = compute_GAB_equilibrium_moisture_am(relative_humidities_eq)
moistures_cryst_eq = compute_GAB_equilibrium_moisture_cryst(relative_humidities_eq)