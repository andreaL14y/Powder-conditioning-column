# IMPORTS
import numpy as np
import math


# constant = k_GP * surface_area * pressure/pressure
###################################### MAIN EQUATIONS (1-4) ############################################################
def compute_moisture_particle_vector(moisture_particle_vector, alpha, N, relative_humidity_vector, dt, constant_vector,
                                     verbose=False):
    change_moisture_x = constant_vector * (
                relative_humidity_vector - compute_equilibrium_moisture_vector(alpha, moisture_particle_vector, N))  # dX/dt
    # print(change_moisture_x[time])
    moisture_difference_x = change_moisture_x * dt  # dX
    moisture_particle_current = moisture_particle_vector + moisture_difference_x
    return moisture_particle_current


def compute_moisture_gas_vector(
        moisture_particle_vector, moisture_gas_vector, alpha, N, relative_humidity_vector, dt, constant_vector, velocity,
                         diffusivity, gradient_vector, laplacian_vector, density_gas, density_particle, porosity,
                         verbose = False):

    change_moisture_diffusion = diffusivity * density_gas * (1 - porosity) * laplacian_vector
    change_moisture_absorption = - constant_vector * density_particle * porosity * \
                                 (relative_humidity_vector -
                                  compute_equilibrium_moisture_vector(alpha, moisture_particle_vector, N))
    # change_moisture_absorption = 0
    change_moisture = (change_moisture_diffusion + change_moisture_absorption) / (density_gas * (1 - porosity)) - \
                      velocity * gradient_vector
    moisture_gas_current = moisture_gas_vector + change_moisture * dt

    return moisture_gas_current


def compute_temperature_particle_vector(
        temp_particle_vector, constant_vector, dt, conductivity, laplacian_vector, density, alpha, moisture_vector,
        relative_humidity_vector, N, heat_of_vaporization, heat_transfer_coefficient, specific_surface, temp_gas_vector,
        heat_capacity, verbose = False):

    conduction = conductivity * laplacian_vector / density
    heat_of_sorption = constant_vector * (
            relative_humidity_vector - compute_equilibrium_moisture_vector(alpha, moisture_vector, N)) * heat_of_vaporization
    heat_transfer = heat_transfer_coefficient * specific_surface * (temp_gas_vector - temp_particle_vector)

    change_temperature = (conduction + heat_of_sorption + heat_transfer) / heat_capacity
    temp_particle_vector += change_temperature * dt

    return temp_particle_vector

def compute_temperature_gas_vector(
        temp_particle_vector, constant_vector, dt, conductivity_gas, laplacian_vector, density_gas, alpha, moisture_vector, N, heat_capacity_vapor,
        relative_humidity_vector, heat_transfer_coefficient, specific_surface, temp_gas_vector, heat_capacity_wet_gas, velocity,
        temp_gradient_vector, porosity, density_particle, verbose=False):

    conduction = conductivity_gas * (1 - porosity) * laplacian_vector

    heat_of_sorption = density_particle * porosity * constant_vector * \
                       (relative_humidity_vector - compute_equilibrium_moisture_vector(alpha, moisture_vector, N)) * \
                       heat_capacity_vapor * (temp_gas_vector - temp_particle_vector)

    heat_transfer = -heat_transfer_coefficient * density_particle * porosity * specific_surface * \
                    (temp_gas_vector - temp_particle_vector)

    change_temperature = (conduction + heat_of_sorption + heat_transfer) / \
                         (density_gas * (1 - porosity) * heat_capacity_wet_gas) - velocity * temp_gradient_vector

    temp_gas_vector += change_temperature * dt
    return temp_gas_vector


######################################### ONE-TIME USE #################################################################
def volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute):
    volumetric_flow_rate = volumetric_flow_rate_liters_per_minute / (60 * 10 ** 3)  # from l/min -> cubic meters per s
    return volumetric_flow_rate


def compute_velocity(volumetric_flow_rate_liters_per_minute, length, diameter, volume_fraction_powder):
    volumetric_flow_rate = volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute)
    area_column = (diameter / 2) ** 2 * np.pi
    fraction_gas = 1 - volume_fraction_powder
    # volume_gas_in_tube = area_column * length * fraction_gas
    velocity = volumetric_flow_rate / (area_column * fraction_gas)  # only area with gas, not powder
    return velocity


def spec_surface_area(particle_diameter, particle_density):  # CORRECT :)
    # S_P in m^2/kg - assuming spherical particles S_P=surface_area/volume*density
    r = particle_diameter / 2
    SSA = 3 / (r * particle_density)
    return SSA


def compute_initial_moisture_particle(alpha, N, relative_humidity):
    # moisture_particle = (-np.log(-(relative_humidity - 1)) / alpha) ** (1 / N)
    # print('MP:',moisture_particle)
    moisture_particle = relative_humidity/alpha
    return moisture_particle


def compute_heat_transfer_coefficient(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
        particle_diameter, Mw, superficial_velocity, molar_concentration, gas_heat_capacity, conductivity_gas):
    reynolds_number = compute_mass_transfer_coefficient_vector(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
        particle_diameter, Mw, superficial_velocity, molar_concentration)[2]
    j_m = 0.61 * reynolds_number ** -0.41
    h_GP = (j_m * gas_density * gas_heat_capacity * superficial_velocity) / (
            (gas_heat_capacity * gas_viscosity / conductivity_gas) ** (2 / 3))
    return h_GP


######################################### RECURRENT ####################################################################
def compute_equilibrium_moisture_vector(alpha, moisture_particle_vector, N):
    indices = len(moisture_particle_vector)
    f_of_x = np.zeros(indices)
    for i in range(indices):
        if moisture_particle_vector[i] < 1/alpha:
            f_of_x[i] = alpha * moisture_particle_vector[i]
        else:
            f_of_x[i] = 1 - np.exp(-alpha * moisture_particle_vector[i] ** N)
    return f_of_x


def compute_p_saturated_vector(A, B, temp_kelvin_vector, C):  # Double-checked and clear! :)
    temp_celsius = temp_kelvin_vector - 273.15
    p_saturated = 10 ** (A - B / (temp_celsius + C))
    p_saturated_pascal_vector = p_saturated * 133.322  # Torr to Pascal
    return p_saturated_pascal_vector


def compute_molar_concentration_vector(relative_humidity_vector, pressure_saturated_vector, R, temp_vector):
    molar_concentration = relative_humidity_vector * pressure_saturated_vector / (R * temp_vector)
    return molar_concentration


def compute_partial_pressure_moisture_vector(molar_concentration_vector, R_gas_constant, temperature_vector):  # c = molar_concentration, ideal gas law
    partial_pressure_moisture = molar_concentration_vector * R_gas_constant * temperature_vector
    return partial_pressure_moisture


def compute_mass_transfer_coefficient_vector(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
        particle_diameter, Mw, superficial_velocity, molar_concentration_vector):
    SSA = spec_surface_area(particle_diameter, particle_density)
    superficial_mass_velocity = (4 * gas_density * flow_rate) / (np.pi * column_diameter ** 2)  # G_0
    particle_surface_area = porosity_powder * particle_density * SSA  # a

    reynolds_number = superficial_mass_velocity / (particle_surface_area * gas_viscosity)  # Re
    denominator = gas_viscosity / (gas_density * moisture_diffusivity)
    j_m = 0.61 * reynolds_number ** -0.41

    k_gp_vector = j_m * (molar_concentration_vector * Mw) * superficial_velocity / (denominator ** (2 / 3))  # kg/(s*m2)
    return superficial_mass_velocity, particle_surface_area, reynolds_number, k_gp_vector


def compute_Y_from_RH(
        molar_mass_dry_air, molar_mass_moisture, pressure_ambient, relative_humidity, pressure_saturated):
    Y_initial = 1 / (1 + molar_mass_dry_air / molar_mass_moisture * (
            pressure_ambient / relative_humidity - pressure_saturated) / pressure_saturated)
    if Y_initial < 0:
        Y_initial = 0
        print('Computing Y from RH gives negative value. Y now set to 0. Check!')
    return Y_initial


def compute_relative_humidity_from_Y_vector(
        molar_mass_dry_air, molar_mass_moisture, pressure_ambient, Y_current_vector, pressure_saturated_vector):
    denominator = Y_current_vector * pressure_saturated_vector - molar_mass_moisture / molar_mass_dry_air * \
                  pressure_saturated_vector * (Y_current_vector - 1)
    relative_humidity = Y_current_vector * pressure_ambient / denominator
    # np.where(relative_humidity < 0)[0] = 0
    relative_humidity[relative_humidity < 0] = 0
    return relative_humidity


def compute_gradient_vector(vector, space_step):  # Should work for temperature and moisture, where vector
    # is an array input e.g. moisture gas and the index is the looping variable x
    length = len(vector)
    gradient = np.zeros(length)
    gradient[0] = 0
    gradient[length - 1] = 0

    for i in range(1, length-1):
        gradient[i] = (vector[i] - vector[i-1]) / space_step     # TODO: changed to minus. Right or not?
    return gradient


def compute_laplacian_vector(vector, space_step):
    length = len(vector)
    laplacian = np.zeros(length)
    laplacian[0] = 0
    laplacian[length - 1] = 0

    for i in range(1, length - 1):
        laplacian[i] = (vector[i - 1] - 2 * vector[i] + vector[i + 1]) / (space_step ** 2)
    return laplacian
