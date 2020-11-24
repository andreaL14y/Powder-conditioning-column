# IMPORTS
import numpy as np
import math


# constant = k_GP * surface_area * pressure/pressure
###################################### MAIN EQUATIONS (1-4) ############################################################
def compute_moisture_particle(moisture_particle, alpha, N, relative_humidity, dt, constant):
    change_moisture = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle, N))
    moisture_particle_current = moisture_particle + change_moisture * dt
    return moisture_particle_current


def compute_temperature_particle(
        temp_particle, constant, dt, conductivity, laplacian, density, alpha, moisture, relative_humidity, N,
        heat_of_vaporization, heat_transfer_coefficient, specific_surface, temp_gas, heat_capacity, x): #TODO: Xis just a helper variable?
    conduction = conductivity * laplacian / density
    heat_of_sorption = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture, N)) * \
                       heat_of_vaporization
    heat_transfer = heat_transfer_coefficient * specific_surface * (temp_gas - temp_particle)

    # if x == 0:
    #     print('Temp particle: ', temp_particle)
    #     print('Heat of sorption: ', heat_of_sorption)
    #     print('Heat transfer particle: ', heat_transfer, '\n')

    # change_temperature = (conduction + heat_of_sorption + heat_transfer)/heat_capacity

    change_temperature = (heat_of_sorption + heat_transfer) / heat_capacity
    temp_particle += change_temperature * dt
    return temp_particle


def compute_moisture_gas(moisture_particle, moisture_gas, alpha, N, relative_humidity, dt, constant, velocity, diffusivity,
                         gradient_moisture, laplacian, density_gas, density_particle, porosity, x):

    change_moisture_diffusion = diffusivity * density_gas * (1 - porosity) * laplacian
    change_moisture_absorption = - constant * density_particle * porosity * \
                                 (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle, N))
    # change_moisture = (change_moisture_diffusion + change_moisture_absorption) / (density_gas * (1 - porosity)) - \
    #                   velocity * gradient_moisture

    change_moisture = (change_moisture_absorption) / (density_gas * (1 - porosity)) - \
                      velocity * gradient_moisture

    moisture_gas_current = moisture_gas + change_moisture * dt

    # if x == 0:
    #     print('Moisture gas: ', moisture_gas_current)
        # print('Moisture chchange_moisture)
    return moisture_gas_current


def compute_temperature_gas(
        temp_particle, constant, dt, conductivity_gas, laplacian, density_gas, alpha, moisture, N, heat_capacity_vapor,
        relative_humidity, heat_transfer_coefficient, specific_surface, temp_gas, heat_capacity_wet_gas, velocity,
        temp_gradient, porosity, density_particle, x):
    conduction = conductivity_gas * (1 - porosity) * laplacian
    heat_of_sorption = density_particle * porosity * constant * \
                       (relative_humidity - compute_equilibrium_moisture(alpha, moisture, N)) * heat_capacity_vapor * \
                       (temp_gas - temp_particle)
    # print('Heat of sorption: ', heat_of_sorption)

    heat_transfer = -heat_transfer_coefficient * density_particle * porosity * specific_surface * \
                    (temp_gas - temp_particle)
    if x == 0:
        print('Heat transfer gas: ', heat_transfer/(porosity * density_particle))
        print('Temp diff: ', temp_gas - temp_particle)

    # change_temperature = (heat_of_sorption + heat_transfer) / (density_gas * (1 - porosity) * heat_capacity_wet_gas) - \
    #                      velocity * temp_gradient
    change_temperature = (heat_transfer) / (density_gas * (1 - porosity) * heat_capacity_wet_gas)
    temp_gas += change_temperature * dt
    return temp_gas


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


def spec_surface_area(particle_diameter, particle_density):
    # S_P in m^2/kg - assuming spherical particles S_P=surface_area/volume*density
    r = particle_diameter / 2
    SSA = (4 * np.pi * r ** 2) / (4 / 3 * np.pi * r ** 3 * particle_density)
    return SSA


def compute_initial_moisture_particle(alpha, N, relative_humidity):
    moisture_particle = (-np.log(-(relative_humidity - 1)) / alpha) ** (1 / N)
    return moisture_particle


######################################### RECURRENT ####################################################################
def compute_equilibrium_moisture(alpha, moisture_particle, N):
    f_of_x = 1 - np.exp(-alpha * moisture_particle ** N)
    return f_of_x


def compute_mass_transfer_coefficient(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
        particle_diameter, Mw, superficial_velocity, molar_concentration):
    SSA = spec_surface_area(particle_diameter, particle_density)
    superficial_mass_velocity = (4 * gas_density * flow_rate) / (np.pi * column_diameter ** 2)  # G_0
    particle_surface_area = porosity_powder * particle_density * SSA  # a
    reynolds_number = superficial_mass_velocity / (particle_surface_area * gas_viscosity)  # Re
    denominator = gas_viscosity / (gas_density * moisture_diffusivity)
    j_m = 0.61 * reynolds_number ** -0.41

    k_gp = j_m * (molar_concentration * Mw) * superficial_velocity / (denominator ** (2 / 3))  # kg/(s*m2)
    return superficial_mass_velocity, particle_surface_area, reynolds_number, k_gp


def compute_heat_transfer_coefficient(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
        particle_diameter, Mw, superficial_velocity, molar_concentration, gas_heat_capacity):

    superficial_mass_velocity, particle_surface_area, reynolds_number, k_gp = compute_mass_transfer_coefficient(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
        particle_diameter, Mw, superficial_velocity, molar_concentration)
    j_m = 0.61 * reynolds_number ** -0.41
    h_GP = (j_m * gas_density * gas_heat_capacity * superficial_velocity) / (
                (gas_heat_capacity * gas_viscosity / k_gp) ** (2 / 3))
    return h_GP


def compute_p_saturated(A, B, temp_kelvin, C):  # Double-checked and clear! :)
    temp_celsius = temp_kelvin - 273.15
    p_saturated = 10 ** (A - B / (temp_celsius + C))
    p_saturated_pascal = p_saturated * 133.322  # Torr to Pascal
    return p_saturated_pascal


def compute_partial_pressure_moisture(molar_concentration, R_gas_constant, temperature):  # c = molar_concentration
    partial_pressure_moisture = molar_concentration * R_gas_constant * temperature
    return partial_pressure_moisture


def compute_relative_humidity(partial_pressure_moisture, pressure_saturated):
    relative_humidity = partial_pressure_moisture / pressure_saturated
    # print(partial_pressure_moisture, pressure_saturated)
    return relative_humidity

# def compute_partial_pressure_moisture(molar_concentration, R_gas_constant, temp): # c = molar_concentration
#     partial_pressure_moisture = R_gas_constant * temp * molar_concentration
#     return partial_pressure_moisture


# def compute_relative_humidity(partial_p_moisture, p_saturated):
#     relative_humidity = partial_p_moisture/(p_saturated)
#     if relative_humidity > 1:
#         print('Relative humidity larger than 1, error!')
#     return relative_humidity

# def compute_relative_humidity():
#     relative_humidity = 0.5
#     if relative_humidity > 1:
#         print('Relative humidity larger than 1, error!')
#     return relative_humidity


# def compute_molar_concentration(relative_humidity, pressure_saturated, R, temp):
#    print(relative_humidity)
#    molar_concentration = 0.01 * relative_humidity * pressure_saturated / (R * temp)
#    return molar_concentration
