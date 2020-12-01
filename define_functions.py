# IMPORTS
import numpy as np
import math


# constant = k_GP * surface_area * pressure/pressure
###################################### MAIN EQUATIONS (1-4) ############################################################
def compute_moisture_particle(moisture_particle, alpha, N, relative_humidity, dt, constant, moisture_gas,
                              density_particle, porosity, density_gas):  # moisture_gas added

    change_moisture_x = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle, N))      # dX/dt
    moisture_difference_x = change_moisture_x * dt                                                                      # dX
    moisture_particle_current = moisture_particle + moisture_difference_x

    # moisture_difference_y = change_moisture_y * dt
    moisture_difference_y_absorption = - change_moisture_x * density_particle * porosity/(density_gas * (1-porosity))      # dY/dt
    # moisture_difference_y_diffusion = diffusivity * density_gas * (1 - porosity) * laplacian
    # moisture_difference_y = moisture_difference_y_absorption + moisture_difference_y_diffusion
    moisture_gas_current = moisture_gas + moisture_difference_y_absorption * dt
    return moisture_particle_current, moisture_gas_current, moisture_difference_y_absorption*dt   # return only particle


# def compute_moisture_particle(moisture_particle, alpha, N, relative_humidity, dt, constant):  # old version
#     change_moisture_x = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle, N))      # dX/dt
#     moisture_difference_x = change_moisture_x * dt                                                                      # dX
#     moisture_particle_current = moisture_particle + moisture_difference_x
#
#     return moisture_particle_current


def compute_temperature_particle(
        temp_particle, constant, dt, conductivity, laplacian, density, alpha, moisture, relative_humidity, N,
        heat_of_vaporization, heat_transfer_coefficient, specific_surface, temp_gas, heat_capacity, x):
    conduction = conductivity * laplacian / density
    heat_of_sorption = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture, N)) * \
                       heat_of_vaporization
    heat_transfer = heat_transfer_coefficient * specific_surface * (temp_gas - temp_particle)

    # change_temperature = (conduction + heat_of_sorption + heat_transfer)/heat_capacity
    # if x == 0:
    #     print('\n\nTemp particle: ', temp_particle)
    change_temperature = (heat_of_sorption + heat_transfer) / heat_capacity
    temp_particle += change_temperature * dt
    # if x == 0:
    #     print('Temp particle: ', temp_particle)
    #     print('Heat of sorption: ', heat_of_sorption)
    #     print('Heat transfer particle: ', heat_transfer)
    #     print('Temp change particle: ', change_temperature)
    return temp_particle


def compute_relative_humidity(moisture_particle, relative_humidity, alpha, N, dt, constant, velocity, diffusivity,
                         gradient_moisture, laplacian, density_gas, density_particle, porosity, x):

    change_moisture_diffusion = diffusivity * density_gas * (1 - porosity) * laplacian
    change_moisture_absorption = - constant * density_particle * porosity * \
                                 (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle, N))
    if x == 0:
        print('Moisture diffusion: ', change_moisture_diffusion)
        print('Moisture absorption: ', change_moisture_absorption)
        # print('u * nablaY: ', velocity * gradient_moisture)
        # print('Equilibrium: ', compute_equilibrium_moisture(alpha, moisture_particle, N))
        # print('Moisture particle in eq. fn: ', moisture_particle)

    change_moisture = (change_moisture_diffusion + change_moisture_absorption) / (density_gas * (1 - porosity)) - \
                      velocity * gradient_moisture

    relative_humidity_current = relative_humidity + change_moisture * dt
    # if x == 0:
    #     print('Moisture gas: ', moisture_gas_current)
        # print('Moisture chchange_moisture)
    return relative_humidity_current


# def compute_moisture_gas(moisture_particle, moisture_gas, alpha, N, relative_humidity, dt, constant, velocity, diffusivity,
#                          gradient_moisture, laplacian, density_gas, density_particle, porosity, x):
#
#     change_moisture_diffusion = diffusivity * density_gas * (1 - porosity) * laplacian
#     change_moisture_absorption = - constant * density_particle * porosity * \
#                                  (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle, N))
#     if x == 0:
#         print('Moisture diffusion: ', change_moisture_diffusion)
#         print('Moisture absorption: ', change_moisture_absorption)
#         # print('u * nablaY: ', velocity * gradient_moisture)
#         # print('Equilibrium: ', compute_equilibrium_moisture(alpha, moisture_particle, N))
#         # print('Moisture particle in eq. fn: ', moisture_particle)
#     # change_moisture = (change_moisture_diffusion + change_moisture_absorption) / (density_gas * (1 - porosity)) - \
#     #                   velocity * gradient_moisture
#
#     change_moisture = (change_moisture_diffusion + change_moisture_absorption) / (density_gas * (1 - porosity)) - \
#                       velocity * gradient_moisture
#
#     moisture_gas_current = moisture_gas + change_moisture * dt
#     # if x == 0:
#     #     print('Moisture gas: ', moisture_gas_current)
#         # print('Moisture chchange_moisture)
#     return moisture_gas_current


def compute_temperature_gas(
        temp_particle, constant, dt, conductivity_gas, laplacian, density_gas, alpha, moisture, N, heat_capacity_vapor,
        relative_humidity, heat_transfer_coefficient, specific_surface, temp_gas, heat_capacity_wet_gas, velocity,
        temp_gradient, porosity, density_particle, x):

    conduction = conductivity_gas * (1 - porosity) * laplacian
    heat_of_sorption = density_particle * porosity * constant * \
                       (relative_humidity - compute_equilibrium_moisture(alpha, moisture, N)) * heat_capacity_vapor * \
                       (temp_gas - temp_particle)
    heat_transfer = -heat_transfer_coefficient * density_particle * porosity * specific_surface * \
                    (temp_gas - temp_particle)
    # if x == 0:
    #     print('RH for gas temp: ', relative_humidity)
    #     print('\nTemp gas: ', temp_gas)
    #     print('Heat transfer gas: ', heat_transfer)
    #     print('Temp diff: ', temp_gas - temp_particle)
        # print('h_GP: ', heat_transfer_coefficient)

    # change_temperature = (conduction + heat_of_sorption + heat_transfer) / (density_gas * (1 - porosity) * heat_capacity_wet_gas) - \
    #                      velocity * temp_gradient
    change_temperature = (heat_transfer + heat_of_sorption) / (density_gas * (1 - porosity) * heat_capacity_wet_gas)    # TODO: simplified
    temp_gas += change_temperature * dt
    # if x == 0:
    #     print('Temp change gas/s: ', change_temperature)
    #     print('Temp gas: ', temp_gas)
    return temp_gas

# def compute_relative_humidity(
#         moisture_particle, alpha, N, moisture_difference, dt, constant, velocity, diffusivity,
#                          gradient_moisture, laplacian, density_gas, density_particle, porosity, x):
#
#     change_moisture_diffusion = diffusivity * density_gas * (1 - porosity) * laplacian
#
#     LH = change_moisture_diffusion - density_gas * (1-porosity) * (moisture_difference/dt + velocity * gradient_moisture)
#     # k_GP * surface_area * pressure / pressure
#     relative_humidity = LH / (constant * density_particle * porosity) + compute_equilibrium_moisture(alpha, moisture_particle, N)
#     if x == 0:
#         print('Eq. equation: ', compute_equilibrium_moisture(alpha, moisture_particle, N))
#     #     print('RH: ', relative_humidity)
#         # print('Moisture absorption: ', change_moisture_absorption)
#
#     return relative_humidity


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


def spec_surface_area(particle_diameter, particle_density):                         # CORRECT :)
    # S_P in m^2/kg - assuming spherical particles S_P=surface_area/volume*density
    r = particle_diameter / 2
    SSA = 3 / (r * particle_density)
    return SSA


def compute_initial_moisture_particle(alpha, N, relative_humidity):
    moisture_particle = (-np.log(-(relative_humidity - 1)) / alpha) ** (1 / N)
    return moisture_particle


######################################### RECURRENT ####################################################################
def compute_equilibrium_moisture(alpha, moisture_particle, N):
    if moisture_particle < alpha: 
        f_of_x = alpha * moisture_particle
    else: 
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
        particle_diameter, Mw, superficial_velocity, molar_concentration, gas_heat_capacity, conductivity_gas):

    reynolds_number = compute_mass_transfer_coefficient(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
        particle_diameter, Mw, superficial_velocity, molar_concentration)[2]
    j_m = 0.61 * reynolds_number ** -0.41
    h_GP = (j_m * gas_density * gas_heat_capacity * superficial_velocity) / (
                (gas_heat_capacity * gas_viscosity / conductivity_gas) ** (2 / 3))
    return h_GP


def compute_p_saturated(A, B, temp_kelvin, C):  # Double-checked and clear! :)
    temp_celsius = temp_kelvin - 273.15
    p_saturated = 10 ** (A - B / (temp_celsius + C))
    p_saturated_pascal = p_saturated * 133.322  # Torr to Pascal
    return p_saturated_pascal


def compute_partial_pressure_moisture(molar_concentration, R_gas_constant, temperature):  # c = molar_concentration
    partial_pressure_moisture = molar_concentration * R_gas_constant * temperature
    return partial_pressure_moisture

def compute_partial_pressure_moisture2(molar_mass, R_gas_constant, temperature, gas_density):
    pp = gas_density * R_gas_constant/molar_mass * temperature
    return pp 


def compute_gradient(vector, index, space_step):    # Should work for temperature and moisture, where vector
                                        #is an array input e.g. moisture gas and the index is the looping variable x
    length = len(vector)
    if index > (length-2):
        grad = 0 
    else: 
        grad = (vector[index+1] - vector[index])/space_step
    return grad

def compute_laplacian(vector, index, space_step):
    length = len(vector)
    if index > (length-2):
        laplacian = 0
    elif index == 0: 
        laplacian = 0
    else: 
        laplacian = (vector[index-1] + 2*vector[index] + vector[index+1])/(space_step**2)
    return laplacian


def compute_molar_concentration(relative_humidity, pressure_saturated, R, temp):
   molar_concentration = relative_humidity * pressure_saturated / (R * temp)        #TODO: erased *0.01
   return molar_concentration


# def compute_partial_pressure_moisture(molar_concentration, R_gas_constant, temp): # c = molar_concentration
#     partial_pressure_moisture = R_gas_constant * temp * molar_concentration
#     return partial_pressure_moisture


# def compute_relative_humidity(partial_p_moisture, p_saturated):
#     relative_humidity = partial_p_moisture/(p_saturated)
#     if relative_humidity > 1:
#         print('Relative humidity larger than 1, error!')
#     return relative_humidity
