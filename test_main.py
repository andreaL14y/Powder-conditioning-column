import numpy as np
from define_functions import*

########################################### PARAMETERS #################################################################
porosity_powder = 0.6
N = 1                                                                           # Parameter, assume 1
alpha_parameter = 75                                                            # parameter, 10 < alpha < 100
gas_density = 1
particle_density = 1500
moisture_density = 4.85 * 10**(-3)                                              # water vapor density at room temp
particle_diameter= 0.00001                                                      # m
heat_of_vaporization = 1000                                                     # delta_H
gas_viscosity = 10 ** -5                                                        # mu_G in kg/(m*s)
moisture_diffusivity = 10 ** -5                                                 # D_G
R_gas_constant = 8.314                                                          # J/(K*mol) ideal gas constant
molar_mass_moisture = 18/1000                                                   # kg/mol for water vapor
molar_mass_dry_air = 28.97/1000                                                 # kg/mol
bed_length = 0.2                                                                # m
column_diameter = 0.1                                                           # m
cross_sectional_area = np.pi * (column_diameter/2)**2
volume_total = cross_sectional_area * bed_length

# molar_concentration_moisture = moisture_density / molar_mass_moisture           # c at room temperature , moles/m3

specific_surface_area = spec_surface_area(particle_diameter, particle_density)  # m2/kg, SSA
volumetric_flow_rate_liters_per_minute = 1                                      # l/min
flow_rate = volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute)  # m3/s
superficial_velocity = flow_rate/(np.pi*(column_diameter/2)**2)                 # superficial velocity U in m/s

# Heat capacities (C's) and conductivities (lambdas)
moisture_vapor_heat_capacity = 2000         # J/(kg*K), C_P_V, heat capacity gas water vapor)
moisture_liquid_heat_capacity = 4000        # J/(kg*K) C?? Heat capacity liquid water
particle_heat_capacity = 1000               # C_P_P AND C_P_WP
gas_heat_capacity = 1000                    # C_P_WG

conductivity_particle = 0.1                 # W/(m*K)
conductivity_gas = 0.01                     # W/(m*K)

# Antoine constants for water
A = 8.07131
B = 1730.630
C = 233.426

################################## INITIAL CONDITIONS ##################################################################
temp_initial = 293.15                       # K, room temperature 20 degrees
relative_humidity_bed_initial = 0.2         # starting condition
relative_humidity_gas_initial = 0.5         # humidity of flowing gas
pressure_ambient = 101325                   # atmospheric pressure, Pa
# dt = 0.01                                   # for now
boiling_temp = 273.15 + 100                 # for water

pressure_saturated_initial = compute_p_saturated(A, B, temp_initial, C)
partial_pressure_moisture_initial = pressure_saturated_initial * relative_humidity_gas_initial

molar_concentration_moisture_initial = partial_pressure_moisture_initial/(R_gas_constant * temp_initial) # TODO: Delete? What did we do here? Thats old right? we have c defined above already
# print('Johans c: ', molar_concentration_moisture_initial)
# print('Our c: ', molar_concentration_moisture)


moisture_particle_initial = compute_initial_moisture_particle(alpha_parameter, N, relative_humidity_bed_initial)

k_GP_initial = compute_mass_transfer_coefficient(
    moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
    particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture_initial)[3]

constant_initial = k_GP_initial * specific_surface_area * pressure_saturated_initial / pressure_ambient # just some simplification

heat_transfer_coefficient_initial = compute_heat_transfer_coefficient(
    moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
    particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture_initial, gas_heat_capacity)

# initial_moisture_particle = compute_initial_moisture_particle(alpha_parameter, N, relative_humidity_bed_initial) #TODO: delete, same as moisture_particle_initial above

# print('k_GP: ', k_GP_initial)
# print('h_GP: ', heat_transfer_coefficient_initial)

#################################### TEMPORARY SIMPLIFICATIONS #########################################################
laplacian = 1                   # TODO: compute (later, not present in simplification)
temp_gradient = 1               # TODO: compute (but how?)
gradient_moisture = 1           # TODO: compute (later, not present in simplification)
laplacian_moisture = 1          # TODO: compute (later, not present in simplification)

########################################### TESTING ####################################################################
moisture_converged = 0.00626672
test = compute_equilibrium_moisture(alpha_parameter, moisture_converged, N)
# print(test)

########################################### UNUSED #####################################################################
# volume_gas = (1-porosity_powder) * volume_total
# mass_gas = gas_density * volume_gas
# moles_gas = mass_gas/molar_mass_dry_air                                        # n in ideal gas law
# molar_concentration_dry_air = moles_gas/volume_gas

# print(molar_concentration_dry_air, 'molar concentration')
# molar_volume_gas = volume_gas/moles_gas


# partial_pressure = compute_partial_pressure_moisture(molar_concentration_dry_air, R_gas_constant, initial_temp)
# print('Saturated pressure: ', pressure_saturated)
# print('Partial pressure: ', partial_pressure)