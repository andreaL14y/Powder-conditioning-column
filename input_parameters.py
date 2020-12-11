import numpy as np
from define_functions import*

########################################### PARAMETERS #################################################################
# Constants
R_gas_constant = 8.314                                                          # ideal gas constant, J/(K*mol)
kelvin = 273.15                                                                 # Conversion to degrees
A = 8.07131                                                                     # Antoine constant for water
B = 1730.630                                                                    # Antoine constant for water
C = 233.426                                                                     # Antoine constant for water

# Material specific for gas and powder
porosity_powder = 0.6
N = 1                                                                           # Parameter, assume 1
alpha_parameter = 25                                                            # parameter, 10 < alpha < 100
gas_density = 1
particle_density = 1500
particle_diameter = 0.00001                                                     # m
heat_of_vaporization = 1000*1000                                                # delta_H, J/kg
gas_viscosity = 10 ** -5                                                        # mu_G, kg/(m*s)
moisture_diffusivity = 10 ** -5                                                 # D_G, m^2/s
molar_mass_moisture = 18/1000                                                   # kg/mol for water vapor
molar_mass_dry_air = 28.97/1000                                                 # kg/mol
specific_surface_area = spec_surface_area(particle_diameter, particle_density)  # m2/kg, SSA

moisture_vapor_heat_capacity = 2000         # J/(kg*K), C_P_V, heat capacity gas water vapor)
moisture_liquid_heat_capacity = 4000        # J/(kg*K) C?? Heat capacity liquid water
particle_heat_capacity = 1000               # C_P_P AND C_P_WP
gas_heat_capacity = 1000                    # C_P_WG
conductivity_particle = 0.1                 # W/(m*K)
conductivity_gas = 0.01                     # W/(m*K)
boiling_temp = kelvin + 100                 # for water

# Cylinder and flow specific
bed_length = 0.2                                                                        # m
column_diameter = 0.1                                                                   # m
cross_sectional_area = np.pi * (column_diameter/2)**2                                   # m^2
volumetric_flow_rate_liters_per_minute = 1                                              # l/min
flow_rate = volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute)  # m3/s
superficial_velocity = flow_rate/(np.pi*(column_diameter/2)**2)                         # superficial velocity U in m/s
gas_velocity = compute_velocity(volumetric_flow_rate_liters_per_minute, bed_length, column_diameter, porosity_powder)

temp_initial = kelvin + 20                  # K, room temperature 20 degrees
relative_humidity_bed_initial = 0.2         # starting condition
relative_humidity_gas_initial = 0.9         # humidity of flowing gas
pressure_ambient = 101325                   # atmospheric pressure, Pa

# Unused?
volume_total = cross_sectional_area * bed_length

################################## INITIAL CONDITIONS ##################################################################
pressure_saturated_initial = compute_p_saturated(A, B, temp_initial, C)
partial_pressure_moisture_initial = pressure_saturated_initial * relative_humidity_gas_initial
molar_concentration_moisture_initial = compute_molar_concentration(
    relative_humidity_gas_initial, pressure_saturated_initial, R_gas_constant, temp_initial)

moisture_gas_initial_bed = compute_Y_from_RH(molar_mass_dry_air, molar_mass_moisture, pressure_ambient,
                                         relative_humidity_bed_initial, pressure_saturated_initial)

moisture_gas_initial_in = compute_Y_from_RH(molar_mass_dry_air, molar_mass_moisture, pressure_ambient,
                                         relative_humidity_gas_initial, pressure_saturated_initial)

moisture_particle_initial = compute_initial_moisture_particle(alpha_parameter, N, relative_humidity_bed_initial)
moisture_particle_saturated = compute_initial_moisture_particle(alpha_parameter, N, relative_humidity_gas_initial)

k_GP_initial = compute_mass_transfer_coefficient(
    moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
    particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture_initial)[3]

constant_initial = k_GP_initial * specific_surface_area * pressure_saturated_initial / pressure_ambient # just some simplification

heat_transfer_coefficient_initial = compute_heat_transfer_coefficient(
    moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, flow_rate,
    particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture_initial,
    gas_heat_capacity, conductivity_gas)

laplacian_initial = 0
temp_gradient_initial = 0
gradient_moisture_initial = 0
laplacian_moisture_initial = 0

########################################### TESTING ####################################################################
# alpha_test = 75
# RH = 0.9
# gas_m = compute_Y_from_RH(molar_mass_dry_air, molar_mass_moisture, pressure_ambient, 0.9, pressure_saturated_initial)
# x_sat = compute_initial_moisture_particle(alpha_test, 1, RH)
# print(x_sat)
# print(1/alpha_test)

########################################### UNUSED #####################################################################
# moisture_density = 4.85 * 10**(-3)                                              # water vapor density at room temp
# molar_concentration_moisture = moisture_density / molar_mass_moisture           # c at room temperature , moles/m3
# molar_concentration_moisture_initial = partial_pressure_moisture_initial/(R_gas_constant * temp_initial)


# mass_gas = gas_density * volume_gas
# moles_gas = mass_gas/molar_mass_dry_air                                        # n in ideal gas law
# molar_concentration_dry_air = moles_gas/volume_gas
