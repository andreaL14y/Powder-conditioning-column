import numpy as np

########################################### PARAMETERS #################################################################
# Constants
r_gas_constant = 8.314                                                          # ideal gas constant, J/(K*mol)
h_planck = 6.62607004 * 10 **(-34)                                              # Planck's constant, m2 * kg / s
k_boltzmann = 1.38064852 * 10 **(-23)                                           # Boltzmann constant, m2 * kg / (s^2 * K)
kelvin = 273.15                                                                 # Conversion to degrees

# Water
A = 8.07131                                                                     # Antoine constant for water
B = 1730.630                                                                    # Antoine constant for water
C = 233.426                                                                     # Antoine constant for water
glass_temp_water_1 = 136                                                        # glass transition water, Kelvin
glass_temp_water_2 = 165                                                        # glass transition water alt 2, Kelvin
gas_viscosity = 10 ** -5                                                        # mu_G, kg/(m*s)
moisture_diffusivity = 10 ** -5                                                 # D_G, m^2/s
molar_mass_moisture = 18/1000                                                   # kg/mol for water vapor
molar_mass_dry_air = 28.97/1000                                                 # kg/mol

moisture_vapor_heat_capacity = 2000                                             # J/(kg*K), C_PV, heat cap water vapor
moisture_liquid_heat_capacity = 4000                                            # J/(kg*K), C_M, heat cap liquid water
boiling_temp = kelvin + 100                                                     # for water
gas_heat_capacity = 1000                                                        # C_P_WG, heat cap gas
conductivity_gas = 0.01                                                         # W/(m*K), lambda, conductivity gas
gas_density = 1                                                                 # kg/m^3


# Lactose
glass_temp_lactose = 101 + kelvin                                               # Kelvin
particle_density = 1500                                                         # kg/m^3
particle_diameter = 0.000002                                                    # m, from page 46 in thesis
heat_of_sorption = (185.55 * 1000 / 342.3) * 1000                               # delta_H, J/kg, enthalpy of sorption, kJ/mol * mol/g * g/kg
heat_of_crystallization = 43.1 * 1000                                           # delta_H, J/kg, enthalpy of crystallization, kJ/mol * mol/g * g/kg
particle_heat_capacity = 417.6 * 1000/342.3                                     # C_P,P & C_P,WP, heat cap particle, J/(kg * K)
conductivity_particle = 0.1                                                     # W/(m*K), lambda, conductivity p


# Material specific for gas and powder
porosity_powder = 0.5
N = 1                                                                           # Parameter, assume 1
alpha_parameter = 25                                                            # parameter, 10 < alpha < 100


# Variables
amorphous_material_initial = 0.16                                               # 16 % am left in batch L6.5 before cond, page 82 thesis

bed_length = 0.032                                                               # m
bed_length = 0.024                                                               # m
# bed_length = 0.07                                                               # m
column_diameter = 0.024                                                          # m
cross_sectional_area = np.pi * (column_diameter/2)**2                            # m^2
volumetric_flow_rate_liters_per_minute = 0.5                                     # l/min
volumetric_flow_rate_liters_per_minute = 0.06                                     # l/min
rotation_time_interval = 19 * 60                                                 # s after which to rotate
rotation_time_interval = 58 * 60                                                 # s after which to rotate

temp_initial = kelvin + 35                                                       # K, room temperature 24 degrees
temp_initial = kelvin + 24                                                       # K, room temperature 24 degrees
temp_walls = temp_initial                                                         # At cylinder walls, cooling
relative_humidity_bed_initial = 0.3                                              # humidity in bed, starting condition
relative_humidity_gas_inlet = 0.45                                                # humidity of flowing gas
relative_humidity_gas_end = 0.3                                                  # humidity at end of cylinder
pressure_ambient = 101325                                                        # atmospheric pressure, Pa

mass_powder = (1 - porosity_powder) * particle_density * \
              bed_length * (column_diameter/2)**2 * 1000                         # weight converted to g
print(mass_powder, 'g')
