import numpy as np
from define_functions import*

########################################### PARAMETERS #################################################################
porosity_powder = 0.6
N = 1                                                                           # Parameter, assume 1
alpha_parameter = 75                                                            # parameter, 10 < alpha < 100
gas_density = 1
particle_density = 1500
particle_diameter= 0.00001                                                      # m
heat_of_vaporization = 1000                                                     # delta_H
gas_viscosity = 10 ** -5                                                        # mu_G in kg/(m*s)
moisture_diffusivity = 10 ** -5                                                 # D_G
R_gas_constant = 8.314                                                          # J/(K*mol) ideal gas constant
molar_mass_moisture = 18/1000                                                   # kg/mol for water
# molar_mass_moisture = molar_mass_moisture_kg(molar_mass_moisture)             # converted in kg/mol
bed_length = 0.2                                                                # m
column_diameter = 0.1                                                           # m
specific_surface_area = spec_surface_area(particle_diameter, particle_density)  # 1/m, SSA
volumetric_flow_rate_liters_per_minute = 1                                      # l/min
flow_rate = volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute)  # m3/s
gas_velocity = compute_velocity(volumetric_flow_rate_liters_per_minute, bed_length, column_diameter, porosity_powder) # u in m/s TODO: u compared to U?
superficial_velocity = flow_rate/(np.pi*(column_diameter/2)**2)                 # superficial velocity U in m/s

# Heat capacities (C's) and conductivities (lambdas)
moisture_vapor_heat_capacity = 2000         # J/(kg*K), C_P_V, heat capacity gas water (steam)
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
initial_temp = 293.15               # K, room temperature 20 degrees
initial_relative_humidity = 0.2     # starting condition
relative_humidity_gas = 0.5         # humidity of flowing gas
pressure_ambient = 101325           # atmospheric pressure, Pa
dt = 0.01                           # for now

######################################## VARIABLES #####################################################################
pressure_saturated = compute_p_saturated(A, B, initial_temp, C)
relative_humidity = compute_relative_humidity()
molar_concentration = compute_molar_concentration(relative_humidity, pressure_saturated, R_gas_constant, initial_temp) # mol/m^3 corresp. to c in our equation
k_GP = mass_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density,
                                 particle_density, flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity,
                                 molar_concentration)[3] # not yet the true value since p_sat is needed for the computation of c

constant = k_GP * specific_surface_area * pressure_saturated / pressure_ambient # just some simplification

#################################### TEMPORARY SIMPLIFICATIONS #########################################################
laplacian = 1                   # TODO: compute (later, not present in simplification)
temp_gradient = 1               # TODO: compute (but how?)
heat_transfer_coefficient = 1   # TODO: what is h_GP?
