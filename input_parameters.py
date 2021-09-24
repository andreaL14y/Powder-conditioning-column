import numpy as np
import time
import scipy.optimize
from scipy.integrate import odeint
# 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

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

########################################### PARAMETERS #################################################################
# Constants
r_gas_constant = 8.314                                                          # ideal gas constant, J/(K*mol)
h_planck = 6.62607004 * 10 **(-34)                                              # Planck's constant, m2 * kg / s
k_boltzmann = 1.38064852 * 10 **(-23)                                           # Boltzmann constant, m2 * kg / (s^2 * K)
kelvin = 273.15                                                                 # Conversion to degrees
pressure_ambient = 101325                                                       # atmospheric pressure, Pa

# Water, liquid
A = 8.07131                                                                     # Antoine constant for water
B = 1730.630                                                                    # Antoine constant for water
C = 233.426                                                                     # Antoine constant for water
# glass_temp_water_1 = 136                                                      # glass transition water, Kelvin
glass_temp_water_1 = 138.15                                                     # glass transition water, Kelvin
molar_mass_moisture = 18/1000                                                   # kg/mol for water
density_water = 1000                                                            # kg/m3
heat_capacity_water = 4200                                                      # J/(kg*K), C_M, heat cap liquid water
boiling_temp = kelvin + 100                                                     # for water

# Water vapor
molar_mass_dry_air = 28.97/1000                                                 # kg/mol
gas_viscosity = 10 ** -5                                                        # mu_G, kg/(m*s)
conductivity_gas = 0.01                                                         # J/(s*m*K), lambda, conductivity gas
heat_capacity_air = 1000                                                        # J/(kg*K), C_P_WG, heat cap air
heat_capacity_vapor = 2000                                                      # J/(kg*K), C_PV, heat cap water vapor
density_gas = 1                                                                 # kg/m^3
moisture_diffusivity = 2.42 * 10 ** (-5)                                        # D_G, m^2/s
heat_of_evaporation_water = 2.5 * 10**(6)                                       # latent heat of evaporation water J/kg
heat_of_water_binding = heat_of_evaporation_water * 1.3                         # latent heat of evaporation water J/kg
heat_binding_diffs = heat_of_water_binding - heat_of_evaporation_water          # difference between binding lactose/evap, J/kg. Positive

# Lactose
glass_temp_lactose = 101 + kelvin                                               # Kelvin
glass_temp_lactose = 114 + kelvin                                               # Kelvin TODO: this value works much better! Gives broader and lower integral
density_particle = 1540                                                         # kg/m^3
heat_capacity_particle = 1252                                                   # ca 1220 C_ps heat cap particle, J/(kg * K)
conductivity_particle = 0.1                                                     # W/(m*K), lambda, conductivity p

# Powder
particle_diameter = 0.000002                                                    # m, from page 46 in thesis
density_powder = 218                                                            # kg/m3
porosity_powder = 1 - density_powder/(density_particle)                         # 0.85
heat_of_crystallization = 43.1 * 1000                                           # delta_H, J/kg, enthalpy of crystallization
heat_of_crystallization = 32 * 1000                                             # delta_H, J/kg, enthalpy of crystallization
specific_surface_area = 6 * (1- porosity_powder)/(particle_diameter * density_particle)       # m2/kg, SSA
diffusivity_eff = moisture_diffusivity * porosity_powder * 0.5                  # m2/s
diffusivity_eff = moisture_diffusivity * porosity_powder * 0.75                 # TODO: new
conductivity_tot = porosity_powder * conductivity_gas + (1 - porosity_powder) * conductivity_particle
conductivity_tot = 0.17

# Sample
weight = 0.00015                                                                # kg, 150 mg
volume = weight/density_powder                                                  # m3
out_diam_TAM_cyl = 0.013                                                        # m, measured as 11 mm
in_diam_TAM_cyl = 0.011                                                         # m, measured as 11 mm
wall_thickness_TAM_cyl = (out_diam_TAM_cyl- in_diam_TAM_cyl)/2                  # m, measured as 11 mm
# radius_vial = out_diam_TAM_cyl / 2
in_area_TAM_cyl = np.pi * (in_diam_TAM_cyl / 2) ** 2                            # m2
total_height_powder = volume / in_area_TAM_cyl                                  # how high is powder filled

amorphous_material_initial = 0.16                                               # 16 % am left in batch L6.5 before cond, page 82 thesis
# amorphous_material_initial = 0.01                                             # 16 % am left in batch L6.5 before cond, page 82 thesis
height_TAM_cyl = 0.035
height_air_TAM_cyl = 0.019
height_air_TAM_cyl = height_air_TAM_cyl - total_height_powder
volume_ampoule = in_area_TAM_cyl * height_air_TAM_cyl
fraction_powder = volume/volume_ampoule
conductivity_glass = 0.8                # W/( m K )

in_diam_ampoule = 0.003                                                         # m, measured 3 mm
area_ampoule = np.pi * (in_diam_ampoule / 2) ** 2                               # m2
temp_initial_celsius = 25
temp_initial = kelvin + temp_initial_celsius                                    # K, room temperature 24 degrees
relative_humidity_bed_initial = 0.2                                             # humidity in bed, starting condition
relative_humidity_gas_inlet = 0.58                                              # humidity of flowing gas
relative_humidity_gas_inlet = 0.5757                                            # humidity of flowing gas

contact_area_powder = in_area_TAM_cyl + np.pi * in_diam_TAM_cyl * total_height_powder
