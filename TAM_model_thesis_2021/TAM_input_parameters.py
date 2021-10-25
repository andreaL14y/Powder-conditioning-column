import numpy as np
import time
import scipy.optimize
from scipy.integrate import odeint, simps
from scipy.signal import savgol_filter

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

############################################# COLORS ###################################################################
c_palette = np.array([[0, 18, 25],     [0, 95, 115],    [10, 147, 150],
                      [148, 210, 189], [233, 216, 166], [251, 248, 239],
                      [238, 155, 0],   [202, 103, 2],   [187, 62, 3],
                      [174, 32, 18],   [155, 34, 38]])
c_palette       = c_palette/255
fig_size        = (10, 6.5)
grid_color      = c_palette[1]
grid_alpha      = 0.4
grid_alpha_min  = 0.1
tabs            = 40

########################################### EQUATION SPECIFIC ##########################################################
# Constants
r_gas_constant      = 8.314                             # ideal gas constant, J/(K*mol)
kelvin              = 273.15                            # Conversion to degrees
pressure_ambient    = 101325                            # atmospheric pressure, Pa

M0_cr               = 1.68 * 10 ** (-4)                 # GAB as fitted by Bronlund
c_cr                = 8.8                               # GAB as fitted by Bronlund
f_cr                = 0.878                             # GAB as fitted by Bronlund
h_cr                = 30                                # GAB as fitted by Bronlund

M0_am               = 0.0488                            # GAB as fitted by Bronlund
c_am                = 3.23                              # GAB as fitted by Bronlund
f_am                = 1.16                              # GAB as fitted by Bronlund

M0_am2              = 0.0491                            # as fitted by Roos & Karel, not used in this work
c_am2               = 4.33                              # as fitted by Roos & Karel, not used in this work
f_am2               = 1.18                              # as fitted by Roos & Karel, not used in this work

n_A                 = 3                                 # Avrami-Bronlund kinetics parameter
avrami_exponent     = (n_A - 1) / n_A                   # Avrami-Bronlund kinetics parameter
c_1                 = 3.54 * 10**4                      # Avrami-Bronlund kinetics parameter
c_2                 = 108.4                             # Avrami-Bronlund kinetics parameter
c_3                 = 3 * 10**27                        # Avrami-Bronlund kinetics parameter

k_param = 6.7                                           # Gordon & Taylor parameter
########################################### MATERIAL SPECIFIC ##########################################################
# Water, liquid
A                   = 8.07131                                               # Antoine constant for water
B                   = 1730.630                                              # Antoine constant for water
C                   = 233.426                                               # Antoine constant for water
glass_temp_water    = 138.15                                                # glass transition water, Kelvin
density_water       = 1000                                                  # kg/m3
heat_capacity_water = 4200                                                  # J/(kg*K), C_M, heat cap liquid water

# Water vapor
molar_mass_dry_air  = 28.97/1000                                            # kg/mol
conductivity_gas    = 0.01                                                  # J/(s*m*K), lambda, conductivity gas
heat_capacity_air   = 1000                                                  # J/(kg*K), C_P_WG, heat cap air
heat_capacity_vapor = 2000                                                  # J/(kg*K), C_PV, heat cap water vapor
density_gas         = 1                                                     # kg/m^3
moisture_diffusivity        = 2.42 * 10 ** (-5)                             # D_G, m^2/s
heat_of_evaporation_water   = 2.5 * 10**(6)                                 # latent heat of evaporation water J/kg
heat_of_water_binding       = heat_of_evaporation_water * 1.0025            # latent heat of evaporation water J/kg, 1.005 for exp 29
heat_binding_diffs          = heat_of_water_binding - heat_of_evaporation_water     # differefence heat of binding to lactose va evaporation, J/kg

# Lactose
glass_temp_lactose  = 114 + kelvin                                          # Kelvin
density_particle    = 1540                                                  # kg/m^3
heat_capacity_particle  = 1252                                              # ca 1220 C_ps heat cap particle, J/(kg * K)
conductivity_particle   = 0.1                                               # W/(m*K), lambda, conductivity p

# Powder
density_powder      = 218                                                   # kg/, as measured at AZ
porosity_powder     = 1 - density_powder/(density_particle)                 # ca 0.85
heat_of_crystallization = 32 * 1000                                         # delta_H, J/kg, enthalpy of crystallization
diffusivity_eff     = moisture_diffusivity * porosity_powder * 0.7          # turtousity and constrictivity
conductivity_tot    = 0.17

# Sample
out_diam_TAM_cyl    = 0.013                                                 # m, measured as 13 mm
in_diam_TAM_cyl     = 0.011                                                 # m, measured as 11 mm
wall_thickness_TAM_cyl  = (out_diam_TAM_cyl- in_diam_TAM_cyl)/2
in_area_TAM_cyl     = np.pi * (in_diam_TAM_cyl / 2) ** 2                    # m2

height_TAM_cyl      = 0.035
height_air_TAM_cyl  = 0.019

in_diam_ampoule         = 0.003                                             # m, measured 3 mm
area_ampoule            = np.pi * (in_diam_ampoule / 2) ** 2                # m2
relative_humidity_bed_initial   = 0.2                                       # humidity in bed, starting condition

####################################### EXPERIMENT SPECIFIC ############################################################
# Choose number of discretized sections to compute; increase for finer, but also slower, results
n_space_steps           = 3

# Choose which experiments to run
exp_list = np.array([13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 28])

# Choose to compute with or without hindered diffusion in amorphous matrix
hindered_am             = True
# hindered_am             = False

temp_initial_celsius    = 25
temp_initial            = temp_initial_celsius + kelvin

# Setup conditions according to experiment
def create_conditions(exp):
    if exp == 13:
        relative_humidity_gas_inlet = 0.5289    # humidity of salt solution
        weight = 154                            # measured weight of sample
        batch = 2
        hours = 5.5                             # adjusted as to include full crystallization peak
        amorphous_material_initial = 0.21       # different for same batch depending on time between runs

    elif exp == 14 or exp==21 or exp==25 or exp==26:
        relative_humidity_gas_inlet = 0.5757
        batch = 2
        if exp==14:
            weight = 153
            amorphous_material_initial = 0.21
        elif exp == 21:
            amorphous_material_initial = 0.2
            weight = 156
        elif exp == 25:
            weight = 151
            amorphous_material_initial = 0.19
        elif exp == 26:
            weight = 150
            amorphous_material_initial = 0.19
        hours = 3.5

    elif exp == 15:
        relative_humidity_gas_inlet = 0.753
        weight = 155
        amorphous_material_initial = 0.21
        hours = 2.5
        batch = 2

    elif exp == 16:
        batch = 1
        relative_humidity_gas_inlet = 0.5289
        weight = 157
        amorphous_material_initial = 0.14
        hours = 4

    elif exp == 17:
        batch = 1
        relative_humidity_gas_inlet = 0.5757
        weight = 158
        amorphous_material_initial = 0.14
        hours = 3

    elif exp == 18:
        batch = 1
        relative_humidity_gas_inlet = 0.753
        weight = 157
        amorphous_material_initial = 0.14
        hours = 3

    elif exp == 19:
        batch = 2
        relative_humidity_gas_inlet = 0.3817
        weight = 307
        amorphous_material_initial = 0.2
        hours = 15

    elif exp == 20:
        batch = 2
        relative_humidity_gas_inlet = 0.5289
        weight = 305
        amorphous_material_initial = 0.2
        hours = 8.5

    elif exp == 22:
        batch = 1
        relative_humidity_gas_inlet = 0.3817
        weight = 307
        amorphous_material_initial = 0.13
        hours = 15

    elif exp == 23:
        batch = 1
        relative_humidity_gas_inlet = 0.5289
        weight = 301
        amorphous_material_initial = 0.13
        hours = 10

    elif exp == 24:
        batch = 1
        relative_humidity_gas_inlet = 0.3278
        weight = 307
        amorphous_material_initial = 0.13
        hours = 20

    elif exp == 28:
        batch = 2
        relative_humidity_gas_inlet = 0.5289
        weight = 307
        amorphous_material_initial = 0.0009             # uncertain value
        hours = 2

    else:
        print('No experiment defined.')
        particle_diameter = 2 * 10 ** (-6)
        relative_humidity_gas_inlet = 0.58
        amorphous_material_initial = 0.16                           # 16 % am left in batch L6.5 before cond, page 82 thesis

    weight *= 10 **(-6)                             # mg to kg
    if batch == 1:
        particle_diameter   = 2.75 * 10**(-6)       # measured for batch 1
    elif batch == 2:
        particle_diameter   = 1.89 * 10**(-6)       # measured for batch 2

    resolution = int(hours * 800)                   # number of outputs, can be adjusted, but this seems to work well
    if exp == 28:
        resolution = int(hours * 200)               # fewer points needed when nothing happens
    return batch, amorphous_material_initial, relative_humidity_gas_inlet, particle_diameter, weight, hours, resolution

################################ HINDERED DIFFUSION AMORPHOUS PARTICLES ################################################
am_layer_thickness  = 0.25*10**(-6)
am_diffusivity      = 9.12*10**(-17)