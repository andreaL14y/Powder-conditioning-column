import numpy as np
from define_functions import*

################# PARAMETERS #######################
porosity_powder = 0.6
N = 1                # Parameter, assume 1
alpha_parameter = 75 # parameter, 10 < alpha < 100

gas_density = 1
particle_density = 1500
heat_of_vaporization = 1000 # delta_H
molar_mass_moisture = 18 # g/mol for water

# Velocity
bed_length = 0.2 # m
column_diameter = 0.1
volumetric_flow_rate = 1 # l/min
velocity = compute_velocity(volumetric_flow_rate, bed_length, column_diameter, porosity_powder) # m/s
superficial_velocity = volumetric_flow_rate/(np.pi*(column_diameter/2)**2) # U

# C_P_P, C_P_WP, C_P_WG, C_P_V and lambdas
moisture_vapor_heat_capacity = 2000 # J/(kg*K), C_P_V
moisture_liquid_heat_capacity = 4000 # J/(kg*K) C_P_WP
particle_heat_capacity = 1000 # C_P_P
gas_heat_capacity = 1000 # C_P_WG

conductivity_particle = 0.1 # W/(m*K)
conductivity_gas = 0.01 # W/(m*K)

# Antoine constants for water
A = 8.07131
B = 1730.630
C = 233.426

# Made up by me
k_GP = 1                # TODO: compute
specific_surface_area = 1
pressure_saturated = 1  # TODO: compute, knowing p_sat = 101325 Pa at T = T_boil, compute_p_saturated(A, B, temp_kelvin, C)
pressure = 1            # TODO: what is it actually?
constant = k_GP*specific_surface_area*pressure_saturated/pressure # just some simplification
relative_humidity = 0.75 # Guessed parameter
laplacian = 1           # TODO: compute (later, not present in simplification)
heat_transfer_coefficient = 1 # TODO: what it h_GP?
dt = 0.01

# Initial conditions
initial_temp = 300 # K, room temperature-ish
moisture = 0

# Just for testing
moisture_particle_next = 0
moisture_particle_old = 1

current_temp = initial_temp
i = 0
while moisture_particle_next != moisture_particle_old:
    i += 1
    moisture_particle_old = moisture_particle_next
    moisture_particle_next = compute_moisture_particle(moisture_particle_old, alpha_parameter, N, relative_humidity, dt, constant)
    # print(moisture_particle_next)

    current_temp = compute_temperature_particle(current_temp, constant, dt, conductivity_particle, laplacian, particle_density,
                                                alpha_parameter, moisture, relative_humidity, N, heat_of_vaporization,
                                                heat_transfer_coefficient, specific_surface_area, initial_temp,
                                                moisture_liquid_heat_capacity)
    # print(current_temp)
print(i)
print(i*dt)

print(compute_p_saturated(A, B, initial_temp, C)) # pressure in Pascal
