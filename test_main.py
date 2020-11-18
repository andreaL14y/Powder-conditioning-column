import numpy as np
from define_functions import*

################# PARAMETERS #######################
porosity_powder = 0.6
N = 1                # Parameter, assume 1
alpha_parameter = 75 # parameter, 10 < alpha < 100

gas_density = 1
particle_density = 1500
particle_diameter= 0.00001 #m
heat_of_vaporization = 1000 # delta_H
molar_mass_moisture = 18 # g/mol for water
Mw = molar_mass_moisture_kg(molar_mass_moisture) #converted in kg/mol

# Velocity
bed_length = 0.2 # m
column_diameter = 0.1
specific_surface_area = spec_surface_area(particle_diameter, particle_density) #SSA

volumetric_flow_rate_liters_per_minute = 1 # l/min
velocity = compute_velocity(volumetric_flow_rate_liters_per_minute, bed_length, column_diameter, porosity_powder) # m/s
flow_rate=volumetric_flow_rate_m3_per_second(volumetric_flow_rate_liters_per_minute) # TODO: in m^3/s -> That is Q? Was the value 0.005 wrong?

superficial_velocity = flow_rate/(np.pi*(column_diameter/2)**2) # U in m/s TODO: I exchanged flow_rate_per_minute with flow_rate(per second)
relative_velocity_gp = superficial_velocity/(np.pi*(column_diameter/2)**2) # :TODO check! m/s U in K_GP equation
print(superficial_velocity, relative_velocity_gp)

# C_P_P, C_P_WP, C_P_WG, C_P_V and lambdas
moisture_vapor_heat_capacity = 2000 # J/(kg*K), C_P_V
moisture_liquid_heat_capacity = 4000 # J/(kg*K) C_P_WP
particle_heat_capacity = 1000 # C_P_P
gas_heat_capacity = 1000 # C_P_WG

conductivity_particle = 0.1 # W/(m*K)
conductivity_gas = 0.01 # W/(m*K)

gas_viscosity = 10 ** -5 #mu_G in kg/(m*s)
moisture_diffusivity = 10 ** -5 #D_G

# Antoine constants for water
A = 8.07131
B = 1730.630
C = 233.426

R = 8.314 # J/K*mol ideal gas constant (fixed)

# Made up by me
temperature = 350 # TODO: seems like the molar concentration changes with the temperature i.e. is dependent of time and position in our model
pressure_saturated = 1  # TODO: compute, knowing p_sat = 101325 Pa at T = T_boil, compute_p_saturated(A, B, temp_kelvin, C)
pressure = 1            # TODO: what is it actually?
relative_humidity = 0.75 # Guessed parameter
molar_concentration= 0.01 * relative_humidity * pressure_saturated/R * temperature #mol/m^3 corresp. to c in our equation # TODO: check!!

k_GP = mass_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, 
                              particle_density, flow_rate, particle_diameter, Mw, relative_velocity_gp,
                                 molar_concentration)[3] #not yet the true value since p_sat is needed for the computation of c
k_GP = mass_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density,
                              particle_density, flow_rate, particle_diameter, Mw, superficial_velocity,
                                 molar_concentration)[3] #not yet the true value since p_sat is needed for the computation of c
print('k_GP is:', k_GP)
constant = k_GP*specific_surface_area*pressure_saturated/pressure # just some simplification
laplacian = 1           # TODO: compute (later, not present in simplification)
temp_gradient = 1       # TODO: compute (but how?)
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
    print(current_temp)
print(i)
print(i*dt)

print('P_sat: ', compute_p_saturated(A, B, initial_temp, C)) # pressure in Pascal

print('SSA: ', specific_surface_area)

# G_0, a, Re, k_gp=mass_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density,
#                               particle_density, flow_rate, particle_diameter, Mw, relative_velocity_gp, molar_concentration)
