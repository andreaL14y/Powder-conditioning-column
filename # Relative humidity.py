# Relative humidity
import math
import numpy as np
from test_main import*

# Relative humidity

def relative_humidity(vapor_pressure, temp_kelvin):
    p_saturated=compute_p_saturated(A, B, temp_kelvin, C)
    RH = vapor_pressure/p_saturated
    return RH

def vapor_pressure_moisture(RH, temp_kelvin):
    p_saturated = compute_p_saturated(A, B, temp_kelvin, C)
    vapor_pressure_moisture = RH * p_saturated
    return vapor_pressure_moisture

#### initial conditions
RH_init=0.2
def X_P_init(alpha, N, RH_init):
    X_P_init=(-np.log(-(RH_init -1))/alpha)**(1/N)
    return X_P_init

#### h_GP
def heat_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, 
                              particle_density, flow_rate, particle_diameter, Mw, superficial_velocity, molar_concentration,
                              gas_heat_capacity):
    superficial_mass_velocity, particle_surface_area, reynolds_number, k_gp= compute_mass_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density,
                                                                                                               particle_density, flow_rate, particle_diameter, Mw, superficial_velocity, molar_concentration)
    j_m=0.61* reynolds_number ** -0.41
    h_GP=(j_m*gas_density*gas_heat_capacity*superficial_velocity)/((gas_heat_capacity*gas_viscosity/k_gp)**(2/3))
    return h_GP

print(heat_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density,
                                particle_density, flow_rate, particle_diameter, Mw, superficial_velocity, molar_concentration_moisture,
                                gas_heat_capacity))

