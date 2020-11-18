# IMPORTS
import numpy as np
import math 

def volumetric_flow_rate_meters_per_second(volumetric_flow_rate_liters_per_minute):
    volumetric_flow_rate = volumetric_flow_rate_liters_per_minute / (60 * 10 ** 3)  # from l/min -> cubic meters per s
    return volumetric_flow_rate

def compute_velocity(volumetric_flow_rate_liters_per_minute, length, diameter, volume_fraction_powder):
    volumetric_flow_rate = volumetric_flow_rate_meters_per_second(volumetric_flow_rate_liters_per_minute)
    area_column = (diameter / 2) ** 2 * np.pi
    fraction_gas = 1 - volume_fraction_powder
    # volume_gas_in_tube = area_column * length * fraction_gas
    velocity = volumetric_flow_rate / (area_column * fraction_gas)  # only area with gas, not powder
    return velocity


def compute_equilibrium_moisture(alpha, moisture_particle, N):
    f_of_x = 1 - np.exp(-alpha * moisture_particle ** N)
    return f_of_x


def compute_moisture_particle(moisture_particle, alpha, N, relative_humidity, dt, constant):
    change_moisture = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle, N))
    moisture_particle_current = moisture_particle + change_moisture * dt
    return moisture_particle_current


def compute_temperature_particle(temperature, constant, dt, conductivity, laplacian, density, alpha, moisture,
                                 relative_humidity, N, heat_of_vaporization, heat_transfer_coefficient, specific_surface,
                                 temperature_gas, heat_capacity):
    conduction = conductivity * laplacian / density
    heat_of_sorption = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture, N)) * heat_of_vaporization
    heat_transfer = heat_transfer_coefficient * specific_surface * (temperature_gas-temperature)

    # change_temperature = (conduction + heat_of_sorption + heat_transfer)/heat_capacity
    change_temperature = (heat_of_sorption + heat_transfer) / heat_capacity
    temperature += change_temperature * dt
    return temperature


def compute_p_saturated(A, B, temp_kelvin, C):
    temp_celsius = temp_kelvin - 273.15
    p_saturated = np.exp(A-B/(temp_celsius+C))
    p_saturated_pascal = p_saturated * 133.322
    return p_saturated_pascal


def spec_surface_area(particle_diameter, particle_density): #S_P in m^2/kg - assuming spherical particles S_P=surface_area/volume*density
    r=particle_diameter/2
    SSA=(4*np.pi*r**2)/(4/3*np.pi*r**3*particle_density)
    return SSA


def molar_mass_moisture_kg(molar_mass_moisture):
    Mw=molar_mass_moisture/1000
    return Mw
    

def mass_transfer_coefficient(moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, 
                              particle_density, flow_rate, particle_diameter, Mw, relative_velocity_gp, molar_concentration):
    SSA=spec_surface_area(particle_diameter, particle_density)
    superficial_mass_velocity=(4*particle_density*flow_rate)/(np.pi*column_diameter**2) #G_0
    particle_surface_area=porosity_powder*particle_density*SSA #a
    reynolds_number=superficial_mass_velocity/particle_surface_area*gas_viscosity #Re
    denominator=gas_viscosity/(gas_density*moisture_diffusivity)
    j_m=0.61* reynolds_number ** -0.41
    k_gp=j_m*molar_concentration*relative_velocity_gp*Mw/denominator ** (2/3)
    return superficial_mass_velocity, particle_surface_area, reynolds_number, k_gp
