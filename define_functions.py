import numpy as np


def compute_velocity(volumetric_flow_rate, length, diameter, volume_fraction_powder):
    volumetric_flow_rate = volumetric_flow_rate / (60 * 10 ** 3)  # from l/min -> cubic meters per s
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
