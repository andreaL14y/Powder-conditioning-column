# IMPORTS
import math 

# Variables

flow_rate=0.005 #Q in m^3/s
particle_density=1500 #rho_P in kg/m^3
gas_density=1 #rho_G in kg/m^3
particle_porosity=0.6 #alpha_P
diameter=0.1 #D in m
gas_viscosity=10 ** -5 #mu_G in kg/m*s
moisture_diffusivity=10 ** -5 #D_G

# Unknowns (just set to some values)
superficial_surface_area=0.01 #S_P in m^2/kg
cU=0.5 #Unknown


def mass_transfer_coefficient():
    superficial_mass_velocity=(4*particle_density*flow_rate)/(math.pi*diameter**2) #G_0
    particle_surface_area=particle_porosity*particle_density*superficial_surface_area #a
    reynolds_number=superficial_mass_velocity/particle_surface_area*gas_viscosity #Re
    denominator=gas_viscosity/(gas_density*moisture_diffusivity)
    j_m=0.61* reynolds_number ** -0.41
    k_gp=j_m*cU//denominator ** (2/3)
    return superficial_mass_velocity, particle_surface_area, reynolds_number, k_gp

# Variables computed in the program
G_0, a, Re, k_gp=mass_transfer_coefficient()
print(k_gp)