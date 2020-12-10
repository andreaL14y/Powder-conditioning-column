import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from vectorized_define_functions import*
from define_functions import*
from test_main import*

def tank(h,t):
    c1 = 0.13
    c2 = 0.2
    Ac = 2
    qout1 = c1 * h[0]**0.5
    qout2 = c2 * h[1]**0.5

    qin = 1.5

    dhdt1 = (qin - qout1)/Ac
    dhdt2 = (qin - qout2)/Ac

    return [dhdt1, dhdt2]

def moisture(moisture_matrix, t, alpha, N, diffusivity, density_gas, density_particle, porosity, velocity,
             heat_transfer_coefficient, heat_capacity_vapor, heat_capacity_wet_gas, heat_capacity_particle, specific_surface_area,
             conductivity_gas, conductivity_particle, bed_length):

    # print(moisture_matrix)
    moisture_gas_vector = moisture_matrix[0]
    moisture_particle_vector = moisture_matrix[1]
    temp_gas_vector = moisture_matrix[2]
    temp_particle_vector = moisture_matrix[3]

    ##################################### UPDATE PARAMETERS ############################################################
    pressure_saturated = compute_p_saturated_vector(A, B, temp_gas_vector, C)

    relative_humidity = compute_relative_humidity_from_Y(
        molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas_vector, pressure_saturated)
    # print('RH :', relative_humidity)

    molar_concentration_moisture = compute_molar_concentration_vector(
        relative_humidity, pressure_saturated, R_gas_constant, temp_gas_vector)

    mass_transfer_coefficient = compute_mass_transfer_coefficient(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density,
        flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture)[3]

    constant = mass_transfer_coefficient * specific_surface_area * pressure_saturated / pressure_ambient

    laplacian_moisture_gas = compute_laplacian([moisture_gas_initial_in, moisture_gas_vector, moisture_gas_initial_bed], 1, 1)
    gradient_moisture_gas = compute_gradient([moisture_gas_initial_in, moisture_gas_vector, moisture_gas_initial_bed], 1, 1)

    laplacian_temp_gas = compute_laplacian([temp_initial, temp_gas_vector, temp_initial], index=1, space_step=bed_length)
    gradient_temp_gas = compute_gradient([temp_initial, temp_gas_vector, temp_initial], index=1, space_step=bed_length)

    laplacian_temp_particle = compute_laplacian([temp_initial, temp_particle_vector, temp_initial], 1, 1)
    gradient_temp_particle = compute_gradient([temp_initial, temp_particle_vector, temp_initial], 1, 1)

    ##################################### UPDATE MOISTURE ##############################################################
    change_moisture_diffusion_gas = diffusivity * density_gas * (1 - porosity) * laplacian_moisture_gas

    change_moisture_absorption_gas = - constant * density_particle * porosity * \
                                 (relative_humidity -
                                  compute_equilibrium_moisture(alpha, moisture_particle_vector, N))
    # change_moisture_absorption_gas = 0

    # print('Change moisture diffusion: ', change_moisture_diffusion_gas)
    # print('Laplacian moisture gas: ', change_moisture_diffusion_gas)
    # print('Change moisture gradient: ', -velocity * gradient_moisture_gas)
    change_m_gas = (change_moisture_diffusion_gas + change_moisture_absorption_gas) / (density_gas * (1 - porosity)) - \
                      velocity * gradient_moisture_gas

    change_m_particle = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle_vector, N))

    ##################################### UPDATE TEMP ##################################################################
    conduction_gas = conductivity_gas * (1 - porosity) * laplacian_temp_gas

    heat_of_sorption = density_particle * porosity * constant * \
                       (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle_vector, N)) * \
                       heat_capacity_vapor * (temp_gas_vector - temp_particle_vector)

    heat_transfer = -heat_transfer_coefficient * density_particle * porosity * specific_surface_area * \
                    (temp_gas_vector - temp_particle_vector)

    change_temp_gas = (conduction_gas + heat_of_sorption + heat_transfer) / \
                      (density_gas * (1 - porosity) * heat_capacity_wet_gas) - velocity * gradient_temp_gas


    conduction_particle = conductivity_particle * laplacian_temp_particle / particle_density
    heat_of_sorption = constant * (relative_humidity - compute_equilibrium_moisture(alpha, moisture_particle_vector, N)) * \
                       heat_of_vaporization

    heat_transfer = heat_transfer_coefficient * specific_surface_area * (temp_gas_vector - temp_particle_vector)

    change_temp_particle = (conduction_particle + heat_of_sorption + heat_transfer) / heat_capacity_particle

    return [change_m_gas, change_m_particle, change_temp_gas, change_temp_particle]

# moisture_gas_initial_in = np.zeros((5, 1)) + moisture_gas_initial_in
# moisture_particle_initial = np.zeros((5, 1)) + moisture_particle_initial

# moisture0 = np.column_stack((moisture_gas_initial_in, moisture_particle_initial))
system0 = np.array([moisture_gas_initial_bed, moisture_particle_initial, temp_initial, temp_initial])
moisture0 = np.array([moisture_gas_initial_bed, moisture_particle_initial])
t = np.linspace(0, 50000000, 10000)

system = odeint(moisture, system0, t, args=(
    alpha_parameter, N, moisture_diffusivity, gas_density, particle_density, porosity_powder,
    gas_velocity, heat_transfer_coefficient_initial, moisture_vapor_heat_capacity, gas_heat_capacity,
    particle_heat_capacity, specific_surface_area, conductivity_gas, conductivity_particle, bed_length))


seconds = t[-1]
hours = seconds/3600
# print('\nMoisture change: ', moisture_change)
# print('Max moisture particles: ', moisture_particle_saturated)
# plt.plot(t, moisture_change[:,1])
# plt.legend(['gas', 'particle'])
# plt.show()

fig, (ax1, ax2) = plt.subplots(2, 1)
# ax1.plot(t, system[:, 0], label = 'Gas moisture')
# ax2.plot(t, system[:, 1], label = 'Particle moisture')
ax1.plot(t, system[:, 2], label = 'Gas temp')
ax2.plot(t, system[:, 3], label = 'Particle temp')
# ax1.set_ylim(293.1, 293.2)
# ax2.set_ylim(292, 294)

# ax[2].set_ylim(0, 500)

ax1.legend(loc="best")
ax2.legend(loc="best")

plt.title(f'Time: {hours} h')
plt.show()