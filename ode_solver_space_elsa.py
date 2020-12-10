import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from vectorized_define_functions import*
from define_functions import*
from test_main import*

def moisture(moisture_matrix, t, space_step, n_space_steps, alpha, N, diffusivity, density_gas, density_particle,
             porosity, velocity, heat_transfer_coefficient, heat_capacity_vapor, heat_capacity_wet_gas,
             heat_capacity_particle, specific_surface_area, conductivity_gas, conductivity_particle, bed_length):

    moisture_gas_vector = moisture_matrix[0:n_space_steps]
    moisture_particle_vector = moisture_matrix[n_space_steps:n_space_steps*2]
    temp_gas_vector = moisture_matrix[n_space_steps*2:n_space_steps*3]
    temp_particle_vector = moisture_matrix[n_space_steps*3:n_space_steps*4]

    ##################################### UPDATE PARAMETERS ############################################################
    equilibrium_state = compute_equilibrium_moisture_vector(alpha, moisture_particle_vector, N)

    pressure_saturated = compute_p_saturated_vector(A, B, temp_gas_vector, C)

    relative_humidity = compute_relative_humidity_from_Y_vector(
        molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas_vector, pressure_saturated)

    molar_concentration_moisture = compute_molar_concentration_vector(
        relative_humidity, pressure_saturated, R_gas_constant, temp_gas_vector)

    mass_transfer_coefficient = compute_mass_transfer_coefficient(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density,
        flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture)[3]

    constant = mass_transfer_coefficient * specific_surface_area * pressure_saturated / pressure_ambient

    laplacian_moisture_gas = compute_laplacian_moisture_vector(moisture_gas_vector, space_step=space_step)
    gradient_moisture_gas = compute_gradient_moisture_vector(moisture_gas_vector, space_step=space_step)

    laplacian_temp_gas = compute_laplacian_temp_vector(temp_gas_vector, space_step=space_step)
    gradient_temp_gas = compute_gradient_temp_vector(temp_gas_vector, space_step=space_step)

    laplacian_temp_particle = compute_laplacian_temp_vector(temp_particle_vector, space_step=space_step)
    gradient_temp_particle = compute_gradient_temp_vector(temp_particle_vector, space_step=space_step)

    ##################################### UPDATE MOISTURE ##############################################################
    change_moisture_diffusion_gas = diffusivity * density_gas * (1 - porosity) * laplacian_moisture_gas

    change_moisture_absorption_gas = - constant * density_particle * porosity * \
                                 (relative_humidity - equilibrium_state)
    change_m_gas = (change_moisture_diffusion_gas + change_moisture_absorption_gas) / (density_gas * (1 - porosity)) - \
                      velocity * gradient_moisture_gas

    change_m_particle = constant * (relative_humidity - equilibrium_state)

    ##################################### UPDATE TEMP ##################################################################
    conduction_gas = conductivity_gas * (1 - porosity) * laplacian_temp_gas

    heat_of_sorption = density_particle * porosity * constant * (relative_humidity - equilibrium_state) * \
                       heat_capacity_vapor * (temp_gas_vector - temp_particle_vector)

    heat_transfer = -heat_transfer_coefficient * density_particle * porosity * specific_surface_area * \
                    (temp_gas_vector - temp_particle_vector)

    change_temp_gas = (conduction_gas + heat_of_sorption + heat_transfer) / \
                      (density_gas * (1 - porosity) * heat_capacity_wet_gas) - velocity * gradient_temp_gas


    conduction_particle = conductivity_particle * laplacian_temp_particle / particle_density
    heat_of_sorption = constant * (relative_humidity - equilibrium_state) * heat_of_vaporization

    heat_transfer = heat_transfer_coefficient * specific_surface_area * (temp_gas_vector - temp_particle_vector)

    change_temp_particle = (conduction_particle + heat_of_sorption + heat_transfer) / heat_capacity_particle

    # change_temp_gas = np.zeros(n_space_steps)
    # change_temp_particle = np.zeros(n_space_steps)
    # test = (np.concatenate([change_m_gas, change_m_particle, change_temp_gas, change_temp_particle]))
    return np.concatenate([change_m_gas, change_m_particle, change_temp_gas, change_temp_particle])

space_steps = 5
length = np. linspace(0, bed_length, space_steps)
space_step = bed_length/space_steps
t = np.linspace(0, 1000000, 500)

moisture_gas_initial = np.zeros(space_steps) + moisture_gas_initial_bed
moisture_particle_initial = np.zeros(space_steps) + moisture_particle_initial
# all_moisture = (np.concatenate([moisture_gas_initial, moisture_particle_initial]))

temp_gas_initial = np.zeros(space_steps) + temp_initial
temp_particle_initial = np.zeros(space_steps) + temp_initial

# system0 = np.array([moisture_gas_initial, moisture_particle_initial, temp_gas_initial, temp_particle_initial])
system0 = np.concatenate([moisture_gas_initial, moisture_particle_initial, temp_gas_initial, temp_particle_initial])

system = odeint(moisture, system0, t, args=(
    space_step, space_steps, alpha_parameter, N, moisture_diffusivity, gas_density, particle_density, porosity_powder,
    gas_velocity, heat_transfer_coefficient_initial, moisture_vapor_heat_capacity, gas_heat_capacity,
    particle_heat_capacity, specific_surface_area, conductivity_gas, conductivity_particle, bed_length))

# print(system)
seconds = t[-1]
hours = seconds/3600

fig, ax = plt.subplots(5, 4)
moisture_color = 'navy'
temp_color = 'orangered'
# ax1.plot(t, system[:, 0], label = 'Gas moisture')
# ax2.plot(t, system[:, 1], label = 'Particle moisture')
ax[0, 0].plot(t, system[:, 0], label = 'Gas moisture section 1', c = moisture_color)
ax[1, 0].plot(t, system[:, 1], label = 'Gas moisture section 2', c = moisture_color)
ax[2, 0].plot(t, system[:, 2], label = 'Gas moisture section 3', c = moisture_color)
ax[3, 0].plot(t, system[:, 3], label = 'Gas moisture section 4', c = moisture_color)
ax[4, 0].plot(t, system[:, 4], label = 'Gas moisture section 5', c = moisture_color)

ax[0, 1].plot(t, system[:, 5], label = 'Particle moisture section 1', c = moisture_color)
ax[1, 1].plot(t, system[:, 6], label = 'Particle moisture section 2', c = moisture_color)
ax[2, 1].plot(t, system[:, 7], label = 'Particle moisture section 3', c = moisture_color)
ax[3, 1].plot(t, system[:, 8], label = 'Particle moisture section 4', c = moisture_color)
ax[4, 1].plot(t, system[:, 9], label = 'Particle moisture section 5', c = moisture_color)

ax[0, 2].plot(t, system[:, 10], label = 'Gas temp section 1', c = temp_color)
ax[1, 2].plot(t, system[:, 11], label = 'Gas temp section 2', c = temp_color)
ax[2, 2].plot(t, system[:, 12], label = 'Gas temp section 3', c = temp_color)
ax[3, 2].plot(t, system[:, 13], label = 'Gas temp section 4', c = temp_color)
ax[4, 2].plot(t, system[:, 14], label = 'Gas temp section 5', c = temp_color)

ax[0, 3].plot(t, system[:, 15], label = 'Particle temp section 1', c = temp_color)
ax[1, 3].plot(t, system[:, 16], label = 'Particle temp section 2', c = temp_color)
ax[2, 3].plot(t, system[:, 17], label = 'Particle temp section 3', c = temp_color)
ax[3, 3].plot(t, system[:, 18], label = 'Particle temp section 4', c = temp_color)
ax[4, 3].plot(t, system[:, 19], label = 'Particle temp section 5', c = temp_color)

# ax[0, 2].set_ylim(293, 294)

# ax[1].legend(loc="best")
# ax[2].legend(loc="best")

# plt.title(f'Time: {hours} h')
plt.show()