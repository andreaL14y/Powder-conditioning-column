import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from vectorized_define_functions import *
from test_main import *

def moisture(moisture_matrix, t, space_step, n_space_steps, alpha, N, diffusivity, density_gas, density_particle,
             porosity, velocity, heat_transfer_coefficient, heat_capacity_vapor, heat_capacity_wet_gas,
             heat_capacity_particle, specific_surface_area, conductivity_gas, conductivity_particle, bed_length):
    moisture_gas_vector = moisture_matrix[0:n_space_steps]
    moisture_particle_vector = moisture_matrix[n_space_steps:n_space_steps * 2]
    temp_gas_vector = moisture_matrix[n_space_steps * 2:n_space_steps * 3]
    temp_particle_vector = moisture_matrix[n_space_steps * 3:n_space_steps * 4]

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


space_steps = 6
length = np.linspace(0, bed_length, space_steps)
space_step = bed_length / space_steps
t = np.linspace(0, 500000, 500)

moisture_gas_initial_all = np.zeros(space_steps) + moisture_gas_initial_bed
moisture_particle_initial_all = np.zeros(space_steps) + moisture_particle_initial

temp_gas_initial = np.zeros(space_steps) + temp_initial
temp_particle_initial = np.zeros(space_steps) + temp_initial

system0 = np.concatenate(
    [moisture_gas_initial_all, moisture_particle_initial_all, temp_gas_initial, temp_particle_initial])

system = odeint(moisture, system0, t, args=(
    space_step, space_steps, alpha_parameter, N, moisture_diffusivity, gas_density, particle_density, porosity_powder,
    gas_velocity, heat_transfer_coefficient_initial, moisture_vapor_heat_capacity, gas_heat_capacity,
    particle_heat_capacity, specific_surface_area, conductivity_gas, conductivity_particle, bed_length))

######################################### PLOTTING #####################################################################
# Convert to easier-to-read units
seconds = t[-1]
hours = seconds / 3600
t /= 3600
system[:, space_steps * 2:space_steps * 3] -= kelvin
max_temp_gas = np.max(system[:, space_steps * 2:space_steps * 3])

system[:, space_steps * 3:space_steps * 4] -= kelvin
max_temp_particle = np.max(system[:, (space_steps * 3):(space_steps * 4)])

print(f'Time computed is: {int(hours)} hours')
print('Max temperature in gas is: {:.4f} degrees Celcius'.format(max_temp_gas))
print('Max temperature in particles is: {:.4f} degrees Celcius'.format(max_temp_particle))

fig, ax = plt.subplots(space_steps, 4)
fig.suptitle(f'Moisture & temperature in sections ove time. Total time: {int(hours)} hours.', fontsize=14)
ax[0, 0].set_title('Moisture Gas')
ax[0, 1].set_title('Moisture Particle')
ax[0, 2].set_title('Temp Gas')
ax[0, 3].set_title('Temp Particle')

# Set styles
moisture_color = 'navy'
temp_color = 'orangered'
initial_line = 'dashed'
initial_color = 'gray'
saturated_color = 'lightcoral'

counter = 0
epsilon = 0.001
for feature in range(4):
    if feature == 0 or feature == 1:
        current_color = moisture_color
    else:
        # print('yo')
        current_color = temp_color

    for step in range(space_steps):
        ax[step, feature].plot(t, system[:, counter], c=current_color)

        if feature == 0:
            patch = mpatches.Patch(color=current_color, label=f'Y {step}')
            ax[step, feature].set_ylim(moisture_gas_initial_bed - epsilon, moisture_gas_initial_in + epsilon)
            ax[step, feature].hlines(moisture_gas_initial_bed, 0, t[-1], colors=initial_color, linestyles=initial_line)
            ax[step, feature].hlines(moisture_gas_initial_in, 0, t[-1], colors=saturated_color, linestyles=initial_line)

        elif feature == 1:
            patch = mpatches.Patch(color=current_color, label=f'X {step}')
            ax[step, feature].set_ylim(moisture_particle_initial - epsilon, moisture_particle_saturated + epsilon)
            ax[step, feature].hlines(moisture_particle_initial, 0, t[-1], colors=initial_color, linestyles=initial_line)
            ax[step, feature].hlines(moisture_particle_saturated, 0, t[-1], colors=saturated_color, linestyles=initial_line)

        elif feature == 2:
            patch = mpatches.Patch(color=current_color, label=f'TG {step}')
            ax[step, feature].set_ylim(temp_initial - (kelvin + epsilon), max_temp_gas + 1)
            ax[step, feature].hlines(temp_initial - kelvin, 0, t[-1], colors=initial_color, linestyles=initial_line)
            ax[step, feature].hlines(max_temp_gas, 0, t[-1], colors=saturated_color, linestyles=initial_line)
        else:
            patch = mpatches.Patch(color=current_color, label=f'TP {step}')
            ax[step, feature].set_ylim(temp_initial - (kelvin + epsilon), max_temp_particle + 1)
            ax[step, feature].hlines(temp_initial - kelvin, 0, t[-1], colors=initial_color, linestyles=initial_line)
            ax[step, feature].hlines(max_temp_particle, 0, t[-1], colors=saturated_color, linestyles=initial_line)

        ax[step, feature].legend(handles=[patch], loc="upper left")
        ax[step, feature].grid()
        # ax[step, feature].xaxis.grid()
        counter += 1

plt.show()
