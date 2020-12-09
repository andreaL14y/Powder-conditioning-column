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

def moisture(moisture_matrix, t, temp_gas, alpha, N, diffusivity, density_gas, density_particle, porosity, velocity):
    # print(moisture_matrix)
    moisture_gas_vector = moisture_matrix[0]
    moisture_particle_vector = moisture_matrix[1]

    ##################################### UPDATE PARAMETERS ############################################################
    pressure_saturated = compute_p_saturated_vector(A, B, temp_gas, C)

    relative_humidity = compute_relative_humidity_from_Y(
        molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas_vector, pressure_saturated)
    # print('RH :', relative_humidity)

    molar_concentration_moisture = compute_molar_concentration_vector(
        relative_humidity, pressure_saturated, R_gas_constant, temp_gas)

    mass_transfer_coefficient = compute_mass_transfer_coefficient(
        moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density,
        flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture)[3]

    constant = mass_transfer_coefficient * specific_surface_area * pressure_saturated / pressure_ambient

    laplacian_moisture_gas = compute_laplacian([moisture_gas_initial_in, moisture_gas_vector, moisture_gas_initial_bed], 1, 1)
    gradient_moisture_gas = compute_gradient([moisture_gas_initial_in, moisture_gas_vector, moisture_gas_initial_bed], 1, 1)


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

    print(change_m_gas, change_m_particle)
    return [change_m_gas, change_m_particle]

# moisture_gas_initial_in = np.zeros((5, 1)) + moisture_gas_initial_in
# moisture_particle_initial = np.zeros((5, 1)) + moisture_particle_initial

# moisture0 = np.column_stack((moisture_gas_initial_in, moisture_particle_initial))
moisture0 = np.array([moisture_gas_initial_bed, moisture_particle_initial])
t = np.linspace(0, 5000, 100)

moisture_change = odeint(
    moisture, moisture0, t, args=(temp_initial, alpha_parameter, N, moisture_diffusivity,
             gas_density, particle_density, porosity_powder, gas_velocity))


# print('\nMoisture change: ', moisture_change)
# print('Max moisture particles: ', moisture_particle_saturated)
# plt.plot(t, moisture_change[:,1])
# plt.legend(['gas', 'particle'])
# plt.show()

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(t, moisture_change[:, 0], label = 'Gas moisture')
ax2.plot(t, moisture_change[:, 1], label = 'Particle moisture')

# ax[2].set_ylim(0, 500)

ax1.legend(loc="best")
ax2.legend(loc="best")

plt.title(f'Time: {t[-1]} s')
plt.show()