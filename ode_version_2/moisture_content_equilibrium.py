import numpy as np
import scipy.optimize
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
# from input_parameters import *

N_cryst = 0.9                                                                         # Parameter, assume 1
# alpha_parameter_am = 25                                                            # parameter, 10 < alpha < 100
alpha_parameter_cryst = 1733                                                          # parameter, 10 < alpha < 100

######################################## FUNCTION & DATA ###############################################################
def compute_air_equilibrium(moisture_particle_vector, alpha, N):
    equilibrium_air = (1 - np.exp(-alpha * moisture_particle_vector ** N))
    return equilibrium_air

# def compute_air_equilibrium_am(moisture_particle_vector):
#     eq_m_am = (1 - np.exp(-alpha_parameter_am * moisture_particle_vector ** N_am))
#     return eq_m_am


def compute_GAB_equilibrium_moisture_am(relative_humidity_vector):
    M0 = 0.0488
    c = 3.23
    f = 1.16
    equilibrium_moisture_particle_vector = M0 * c * f * relative_humidity_vector / \
                                           ((1 - f * relative_humidity_vector) * (
                                                       1 - (1 - c) * f * relative_humidity_vector))

    return equilibrium_moisture_particle_vector


def compute_GAB_equilibrium_moisture_cryst(relative_humidity_vector):
    M0 = 1.68 * 10 ** (-4)
    c = 8.8
    f = 0.878
    h = 30
    H = 1 + (1 - f) * (f * relative_humidity_vector) ** h / (f * (1 - relative_humidity_vector))
    H_prime = 1 + (H - 1) * (1 - f * relative_humidity_vector) / (H * (1 - relative_humidity_vector)) * (
                h + (1 - h) * relative_humidity_vector)

    equilibrium_moisture_particle_vector = M0 * c * f * relative_humidity_vector * H * H_prime / \
                                           ((1 - f * relative_humidity_vector) * (
                                                       1 + (c * H - 1) * f * relative_humidity_vector))

    return equilibrium_moisture_particle_vector


def equilibrium_curve_fit(X, alpha, N):
    return 1 - np.exp(-alpha * X**N)

relative_humidities = np.linspace(0, 0.85, 10000)
moistures_am = compute_GAB_equilibrium_moisture_am(relative_humidities)

moistures_cryst = compute_GAB_equilibrium_moisture_cryst(relative_humidities)
eq_air_cryst = compute_air_equilibrium(moistures_cryst, alpha_parameter_cryst, N_cryst)

############################################ FITTING ###################################################################
start_params = (3, 0.7) # start with values near those we expect
params, cv = scipy.optimize.curve_fit(equilibrium_curve_fit, moistures_am, relative_humidities, start_params)
alpha_parameter_am, N_am = params

# print('alpha, N:', params)

eq_air_am = equilibrium_curve_fit(moistures_am, alpha_parameter_am, N_am)

############################################ PLOTTING ##################################################################
# fig, axs = plt.subplots(2, figsize=(15, 15))
# fig.suptitle('Equilibrium GAB vs fitted curves')
# axs[0].plot(relative_humidities, moistures_cryst, label='Moisture content')
# axs[0].plot(eq_air_cryst, moistures_cryst, label='Fitted curve')
# axs[0].legend()
# axs[0].grid()
# axs[0].set_title('Crystalline lactose')
#
# axs[1].plot(relative_humidities, moistures_am, label='Moisture content')
# axs[1].plot(eq_air_am, moistures_am, label='Fitted curve')
# axs[1].legend()
# axs[1].grid()
# axs[1].set_ylim(0, 0.9)
# axs[1].set_title('Amorhpous lactose')
#
# plt.xlabel('RH')
# plt.ylabel('Moisture, g/g lactose')
#
# plt.show()
