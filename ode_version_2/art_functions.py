import numpy as np
from input_parameters import*
from glass_transition_curve import compute_glass_temp_mix, glass_temp_water_1, glass_temp_lactose
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# Constants for lactose
m1 = 21000
c1 = -16.7
k1 = 4.7

m2 = -95
c2 = c1
k2 = k1

def compute_entropy_of_activation_S(moisture):
    entropy_of_activation_S = m2 * c2 * k2 * moisture/((1 - k2 * moisture) * (1 - k2 * moisture + c2 * k2 * moisture))
    return entropy_of_activation_S


def compute_enthalpy_of_activation_H(moisture):
    enthalpy_of_activation_H = m1 * c1 * k1 * moisture/((1 - k1 * moisture) * (1 - k1 * moisture + c1 * k1 * moisture))
    return enthalpy_of_activation_H


def compute_crystal_growth_rate(moisture, temp):
    t_g = compute_glass_temp_mix(1-moisture, glass_temp_lactose, glass_temp_water_1)
    entropy_S = compute_entropy_of_activation_S(moisture)
    enthalpy_H = compute_enthalpy_of_activation_H(moisture)

    r_hand_1 = k_boltzmann / h_planck
    r_hand_2 = - enthalpy_H / (r_gas_constant * temp)
    r_hand_3 = entropy_S / r_gas_constant

    ln_stuff = r_hand_2 + r_hand_3

    k = np.exp(ln_stuff) * r_hand_1 * temp             # 0.07 g/g m, 0.93 g/g lactose
    reaction_rate = moisture * k
    # import sys
    # import warnings
    # warnings.filterwarnings("error")
    # try:
    #     k = np.where(temp < t_g, 0, np.exp(r_hand_1 + r_hand_2)/temp)
    #     # k = np.where(temp < t_g, 0, np.exp(r_hand_1 + r_hand_2)/temp)
    # except:
    #     print('\nHey I have issues!')
    #     print( temp < t_g )
    #     sys.exit()

    # if temp[temp > t_g] != []:
    #     print(temp[temp > t_g])

    # if temp > t_g:
    #     entropy_S = compute_entropy_of_activation_S(moisture)
    #     enthalpy_H = compute_enthalpy_of_activation_H(moisture)
    #     r_hand_1 = np.log(k_boltzmann/h_planck) + entropy_S/r_gas_constant
    #     r_hand_2 = -(enthalpy_H/r_gas_constant)/temp
    #     k = np.exp(r_hand_1 + r_hand_2)/temp
    # else:
    #     k = 0
    return k, reaction_rate


def compute_moisture_content(dt, k, moisture_content):
    moisture_content_change = moisture_content * k
    moisture_content += moisture_content_change * dt
    return moisture_content


# k_test, reaction_rate = compute_crystal_growth_rate(0.088, 15 + kelvin)
# print(f'k is: {k_test:.4f}')
# k_test, reaction_rate = compute_crystal_growth_rate(0.07, 25 + kelvin)
# print(f'k is: {k_test:.4f}')
# k_test, reaction_rate = compute_crystal_growth_rate(0.036, 40 + kelvin)
# print(f'k is: {k_test:.4f}')
# # print(f'Reaction rate is: {reaction_rate:.4f}/s')
# k_test = np.exp(-9) * (60 + kelvin)
# print('ln(k/T) = - 9 gives k =', k_test)
#
# should_be = np.log( 0.0007/(60 + kelvin) )
# print('ln(k/T) should be:', should_be)
#
# temps = np.array([15, 25, 40], dtype='float64')
# initial_moistures = np.array([0.088, 0.07, 0.036])
# temps += kelvin
# max_times = np.array([1000, 600, 350])
# max_time = 1400         # seconds
# colors = ['blue', 'green', 'red']
# for setup in range(3):
#     # squared = 0
#     max_time = max_times[setup]
#     moistures = np.zeros(max_time)
#     moistures[0] = initial_moistures[setup]
#     temp = temps[setup]
#     for time in range(max_time-1):
#         k = -compute_crystal_growth_rate(moistures[time], temp)[0]
#         moistures[time + 1] = compute_moisture_content(1, k, moistures[time])
#     min_m = min(moistures)
#     max_m = max(moistures)
#     moistures = (moistures - min_m)/(max_m - min_m)
#
#     times = np.arange(0, max_time)
#     plt.plot(times, moistures, label=f'Temp {temp-kelvin:.2f}', color=colors[setup])
# plt.xlabel('time, min')
# plt.ylabel('Moisture content')
# plt.legend()

# # PLOTTING
initial_moisture = 0.0125
# fig, axs = plt.subplots(4, figsize=(15,15))
# fig.suptitle(f'At temp {temp:.2f} K = {temp-kelvin:.2f} deg C')
#
# m = np.linspace(initial_moisture, 0.12, 150)
# time = np.linspace(0, 1000, 1000)
# # m = np.linspace(0.01, 0.12, 150)
#
# k, r_r = compute_crystal_growth_rate(m, temp)
# m_max_index = np.argmax(k)
# h = compute_enthalpy_of_activation_H(m)
# s = compute_entropy_of_activation_S(m)
#
# axs[0].plot(m, h, label='Enthalpy H')
# axs[1].plot(m, s, label='Entropy S')
# axs[2].plot(m, k, label='k')
# axs[3].plot(m, r_r, label='Reaction rate')
#
# axs[0].set(ylim=[0, 450000])
# axs[1].set(ylim=[-1000, 0], label='Entropy S')
# axs[2].vlines(m[m_max_index], 0, 0.07, color='darkred', linestyles='dashed')
#
# axs[3].set(ylabel='dX/dt = k*X')
# print(f'm_max is: {m[m_max_index]:.2f}')
# axs[3].vlines(m[m_max_index], 0, 0.0035, color='darkred', linestyles='dashed')
#
# axs[0].legend()
# axs[1].legend()
# axs[2].legend()
# axs[3].legend()
#
# axs[0].grid()
# axs[1].grid()
# axs[2].grid()
# axs[3].grid()
#
# plt.xlabel('Moisture content X')


# plt.show()


