import numpy as np
import scipy.optimize
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def monoExp(data, starting_point, k_parameter, rest):
    return starting_point * np.exp(-k_parameter * data) + rest


########################################### FROM FIG 4 #################################################################
time = np.array([0, 3, 6, 9])
multiply = 1 #* 60 * 2
time *= multiply

enthalpy_2030 = np.array([4.6, 2.1, 1.2, 1.1])
enthalpy_2080 = np.array([4.6, 2,   1.2, 1.0])
enthalpy_6030 = np.array([4.6, 1.5, 0.9, 0.0])
enthalpy_6080 = np.array([4.6, 1,   0.1, 0.0])
data = np.array([enthalpy_2030, enthalpy_2080, enthalpy_6030, enthalpy_6080])

############################################# FITTING ##################################################################
x2 = np.linspace(0, 20 * multiply, 200)
parameters = np.zeros((4, 3))

colors = ['green', 'black', 'red', 'hotpink']
for set in range(4):
    initial_parameters = (4.6, .005, .5)  # start with values near those we expect
    enthalpy = data[set]
    parameters[set], cv = scipy.optimize.curve_fit(monoExp, time, enthalpy, initial_parameters)
    if parameters[set][2] < 0:
        parameters[set][2] = 0

    starting_point, k_parameter, rest = parameters[set]
    y2 = monoExp(x2, starting_point, k_parameter, rest)
#     plt.plot(x2, y2, '-', color=colors[set], label="fitted")
#     plt.plot(time, enthalpy, 'X', color=colors[set], label="Original data")
#
#
# plt.hlines(0, 0, 20 * multiply, linestyles='--', color='gray')
#
# plt.title("Extrapolated Exponential Curve")
# print(f'\nParameters are: \nT 20, RH 30: \nStarting point: {parameters[0, 0]:.2f}, k: {parameters[0, 1]:.2f} and rest: {parameters[0, 2]:.2f}'
#       f'\nT 20, RH 80: \nStarting point: {parameters[1, 0]:.2f}, k: {parameters[1, 1]:.2f} and rest: {parameters[1, 2]:.2f}'
#       f'\nT 60, RH 30: \nStarting point: {parameters[2, 0]:.2f}, k: {parameters[2, 1]:.2f} and rest: {parameters[2, 2]:.2f}'
#       f'\nT 60, RH 60: \nStarting point: {parameters[3, 0]:.2f}, k: {parameters[3, 1]:.2f} and rest: {parameters[3, 2]:.2f}\n')
# plt.legend()
# plt.show()
