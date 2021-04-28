import numpy as np
import scipy.optimize
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def monoExp(data, starting_point, k_parameter, rest):
    return starting_point * np.exp(-k_parameter * data) + rest


########################################### FROM FIG 4 #################################################################
time = np.array([0, 3, 6, 9])
multiply = 1 #7 * 24 #* 60 * 60
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
    initial_parameters = (4.6, .05, .5)  # start with values near those we expect
    enthalpy = data[set]
    parameters[set], cv = scipy.optimize.curve_fit(monoExp, time, enthalpy, initial_parameters)
    if parameters[set][2] < 0:
        parameters[set][2] = 0

    starting_point, k_parameter, rest = parameters[set]
    y2 = monoExp(x2, starting_point, k_parameter, rest)
    plt.plot(x2, y2, '-', color=colors[set], label="fitted")
    plt.plot(time, enthalpy, 'X', color=colors[set], label="Original data")


plt.hlines(0, 0, 20 * multiply, linestyles='--', color='gray')

plt.title("Extrapolated Exponential Curve")
print(f'Parameters are start: {starting_point:.2f}, k: {k_parameter:.2f} and rest: {rest:.2f}')
plt.legend()
plt.show()
