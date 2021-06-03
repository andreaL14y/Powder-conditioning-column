import numpy as np
import scipy.optimize
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

########################################################################################################################
def compute_k_vector_old(temperature_vector, humidity_vector, glass_transition_temp_vector):
    height_steps, length_steps = temperature_vector.shape
    parameters_vector = np.zeros((height_steps, length_steps, 3))
    for h in range(height_steps):
        for l in range(length_steps):
            temp = temperature_vector[h, l]
            rh = humidity_vector[h, l]
            if temp < glass_transition_temp_vector[h, l]:
                starting_point, k_parameter, rest = 0, 0, 0
            elif temp < 40 + kelvin:
                if rh < 55:
                    starting_point, k_parameter, rest = crystallization_parameters[0]
                else:
                    starting_point, k_parameter, rest = crystallization_parameters[1]
            else:
                if rh < 55:
                    starting_point, k_parameter, rest = crystallization_parameters[2]
                else:
                    starting_point, k_parameter, rest = crystallization_parameters[3]
            k_parameter /= 60 * 60 #* 7 * 24                                     # Unit transformed from 1/h to 1/s

            parameters_vector[h, l, :] = starting_point, k_parameter, rest
    return parameters_vector
########################################################################################################################

def crystallization_speed_curves(data, starting_point, k_parameter, rest):
    return starting_point * np.exp(-k_parameter * data) + rest

def glass_transition_curve(weight_fraction_lactose, param_1, param_2):
    return param_1 + 1/(np.exp(-param_2 * weight_fraction_lactose))


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
crystallization_parameters = np.zeros((4, 3))

colors = ['green', 'black', 'red', 'hotpink']
for set in range(4):
    initial_parameters = (4.6, .005, .5)  # start with values near those we expect
    enthalpy = data[set]
    crystallization_parameters[set], cv = scipy.optimize.curve_fit(crystallization_speed_curves, time, enthalpy, initial_parameters)
    if crystallization_parameters[set][2] < 0:
        crystallization_parameters[set][2] = 0

    crystallization_param_1, crystallization_param_2, crystallization_param_3 = crystallization_parameters[set]
    y2 = crystallization_speed_curves(x2, crystallization_param_1, crystallization_param_2, crystallization_param_3)
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

weight_fraction_lactose = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
temperature = np.array([-135, -127, -112, -90, -50, 100])

# perform the fit
p0 = (-135, -1) # start with values near those we expect
params, cv = scipy.optimize.curve_fit(glass_transition_curve, weight_fraction_lactose, temperature, p0)
glass_transition_param_1, glass_transition_param_2 = params

