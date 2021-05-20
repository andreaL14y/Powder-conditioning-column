import numpy as np
import scipy.optimize
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

############################################### DATA ###################################################################
kelvin = 273.15
weight_fraction_lactose = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
temperature = np.array([-135, -127, -112, -90, -50, 101]) + kelvin
glass_temp_water_1 = 136            # Kelvin
glass_temp_water_2 = 165
glass_temp_lactose = 101 + kelvin   # Kelvin
room_temp = 20 + kelvin

############################################ FUNCTIONS #################################################################
def monoExp(x, m, t):
    return m + 1 / (np.exp(-t * x))


def glass_temp_mix(w1, t_g1, t_g2):
  w2 = 1 - w1
  k_param = 6.7                                                                 # Glass Transitions and Crystallization in Milk Powders
  glass_temp_mix = (w1 * t_g1 + w2 * t_g2 * k_param)/(w1 + w2 * k_param)
  return(glass_temp_mix)

############################################# FITTING ##################################################################
initial_params = (-135, -1) # start with values near those we expect
params, cv = scipy.optimize.curve_fit(monoExp, weight_fraction_lactose, temperature, initial_params)
m, t = params

######################################## GORDON & TAYLOR ###############################################################
weight_fractions_lactose = np.linspace(0, 1, 100)
glass_temps_1 = glass_temp_mix(weight_fractions_lactose, glass_temp_lactose, glass_temp_water_1)
glass_temps_2 = glass_temp_mix(weight_fractions_lactose, glass_temp_lactose, glass_temp_water_2)


############################################# PLOT #####################################################################
# fig = plt.figure(figsize=(8, 12), dpi= 100, facecolor='w', edgecolor='k')
# plt.plot(weight_fractions_lactose, glass_temps_1, color='crimson', label='T_G water 136 K')
# plt.plot(weight_fractions_lactose, glass_temps_2, color='green', label='T_G water 165 K')
# plt.plot(weight_fractions_lactose, monoExp(weight_fractions_lactose, m, t), '--', color='orange', label="Fitted to data as before")
#
# plt.hlines(room_temp, 0, 1, linestyles='dashed', label='Room temp')
# ax = fig.add_subplot(1, 1, 1)
#
# # Major ticks every 20, minor ticks every 5
# major_ticks = np.linspace(0, 1, 11)
# minor_ticks = np.linspace(0, 1, 101)
# major_ticks_y = np.arange(130, 400, 30)
# minor_ticks_y = np.arange(130, 400, 10)
#
# ax.set_xticks(major_ticks)
# ax.set_xticks(minor_ticks, minor=True)
# ax.set_yticks(major_ticks_y)
# ax.set_yticks(minor_ticks_y, minor=True)
#
# # Or if you want different settings for the grids:
# ax.grid(which='minor', alpha=0.2)
# ax.grid(which='major', alpha=0.5)
#
# #plt.grid()
# plt.xlim(0.8, 1)
# plt.ylim(220, 400)
# plt.ylabel('Temp K')
# plt.xlabel('Weight fraction lactose')
# plt.legend()
# plt.show()


