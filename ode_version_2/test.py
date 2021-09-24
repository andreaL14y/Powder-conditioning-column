from initial_conditions import *
import time
import numpy as np
kelvin = 273.15

a = np.array([1, 2, 3])
b = np.array([0, 1, 2])
print(a*b)

def f(x):
    A = np.array([1, 2, 3])
    # A + x
    diff = np.mean(A + x, axis=0)
    return diff



# T = np.linspace(18, 40, 50)
# wa = np.linspace(relative_humidity_bed_initial, relative_humidity_gas_inlet, 50)
#
#
# H1 = compute_H_from_aw_temp(wa, temp_initial_celsius)
# densities1 = compute_density_air(temp_initial_celsius, H1)
#
# H2 = compute_H_from_aw_temp(wa, temp_initial_celsius+10)
# densities2 = compute_density_air(temp_initial_celsius+10, H2)
#
# # plt.plot(H1, densities2, label=f'{temp_initial_celsius}')
# # plt.plot(H2, densities2, label=f'{temp_initial_celsius+10}')
# plt.plot(wa, H1, label=f'{temp_initial_celsius}')
# plt.plot(wa, H2, label=f'{temp_initial_celsius+10}')
# plt.ylabel('H')
# plt.xlabel('wa')
# plt.legend()
# plt.show()

# def compute_H_and_M_iteratively(water_activity, total_water, temp_celsius, amount_am):
#     pressure_water = compute_pressure_water(temp_celsius)
#     m_void = 18 * water_activity * pressure_water/(29 * (pressure_ambient - water_activity * pressure_water))
#
#     m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
#     m_am = compute_GAB_equilibrium_moisture_am(water_activity)
#     m_particle_tot = m_am * amount_am + m_cryst * (1 - amount_am)
#
#     W = m_void * porosity_powder * density_gas + m_particle_tot * (1 - porosity_powder) * density_particle
#     diff = (W - total_water)**2
#     error = diff.mean(axis=0)
#     return error
#
# cons = [{"type": "ineq", 'fun': lambda x: x}, {"type": "ineq", 'fun': lambda x: 1 - x}] # water_activity, total_water, temp_celsius, amount_am
#
# aw_sur_vector = relative_humidity_gas_inlet
# aw_sur_vector = 0.8
# total_water = system_conc_initial
# temp_celsius_prev = temp_initial_celsius #- 30
# amount_am = amorphous_material_initial
#
# comp_time = 0
# table_time = 0
# for i in range(2):
#     am_amorph = np.zeros(2) + amorphous_material_initial
#     total_water_avg = np.zeros(2) + total_water
#
#     # table_start = time.time()
#     # table_H = compute_H_and_M_gas_tables(am_amorph, total_water_avg)
#     # table_time += time.time() - table_start
#
#     start = time.time()
#     optimized_aw = scipy.optimize.minimize(compute_H_and_M_iteratively, aw_sur_vector, args=(total_water, temp_celsius_prev, amount_am),
#                                            constraints=cons, method='SLSQP', options={'ftol': 1e-6})
#     print(optimized_aw.x)
#     comp_time += time.time() - start
#
# print('Optimization time:', comp_time)
# print('Table time:', table_time)



