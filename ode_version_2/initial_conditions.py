from functions import *

def compute_density_air(temp_celsius, H):
    top = kelvin
    bottom = 22.4 * (kelvin + temp_celsius) * (1/29 + H/18)
    density = top/bottom
    return density

p_saturated                     = compute_p_saturated_vector(temp_initial)

# Initial surroundings, RH 58
molar_conc_sur                  = compute_molar_concentration_vector(relative_humidity_gas_inlet, p_saturated, temp_initial)
vapor_pressure_sur              = compute_partial_pressure_moisture_vector(molar_conc_sur, temp_initial)
water_activity_sur              = vapor_pressure_sur / p_saturated
m_gas_sur                       = 18 * vapor_pressure_sur/(29 * (pressure_ambient - vapor_pressure_sur))
m_gas_conc_sat                  = m_gas_sur * porosity_powder * density_gas

# Initial system, RH 20
# Gas
molar_concentration_void_initial= compute_molar_concentration_vector(relative_humidity_bed_initial, p_saturated, temp_initial)
vapor_pressure_void_initial     = compute_partial_pressure_moisture_vector(molar_concentration_void_initial, temp_initial)
water_activity_void_initial     = vapor_pressure_void_initial / p_saturated
m_void_initial                   = 18 * vapor_pressure_void_initial / (29 * (pressure_ambient - vapor_pressure_void_initial))
m_void_conc_initial              = m_void_initial * porosity_powder * density_gas

# Particle
m_particle_cryst_initial          = compute_GAB_equilibrium_moisture_cryst(water_activity_void_initial)
m_particle_cryst_conc_initial     = m_particle_cryst_initial * (1 - porosity_powder) * density_particle
m_particle_am_initial             = compute_GAB_equilibrium_moisture_am(water_activity_void_initial)
m_particle_am_conc_initial        = m_particle_am_initial * (1 - porosity_powder) * density_particle
m_particle_tot_initial          = m_particle_cryst_initial * (1 - amorphous_material_initial) + m_particle_am_initial * amorphous_material_initial
m_particle_tot_conc_initial     = m_particle_cryst_conc_initial * (1 - amorphous_material_initial) + m_particle_am_conc_initial * amorphous_material_initial

system_conc_initial             = m_particle_cryst_conc_initial * (1 - amorphous_material_initial) + m_particle_am_conc_initial * amorphous_material_initial + m_void_conc_initial

# Energy
enthalpy_vapor_initial          = heat_capacity_vapor * temp_initial_celsius + heat_of_evaporation_water
enthalpy_gas_initial            = heat_capacity_air * temp_initial_celsius + m_void_initial * enthalpy_vapor_initial
enthalpy_surrounding_initial    = heat_capacity_air * temp_initial_celsius + m_gas_sur * enthalpy_vapor_initial
enthalpy_particle_initial       = heat_capacity_particle * temp_initial_celsius + m_particle_tot_initial * heat_capacity_water * temp_initial_celsius
enthalpy_powder_initial         = enthalpy_gas_initial * porosity_powder * density_gas + enthalpy_particle_initial * (1 - porosity_powder) * density_particle

hot = 100
enthalpy_vapor_hot          = heat_capacity_vapor * hot + heat_of_evaporation_water
enthalpy_gas_hot            = heat_capacity_air * hot + m_void_initial * enthalpy_vapor_hot
enthalpy_particle_hot       = heat_capacity_particle * hot + m_particle_tot_initial * heat_capacity_water * hot
enthalpy_powder_hot         = enthalpy_gas_hot * porosity_powder * density_gas + enthalpy_particle_hot * (1 - porosity_powder) * density_particle

# Saturated system
m_particle_cryst_sat            = compute_GAB_equilibrium_moisture_cryst(relative_humidity_gas_inlet)
m_particle_cryst_conc_sat       = m_particle_cryst_sat * (1 - porosity_powder) * density_particle
m_particle_am_sat               = compute_GAB_equilibrium_moisture_am(relative_humidity_gas_inlet)
m_particle_am_conc_sat          = m_particle_am_sat * (1 - porosity_powder) * density_particle
m_particle_tot_conc_sat         = m_particle_cryst_conc_sat * (1-amorphous_material_initial) + m_particle_am_conc_sat * amorphous_material_initial
m_particle_tot_sat              = m_particle_cryst_sat * (1-amorphous_material_initial) + m_particle_am_sat * amorphous_material_initial
system_conc_sat                 = m_particle_cryst_conc_sat * (1-amorphous_material_initial) + m_particle_am_conc_sat * amorphous_material_initial + m_gas_conc_sat
m_gas_sat                       = m_gas_conc_sat/(porosity_powder * density_gas)

density_sat = compute_density_air(temp_initial_celsius, m_gas_sat)
gas_density_initial = compute_density_air(temp_initial_celsius, m_void_initial)
gas_density_sur = compute_density_air(temp_initial_celsius, m_gas_sur)

gas_fraction = m_void_conc_initial / (m_void_conc_initial + m_particle_tot_conc_initial)
gas_fraction_initial = gas_fraction
gas_fraction_2 = m_gas_conc_sat/(m_gas_conc_sat + m_particle_tot_conc_sat)
gas_fraction_min = gas_fraction_2

gas_fraction_3 = m_gas_conc_sat/(m_gas_conc_sat + m_particle_cryst_conc_sat)
gas_fraction_max_tot = gas_fraction_3


####################################### CREATE TABLES ##################################################################
# H = np.linspace(m_void_initial, .01999, 10000)
# print(m_gas_sat)
# print(m_void_initial)
# n_H = np.size(H)
# a_w = compute_water_activity_from_m_void(H, p_saturated)
# M_am = compute_GAB_equilibrium_moisture_am(a_w)
# M_cr = compute_GAB_equilibrium_moisture_cryst(a_w)
#
# am_amorph = np.linspace(0, amorphous_material_initial, 100000)
# n_am_amorph = np.size(am_amorph)
# W = np.zeros([n_am_amorph, n_H])
# print('Start creating W')
# for am in range(n_am_amorph):
#     if am % 1000 == 0:
#         print('Currently at iteration', am, 'of', n_am_amorph)
#     M_tot = am_amorph[am] * M_am + (1-am_amorph[am]) * M_cr
#     W[am, :] = M_tot * density_particle * (1-porosity_powder) + H * density_gas * porosity_powder
# print('Done creating W')
# print('W is:')
# print(W[0:10, 0:5])
# np.save('X_table', am_amorph)
# np.save('H_table', H)
# np.save('W[X, H]_table', W)

############################### UPDATED #####################################################
# water_activities = np.linspace(0, 1, 1000)
# temps_celsius = np.linspace(0, 60, 1000)
# am_amorphs = np.linspace(0, amorphous_material_initial, 1000)
#
# n_was = np.size(water_activities)
# n_temps = np.size(temps_celsius)
# n_ams = np.size(am_amorphs)
#
# M_ams = compute_GAB_equilibrium_moisture_am(water_activities)
# M_crs = compute_GAB_equilibrium_moisture_cryst(water_activities)
#
# m_voids = np.zeros([n_was, n_temps])
# W = np.zeros([n_was, n_temps, n_ams])
# for t in range(np.size(temps_celsius)):
#     if t % 100 == 0:
#         print('t is', t)
#     temp_c = temps_celsius[t]
#     m_voids = compute_H_from_aw_temp(water_activities, temp_c)
#
#     for a in range(n_ams):
#         am_am = am_amorphs[a]
#         M_tot = am_am * M_ams + (1 - am_am) * M_crs
#         W[:, t, a] = m_voids * density_gas * porosity_powder + M_tot * density_particle * (1-porosity_powder)
#
# print('Done creating W')
# print('W is:')
# print(W[0:10, 0:5, 0])
# np.save('ams_table', am_amorphs)
# np.save('was_table', water_activities)
# np.save('temps_table', temps_celsius)
# np.save('W[wa, T, am]_table', W)

# ##################### TEMPERATURE ################################################
# H = np.linspace(m_void_initial, .01999, 1000)
# M = np.linspace(m_particle_cryst_initial, m_particle_am_sat * amorphous_material_initial, 1000)
# T = np.linspace(10, 100, 1000)
#
# n_H = np.size(H)
# n_M = np.size(M)
# n_T = np.size(T)
#
# Q = np.zeros([n_H, n_M, n_T])
# print('Start creating T')
# for h in range(n_H):
#     H_current = H[h]
#
#     if h % 100 == 0:
#         print('Currently at iteration', h, 'of', n_H)
#
#     for m in range(n_M):
#         M_current = M[m]
#         temp_celcius = T
#         enthalpy_vapor = heat_capacity_vapor * temp_celcius + heat_of_evaporation_water
#         enthalpy_gas = heat_capacity_air * temp_celcius + H_current * enthalpy_vapor
#         enthalpy_solid = temp_celcius * (heat_capacity_particle + M_current * particle_wet_heat_capacity)
#         Q[h, m, :] = enthalpy_gas * porosity_powder * density_gas + enthalpy_solid * (1 - porosity_powder) * density_particle
#
#
# print('Done creating Q')
# print('Q is:')
# print(Q[0:10, 0:5, 0])
# print(Q[0:10, 0:5, 1])
# np.save('HT_table', H)
# np.save('MT_table', M)
# np.save('T_table', T)
# np.save('Q[H, M, T]_table', Q)