from functions import *

p_saturated                     = compute_p_saturated_vector(temp_initial)

# Initial conditions surroundings
molar_conc_sur                  = compute_molar_concentration_vector(relative_humidity_gas_inlet, p_saturated, temp_initial)
vapor_pressure_sur              = compute_partial_pressure_moisture_vector(molar_conc_sur, temp_initial)
water_activity_sur              = vapor_pressure_sur / p_saturated
m_gas_sur                       = 18 * vapor_pressure_sur/(29 * (pressure_ambient - vapor_pressure_sur))
m_gas_conc_sat                  = m_gas_sur * porosity_powder * gas_density

# Initial gas moisture at 0.2 RH
molar_concentration_void_initial= compute_molar_concentration_vector(relative_humidity_bed_initial, p_saturated, temp_initial)
vapor_pressure_void_initial     = compute_partial_pressure_moisture_vector(molar_concentration_void_initial, temp_initial)
water_activity_void_initial     = vapor_pressure_void_initial / p_saturated
m_gas_initial                   = 18 * vapor_pressure_void_initial/(29 * (pressure_ambient - vapor_pressure_void_initial))
m_gas_conc_initial              = m_gas_initial * porosity_powder * gas_density

# Initial moisture cryst powder at start
m_powder_cryst_initial          = compute_GAB_equilibrium_moisture_cryst(water_activity_void_initial)
m_powder_cryst_conc_initial     = m_powder_cryst_initial * (1 - porosity_powder) * particle_density
m_powder_am_initial             = compute_GAB_equilibrium_moisture_am(water_activity_void_initial)
m_powder_am_conc_initial        = m_powder_am_initial * (1 - porosity_powder) * particle_density
m_particle_tot_conc_initial     = m_powder_cryst_conc_initial * (1-amorphous_material_initial) + m_powder_am_conc_initial * amorphous_material_initial


# Saturated system
m_particle_cryst_sat            = compute_GAB_equilibrium_moisture_cryst(relative_humidity_gas_inlet)
m_particle_cryst_conc_sat       = m_particle_cryst_sat * (1 - porosity_powder) * particle_density
m_particle_am_sat               = compute_GAB_equilibrium_moisture_am(relative_humidity_gas_inlet)
m_particle_am_conc_sat          = m_particle_am_sat * (1 - porosity_powder) * particle_density
m_particle_tot_conc_sat         = m_particle_cryst_conc_sat * (1-amorphous_material_initial) + m_particle_am_conc_sat * amorphous_material_initial
system_conc_sat                 = m_particle_cryst_conc_sat * (1-amorphous_material_initial) + \
                                m_particle_am_conc_sat * amorphous_material_initial + m_gas_conc_sat
m_gas_sat                       = m_gas_conc_sat/(porosity_powder * gas_density)

# print(system_conc_sat)
# print(m_particle_tot_conc_initial)
gas_fraction = m_gas_conc_initial/(m_gas_conc_initial + m_particle_tot_conc_initial)
gas_fraction_initial = gas_fraction
gas_fraction_2 = m_gas_conc_sat/(m_gas_conc_sat + m_particle_tot_conc_sat)
gas_fraction_min = gas_fraction_2

gas_fraction_3 = m_gas_conc_sat/(m_gas_conc_sat + m_particle_cryst_conc_sat)
gas_fraction_max_tot = gas_fraction_3


####################################### CREATE TABLES ##################################################################
# H = np.linspace(m_gas_initial, .01999, 100000)
# # print(m_gas_sat)
# # print(m_gas_initial)
# n_H = np.size(H)
# a_w = compute_water_activity_from_m_void(H, p_saturated)
# # print(a_w)
# # print(np.where(a_w >= 0.99))
# M_am = compute_GAB_equilibrium_moisture_am(a_w)
# M_cr = compute_GAB_equilibrium_moisture_cryst(a_w)
#
# am_amorph = np.linspace(0, amorphous_material_initial, 10000)
# n_am_amorph = np.size(am_amorph)
# W = np.zeros([n_am_amorph, n_H])
# print('Start creating W')
# for am in range(n_am_amorph):
#     if am % 1000 == 0:
#         print('Currently at iteration', am, 'of', n_am_amorph)
#     M_tot = am_amorph[am] * M_am + (1-am_amorph[am]) * M_cr
#     W[am, :] = M_tot * particle_density * (1-porosity_powder) + H * gas_density * porosity_powder
# print('Done creating W')
# print('W is:')
# print(W[0:10, 0:5])
# np.save('X_table', am_amorph)
# np.save('H_table', H)
# np.save('W[X, H]_table', W)


