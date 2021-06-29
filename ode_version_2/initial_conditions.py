from functions import *

p_saturated                     = compute_p_saturated_vector(temp_initial)

# Initial conditions surroundings
molar_conc_sur                  = compute_molar_concentration_vector(relative_humidity_gas_inlet, p_saturated, temp_initial)
vapor_pressure_sur              = compute_partial_pressure_moisture_vector(molar_conc_sur, temp_initial)
water_activity_sur              = vapor_pressure_sur / p_saturated
m_gas_sur                       = 18 * vapor_pressure_sur/(29 * (pressure_ambient - vapor_pressure_sur))

# Initial gas moisture at 0.2 RH
molar_concentration_void_initial= compute_molar_concentration_vector(relative_humidity_bed_initial, p_saturated, temp_initial)
vapor_pressure_void_initial     = compute_partial_pressure_moisture_vector(molar_concentration_void_initial, temp_initial)
water_activity_void_initial     = vapor_pressure_void_initial / p_saturated
m_gas_initial                   = 18 * vapor_pressure_void_initial/(29 * (pressure_ambient - vapor_pressure_void_initial))
m_gas_conc_initial              = m_gas_initial * porosity_powder * gas_density
m_gas_conc_sat                  = m_gas_sur * porosity_powder * gas_density

# Initial moisture cryst powder at start
m_powder_cryst_initial          = compute_GAB_equilibrium_moisture_cryst(water_activity_void_initial)
m_powder_cryst_conc_initial     = m_powder_cryst_initial * (1 - porosity_powder) * particle_density
m_powder_am_initial             = compute_GAB_equilibrium_moisture_am(water_activity_void_initial)
m_powder_am_conc_initial        = m_powder_am_initial * (1 - porosity_powder) * particle_density

# Saturated system
m_particle_cryst_sat            = compute_GAB_equilibrium_moisture_cryst(relative_humidity_gas_inlet)
m_particle_cryst_conc_sat       = m_particle_cryst_sat * (1 - porosity_powder) * particle_density
m_particle_am_sat               = compute_GAB_equilibrium_moisture_am(relative_humidity_gas_inlet)
m_particle_am_conc_sat          = m_particle_am_sat * (1 - porosity_powder) * particle_density

system_conc_sat                 = m_particle_cryst_conc_sat * (1-amorphous_material_initial) + \
                                m_particle_am_conc_sat * amorphous_material_initial + m_gas_conc_sat
