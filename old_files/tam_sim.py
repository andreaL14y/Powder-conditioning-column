from v_ode_functions_2 import*
from glass_transition_curve import compute_glass_temp_mix, glass_temp_lactose, glass_temp_water_1
counter = 0
from old_files.art_functions import compute_crystal_growth_rate
from plot_functions import compute_air_equilibrium
from moisture_sorption_constant import compute_moisture_change

###################################### MAIN EQUATIONS (1-4) ############################################################
def conditioning_column(moisture_matrix, t, temp_gas, relative_humidity):
    global counter
    counter += 1

    ################################### CREATE ALL MATRICES ############################################################
    moisture_particle_cryst, moisture_particle_am, temp_particle, amorphous_material = moisture_matrix
    glass_transition_temp = compute_glass_temp_mix(1 - moisture_particle_am, glass_temp_lactose, glass_temp_water_1)
    diff_temp = temp_particle - glass_transition_temp

    ##################################### UPDATE PARAMETERS ############################################################
    pressure_saturated = compute_p_saturated_vector(temp_gas)

    molar_concentration_moisture = compute_molar_concentration_vector(relative_humidity, pressure_saturated, temp_gas)
    # mass_transfer_coefficient = compute_mass_flux_diffusion_vector(molar_concentration_moisture)
    # print("New way:", mass_transfer_coefficient)
    mass_transfer_coefficient = compute_mass_transfer_coefficient_vector(molar_concentration_moisture)[3]
    # print("Old way:", mass_transfer_coefficient)
    mass_transfer_coefficient /= 200                                                                                    # TODO: CHANGED

    constant = mass_transfer_coefficient * specific_surface_area * pressure_saturated / pressure_ambient

    # constant /= 130

    # print(mass_transfer_coefficient, specific_surface_area, pressure_saturated)
    # print("Pressure ambient", pressure_ambient)
    # print("Constant", constant)

    equilibrium_state_cryst_air = compute_air_equilibrium(moisture_particle_cryst, alpha_parameter_cryst, N_cryst)
    equilibrium_state_am_air = compute_air_equilibrium(moisture_particle_am, alpha_parameter_am, N_am)

    ############################## UPDATE AMORPHOUS MATERIAL ###########################################################
    change_amorphous_material = compute_amorphicity(moisture_particle_am, temp_particle, amorphous_material,
                                                    glass_transition_temp)     # will be negative

    ##################################### UPDATE MOISTURE ##############################################################
    total_equilibrium_term = ((1 - amorphous_material) * (relative_humidity - equilibrium_state_cryst_air) +
                              amorphous_material * (relative_humidity - equilibrium_state_am_air))

    #TODO: stop sorbing after start of crystallization
    #TODO: include crystallization moisture?

    change_m_particle_cryst = constant * (relative_humidity - equilibrium_state_cryst_air)  # sorbing from gas; eq_state_cryst is the air moisture equilibrium for current moisture powder

    # change_m_particle_am_sorption = constant * (relative_humidity - equilibrium_state_am_air)  # sorbing from gas; eq_state_am is the powder moisture equilibrium
    # print("Old change:", change_m_particle_am_sorption)
    change_m_particle_am_sorption = compute_moisture_change(relative_humidity, porosity_powder, moisture_diffusivity, 1)
    print("New change:", change_m_particle_am_sorption)

    # change_m_particle_am_crystallization = change_amorphous_material * constant * (moisture_particle_am - moisture_particle_cryst)  # negative, losing moisture
    change_m_particle_am_crystallization = - (moisture_particle_am - moisture_particle_cryst)  # negative, losing moisture

    change_m_particle_am = change_m_particle_am_sorption #+ change_m_particle_am_crystallization * change_amorphous_material

    ##################################### UPDATE TEMP ##################################################################

    heat_of_sorption_particle = heat_of_sorption * constant * total_equilibrium_term
    heat_transfer_particle  = heat_transfer_coefficient * specific_surface_area * (temp_gas-temp_particle)            # TODO: ORIGINAL
    # heat_transfer_particle  = heat_transfer_coefficient * constant * (temp_gas-temp_particle)
    # heat_transfer_particle  = 0

    heat_of_cryst_particle  = heat_of_crystallization * constant * (-change_amorphous_material)     # kind of copied from heat of sorption

    change_temp_particle    = (heat_of_sorption_particle + heat_transfer_particle + heat_of_cryst_particle) / heat_capacity_particle

    ##################################### CONCATENATE ##################################################################
    moisture_matrix_updated = change_m_particle_cryst, change_m_particle_am, change_temp_particle, change_amorphous_material
    return moisture_matrix_updated


def compute_amorphicity(
        moisture_particle_am_vector, temp_particle_vector, amorphous_material_vector, glass_transition_temp_vector):

    # Only change if above Tg
    k_parameter = compute_k_vector(temp_particle_vector, moisture_particle_am_vector, glass_transition_temp_vector)

    change_amorphous_material = - k_parameter * amorphous_material_vector
    return change_amorphous_material


def compute_k_vector(temperature_vector, moisture_content_vector, glass_transition_temp_vector):
    k_param = compute_crystal_growth_rate(moisture_content_vector, temperature_vector)
    diff_temp = temperature_vector - glass_transition_temp_vector

    sigma = 1/(1 + np.exp(-5 * diff_temp + 3))
    sigma = 1/(1 + np.exp(-0.3 * diff_temp + 4))                                                                        # TODO: CHANGED

    # (1)/(1+ℯ^(-0.3 x+4))
    k_param *= sigma
    # if diff_temp < 30:
    #     k_param = 0
    k_param /= 30                                                                                                       # TODO: CHANGED
    return k_param

def normalize_data(data):
    min = np.min(data)
    max = np.max(data)
    normalized_data = (data - min)/(max - min)
    return normalized_data

   # if counter % 100 == 0:
    #     # print('Moisture in particle', moisture_particle_cryst)
    #     # print('Sat moisture in particle from RH', m_test_sat)
    #     # print('Saturation eq state cryst', saturation_test)
    #     print('Diff eq am', equilibrium_state_am - moisture_particle_am)
    #     print('Diff eq cryst', relative_humidity - equilibrium_state_cryst)
    #     print('Change m am', change_m_particle_am)
    #     print('Change m cr', change_m_particle_cryst)
    #     # print('Change am-cr', change_amorphous_material)
    #     print(constant)
    #     print(mass_transfer_coefficient, specific_surface_area, pressure_saturated/pressure_ambient)
    #     # = mass_transfer_coefficient * specific_surface_area * pressure_saturated / pressure_ambient)
    #     print('\n')
    # if change_amorphous_material < 0:
    #     print('Change am-cr', change_amorphous_material)