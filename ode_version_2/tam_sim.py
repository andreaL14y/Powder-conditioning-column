from v_ode_functions_2 import*
from glass_transition_curve import compute_glass_temp_mix, glass_temp_lactose, glass_temp_water_1
counter = 0
from art_functions import compute_entropy_of_activation_S, compute_enthalpy_of_activation_H, compute_crystal_growth_rate

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
    mass_transfer_coefficient = compute_mass_transfer_coefficient_vector(molar_concentration_moisture)[3]
    constant = mass_transfer_coefficient * specific_surface_area * pressure_saturated / pressure_ambient

    m_test_sat = compute_moisture_particle_from_RH_cryst(relative_humidity)
    saturation_test = compute_air_equilibrium_cryst(moisture_cryst_particle_saturated)
    # equilibrium_state_cryst = compute_equilibrium_moisture_vector(moisture_particle_cryst)
    equilibrium_state_cryst = compute_GAB_equilibrium_moisture_cryst(relative_humidity)      #TODO: new eq function
    equilibrium_state_am = compute_GAB_equilibrium_moisture_am(relative_humidity)
    y = compute_Y_from_RH(0.45, pressure_saturated)

    # moisture_estimate = np.where(moistures_cryst_eq > moisture_particle_cryst)
    moisture_estimate_cryst = np.argwhere(moistures_cryst_eq > moisture_particle_cryst)[0]
    equilibrium_state_cryst_air = relative_humidities_eq[moisture_estimate_cryst]

    moisture_estimate_am = np.argwhere(moistures_am_eq > moisture_particle_am)[0]
    equilibrium_state_am_air = relative_humidities_eq[moisture_estimate_am]

    print(moisture_estimate_cryst, moistures_cryst_eq[moisture_estimate_cryst], equilibrium_state_cryst_air)
    print(moisture_estimate_am, moistures_am_eq[moisture_estimate_am], equilibrium_state_am_air)


    ############################## UPDATE AMORPHOUS MATERIAL ###########################################################
    print(' ')
    print(diff_temp)
    change_amorphous_material = compute_amorphicity(moisture_particle_am, temp_particle, amorphous_material, glass_transition_temp)     # will be negative
    # print('Moisture particle', moisture_particle_am)
    # print('Temp particle', temp_particle)
    # print('Am amount', amorphous_material)
    # print('Change am', change_amorphous_material)


    ##################################### UPDATE MOISTURE ##############################################################
    total_equilibrium_term = ((1 - amorphous_material) * (relative_humidity - equilibrium_state_cryst_air) +
                              amorphous_material * (relative_humidity - equilibrium_state_am_air))

    #TODO: stop sorbing after start of crystallization
    #TODO: include crystallization moisture?

    gas_contact_fraction = 0.02         # TODO: MADE UP!
    change_m_particle_cryst = gas_contact_fraction * constant * (relative_humidity - equilibrium_state_cryst_air)  # sorbing from gas; eq_state_cryst is the air moisture equilibrium for current moisture powder

    change_m_particle_am_sorption = gas_contact_fraction * constant * (relative_humidity - equilibrium_state_am_air)  # sorbing from gas; eq_state_am is the powder moisture equilibrium
    change_m_particle_am_crystallization = change_amorphous_material * constant * (moisture_particle_am - moisture_particle_cryst)  # negative, losing moisture
    # change_m_particle_am_crystallization = - (moisture_particle_am - moisture_particle_cryst)  # negative, losing moisture

    exponent = np.exp(-5 * diff_temp + 1)
    sigma = 1/(1 + exponent)
    sigma_rev = exponent/(1 + exponent)

    change_m_particle_am_crystallization *= sigma
    change_m_particle_am_sorption *= sigma_rev * amorphous_material

    # if diff_temp > 0:
    #     change_m_particle_am_sorption = 0
    # else:
    #     change_m_particle_am_crystallization = 0

    change_m_particle_am = change_m_particle_am_sorption + change_m_particle_am_crystallization

    ##################################### UPDATE TEMP ##################################################################

    heat_of_sorption_particle = constant * total_equilibrium_term * heat_of_sorption
    heat_transfer_particle  = heat_transfer_coefficient * specific_surface_area * (temp_gas-temp_particle)

    heat_of_cryst_particle  = heat_of_crystallization * constant * (-change_amorphous_material)     # kind of copied from heat of sorption

    change_temp_particle    = (heat_of_sorption_particle + heat_transfer_particle + heat_of_cryst_particle) / particle_heat_capacity

    ##################################### CONCATENATE ##################################################################
    moisture_matrix_updated = change_m_particle_cryst, change_m_particle_am, change_temp_particle, change_amorphous_material
    return moisture_matrix_updated


def compute_amorphicity(
        moisture_particle_am_vector, temp_particle_vector, amorphous_material_vector, glass_transition_temp_vector):

    # Only change if above Tg
    k_parameter = compute_k_vector(temp_particle_vector, moisture_particle_am_vector, glass_transition_temp_vector)
    print('k is',k_parameter)
    change_amorphous_material = - k_parameter * amorphous_material_vector
    return change_amorphous_material


def compute_k_vector(temperature_vector, moisture_content_vector, glass_transition_temp_vector):
    k_param = compute_crystal_growth_rate(moisture_content_vector, temperature_vector)
    diff_temp = temperature_vector - glass_transition_temp_vector

    # k_vector *= 1/(1 + np.exp(-3 * diff_temp - 3))          # Sigmoid to ramp up, TODO: change this, nothing <0
    sigma = 1/(1 + np.exp(-5 * diff_temp + 3))
    k_param *= sigma
    if diff_temp < 0:
        k_param = 0
    # elif diff_temp > 10:
    #     print(diff_temp, 'is diff_temp')
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