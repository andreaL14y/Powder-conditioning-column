from input_parameters import *

def compute_glass_temp_mix(w1_lactose, t_g_lactose, t_g2):
  w2 = 1 - w1_lactose
  k_param = 6.7                                                                 # Glass Transitions and Crystallization in Milk Powders
  glass_temp_mix = (w1_lactose * t_g_lactose + w2 * t_g2 * k_param)/(w1_lactose + w2 * k_param)
  return(glass_temp_mix)


def compute_heat_of_sorption(water_activity, temp):

    heat_of_sorption = - np.log(water_activity) * r_gas_constant * temp     # J/mol
    heat_of_sorption /= molar_mass_moisture                                 # J/mol / kg/mol = J/kg
    return heat_of_sorption


def compute_kinetics_avrami(am_amorph, temp_diff, verbose=False):
    if am_amorph >= 1:
        am_amorph = 0.9999
    elif am_amorph < 0:
        am_amorph = 0

    Y = am_amorph
    crystallinity = 1-Y

    n_A = 3
    c_1 = 3.54 * 10**4
    c_2 = 108.4
    c_3 = 3 * 10**27
    R = r_gas_constant

    tabs = 40
    if temp_diff < 0 or am_amorph == 0:
        change_crystallinity = 0
    else:
        K = c_3 * (np.exp( -c_1 / (R * (c_2 + temp_diff)) ) )**3
        change_amorphicity_1 = n_A * K * am_amorph
        change_amorphicity_2 = (-np.log(am_amorph)/K)**((n_A-1)/n_A)
        change_crystallinity = change_amorphicity_1 * change_amorphicity_2

        if verbose:
            print('\nIn compute_kinetics_avrami')
            print('K', K)
            print('Am amorph:'.ljust(tabs), am_amorph)
            print('Crystallinity:'.ljust(tabs), crystallinity)
            print('Change p1:'.ljust(tabs), change_amorphicity_1)
            print('Change p2:'.ljust(tabs), change_amorphicity_2)
            print('Change total:'.ljust(tabs), change_crystallinity)  # positive
            print('Leaving function\n')

        if change_crystallinity < 0:
            print('Change crystallinity < 0! Das nix gut diese!!!')

    change_amorphicity = - change_crystallinity
    if change_amorphicity > 0:
        change_amorphicity = 0

    if change_amorphicity < -1:
        change_amorphicity = -1
    return change_amorphicity


def compute_diffusion_laplacian(m_void, m_surrounding):
    m_diff = (m_surrounding - m_void) * gas_density * porosity_powder                                # kg/m3
    step = 0.001                                                    # 1mm step...?
    step = 0.006
    laplacian = m_diff/step**2                                     # kg/m5
    return laplacian


def normalize_data(data):
    min = np.min(data)
    max = np.max(data)
    normalized_data = (data - min)/(max - min)
    return normalized_data

############################### PRESSURES & CONCENTRATIONS #############################################################
def compute_p_saturated_vector(temp_kelvin_vector):
    temp_celsius = temp_kelvin_vector - kelvin
    p_saturated = 10 ** (A - B / (temp_celsius + C))
    p_saturated_pascal_vector = p_saturated * 133.322                                           # Torr to Pascal
    return p_saturated_pascal_vector


def compute_molar_concentration_vector(relative_humidity_vector, pressure_saturated_vector, temp_vector):
    molar_concentration = relative_humidity_vector * pressure_saturated_vector / (r_gas_constant * temp_vector)
    return molar_concentration


def compute_partial_pressure_moisture_vector(molar_concentration_vector, temperature_vector):           # vapor pressue
    partial_pressure_moisture = molar_concentration_vector * r_gas_constant * temperature_vector
    return partial_pressure_moisture


############################### MOISTURE ISOTHERMS AND EQUILIBRIUMS ####################################################
def compute_GAB_equilibrium_moisture_am(relative_humidity_vector):
    M0 = 0.0488
    c = 3.23
    f = 1.16
    equilibrium_moisture_particle_vector = M0 * c * f * relative_humidity_vector / \
                                           ((1 - f * relative_humidity_vector) * (
                                                       1 - (1 - c) * f * relative_humidity_vector))
    # f(x) = moisture_particle(relative_humidity) = kg/kg
    return equilibrium_moisture_particle_vector


def compute_GAB_equilibrium_moisture_cryst(relative_humidity_vector):
    M0 = 1.68 * 10 ** (-4)
    c = 8.8
    f = 0.878
    h = 30
    H = 1 + (1 - f) * (f * relative_humidity_vector) ** h / (f * (1 - relative_humidity_vector))
    H_prime = 1 + (H - 1) * (1 - f * relative_humidity_vector) / (H * (1 - relative_humidity_vector)) * (
                h + (1 - h) * relative_humidity_vector)

    equilibrium_moisture_particle_vector = M0 * c * f * relative_humidity_vector * H * H_prime / \
                                           ((1 - f * relative_humidity_vector) * (
                                                       1 + (c * H - 1) * f * relative_humidity_vector))
    return equilibrium_moisture_particle_vector


def compute_derivative_m_GAB(water_activity, cryst=True):
    if cryst:
        M0 = 1.68 * 10 ** (-4)
        c = 8.8
        f = 0.878
        h = 30
    else:
        M0 = 0.0488
        c = 3.23
        f = 1.16

    top = c * f * M0 * (c * f**2 * water_activity**2 - f**2 * water_activity**2 +1)
    bottom = (f * water_activity - 1)**2 * (c * f * water_activity - f * water_activity + 1)**2

    derivative = top/bottom
    return derivative


def compute_water_activity_from_m_void(m_void, p_saturated):
    vapor_pressure = 29 * m_void * pressure_ambient / (18 + 29 * m_void)
    water_activity = vapor_pressure / p_saturated
    return water_activity


def compute_moisture_change_am(water_activity, porosity, diffusion, laplacian):
    derivative = compute_derivative_m_GAB(water_activity, cryst=False)

    top = diffusion * laplacian
    bottom = (1 - porosity) * derivative + porosity
    moisture_change = top/bottom
    return moisture_change


############################### CURVE FIT ACCORDING TO AZ ##############################################################
def compute_moisture_particle_from_RH(relative_humidity, alpha, N):
    moisture_particle = (-np.log(1-relative_humidity)/alpha)**(1/N)
    return moisture_particle


def compute_air_equilibrium(moisture_particle_vector, alpha, N):
    equilibrium_air = (1 - np.exp(-alpha * moisture_particle_vector ** N))
    return equilibrium_air


def equilibrium_curve_fit(X, alpha, N):
    return 1 - np.exp(-alpha * X**N)


###################### COMPUTE a AND N AMORHPOUS LACTOSE FROM GAB ######################################################
relative_humidities     = np.linspace(0, 0.85, 10000)
moistures_am            = compute_GAB_equilibrium_moisture_am(relative_humidities)
moistures_am_der        = compute_derivative_m_GAB(relative_humidities, cryst=False)

moistures_cryst         = compute_GAB_equilibrium_moisture_cryst(relative_humidities)
eq_air_cryst            = compute_air_equilibrium(moistures_cryst, alpha_parameter_cryst, N_cryst)

start_params            = (3, 0.7) # start with values near those we expect
params, cv              = scipy.optimize.curve_fit(equilibrium_curve_fit, moistures_am, relative_humidities, start_params)
alpha_parameter_am, N_am = params