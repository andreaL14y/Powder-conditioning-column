from input_parameters import *
Xs = np.load('X_table.npy')
Hs = np.load('H_table.npy')
Ws = np.load('W[X, H]_table.npy')
tabs = 50

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
    n_steps = am_amorph.size
    am_amorph[am_amorph >= 1] = 0.9999
    am_amorph[am_amorph < 0] = 0

    crystallinity = 1-am_amorph
    Y = crystallinity
    # Y = am_amorph
    # crystallinity = 1-Y

    n_A = 3
    c_1 = 3.54 * 10**4
    c_2 = 108.4
    c_3 = 3 * 10**27
    R = r_gas_constant

    tabs = 40
    # if temp_diff < 0 or am_amorph == 0:
    #     change_crystallinity = 0
    # else:
    #     K = c_3 * (np.exp( -c_1 / (R * (c_2 + temp_diff)) ) )**3
    #     change_amorphicity_1 = n_A * K * am_amorph
    #     change_amorphicity_2 = (-np.log(am_amorph)/K)**((n_A-1)/n_A)
    #     change_crystallinity = change_amorphicity_1 * change_amorphicity_2
    #
    #     if verbose:
    #         print('\nIn compute_kinetics_avrami')
    #         print('K', K)
    #         print('Am amorph:'.ljust(tabs), am_amorph)
    #         print('Crystallinity:'.ljust(tabs), crystallinity)
    #         print('Change p1:'.ljust(tabs), change_amorphicity_1)
    #         print('Change p2:'.ljust(tabs), change_amorphicity_2)
    #         print('Change total:'.ljust(tabs), change_crystallinity)  # positive
    #         print('Leaving function\n')

    # change_crystallinity = np.zeros(n_steps)
    change_amorphicity = np.zeros(n_steps)
    for n in range(n_steps):
        am_am = am_amorph[n]
        am_cr = crystallinity[n]
        if am_am == 0:
            change_amorphicity[n] = 0
        else:
            K = c_3 * (np.exp(-c_1 / (R * (c_2 + temp_diff[n])))) ** 3
            change_amorphicity_1 = n_A * K * am_cr
            change_amorphicity_2 = (-np.log(am_cr) / K) ** ((n_A - 1) / n_A)
            change_amorphicity[n] = change_amorphicity_1 * change_amorphicity_2

    if verbose:
        print('\nIn compute_kinetics_avrami')
        print('K', K)
        print('Am amorph:'.ljust(tabs), am_amorph)
        print('Crystallinity:'.ljust(tabs), crystallinity)
        print('Change p1:'.ljust(tabs), change_amorphicity_1)
        print('Change p2:'.ljust(tabs), change_amorphicity_2)
        print('Change total:'.ljust(tabs), change_amorphicity)  # positive
        print('Leaving function\n')

    # change_amorphicity = - change_crystallinity
    # print(change_amorphicity)
    change_amorphicity = -change_amorphicity
    change_amorphicity[change_amorphicity > 0] = 0
    change_amorphicity[change_amorphicity < -1] = -1
    # if change_amorphicity > 0:
    #     change_amorphicity = 0
    #
    # if change_amorphicity < -1:
    #     change_amorphicity = -1
    return change_amorphicity


def compute_diffusion_laplacian(m_void, m_surrounding, step_length):
    n_steps = m_void.size

    m_matrix = np.insert(m_void, 0, m_surrounding)
    m_matrix = np.append(m_matrix, m_void[-1])
    m_matrix *= gas_density * porosity_powder

    top = m_matrix[0:n_steps] - 2 * m_matrix[1:1+n_steps] + m_matrix[2:]
    laplacian = top/step_length**2

    return laplacian
# m_void = np.array([2, 2, 2, 2, 2])
# m_surrounding = 4
# compute_diffusion_laplacian(m_void, m_surrounding, 1)

def normalize_data(data, zero_one = True):
    min = np.min(data)
    max = np.max(data)
    if min != max:
        if zero_one:
            normalized_data = (data - min)/(max - min)
        else:
            normalized_data = (1 - min) * (data - min) / (max - min) + min
    else:
        normalized_data = data/max
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
    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector > 1, 0.9, equilibrium_moisture_particle_vector)
    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector < 0, 0.9, equilibrium_moisture_particle_vector)
    return equilibrium_moisture_particle_vector


# a = compute_GAB_equilibrium_moisture_am(0.85)
# print(a)
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
    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector > 1, 0.9,
                                                    equilibrium_moisture_particle_vector)
    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector < 0, 0,
                                                    equilibrium_moisture_particle_vector)
    return equilibrium_moisture_particle_vector
# a = compute_GAB_equilibrium_moisture_cryst(1)
# print(a)

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
    water_activity = np.where(water_activity > 1, 0.99, water_activity)
    return water_activity


def compute_moisture_change_am(water_activity, porosity, diffusion, laplacian):
    derivative = compute_derivative_m_GAB(water_activity, cryst=False)

    top = diffusion * laplacian
    bottom = (1 - porosity) * derivative + porosity
    moisture_change = top/bottom
    return moisture_change


def compute_H_and_M_gas_tables(am_amorph, total_weight):

    size = np.size(am_amorph)
    m_void = np.zeros(size)
    for n in range(size):
        x = am_amorph[n]
        w = total_weight[n]
        # print(x)
        X = np.where(Xs >= x)       #[0][0]
        if X[0] != []:
            X = X[0][0]
            ws = Ws[X, :]
        else:
            # print('Im in else, and X is:', X)
            ws = Ws[-1, :]

        h_index = np.where(ws >= w)#[0][0]
        if h_index[0] != []:
            h_index = h_index[0][0]
        else:
            print('Im in else, and w is:'.ljust(tabs), w)
            print('The last value in ws is:'.ljust(tabs), ws[-1])
            print('am_amorph is:', am_amorph)

        m_void[n] = Hs[h_index]
    return m_void
# a = compute_H_and_M_gas_tables(0.10, 4.56)

def compute_H_and_M_gas_fractions(total_water, max_diff, m_void, p_saturated, amount_am, m_gas_sur, n_space_steps,
                                  gas_fraction_min, gas_fraction_max, gas_fraction_initial):
    counter_2 = 0
    # Would be if everything equilibrium with m_void
    water_activity = compute_water_activity_from_m_void(m_void, p_saturated)
    m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_am = compute_GAB_equilibrium_moisture_am(water_activity)
    m_particle_tot = m_am * amount_am + m_cryst * (1 - amount_am)

    W = m_void * porosity_powder * gas_density + m_particle_tot * (1 - porosity_powder) * particle_density
    diff = W - total_water

    tabs = 50
    gas_fraction = np.zeros(n_space_steps) + gas_fraction_initial

    while np.any(abs(diff) > max_diff):
        counter_2 += 1
        divider = 5000
        # if counter_2 % 1 == 0:
        if counter_2 % divider == 0:
            print('\nCounter in while:', counter_2)
            # print('m_void:'.ljust(tabs), m_void[0])
            # print('total_water:'.ljust(tabs), total_water)
            # print('w_act:'.ljust(tabs), water_activity[0])
            # print('m_cr:'.ljust(tabs), m_cryst[0])
            # print('m_am:'.ljust(tabs), m_am[0])
            # print('m_p_tot:'.ljust(tabs), m_particle_tot[0])
            # print('W:'.ljust(tabs), W)
            print('Diff is:'.ljust(tabs), diff[0])
            print('Gas fraction now:'.ljust(tabs), gas_fraction)
            print('Gas fraction min:'.ljust(tabs), gas_fraction_min)
            print('Gas fraction max:'.ljust(tabs), gas_fraction_max)
        m_void_conc = gas_fraction * total_water
        m_void = m_void_conc/(porosity_powder * gas_density)

        water_activity = compute_water_activity_from_m_void(m_void, p_saturated)
        m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
        m_am = compute_GAB_equilibrium_moisture_am(water_activity)
        m_particle_tot = m_am * amount_am + m_cryst * (1 - amount_am)
        # m_particle_tot = 1 - m_void

        W = m_void_conc + m_particle_tot * (1 - porosity_powder) * particle_density
        diff = W - total_water

        # max_diff = 0.0001
        gas_to_max = abs(gas_fraction_max - gas_fraction)
        gas_to_min = abs(gas_fraction - gas_fraction_min)
        if counter_2 % divider == 0:
            print('And diff now:'.ljust(tabs), diff[0])
            # print('Before:'.ljust(tabs), gas_fraction[0])
        gas_fraction = np.where(diff > max_diff * 10, gas_fraction - gas_to_min/100, gas_fraction)
        gas_fraction = np.where((diff > max_diff * 2) & (diff < max_diff * 10), gas_fraction - gas_to_min/5000, gas_fraction)
        gas_fraction = np.where((diff > max_diff * 1.5) & (diff < max_diff * 2), gas_fraction - gas_to_min/10000, gas_fraction)
        gas_fraction = np.where((diff > max_diff) & (diff < max_diff * 1.5), gas_fraction - gas_to_min/100000, gas_fraction)
        # gas_fraction = np.where(diff < -max_diff, gas_fraction + gas_to_min/1000, gas_fraction)
        gas_fraction = np.where(diff < -max_diff * 10, gas_fraction + gas_to_max/1000, gas_fraction)
        gas_fraction = np.where((diff < -max_diff * 2) & (max_diff > - max_diff * 10), gas_fraction + gas_to_max/5000, gas_fraction)
        gas_fraction = np.where((diff < -max_diff) & (max_diff > - max_diff * 2), gas_fraction + gas_to_max/100000, gas_fraction)
    return m_void


def compute_H_and_M(total_water, max_diff, m_void, p_saturated, amount_am, m_gas_sur, n_space_steps, gas_fraction):
    counter_2 = 0
    water_activity = compute_water_activity_from_m_void(m_void, p_saturated)
    m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_am = compute_GAB_equilibrium_moisture_am(water_activity)
    m_particle_tot = m_am * amount_am + m_cryst * (1 - amount_am)

    W = m_void * porosity_powder * gas_density + m_particle_tot * (1 - porosity_powder) * particle_density
    diff = W - total_water

    #### NEW IDEA
    m_void = W * gas_fraction/(porosity_powder * gas_density)
    water_activity = compute_water_activity_from_m_void(m_void, p_saturated)
    m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_am = compute_GAB_equilibrium_moisture_am(water_activity)
    m_particle_tot = m_am * amount_am + m_cryst * (1 - amount_am)

    W = m_void * porosity_powder * gas_density + m_particle_tot * (1 - porosity_powder) * particle_density
    diff = W - total_water

    # max_diff = 0.001
    while np.any(abs(diff) > max_diff):

        counter_2 += 1
        if counter_2 % 1000 == 0:
            print('\nIn while loop, counter:', counter_2)
            print('m_void', m_void[0])
            print('water act', water_activity[0])
            print('m_am', m_am[0])
            print('m_cr', m_cryst[0])
            print('M', m_particle_tot[0])
            print('W', W[0])
            print('diff', diff[0])
            print('m_am', amount_am[0])
        m_void = np.where((diff > 1),                             m_void * 0.7, m_void)
        m_void = np.where((diff > 0.1)        & (diff < 1),       m_void * 0.9, m_void)
        m_void = np.where((diff > 0.01)       & (diff < 0.1),     m_void * 0.99, m_void)
        m_void = np.where((diff > max_diff)   & (diff < 0.01),    m_void * 0.999, m_void)

        m_void = np.where((diff < -0.1),                          m_void * 1.05, m_void)
        m_void = np.where((diff < -0.01)      & (diff > -0.1),    m_void * 1.005, m_void)
        m_void = np.where((diff < -max_diff)  & (diff > -0.01),   m_void * 1.0001, m_void)

        # for d in range(n_space_steps):
        #     if diff[d] > 1:
        #         m_void[d] *= 0.7
        #     elif diff[d] > 0.1:
        #         m_void[d] *= 0.9
        #     elif diff[d] > 0.01:
        #         m_void[d] *= 0.99
        #     elif diff[d] > max_diff:
        #         m_void[d] *= 0.999
        #
        #     elif diff[d] < -0.1:
        #         m_void[d] *= 1.05
        #     elif diff[d] < -0.01:
        #         m_void[d] *= 1.005
        #     elif diff[d] < -max_diff:
        #         m_void[d] *= 1.0001

        # Can't hold more moisture than saturation
        m_void = np.where(m_void > m_gas_sur, m_gas_sur, m_void)
        m_void = np.where(m_void < 0, 0.0001, m_void)
        water_activity = compute_water_activity_from_m_void(m_void, p_saturated)
        m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
        m_am = compute_GAB_equilibrium_moisture_am(water_activity)
        m_particle_tot = m_am * amount_am + m_cryst * (1 - amount_am)

        estimated_water = m_void * porosity_powder * gas_density + m_particle_tot * (1 - porosity_powder) * particle_density
        diff = estimated_water - total_water

    return m_void
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