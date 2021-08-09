from input_parameters import *

import_tables = 0
# import_tables = 1
if import_tables == 1:
    print('Begin importing tables')
    import_tables_time = time.time()
    Xs = np.load('X_table.npy')
    Hs = np.load('H_table.npy')
    Ws = np.load('W[X, H]_table.npy')

    Xs_new = np.load('ams_table.npy')
    WAs = np.load('was_table.npy')
    Ts = np.load('temps_table.npy')
    Ws_new = np.load('W[wa, T, am]_table.npy')

    import_tables_time = time.time() - import_tables_time
    print(f'Tables imported in {import_tables_time:.3f} s')

tabs = 50
M0_cr = 1.68 * 10 ** (-4)
c_cr = 8.8
f_cr = 0.878
h_cr = 30

M0_am = 0.0488
c_am = 3.23
f_am = 1.16
# ((1 - f_am * water_activity) * (1 - (1 - c_am) * f_am * water_activity))
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

    n_A = 3
    c_1 = 3.54 * 10**4
    c_2 = 108.4
    c_3 = 3 * 10**27
    R = r_gas_constant

    tabs = 40
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

    # if np.any(abs(change_amorphicity) > 1e-1):
    #     print('\nIn compute_kinetics_avrami')
    #     print(change_amorphicity)

    change_amorphicity = -change_amorphicity
    change_amorphicity[change_amorphicity > 0] = 0
    change_amorphicity[change_amorphicity < -1] = -1

    # if np.any(abs(change_amorphicity) > 1e-5):
    #     print(change_amorphicity)
    #     print('Leaving function\n')
    return change_amorphicity


def compute_laplacian(array, boundary, step_length, verbose=False):
    n_steps = array.size

    matrix = np.insert(array, 0, boundary)                              # Add boundary conditions to array at 0
    matrix = np.append(matrix, array[-1])                               # Add final conditions to array again at end

    top = matrix[0:n_steps] - 2 * matrix[1:1+n_steps] + matrix[2:]
    laplacian = top/(step_length**2)

    if verbose:
        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        print('In Laplacian:')
        print('Array[0]:'.ljust(tabs), f'{array[0]:.1e}')
        print('Array[-1]:'.ljust(tabs), f'{array[-1]:.1e}')
        print('Boundary:'.ljust(tabs), f'{boundary:.1e}')
        print('Top:'.ljust(tabs), top)
        print('Laplacian:'.ljust(tabs), laplacian)
    return laplacian

# m_void = np.array([4, -2, 2, 2])
# m_surrounding = 0
# a = compute_laplacian(m_void, m_surrounding, 1)
# print(a)


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

############################### PRESSURES, TEMPS & CONCENTRATIONS ######################################################
def compute_p_saturated_vector(temp_kelvin_vector):
    temp_celsius = temp_kelvin_vector - kelvin
    p_saturated = 10 ** (A - B / (temp_celsius + C))
    p_saturated_pascal_vector = p_saturated * 133.322                                           # Torr to Pascal
    return p_saturated_pascal_vector


def compute_temp_from_energy(H, M, total_energy, verbose=False):
    if verbose:
        print('\nIn')

    temp_top = total_energy - H * heat_of_evaporation_water * porosity_powder * density_gas
    temp_bottom_gas = porosity_powder * density_gas * (heat_capacity_air + H * heat_capacity_vapor)
    temp_bottom_particle = (1 - porosity_powder) * density_particle * (heat_capacity_particle + M * heat_capacity_water)
    temp_celsius = temp_top/(temp_bottom_gas + temp_bottom_particle)
    temp_kelvin = temp_celsius + kelvin
    return temp_kelvin, temp_celsius


def compute_temp_iteratively(temp_kelvin, H, M, total_energy, verbose=False):
    if verbose:
        print('\nIn')
        print('Temp_kelvin:', temp_kelvin)

    # temp_kelvin = temp_kelvin.flatten()
    temp_celcius = temp_kelvin - kelvin
    enthalpy_vapor = heat_capacity_vapor * temp_celcius + heat_of_evaporation_water
    enthalpy_gas = heat_capacity_air * temp_celcius + H * enthalpy_vapor
    enthalpy_solid = temp_celcius * (heat_capacity_particle + M * heat_capacity_water)

    Q = enthalpy_gas * porosity_powder * density_gas + enthalpy_solid * (1 - porosity_powder) * density_particle
    error = np.abs(total_energy - Q)
    if verbose:
        print('Error before mean:', error)
    error = error.mean(axis=0)

    if verbose:
        print('Error after mean:', error)
        print('Out')
    return error


def compute_pressure_water(temp_celsius):
    exponent = 23.4795 - (3990.56/(temp_celsius + 233.833))
    pressure_water = np.exp(exponent)

    return pressure_water


def compute_molar_concentration_vector(relative_humidity_vector, pressure_saturated_vector, temp_vector):
    molar_concentration = relative_humidity_vector * pressure_saturated_vector / (r_gas_constant * temp_vector)
    return molar_concentration


def compute_partial_pressure_moisture_vector(molar_concentration_vector, temperature_vector):           # vapor pressue
    partial_pressure_moisture = molar_concentration_vector * r_gas_constant * temperature_vector
    return partial_pressure_moisture


############################### MOISTURE ISOTHERMS AND EQUILIBRIUMS ####################################################
def compute_GAB_equilibrium_moisture_am(water_activity, verbose=False):
    equilibrium_moisture_particle_vector = M0_am * c_am * f_am * water_activity / ((1 - f_am * water_activity) * (1 - (1 - c_am) * f_am * water_activity))
    if verbose:
        print('\n')
        print(water_activity)
        print(equilibrium_moisture_particle_vector)
    # f(x) = moisture_particle(relative_humidity) = kg/kg
    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector > 1, 0.9, equilibrium_moisture_particle_vector)
    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector < 0, 0, equilibrium_moisture_particle_vector)
    return equilibrium_moisture_particle_vector


# a = compute_GAB_equilibrium_moisture_am(0.85)
# print(a)
def compute_GAB_equilibrium_moisture_cryst(water_activity):
    # H = 1 + (1 - f_cr) * (f_cr * water_activity) ** h / (f_cr * (1 - water_activity))
    H = 1
    # H_prime = 1 + (H - 1) * (1 - f_cr * water_activity) / (H * (1 - water_activity)) * (h_cr + (1 - h_cr) * water_activity)
    H_prime = 1

    equilibrium_moisture_particle_vector = M0_cr * c_cr * f_cr * water_activity * H * H_prime /\
                                           ((1 - f_cr * water_activity) * (1 + (c_cr * H - 1) * f_cr * water_activity))
    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector > 1, 0.9, equilibrium_moisture_particle_vector)
    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector < 0, 0, equilibrium_moisture_particle_vector)
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


def compute_H_from_aw_temp(water_activity, temp_celsius):
    pressure_water = compute_pressure_water(temp_celsius)
    H = 18 * water_activity * pressure_water/(29 * (pressure_ambient - water_activity * pressure_water))
    return H

def compute_moisture_change_am(water_activity, porosity, diffusion, laplacian):
    derivative = compute_derivative_m_GAB(water_activity, cryst=False)

    top = diffusion * laplacian
    bottom = (1 - porosity) * derivative + porosity
    moisture_change = top/bottom
    return moisture_change


def compute_H_and_M_gas_tables(am_amorph, total_weight):
    size = np.size(am_amorph)
    excess = np.full(size, False, dtype=bool)
    extra_weight = np.zeros(size)

    m_void = np.zeros(size)
    for n in range(size):
        x = am_amorph[n]
        w = total_weight[n]

        # Start with finding the correct amount amorph material. Create the weight table for only that amount.
        X = np.where(Xs >= x)       #[0][0]
        if len(X[0]) != 0:
            X = X[0][0]
            ws = Ws[X, :]
        else:
            print('\nIm in else, and x is:', x)
            print('Largest value in X is:', Xs[-1])
            ws = Ws[-1, :]

        # Find the moisture content void that corresponds to the total weight and correct amount amorph material
        h_index = np.where(ws >= w)#[0][0]
        if len(h_index[0]) != 0:
            h_index = h_index[0][0]
        else:
            # print('\nIm in else, and w is:'.ljust(tabs), w)
            # print('The last value in ws is:'.ljust(tabs), ws[-1])
            # print('am_amorph is:'.ljust(tabs), am_amorph[n])
            # print('Part of cylinder:'.ljust(tabs), n)
            h_index = np.size(Hs) - 1
            excess[n] = True
            #np.size(Qs_array) - 1

        m_void[n] = Hs[h_index]
        if excess[n]:
            extra_weight[n] = total_weight[n] - ws[h_index]
            print('\nTotal water[n]:'.ljust(tabs), total_weight[n])
            print('Max water in current ws:'.ljust(tabs), ws[h_index])
    return m_void, extra_weight


def compute_H_and_M_iteratively(water_activity, total_water, temp_celsius, amount_am):
    m_void = compute_H_from_aw_temp(water_activity, temp_celsius)

    m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_am = compute_GAB_equilibrium_moisture_am(water_activity)
    m_particle_tot = m_am * amount_am + m_cryst * (1 - amount_am)

    W = m_void * porosity_powder * density_gas + m_particle_tot * (1 - porosity_powder) * density_particle
    diff = (W - total_water)**2
    error = diff.mean(axis=0)
    return error


def compute_H_and_M_gas_tables_new(am_amorph, total_weight, temp_celsius):  # W[wa, T, am]
    size = np.size(am_amorph)
    extra_weight = np.zeros(size)
    water_activity = np.zeros(size)

    for n in range(size):
        x = am_amorph[n]
        t = temp_celsius[n]

        # Start with finding the correct amount amorph material. Create the weight table for only that amount.
        x_indeces = np.where(Xs_new >= x)       #[0][0]
        if len(x_indeces[0]) != 0:
            x_index = x_indeces[0][0]

        else:
            print('\nIm in else, and x is:', x)
            print('Largest value in X is:', Xs_new[-1])
            x_index = -1

        # Find the moisture content void that corresponds to the total weight and correct amount amorph material
        t_indeces = np.where(Ts >= t)
        if len(t_indeces[0]) != 0:
            t_index = t_indeces[0][0]
        else:
            t_index = -1

        ws = Ws_new[:, t_index, x_index]

        wa_indeces = np.where(ws >= total_weight[n])
        if len(wa_indeces[0]) != 0:
            wa_index = wa_indeces[0][0]
            water_activity[n] = ws[wa_index]
        else:
            water_activity[n] = ws[-1]
            extra_weight[n] = total_weight[n] - ws[-1]

            print('\nTotal water[n]:'.ljust(tabs), total_weight[n])
            print('Max water in current ws:'.ljust(tabs), ws[-1])

    return water_activity, extra_weight


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


# def compute_temp_tables(Q, H, M):
#
#     size = np.size(Q)
#     temp_celsius = np.zeros(size)
#     for n in range(size):
#         q = Q[n]
#         h = H[n]
#         m = M[n]
#
#         # q = np.where(Qs >= q)       #[0][0]
#         h = np.where(HTs >= h)       #[0][0]
#         m = np.where(MTs >= m)       #[0][0]
#
#         if h[0] != []:
#             h_index = h[0][0]         # first place where Qs >= q; first row and first column
#         else:
#             print('Im in else, and h is:', H[n])
#             print('Max in HTs is:', HTs[-1])
#
#         if m[0] != []:
#             m_index = m[0][0]         # first place where Qs >= q; first row and first column
#         else:
#             print('Im in else, and m is:', M[n])
#             print('Max in MTs is:', MTs[-1])
#
#         Qs_array = Qs[h_index, m_index, :]
#
#         q = np.where(Qs_array >= q)  # [0][0]
#         if q[0] != []:
#             q_index = q[0][0]         # first place where Qs >= q; first row and first column
#         else:
#             print('\nIm in else, and q is:', Q[n])
#             print('Max in Qs is:', Qs_array[-1])
#             q_index = np.size(Qs_array) - 1
#
#         temp_celsius[n] = Ts[q_index]
#
#     temp_kelvin = temp_celsius + kelvin
#     return temp_kelvin, temp_celsius
