from input_parameters import *

tabs = 50
M0_cr = 1.68 * 10 ** (-4)
c_cr = 8.8
f_cr = 0.878
h_cr = 30

M0_am = 0.0488
c_am = 3.23
f_am = 1.16

n_A = 3
# n_A = 2.8
# n_A_inverse = 1/n_A
avrami_exponent = (n_A - 1) / n_A
c_1 = 3.54 * 10**4
c_2 = 108.4
c_3 = 3 * 10**27

one_over_29 = 1/29
one_over_18 = 1/18

tabs = 40

def compute_glass_temp_mix(w1_lactose, t_g_lactose, t_g2):
  w2 = 1 - w1_lactose
  k_param = 6.7                                                                 # Glass Transitions and Crystallization in Milk Powders
  glass_temp_mix = (w1_lactose * t_g_lactose + w2 * t_g2 * k_param)/(w1_lactose + w2 * k_param)
  return(glass_temp_mix)


def compute_enthalpy_humid_air(temp_celsius, m_void):
    enthalpy = heat_capacity_air * temp_celsius + m_void * (heat_capacity_vapor * temp_celsius + heat_of_evaporation_water) # kg/m3
    # enthalpy *= density
    return enthalpy


def compute_enthalpy_vapor(temp_celsius):
    enthalpy = heat_capacity_vapor * temp_celsius + heat_of_evaporation_water # kg/m3
    return enthalpy


def compute_heat_of_sorption(water_activity, temp):
    water_activity = np.where(water_activity == 0, 1e-5, water_activity)    # TODO: APPROXIMATION
    heat_of_sorption = - np.log(water_activity) * r_gas_constant * temp     # J/mol
    heat_of_sorption /= molar_mass_moisture                                 # J/mol / kg/mol = J/kg
    return heat_of_sorption


def compute_kinetics_avrami(am_amorph, temp_diff, verbose=False):
    n_steps = am_amorph.size
    am_amorph[am_amorph >= 1] = 0.9999
    am_amorph[am_amorph < 0] = 0
    crystallinity = 1-am_amorph

    change_amorphicity = np.zeros(n_steps)
    for n in range(n_steps):
        am_am = am_amorph[n]
        am_cr = crystallinity[n]
        if am_am == 0:
            change_amorphicity[n] = 0
        else:
            K_exp = -c_1 / (r_gas_constant * (c_2 + temp_diff[n]))
            K = c_3 * (np.exp(K_exp)) ** 3

            change_amorphicity_1 = n_A * K * am_cr
            # change_amorphicity_2 = (-np.log(am_cr) / K) ** ((n_A - 1) / n_A)
            change_amorphicity_2 = (-np.log(am_cr) / K) ** avrami_exponent
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

    change_amorphicity = -change_amorphicity
    change_amorphicity[change_amorphicity > 0] = 0
    change_amorphicity[change_amorphicity < -1] = -1
    return change_amorphicity


def compute_laplacian(array, boundary, step_length, double_bc=False, inflow=False):
    #     m_sur_bounds = np.array([m_gas_sur, m_void[0]])
    n_steps = array.size

    if double_bc:
        matrix = np.insert(array, 0, boundary[0])   # Add boundary conditions to array at 0
        matrix = np.append(matrix, boundary[1])     # Add second boundary conditions at end
    else:
        matrix = np.insert(array, 0, boundary)                              # Add boundary conditions to array at 0
        matrix = np.append(matrix, array[-1])                               # Add final conditions to array again at end

    top = matrix[0:n_steps] - 2 * matrix[1:1 + n_steps] + matrix[2:]
    laplacian = top / (step_length * step_length)

    if inflow:
        from_sur = (matrix[0] - matrix[1]) / (step_length * step_length)
        return laplacian, from_sur

    return laplacian

# a = np.array([1, 2, 4, 0])
# l = compute_laplacian(a, 0, 1)
# print(l)
#
# a = np.array([1, 2, 4, 0])
# l = compute_laplacian(a, 0, 1, energy=True)
# print(l)
#
# a = np.array([1, 2, 4, 0])
# l = compute_laplacian(a, [0, 0], 1, double_bc=True)
# print(l)

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
def compute_air_temp_c_from_H(H, enthalpy):
    temp_c = (enthalpy - H * heat_of_evaporation_water)/(heat_capacity_air + H * heat_capacity_vapor)
    return temp_c


def compute_p_saturated_vector(temp_kelvin_vector):
    temp_celsius = temp_kelvin_vector - kelvin
    p_saturated = 10 ** (A - B / (temp_celsius + C))
    p_saturated_pascal_vector = p_saturated * 133.322                                           # Torr to Pascal
    return p_saturated_pascal_vector


def compute_pressure_water(temp_celsius):
    exponent = 23.4795 - (3990.56/(temp_celsius + 233.833))
    pressure_water = np.exp(exponent)
    return pressure_water


def compute_temp_from_energy(H, M, total_energy, excess_water_W, verbose=False):
    # if np.any(excess_water != 0):
    #     temp_water = excess_water * heat_capacity_water
    excess_water_M = excess_water_W/(density_particle * (1-porosity_powder))
    M = M + excess_water_M

    temp_top = total_energy - H * heat_of_evaporation_water * porosity_powder * density_gas
    temp_bottom_gas = porosity_powder * density_gas * (heat_capacity_air + H * heat_capacity_vapor)
    temp_bottom_particle = (1 - porosity_powder) * density_particle * (heat_capacity_particle + M * heat_capacity_water)

    temp_celsius = temp_top/(temp_bottom_gas + temp_bottom_particle)
    temp_kelvin = temp_celsius + kelvin
    return temp_kelvin, temp_celsius


def compute_molar_concentration_vector(relative_humidity_vector, pressure_saturated_vector, temp_vector):
    molar_concentration = relative_humidity_vector * pressure_saturated_vector / (r_gas_constant * temp_vector)
    return molar_concentration


def compute_partial_pressure_moisture_vector(molar_concentration_vector, temperature_vector):           # vapor pressue
    partial_pressure_moisture = molar_concentration_vector * r_gas_constant * temperature_vector
    return partial_pressure_moisture


############################### MOISTURE ISOTHERMS AND EQUILIBRIUMS ####################################################
def compute_GAB_equilibrium_moisture_am(water_activity, verbose=False):
    water_activity_max = 0.6            # TODO: used to be 0.6
    water_activity = np.where(water_activity > water_activity_max, water_activity_max, water_activity)

    top = M0_am * c_am * f_am * water_activity
    bottom = (1 - f_am * water_activity) * (1 - (1 - c_am) * f_am * water_activity)
    # equilibrium_moisture_particle_vector = np.where(bottom == 0, 0.9, top/bottom)
    equilibrium_moisture_particle_vector = top / bottom
    if verbose:
        print('\n')
        print(water_activity)
        print(equilibrium_moisture_particle_vector)
    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector < 0, 0, equilibrium_moisture_particle_vector)
    return equilibrium_moisture_particle_vector


def compute_GAB_equilibrium_moisture_cryst(water_activity, H_H = True):
    water_activity_limit = 0.85
    # water_activity = np.where(water_activity > water_activity_max, water_activity_max, water_activity)
    if np.all(water_activity < water_activity_limit) or H_H == True:
        H = 1
        H_prime = 1
    else:
        H = 1 + (1 - f_cr) * (f_cr * water_activity) ** h_cr / (f_cr * (1 - water_activity))
        H_prime = 1 + (H - 1) * (1 - f_cr * water_activity) / (H * (1 - water_activity)) * (h_cr + (1 - h_cr) * water_activity)

    top = M0_cr * c_cr * f_cr * water_activity * H * H_prime
    bottom = (1 - f_cr * water_activity) * (1 + (c_cr * H - 1) * f_cr * water_activity)
    equilibrium_moisture_particle_vector = top/bottom
    return equilibrium_moisture_particle_vector


def compute_water_activity_from_m_void(m_void, temp_celsius):
    vapor_pressure = 29 * m_void * pressure_ambient / (18 + 29 * m_void)
    water_pressure = compute_pressure_water(temp_celsius)
    water_activity = vapor_pressure / water_pressure
    water_activity = np.where(water_activity > 1, 0.99, water_activity)
    return water_activity


# def compute_H_from_water_activity_temp(water_activity, temp_celsius):
def compute_H_from_water_activity_temp(water_activity):
    # pressure_water = compute_pressure_water(temp_celsius)
    pressure_water = compute_pressure_water(temp_initial_celsius)

    H = 18 * water_activity * pressure_water/(29 * (pressure_ambient - water_activity * pressure_water))
    return H


# def compute_density_air(temp_celsius, H):
def compute_density_air(H):
    top = kelvin
    # bottom = 22.4 * (kelvin + temp_celsius) * (one_over_29 + H * one_over_18)
    bottom = 22.4 * (kelvin + temp_initial_celsius) * (one_over_29 + H * one_over_18)
    density = top/bottom
    # density = density_gas
    return density


# def compute_H_and_M_iteratively(water_activity, total_water, temp_celsius, amount_am):
def compute_H_and_M_iteratively(water_activity, total_water, amount_am):
    # excess = np.where(water_activity >= 1, True, False)

    # m_void = compute_H_from_water_activity_temp(water_activity, temp_celsius)
    m_void = compute_H_from_water_activity_temp(water_activity)
    m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_am = compute_GAB_equilibrium_moisture_am(water_activity)
    m_particle_tot = m_am * amount_am + m_cryst * (1 - amount_am)
    # density_gas_vector = compute_density_air(temp_celsius, m_void)
    density_gas_vector = compute_density_air(m_void)

    W = m_void * porosity_powder * density_gas_vector + m_particle_tot * (1 - porosity_powder) * density_particle

    diff = (W - total_water) * (W - total_water)
    # diff = W * W - total_water * total_water
    # diff = np.abs(W - total_water)
    error = diff.mean(axis=0)
    return error

