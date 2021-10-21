from TAM_input_parameters import *
one_over_29 = 1/29              # computing-time-saver
one_over_18 = 1/18

def create_initial_conditions(relative_humidity_gas_inlet, amorphous_material_initial):
    p_saturated = compute_p_saturated_vector(temp_initial)

    # Surroundings, RH 58
    water_activity_sur = relative_humidity_gas_inlet
    m_gas_sur = compute_H_from_water_activity_temp(water_activity_sur)
    partial_pressure_moisture = compute_pressure_water(temp_initial_celsius)
    vapor_pressure_sur = water_activity_sur * partial_pressure_moisture
    m_particle_am_sat = compute_GAB_equilibrium_moisture_am(relative_humidity_gas_inlet)
    m_particle_cryst_sat = compute_GAB_equilibrium_moisture_cryst(relative_humidity_gas_inlet)
    m_particle_tot_sat = m_particle_cryst_sat * (
                1 - amorphous_material_initial) + m_particle_am_sat * amorphous_material_initial

    # Initial system, RH 20
    molar_concentration_void_initial = compute_molar_concentration_vector(relative_humidity_bed_initial, p_saturated,
                                                                          temp_initial)
    vapor_pressure_void_initial = compute_partial_pressure_moisture_vector(molar_concentration_void_initial,
                                                                           temp_initial)
    water_activity_void_initial = vapor_pressure_void_initial / p_saturated
    m_void_initial = compute_H_from_water_activity_temp(water_activity_void_initial)
    m_particle_am_initial = compute_GAB_equilibrium_moisture_am(water_activity_void_initial)
    m_particle_cryst_initial = compute_GAB_equilibrium_moisture_cryst(water_activity_void_initial)
    m_particle_tot_initial = m_particle_cryst_initial * (
                1 - amorphous_material_initial) + m_particle_am_initial * amorphous_material_initial

    gas_density_initial = compute_density_air(m_void_initial)
    gas_density_sur = compute_density_air(m_gas_sur)

    # CONCENTRATIONS

    # Surroundings, RH 58
    m_gas_conc_sat = m_gas_sur * porosity_powder * gas_density_sur  # kg/m3

    # Initial system, RH 20
    m_void_conc_initial = m_void_initial * porosity_powder * gas_density_initial
    m_particle_cryst_conc_initial = m_particle_cryst_initial * (1 - porosity_powder) * density_particle
    m_particle_am_conc_initial = m_particle_am_initial * (1 - porosity_powder) * density_particle
    m_particle_tot_conc_initial = m_particle_cryst_conc_initial * (
                1 - amorphous_material_initial) + m_particle_am_conc_initial * amorphous_material_initial

    # System
    system_conc_initial = m_particle_cryst_conc_initial * (
                1 - amorphous_material_initial) + m_particle_am_conc_initial * amorphous_material_initial + m_void_conc_initial

    # ENERGY
    enthalpy_gas_initial = compute_enthalpy_humid_air(temp_initial_celsius, m_void_initial)  # J/kg
    enthalpy_sur_air_initial = compute_enthalpy_humid_air(temp_initial_celsius, m_void_initial)  # J/kg
    enthalpy_sur_vapor = compute_enthalpy_vapor(temp_initial_celsius)
    energy_sur_vapor = enthalpy_sur_vapor * m_void_initial * gas_density_initial
    enthalpy_particle_initial = heat_capacity_particle * temp_initial_celsius + m_particle_tot_initial * heat_capacity_water * temp_initial_celsius
    enthalpy_powder_initial = enthalpy_gas_initial * porosity_powder * gas_density_initial + enthalpy_particle_initial * (
                1 - porosity_powder) * density_particle

    # Saturated system
    m_particle_cryst_conc_sat = m_particle_cryst_sat * (1 - porosity_powder) * density_particle
    m_particle_am_conc_sat = m_particle_am_sat * (1 - porosity_powder) * density_particle
    m_particle_tot_conc_sat = m_particle_cryst_conc_sat * (
                1 - amorphous_material_initial) + m_particle_am_conc_sat * amorphous_material_initial
    system_conc_sat = m_particle_cryst_conc_sat * (
                1 - amorphous_material_initial) + m_particle_am_conc_sat * amorphous_material_initial + m_gas_conc_sat

    return (m_particle_tot_conc_initial, m_void_conc_initial, enthalpy_powder_initial, water_activity_void_initial,
            m_void_initial, temp_initial_celsius, m_gas_sur, m_particle_cryst_initial, m_particle_am_initial,
            m_particle_tot_initial, gas_density_initial, water_activity_sur, m_particle_cryst_conc_initial)


def compute_glass_temp_mix(w1_lactose, t_g_lactose, t_g2):
  w2 = 1 - w1_lactose
  glass_temp_mix = (w1_lactose * t_g_lactose + w2 * t_g2 * k_param)/(w1_lactose + w2 * k_param)
  return(glass_temp_mix)


def compute_kinetics_avrami(am_amorph, temp_diff, single = False):
    # Y = fraction unacomplished change; fraction amorph wanting to turn (ex 0.2)
    # 1-Y = change already done; fraction amorph turned crystalline (ex 0.8)    d(1-Y)/dt = - dY/dt
    if single == False:
        n_steps = am_amorph.size
        am_amorph[am_amorph >= 1] = 0.9999
        am_amorph[am_amorph < 0] = 0
    else:
        n_steps = 1

    change_amorphicity = np.zeros(n_steps)
    for n in range(n_steps):
        if single == False:
            am_am = am_amorph[n]                                                        # Y
            temp_diff_curr = temp_diff[n]
        else:
            am_am = am_amorph

        if am_am == 0:
            change_amorphicity[n] = 0
        else:
            K_exp = -c_1 / (r_gas_constant * (c_2 + temp_diff_curr))
            K = c_3 * (np.exp(K_exp)) ** 3
            change_amorphicity_1 = n_A * K * am_am
            change_amorphicity_2 = (-np.log(am_am) / K) ** avrami_exponent
            change_amorphicity[n] = change_amorphicity_1 * change_amorphicity_2     # d(1-Y)/dt, positive

    change_amorphicity = -change_amorphicity                                        # negative
    change_amorphicity[change_amorphicity > 0] = 0
    change_amorphicity[change_amorphicity < -1] = -1
    if single:
        return change_amorphicity, K
    else:
        return change_amorphicity


def compute_laplacian(array, boundary, step_length, double_bc=False, inflow=False):
    n_steps = array.size

    # if the BC are different at two sides
    if double_bc:
        matrix = np.insert(array, 0, boundary[0])           # Add boundary conditions to array at 0
        matrix = np.append(matrix, boundary[1])             # Add second boundary conditions at end

    # if they are the same
    else:
        if n_steps == 1:
            matrix = np.array([boundary, array, array])
        else:
            matrix = np.insert(array, 0, boundary)           # Add boundary conditions to array at 0
            matrix = np.append(matrix, array[-1])            # Add final conditions to array again at end

    # if computing what flows into powder bed from air above it; hindered diffusion through bed, but not in air above it
    if inflow:
        top = (matrix[0:n_steps] - 2 * matrix[1:1 + n_steps] + matrix[2:]) * moisture_diffusivity
        top[-1] = matrix[-3] * moisture_diffusivity - matrix[-2] * (moisture_diffusivity + diffusivity_eff) + matrix[-1] * diffusivity_eff

    else:
        top = matrix[0:n_steps] - 2 * matrix[1:1 + n_steps] + matrix[2:]

    laplacian = top / (step_length * step_length)
    return laplacian


############################### PRESSURES, TEMPS & CONCENTRATIONS ######################################################
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


def compute_enthalpy_humid_air(temp_celsius, m_void):
    enthalpy = heat_capacity_air * temp_celsius + m_void * (heat_capacity_vapor * temp_celsius + heat_of_evaporation_water) # kg/m3
    return enthalpy


def compute_enthalpy_vapor(temp_celsius):
    enthalpy = heat_capacity_vapor * temp_celsius + heat_of_evaporation_water # kg/m3
    return enthalpy


############################### MOISTURE ISOTHERMS AND EQUILIBRIUMS ####################################################
def compute_GAB_equilibrium_moisture_am(water_activity, limit_aw = True):
    water_activity_max = 0.6
    if limit_aw == False:
        water_activity_max = 1
    water_activity = np.where(water_activity > water_activity_max, water_activity_max, water_activity)

    top = M0_am * c_am * f_am * water_activity
    bottom = (1 - f_am * water_activity) * (1 - (1 - c_am) * f_am * water_activity)
    equilibrium_moisture_particle_vector = top / bottom

    equilibrium_moisture_particle_vector = np.where(equilibrium_moisture_particle_vector < 0, 0, equilibrium_moisture_particle_vector)
    return equilibrium_moisture_particle_vector


def compute_GAB_equilibrium_moisture_cryst(water_activity, H_H = True):
    water_activity_limit = 0.85

    # choose whether to use simple or complex version; complex takes MUCH longer
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


def compute_H_from_water_activity_temp(water_activity):
    pressure_water = compute_pressure_water(temp_initial_celsius)
    H = 18 * water_activity * pressure_water/(29 * (pressure_ambient - water_activity * pressure_water))
    return H


def compute_density_air(H):
    top = kelvin
    bottom = 22.4 * (kelvin + temp_initial_celsius) * (one_over_29 + H * one_over_18)
    density = top/bottom
    return density


def compute_H_and_M_iteratively(water_activity, total_water, amount_am):
    m_void = compute_H_from_water_activity_temp(water_activity)
    m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_am = compute_GAB_equilibrium_moisture_am(water_activity)
    m_particle_tot = m_am * amount_am + m_cryst * (1 - amount_am)
    density_gas_vector = compute_density_air(m_void)

    # compute the total water based on the given water activity
    W = m_void * porosity_powder * density_gas_vector + m_particle_tot * (1 - porosity_powder) * density_particle

    # compare computed W with known total water, and return the squared error (no square root to save computing time)
    diff = (W - total_water) * (W - total_water)
    error = diff.mean(axis=0)
    return error


###################################### NOT USED IN FINAL VERSION #######################################################
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

def compute_air_temp_c_from_H(H, enthalpy):
    temp_c = (enthalpy - H * heat_of_evaporation_water)/(heat_capacity_air + H * heat_capacity_vapor)
    return temp_c