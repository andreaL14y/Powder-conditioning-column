from scipy.integrate import odeint
from input_parameters import *
from v_ode_functions_2 import compute_partial_pressure_moisture_vector, compute_molar_concentration_vector, compute_p_saturated_vector
from moisture_content_equilibrium import compute_GAB_equilibrium_moisture_cryst, derivative_m_isotherm_cryst, derivative_m_isotherm_am
from glass_transition_curve import compute_glass_temp_mix

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

######################################## FUNCTIONS #####################################################################
def compute_diffusion_laplacian(m_void):
    m_diff = m_gas_sur - m_void                                # fraction
    step = 0.001                                                    # 1mm step...?
    step = 0.01
    laplacian = m_diff * gas_density/step**2                                     # kg/m5
    return laplacian


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
    return change_amorphicity


def compute_water_activity_from_m_void(m_void):
    vapor_pressure = 29 * m_void * total_pressure / (18 + 29 * m_void)
    water_activity = vapor_pressure / p_saturated
    return water_activity


def compute_system(initials, t, cryst=True):
    global counter
    counter += 1

    m_void_conc, amount_am = initials
    m_void = m_void_conc/(gas_density * porosity_powder)                                              # fraction
    water_activity = compute_water_activity_from_m_void(m_void)

    if cryst:
        m_particle = compute_GAB_equilibrium_moisture_cryst(water_activity)
    else:
        m_particle = compute_GAB_equilibrium_moisture_am(water_activity)

    m_particle_conc = m_particle * (1-porosity_powder) * particle_density

    glass_temp = compute_glass_temp_mix(1 - m_particle, glass_temp_lactose, glass_temp_water_1)
    temp_diff = temp_initial - glass_temp
    crystallinity = 1 - amount_am                                   # 0 in first step
    if cryst:
        change_amorph = 0
    else:
        change_amorph = compute_kinetics_avrami(amount_am, temp_diff)

    if counter % 100 == 0 and not cryst:
        print('\nCRYST at time', int(t))
        print('Water activity:', water_activity)
        print('m_particle', m_particle)
        print('m_particle_conc', m_particle_conc)
        print('m_void', m_void)
        print('change am', change_amorph)
    laplacian_conc = compute_diffusion_laplacian(m_void)
    diffusion = laplacian_conc * moisture_diffusivity * porosity_powder * 0.5                # kg/m5 * m2/s = kg/(m3 s)

    if cryst:
        derivative = derivative_m_isotherm_cryst(water_activity)
    else:
        derivative = derivative_m_isotherm_am(water_activity)

    derivative *= particle_density * (1-porosity_powder)                                    # Concentration

    top = diffusion
    bottom = porosity_powder + (1-porosity_powder) * derivative
    m_void_change = top/bottom
    system_change = np.array([m_void_change, change_amorph])
    return system_change


def normalize_data(data):
    min = np.min(data)
    max = np.max(data)
    normalized_data = (data - min)/(max - min)
    return normalized_data


################################ PARAMETERS & INITIAL CONDITIONS #######################################################
p_saturated                     = compute_p_saturated_vector(temp_initial)
total_pressure                  = pressure_ambient

# Initial conditions surroundings
molar_conc_sur                  = compute_molar_concentration_vector(relative_humidity_gas_inlet, p_saturated, temp_initial)
vapor_pressure_sur              = compute_partial_pressure_moisture_vector(molar_conc_sur, temp_initial)
water_activity_sur              = vapor_pressure_sur / p_saturated
m_gas_sur                       = 18 * vapor_pressure_sur/(29 * (total_pressure - vapor_pressure_sur))

# Initial gas moisture at 0.2 RH
molar_concentration_void_initial = compute_molar_concentration_vector(relative_humidity_bed_initial, p_saturated, temp_initial)
vapor_pressure_void_initial     = compute_partial_pressure_moisture_vector(molar_concentration_void_initial, temp_initial)
water_activity_void_initial     = vapor_pressure_void_initial / p_saturated
m_gas_initial            = 18 * vapor_pressure_void_initial/(29 * (total_pressure - vapor_pressure_void_initial))
m_gas_conc_initial       = m_gas_initial * porosity_powder * gas_density
m_gas_conc_sat           = m_gas_sur * porosity_powder * gas_density

# Initial moisture cryst powder at start
m_powder_cryst_initial         = compute_GAB_equilibrium_moisture_cryst(water_activity_void_initial)
m_powder_cryst_conc_initial    = m_powder_cryst_initial * (1 - porosity_powder) * particle_density
m_powder_am_initial         = compute_GAB_equilibrium_moisture_am(water_activity_void_initial)
m_powder_am_conc_initial    = m_powder_am_initial * (1 - porosity_powder) * particle_density

# Saturated system
m_particle_cryst_sat           = compute_GAB_equilibrium_moisture_cryst(relative_humidity_gas_inlet)
m_particle_cryst_conc_sat      = m_particle_cryst_sat * (1 - porosity_powder) * particle_density
m_particle_am_sat           = compute_GAB_equilibrium_moisture_am(relative_humidity_gas_inlet)
m_particle_am_conc_sat      = m_particle_am_sat * (1 - porosity_powder) * particle_density

system_conc_sat             = m_particle_cryst_conc_sat * (1-amorphous_material_initial) + \
                              m_particle_am_conc_sat * amorphous_material_initial + m_gas_conc_sat

# Computing times
hours = 2
seconds = hours * 60 * 60
resolution = 100
time_step = seconds/resolution                      # each step in time is this many seconds, s
# resolution = int(seconds/10)
discrete_time = np.linspace(0, seconds, resolution)

# Amorphous material
initial_amorphicity = 1
initial_amorphicity = amorphous_material_initial
initials = np.array([m_gas_conc_initial, initial_amorphicity])

tabs = 50
print('\n')
print('############################################# CONDITIONS #############################################')
print('Water activity surroundings:'.ljust(tabs), f'{water_activity_sur:.2f}')
print('Moisture in the surrounding as fraction:'.ljust(tabs), f' {m_gas_sur:.5f}')
print('Initial moisture gas as fraction:'.ljust(tabs), f'{m_gas_initial:.5f}')
print('Initial moisture gas as conc:'.ljust(tabs), f'{m_gas_conc_initial:.5f} kg/m3')

print('Initial moisture cryst powder as fraction:'.ljust(tabs), f'{m_powder_cryst_initial:.5f}')
print('Initial moisture cryst powder as conc:'.ljust(tabs), f'{m_powder_cryst_conc_initial:.5f} kg/m3')

print('Initial moisture cryst powder as fraction:'.ljust(tabs), f'{m_powder_am_initial:.5f}')
print('Initial moisture cryst powder as conc:'.ljust(tabs), f'{m_powder_am_conc_initial:.5f} kg/m3')

print('Total water content at beginning:'.ljust(tabs), f'{system_conc_sat:.5f}')
print('Max time:'.ljust(tabs), hours, 'hours')
print('\n')
print('############################################# START COMPUTATION #############################################')

counter = 0
computed_system_cryst = odeint(compute_system, initials, discrete_time, mxstep=500)

counter = 0
computed_system_am = odeint(compute_system, initials, discrete_time, args=(False, ), mxstep=500)

print('############################################# COMPUTATION COMPLETE #############################################')

m_void_conc_cryst   = computed_system_cryst[:, 0]
m_void_conc_am      = computed_system_am[:, 0]
amorphicity         = computed_system_am[:, 1]

m_void_cryst                = m_void_conc_cryst/(gas_density * (1-porosity_powder))
water_activity_void_cryst   = compute_water_activity_from_m_void(m_void_cryst)
m_particle_cryst_vector     = compute_GAB_equilibrium_moisture_cryst(water_activity_void_cryst)

m_void_am_vector                = m_void_conc_am / (gas_density * (1 - porosity_powder))
water_activity_void_am_vector   = compute_water_activity_from_m_void(m_void_am_vector)
# m_particle_am_vector           = compute_GAB_equilibrium_moisture_cryst(water_activity_void_am_vector)
m_particle_am_vector            = compute_GAB_equilibrium_moisture_am(water_activity_void_am_vector)

m_powder_total = m_particle_am_vector * amorphicity + m_particle_cryst_vector * (1-amorphicity)

print('\n')
# print('############################################# COMPUTATION COMPLETE #############################################')
tg_s = compute_glass_temp_mix(1 - m_particle_am_vector, glass_temp_lactose, glass_temp_water_1)
t_tg_diff = temp_initial - tg_s

m_powder_diffs = m_particle_am_vector[1:] - m_particle_am_vector[:-1]
m_powder_diffs = m_powder_total[1:] - m_powder_total[:-1]
m_powder_diffs = np.insert(m_powder_diffs, 0, 0, axis=0)
m_powder_diffs[m_powder_diffs < 0] = 0

am_am_diffs = amorphicity[1:] - amorphicity[:-1]
am_am_diffs = np.insert(am_am_diffs, 0, 0, axis=0)
am_cryst_diffs = - am_am_diffs

heat_of_sorption_vector = m_powder_diffs * (heat_of_sorption/1000) / time_step        # J/kg to J/g
heat_of_sorption_vector = m_powder_total * heat_of_sorption/1000         # J/kg to J/g
heat_of_cryst_vector = am_cryst_diffs * (heat_of_crystallization/1000) / time_step    # J/kg to J/g

print(m_powder_diffs[0:20], '\n')
print(heat_of_sorption_vector[0:20])

total_energy_vector = heat_of_cryst_vector + heat_of_sorption_vector

amorphicity_norm = normalize_data(amorphicity)

########################################## PLOT ########################################################################
fig, axs = plt.subplots(3, 3, figsize=(25, 15))
time = discrete_time/3600

grid_color = 'slategray'
temp_color = 'crimson'
m_color = 'navy'
am_color = 'forestgreen'
sat_color = 'gold'
x_min = -0.1
x_max = hours + 0.1

fig.suptitle(f'Moisture sorption TAM at RH {relative_humidity_gas_inlet} and T {temp_initial-kelvin} C')
axs[0, 0].plot(time, m_void_conc_cryst, label='Crystalline', color='lightsteelblue')
axs[0, 0].plot(time, m_void_conc_am, label='Amorhpous', color=m_color)
axs[0, 0].legend()
axs[0, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
axs[0, 0].set_title('Void moisture concentration')
axs[0, 0].set(xlabel='Time, hours', ylabel='Moisture concentration powder void, kg/m3')

# axs[2, 0].hlines(m_particle_am_sat, -1, hours+1, label='Saturated am', color=sat_color, linestyles='--')
axs[2, 0].hlines(m_particle_cryst_sat, -1, hours+1, label='Saturated cryst', color=sat_color, linestyles='--')
axs[2, 0].plot(time, m_powder_total, label='Total moisture', color=m_color)
axs[2, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
axs[2, 0].set_xlim(x_min, x_max)
axs[2, 0].set_title('Total moisture content powder')
axs[2, 0].legend()
axs[2, 0].set(xlabel='Time, hours', ylabel='Moisture content, kg/kg')

axs[0, 1].hlines(0, -1, hours+1, color=sat_color, linestyles='--')
axs[0, 1].plot(time, t_tg_diff, label='T-Tg', color=temp_color)
axs[0, 1].legend()
axs[0, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
axs[0, 1].set_xlim(x_min, x_max)
axs[0, 1].set_title('Glass transition temp')
axs[0, 1].set(xlabel='Time, hours', ylabel='Opeating temp - glass temp')

axs[1, 0].hlines(m_particle_am_sat, -1, hours+1, label='Saturated moisture content', color=sat_color, linestyles='--')
axs[1, 0].plot(time, m_particle_am_vector, label='Moisture content amorphous', color=m_color)
axs[1, 0].legend()
axs[1, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
axs[1, 0].set_xlim(x_min, x_max)
axs[1, 0].set_title('Moisture content amorphous')
axs[1, 0].set(xlabel='Time, hours', ylabel='Moisture, g/g lactose')


axs[1, 1].plot(time, amorphicity_norm, label='Amorphicity', color=am_color)
axs[1, 1].legend()
axs[1, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
axs[1, 1].set_xlim(x_min, x_max)
axs[1, 1].set_ylim(-0.1, 1.1)
axs[1, 1].set_title('Amount amorhpous material')
axs[1, 1].set(xlabel='Time, hours', ylabel='Amorhpous material, g/g lactose')

axs[0, 2].plot(time, heat_of_sorption_vector, label='Heat of sorption', color=temp_color)
axs[0, 2].legend()
axs[0, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
# axs[1, 1].set_ylim(0, 0.9)
axs[0, 2].set_title('Heat due to sorption')
axs[0, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')

axs[1, 2].plot(time, heat_of_cryst_vector, label='Heat of crystallization', color=temp_color)
axs[1, 2].legend()
axs[1, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
# axs[1, 1].set_ylim(0, 0.9)
axs[1, 2].set_title('Heat due to sorption')
axs[1, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')

axs[2, 2].plot(time, total_energy_vector, label='Total energy', color=temp_color)
axs[2, 2].legend()
axs[2, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
# axs[1, 1].set_ylim(0, 0.9)
axs[2, 2].set_title('Total energy')
axs[2, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')

plt.show()

