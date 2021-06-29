from initial_conditions import *

######################################## FUNCTIONS #####################################################################
def compute_system(initials, t):
    global counter
    counter += 1

    m_void_conc, amount_am = initials
    m_void = m_void_conc/(gas_density * porosity_powder)                                              # fraction
    water_activity = compute_water_activity_from_m_void(m_void, p_saturated)

    m_particle_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_particle_am = compute_GAB_equilibrium_moisture_am(water_activity)

    # if cryst:
    #     m_particle = m_particle_cryst
    # else:
    #     m_particle = m_particle_am
    m_particle_total = amount_am * m_particle_am + (1-amount_am) * m_particle_cryst
    m_particle_conc = m_particle_total * (1-porosity_powder) * particle_density

    glass_temp = compute_glass_temp_mix(1 - m_particle_am, glass_temp_lactose, glass_temp_water_1)
    temp_diff = temp_initial - glass_temp

    change_amorph = compute_kinetics_avrami(amount_am, temp_diff)

    if counter % 100 == 0:
        print('\nAT TIME', int(t))
        print('Water activity:', water_activity)
        # print('m_particle', m_particle)
        print('m_particle_conc', m_particle_conc)
        print('m_void', m_void)
        print('change am', change_amorph)

    laplacian_conc = compute_diffusion_laplacian(m_void, m_gas_sur)#
    diffusion = laplacian_conc * moisture_diffusivity * porosity_powder * 0.5                # kg/m5 * m2/s = kg/(m3 s)

    derivative_cryst = compute_derivative_m_GAB(water_activity)
    derivative_am = compute_derivative_m_GAB(water_activity, cryst=False)

    derivative_total = amount_am * derivative_am + (1-amount_am * derivative_cryst)
    derivative_total *= particle_density * (1-porosity_powder)                                      # kg/m3

    m_void_change_cr = -change_amorph * (m_particle_am - m_particle_cryst_sat)                      # fraction
    m_void_change_cr *= particle_density * (1 - porosity_powder) / (gas_density * porosity_powder)  # kg/m3
    m_void_change_cr /= 100

    top = diffusion + m_void_change_cr
    # bottom = porosity_powder + (1-porosity_powder) * derivative
    bottom = porosity_powder + (1-porosity_powder) * derivative_total
    m_void_change_sorption = top/bottom

    print('\nm void change')
    print('Diffusion'.ljust(tabs), diffusion)
    print('Cryst'.ljust(tabs), m_void_change_cr)

    m_void_change = m_void_change_sorption
    system_change = np.array([m_void_change, change_amorph])
    return system_change


def plot_tam():
    fig, axs = plt.subplots(3, 3, figsize=(25, 15))
    time = discrete_time / 3600

    grid_color = 'slategray'
    temp_color = 'crimson'
    m_color = 'navy'
    am_color = 'forestgreen'
    sat_color = 'gold'
    x_min = -0.1
    x_min = 0
    x_max = hours + 0.1
    # x_max = 0.02
    fig.suptitle(f'Moisture sorption TAM at RH {relative_humidity_gas_inlet} and T {temp_initial-kelvin} C')
    # axs[0, 0].plot(time, m_void_conc_cryst, label='Crystalline', color='lightsteelblue', linestyle='--')
    # axs[0, 0].plot(time, m_void_conc_am, label='Amorhpous', color=m_color, linestyle='--')
    axs[0, 0].plot(time, m_void_conc_total, label='Total', color=m_color)
    axs[0, 0].legend()
    axs[0, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 0].set_xlim(x_min, x_max)
    axs[0, 0].set_title('Void moisture concentration')
    axs[0, 0].set(xlabel='Time, hours', ylabel='Moisture concentration powder void, kg/m3')

    axs[1, 0].hlines(m_particle_am_sat, -1, hours+1, label='Saturated am', color=sat_color, linestyles='--')
    axs[1, 0].plot(time, m_particle_am_vector, label='Amorphous', color=m_color)
    axs[1, 0].legend()
    axs[1, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 0].set_xlim(x_min, x_max)
    axs[1, 0].set_title('Moisture content amorphous')
    axs[1, 0].set(xlabel='Time, hours', ylabel='Moisture, g/g lactose')

    # axs[2, 0].hlines(m_particle_am_sat, -1, hours+1, label='Saturated am', color=sat_color, linestyles='--')
    axs[2, 0].hlines(m_particle_cryst_sat, -1, hours+1, label='Saturated cryst', color=sat_color, linestyles='--')
    axs[2, 0].plot(time, m_powder_total, label='Total moisture', color=m_color)
    # axs[2, 0].plot(time, m_particle_cryst_vector, label='Crystalline', color='lightsteelblue')
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

    axs[1, 1].plot(time, amorphicity_norm, label='Amorphicity', color=am_color)
    axs[1, 1].legend()
    axs[1, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 1].set_xlim(x_min, x_max)
    axs[1, 1].set_ylim(-0.1, 1.1)
    axs[1, 1].set_title('Amount amorhpous material')
    axs[1, 1].set(xlabel='Time, hours', ylabel='Amorhpous material, g/g lactose')

    # axs[2, 1].hlines(m_particle_cryst_sat, -1, hours+1, label='Saturated cryst', color=sat_color, linestyles='--')
    # axs[2, 1].plot(time, m_particle_cryst_vector, label='Crystalline', color=m_color)
    axs[2, 1].plot(time, m_powder_diffs, label='Diff moisture', color=m_color)
    axs[2, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 1].set_xlim(x_min, x_max)
    axs[2, 1].set_title('Total moisture content powder')
    axs[2, 1].set_xlim(x_min, x_max)
    axs[2, 1].legend()
    axs[2, 1].set(xlabel='Time, hours', ylabel='Moisture content, kg/kg')

    axs[0, 2].plot(time, heat_of_sorption_vector, label='Heat of sorption', color=temp_color)
    axs[0, 2].legend()
    axs[0, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    # axs[1, 1].set_ylim(0, 0.9)
    axs[0, 2].set_xlim(x_min, x_max)
    axs[0, 2].set_title('Sorption')
    axs[0, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')

    axs[1, 2].plot(time, heat_of_cryst_vector, label='Crystallization', color=temp_color)
    axs[1, 2].legend()
    axs[1, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    # axs[1, 1].set_ylim(0, 0.9)
    axs[1, 2].set_xlim(x_min, x_max)
    axs[1, 2].set_title('Heat due to crystallization')
    axs[1, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')

    axs[2, 2].plot(time, total_energy_vector, label='Total energy', color=temp_color)
    axs[2, 2].legend()
    axs[2, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    # axs[1, 1].set_ylim(0, 0.9)
    axs[2, 2].set_xlim(x_min, x_max)
    axs[2, 2].set_title('Total energy')
    axs[2, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')

    plt.show()

################################ PARAMETERS & INITIAL CONDITIONS #######################################################

# Computing times
hours = 2
seconds = hours * 60 * 60
resolution = 10000
time_step = seconds/resolution                      # each step in time is this many seconds, s
discrete_time = np.linspace(0, seconds, resolution)

# Amorphous material
# initial_amorphicity = 1
initial_amorphicity = amorphous_material_initial
initials = np.array([m_gas_conc_initial, initial_amorphicity])      # kg/m3, kg/kg

tabs = 50
def print_info():
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
    print('Time step:'.ljust(tabs), time_step, 's')
    print('\n')
print_info()

print('############################################# START COMPUTATION #############################################')
counter = 0
computed_system = odeint(compute_system, initials, discrete_time, mxstep=500)

print('############################################# COMPUTATION COMPLETE #############################################')
m_void_conc_total           = computed_system[:, 0]
amorphicity                 = computed_system[:, 1]

m_void_tot_vector           = m_void_conc_total/(gas_density * porosity_powder)
water_activity_void_vector  = compute_water_activity_from_m_void(m_void_tot_vector, p_saturated)
m_particle_cryst_vector     = compute_GAB_equilibrium_moisture_cryst(water_activity_void_vector)
m_particle_am_vector        = compute_GAB_equilibrium_moisture_am(water_activity_void_vector)

m_powder_total = m_particle_am_vector * amorphicity + m_particle_cryst_vector * (1-amorphicity)

tg_s = compute_glass_temp_mix(1 - m_particle_am_vector, glass_temp_lactose, glass_temp_water_1)
t_tg_diff = temp_initial - tg_s

m_powder_diffs = m_powder_total[1:] - m_powder_total[:-1]                   # fraction, kg water/kg dry
m_powder_diffs = np.insert(m_powder_diffs, 0, 0, axis=0)
m_powder_diffs[m_powder_diffs < 0] = 0

am_am_diffs = amorphicity[1:] - amorphicity[:-1]
am_am_diffs = np.insert(am_am_diffs, 0, 0, axis=0)
am_cryst_diffs = - am_am_diffs

heat_of_sorption = compute_heat_of_sorption(water_activity_void_vector, temp_initial)
heat_of_sorption_vector = m_powder_diffs * (heat_of_sorption/1000) / time_step        # J/kg to J/g

heat_of_cryst_vector = am_cryst_diffs * (heat_of_crystallization/1000) / time_step    # J/kg to J/g

total_energy_vector = heat_of_cryst_vector + heat_of_sorption_vector

amorphicity_norm = normalize_data(amorphicity)

########################################## PLOT ########################################################################
plot_tam()

