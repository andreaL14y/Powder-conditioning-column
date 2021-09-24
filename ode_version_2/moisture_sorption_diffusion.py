from initial_conditions import *
import time


######################################## FUNCTIONS #####################################################################
def compute_system(initials, t, n_space_steps, step_length, verbose=False, time_continue=False):
    global counter
    global while_time
    global water_activity_just_used
    divider = 25000

    total_water         = initials[0:n_space_steps]
    amount_am           = initials[n_space_steps:n_space_steps * 2]
    total_energy        = initials[n_space_steps * 2:n_space_steps * 3]
    water_activity_prev = initials[n_space_steps * 3:n_space_steps * 4]
    m_sur               = initials[n_space_steps * 4:n_space_steps * 5]
    temp_celsius_prev   = initials[n_space_steps * 5:n_space_steps * 6]         #TODO: uncomment
    # energy_sur_air    = initials[n_space_steps * 5:n_space_steps * 6]             #TODO: comment

    temp_celsius_sur = compute_air_temp_c_from_H(m_sur, enthalpy_sur_air_initial)
    water_activity_sur_vector = compute_water_activity_from_m_void(m_sur, temp_initial_celsius)
    amount_am                 = np.where(amount_am > amorphous_material_initial, amorphous_material_initial, amount_am)

    counter += 1
    ########################################### MOISTURE ###############################################################
    while_time_start = time.time()

    optimized_aw = scipy.optimize.minimize(compute_H_and_M_iteratively, water_activity_just_used, method='SLSQP',
                                           args=(total_water, amount_am),
                                           bounds=((0, 1),) * n_space_steps, options={'ftol': 1e-07})   #TODO: decrease!

    water_activity          = optimized_aw.x
    water_activity_change   = water_activity - water_activity_prev
    m_void                  = compute_H_from_water_activity_temp(water_activity)

    while_time += time.time() - while_time_start

    m_particle_cryst    = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_particle_am       = compute_GAB_equilibrium_moisture_am(water_activity)
    m_particle_tot      = m_particle_cryst * (1-amount_am) + m_particle_am * amount_am

    m_particle_cryst_prev = compute_GAB_equilibrium_moisture_cryst(water_activity_prev)
    m_particle_am_prev = compute_GAB_equilibrium_moisture_am(water_activity_prev)
    m_particle_tot_prev = m_particle_cryst_prev * (1 - amount_am) + m_particle_am_prev * amount_am      #TODO: am am should be prev as well
    ########################################### ENERGY #################################################################
    gas_density     = compute_density_air(m_void)
    water_W         = porosity_powder * m_void * gas_density + (1 - porosity_powder) * m_particle_tot * density_particle
    excess_water    = total_water - water_W
    excess_water    = np.where((excess_water > 0) & (water_activity == 1), excess_water, 0)

    temp_kelvin, temp_celsius = compute_temp_from_energy(m_void, m_particle_tot, total_energy, excess_water, verbose=False)
    temp_change = temp_celsius - temp_celsius_prev
    # temp_kelvin, temp_celsius = np.zeros(n_space_steps) + temp_initial, np.zeros(n_space_steps) + temp_initial_celsius
    # temp_change = np.zeros(n_space_steps) + 0

    energy_vapor              =  m_void * compute_enthalpy_vapor(temp_celsius) * gas_density

    ######################################### AIR IN AMPOULE ###########################################################
    gas_density_sur = compute_density_air(m_sur)

    m_sur_bounds = np.array([m_gas_sur, m_void[0]])
    m_sur_laplacian, m_from_sur = compute_laplacian(m_sur, m_sur_bounds, step_length_air, double_bc=True, inflow=True)
    m_sur_change = m_sur_laplacian * moisture_diffusivity                           # kg/m5 * m2/s = kg/(m3 s)

    ################################## RESULTING CHANGES IN Q AND W ####################################################
    laplacian_conc = compute_laplacian(m_void * gas_density, m_sur[-1] * gas_density_sur[-1], step_length)
    m_diffusion = laplacian_conc * diffusivity_eff                          # kg/m5 * m2/s = kg/(m3 s)

    energy_vapor_sur = m_sur[-1] * compute_enthalpy_vapor(temp_initial_celsius) * gas_density_sur[-1]   #TODO:!
    energy_diff_laplacian = compute_laplacian(energy_vapor, energy_vapor_sur, step_length)      # energy_vapor[0]
    change_energy_diffusion = diffusivity_eff * energy_diff_laplacian       # m2/s * (J/kg * kg/m3)/m2 = J/(m3 s)

    temp_bounds = np.array([temp_initial_celsius, temp_initial_celsius])
    temp_laplacian = compute_laplacian(temp_celsius, temp_bounds, step_length, double_bc=True)

    change_energy_conduction = conductivity_tot * temp_laplacian

    ###################################### AMORPHICITY #################################################################
    glass_temp = compute_glass_temp_mix(1 - m_particle_am, glass_temp_lactose, glass_temp_water_1)
    temp_diff = temp_kelvin - glass_temp
    change_amorph = compute_kinetics_avrami(amount_am, temp_diff)
    change_moisture     = m_diffusion

    change_energy_crystallization = change_amorph * density_particle * (
                1 - porosity_powder) * heat_of_crystallization              # Negative, 1/s * kg/m3 * J/kg
    change_energy_crystallization = - change_energy_crystallization         # Positive

    change_energy_evap = change_moisture * heat_of_evaporation_water        # positive during sorption kg/m3 * J/kg = J/m3
    # m_diffs_p = m_particle_tot - m_particle_tot_prev                        # positive during sorption, then neg
    # change_energy_binding = m_diffs_p * (1 - porosity_powder) * density_particle * heat_binding_diffs    # positive during sorption, then neg

    change_energy = change_energy_conduction + change_energy_crystallization + change_energy_diffusion - change_energy_evap #+ change_energy_binding


    if verbose and counter % divider == 0:
        hours_print = int(t / 3600)
        minutes_print = int((t - hours_print * 3600) / 60)
        seconds_print = int(t - (minutes_print * 60 + hours_print * 3600))
        time_percentage = t/seconds
        print('')
        print(f'        {time_percentage*100:.1f} % completed: TIME', hours_print, 'h;', minutes_print, 'min;', seconds_print, 's')
        np.set_printoptions(formatter={'float': '{: .2f}'.format})
        print('T-Tg:'.ljust(tabs), temp_diff)
        print('Temp C:'.ljust(tabs), temp_celsius)
        print('Temp C sur:'.ljust(tabs), temp_celsius_sur)
        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        print('Total energy:'.ljust(tabs), total_energy)
        print('')

        print('m_void:'.ljust(tabs), m_void)
        # print('m_gas_sur:'.ljust(tabs), f'{m_gas_sur:.1e}')
        print('M am:'.ljust(tabs), m_particle_am)
        # print('m_diffs_p:'.ljust(tabs), m_diffs_p)
        print('Excess water:'.ljust(tabs), excess_water)

        np.set_printoptions(formatter={'float': '{: .3f}'.format})
        print('Water act:'.ljust(tabs), water_activity)
        print('Water act sur:'.ljust(tabs), water_activity_sur_vector)
        np.set_printoptions(formatter={'float': '{: .5f}'.format})
        print('Am am:'.ljust(tabs), amount_am)

        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        print('Change am:'.ljust(tabs), change_amorph)
        print('')
        #
        print('Change energy diff:'.ljust(tabs), change_energy_diffusion)
        print('Change energy cond:'.ljust(tabs), change_energy_conduction)
        print('Change energy cryst:'.ljust(tabs), change_energy_crystallization)
        # print('Change energy binding:'.ljust(tabs), change_energy_binding)
        # print('Change energy excess water:'.ljust(tabs), change_energy_excess_m)
        # print('Change energy evap:'.ljust(tabs), f'{change_energy_evap_sur:.1e}')
        print('Change energy tot:'.ljust(tabs), change_energy)

        print('N iterations in optimize:'.ljust(tabs), optimized_aw.nfev)
        print(f'        {time_percentage*100:.1f} % completed: TIME', hours_print, 'h;', minutes_print, 'min;', seconds_print, 's')
        print('')

    water_activity_just_used = water_activity

    system_change       = np.concatenate([change_moisture, change_amorph, change_energy,
                                          water_activity_change, m_sur_change, temp_change  #temp_sur_change, change_energy_sur
                                          ]).astype('float64')

    return system_change


def plot_tam():
    fig, axs = plt.subplots(3, 3, figsize=(25, 15))
    time = discrete_time / 3600

    grid_color = 'slategray'
    temp_color = 'crimson'
    m_color = 'navy'
    m_color_light = 'lightsteelblue'
    am_color = 'forestgreen'
    sat_color = 'gold'
    x_min = 0
    x_max = hours + 0.1
    fig.suptitle(f'Moisture sorption TAM at RH {relative_humidity_gas_inlet}, T {temp_initial-kelvin} C and {amorphous_material_initial * 100} % amorphous material initially')

    ############################################# m_void ###############################################################
    axs[0, 0].hlines(m_gas_sur, -1, hours+1, label='Sat sur', linestyles='--', color=sat_color)
    axs[0, 0].plot(time, m_void_avg[:], label='Avg', color=m_color)
    axs[0, 0].plot(time, m_void_all[:, 0], label='1', color=m_color_light, linestyle='--')
    # axs[0, 0].plot(time, m_void_all[:, int(n_space_steps / 3)], label=f'{int(n_space_steps/3)}', color='cyan', linestyle='--')
    axs[0, 0].plot(time, m_void_all[:, n_space_steps-1], label=f'{n_space_steps}', color='cyan', linestyle='--')
    axs[0, 0].legend()
    axs[0, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 0].set_xlim(x_min, x_max)
    axs[0, 0].set_title('Void moisture concentration')
    axs[0, 0].set(xlabel='Time, hours', ylabel='Moisture conc void, kg/m3')

    ############################################# m_am #################################################################
    axs[1, 0].hlines(m_particle_am_sat, -1, hours+1,                    label='Sat am', color=sat_color, linestyles='--')
    axs[1, 0].plot(time, m_particle_am_avg,                             label='Am', color=m_color)
    axs[1, 0].plot(time, m_particle_am_all[:, 0],                       label='1', color=m_color_light, linestyle='--')
    # axs[1, 0].plot(time, m_particle_am_all[:, int(n_space_steps / 3)],  label=f'{int(n_space_steps/3)}', color='cyan', linestyle='--')
    axs[1, 0].plot(time, m_particle_am_all[:, n_space_steps-1],  label=f'{n_space_steps}', color='cyan', linestyle='--')
    # axs[1, 0].plot(time, m_particle_am_all[:, int(2 * n_space_steps/3)], label=f'{int(2 * n_space_steps/3)}', color='darkcyan', linestyle='--')
    # axs[1, 0].plot(time, m_particle_am_all[:, -1], label=f'{n_space_steps-1}', color=m_color, linestyle='--')
    axs[1, 0].legend()
    axs[1, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 0].set_xlim(x_min, x_max)
    axs[1, 0].set_title('Moisture content amorphous')
    axs[1, 0].set(xlabel='Time, hours', ylabel='Moisture, kg/kg lactose')

    ########################################### Total Water W ##########################################################
    axs[2, 0].hlines(m_particle_cryst_sat, -1, hours+1, color=sat_color, linestyles='--') #label='Sat cryst, kg/m3',
    axs[2, 0].plot(time, total_water_avg,                           label='Avg', color=m_color)
    axs[2, 0].plot(time, total_water_all[:, 0],                     label=f'{1}', color=m_color_light, linestyle='--')
    # axs[2, 0].plot(time, total_water_all[:, int(n_space_steps/3)],  label=f'{int(n_space_steps/3)}', color='cyan', linestyle='--')
    axs[2, 0].plot(time, total_water_all[:, n_space_steps-1],  label=f'{n_space_steps}', color='cyan', linestyle='--')
    # axs[2, 0].plot(time, total_water_all[:, int(2 * n_space_steps/3)], label=f'{int(2 * n_space_steps/3)}', color='darkcyan', linestyle='--')
    # axs[2, 0].plot(time, total_water_all[:, -1], label=f'{n_space_steps-1}', color=m_color, linestyle='--')
    # axs[2, 0].plot(time, m_particle_cryst_avg, label='Cr', color='lightsteelblue')
    axs[2, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 0].set_xlim(x_min, x_max)
    axs[2, 0].set_title('Total moisture system')
    axs[2, 0].legend()
    axs[2, 0].set(xlabel='Time, hours', ylabel='Moisture content, kg/m3')

    ########################################### T-Tg ###################################################################
    axs[0, 1].hlines(0, -1, hours+1, color=sat_color, linestyles='--')
    axs[0, 1].plot(time, t_tg_diff_avg,         label='Avg', color=temp_color)
    axs[0, 1].plot(time, t_tg_diff_all[:, 0],   label='1', color='rosybrown', linestyle='--')
    axs[0, 1].plot(time, t_tg_diff_all[:, n_space_steps-1],   label=f'{n_space_steps}', color='darkgoldenrod', linestyle='--')
    axs[0, 1].legend()
    axs[0, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 1].set_xlim(x_min, x_max)
    axs[0, 1].set_title('Glass transition temp')
    axs[0, 1].set(xlabel='Time, hours', ylabel='Opeating temp - glass temp')

    ###################################### Amorphicity X ###############################################################
    axs[1, 1].plot(time, amorphicity_avg, label='Avg', color=am_color)
    axs[1, 1].plot(time, amorphicity_all[:, 0],                         label='1', color='greenyellow', linestyle='--')
    # axs[1, 1].plot(time, amorphicity_all[:, int(n_space_steps/3)],      label=f'{1+int(n_space_steps/3)}', color='green', linestyle='--')
    # axs[1, 1].plot(time, amorphicity_all[:, 2 * int(n_space_steps/3)],  label=f'{1+2*int(n_space_steps/3)}', color='limegreen', linestyle='--')
    axs[1, 1].plot(time, amorphicity_all[:, n_space_steps-1],           label=f'{n_space_steps}', color='darkgreen', linestyle='--')
    axs[1, 1].legend()
    axs[1, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 1].set_xlim(x_min, x_max)
    ylim_max = np.max(amorphicity_avg + 0.02)
    axs[1, 1].set_ylim(-0.02, ylim_max)
    axs[1, 1].set_title('Amount amorhpous material')
    axs[1, 1].set(xlabel='Time, hours', ylabel='Am material, kg/kg lactose')

    ################################### Water Activity RH ##############################################################
    axs[2, 1].plot(time, water_activity_avg, label='Avg w_a', color=am_color)
    axs[2, 1].plot(time, water_activity_all[:, 0], label='1', color='green', linestyle='--')
    # axs[2, 1].plot(time, water_activity_all[:, int(n_space_steps / 3)], label=f'{1+int(n_space_steps/3)}',
    #                color='lightgreen', linestyle='--')
    # axs[2, 1].plot(time, water_activity_all[:, 2 * int(n_space_steps / 3)], label=f'{1+2*int(n_space_steps/3)}',
    #                color='limegreen', linestyle='--')
    axs[2, 1].plot(time, water_activity_all[:, n_space_steps - 1], label=f'{n_space_steps}', color='darkgreen',
                   linestyle='--')
    axs[2, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 1].set_xlim(x_min, x_max)
    axs[2, 1].set_title('Water activity over time')
    axs[2, 1].set(xlabel='Time, hours', ylabel='wa')
    axs[2, 1].legend()

    ################################### Heat of sorption ###############################################################
    # axs[0, 2].plot(time, heat_of_sorption_vector, label='Heat of sorption', color='yellowgreen', linestyle='--')
    # axs[0, 2].plot(time, heat_of_cryst_vector, label='Crystallization', color='darkgoldenrod', linestyle='--')
    # axs[0, 2].plot(time, heat_of_cryst_vector + heat_of_sorption_vector, label='Total', color='darkcyan') #, linestyle='--')

    axs[0, 2].plot(time, energy_avg, label='Total', color=temp_color)
    axs[0, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 2].set_xlim(x_min, x_max)
    axs[0, 2].set_title('Energy')
    axs[0, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')
    axs[0, 2].legend()

    ################################# Heat of crystallization ##########################################################
    axs[1, 2].plot(time, heat_of_cryst_vector, label='Crystallization', color='darkgoldenrod') #, linestyle='--')
    # axs[0, 2].plot(time, total_energy_vector, label='Total', color=temp_color)
    axs[1, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 2].set_xlim(x_min, x_max)
    axs[1, 2].set_title('Heat generated')
    axs[1, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')
    axs[1, 2].legend()

    ####################################### Temperature ################################################################
    # max_cryst = np.max(heat_of_cryst_vector)
    # axs[2, 2].set_ylim(0, max_cryst)
    axs[2, 2].plot(time, temp_avg_celsius,                      label='Avg temp', color=temp_color)
    axs[2, 2].plot(time, temp_all_celsius[:, 0],                label='1', color='goldenrod', linestyle='--')
    axs[2, 2].plot(time, temp_all_celsius[:, n_space_steps-1],  label=f'{n_space_steps}', color='palegoldenrod', linestyle='--')
    axs[2, 2].hlines(temp_initial_celsius, x_min, x_max, label=f'{temp_initial_celsius}', color='firebrick', linestyle='--')
    axs[2, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 2].set_xlim(x_min, x_max)
    axs[2, 2].set_title('Temperature')
    axs[2, 2].set(xlabel='Time, hours', ylabel='T, deg C')
    axs[2, 2].legend()

    # plt.savefig(f'data_and_plots/wa {relative_humidity_gas_inlet}, T {temp_initial_celsius}, am {amorphous_material_initial}, time {hours}, n {n_space_steps}.pdf')     # 1% am
    plt.savefig(f'data_and_plots/wa {relative_humidity_gas_inlet}, T {temp_initial_celsius}, am {amorphous_material_initial}, time {hours}, n {n_space_steps} 2.pdf')     # 1% am
    plt.show()

################################ PARAMETERS & INITIAL CONDITIONS #######################################################

# Computing times
hours = 2.1
# hours = 10
seconds = hours * 60 * 60
resolution = 2000
# resolution = 4000
time_step = seconds/resolution                          # each step in time is this many seconds, s
discrete_time = np.linspace(0, seconds, resolution)
while_time = 0

# n_space_steps = 6
n_space_steps = 4
step_length = total_height_powder/n_space_steps
step_length_air = height_air_TAM_cyl / n_space_steps

m_tot_conc_vector           = np.zeros(n_space_steps) + m_particle_tot_conc_initial + m_void_conc_initial
initial_amorphicity_vector  = np.zeros(n_space_steps) + amorphous_material_initial
energy_initial_vector       = np.zeros(n_space_steps) + enthalpy_powder_initial
water_activity_initial_vector = np.zeros(n_space_steps) + water_activity_void_initial
m_sur_initial_vector        = np.zeros(n_space_steps) + m_void_initial
temp_initial_vector         = np.zeros(n_space_steps) + temp_initial_celsius

#TODO: change initials
initials = np.concatenate([m_tot_conc_vector, initial_amorphicity_vector, energy_initial_vector,
                           water_activity_initial_vector, m_sur_initial_vector, temp_initial_vector])    #, temp_sur_initial_vector

water_activity_just_used = np.zeros(n_space_steps) + water_activity_void_initial
tabs = 50
def print_info():
    print('\n')
    print('############################################# CONDITIONS #############################################')
    print('Height of powder bed:'.ljust(tabs), f'{total_height_powder:.3f} m; {total_height_powder*1000:.2f} mm')
    print('Porosity of powder:'.ljust(tabs), f'{porosity_powder:.2f}')
    print('Water activity surroundings:'.ljust(tabs), f'{water_activity_sur:.2f}')
    print('Initial moisture gas as fraction:'.ljust(tabs), f'{m_void_initial:.5f}')
    print('Initial moisture gas as conc:'.ljust(tabs), f'{m_void_conc_initial:.5f} kg/m3')
    print('Saturated moisture gas as fraction:'.ljust(tabs), f'{m_gas_sur:.5f}')
    print('Saturated moisture gas as conc:'.ljust(tabs), f'{m_gas_conc_sat:.5f} kg/m3 \n')

    print('Initial moisture cryst powder as fraction:'.ljust(tabs), f'{m_particle_cryst_initial:.5f}')
    print('Initial moisture cryst powder as conc:'.ljust(tabs), f'{m_particle_cryst_conc_initial:.5f} kg/m3')

    print('Initial moisture am powder as fraction:'.ljust(tabs), f'{m_particle_am_initial:.5f}')
    print('Initial moisture am powder as conc:'.ljust(tabs), f'{m_particle_am_conc_initial:.5f} kg/m3')

    print('Initial moisture total powder as fraction:'.ljust(tabs), f'{m_particle_tot_initial:.5f}')
    print('Initial moisture total powder as conc:'.ljust(tabs), f'{m_particle_tot_conc_initial:.5f} kg/m3\n')

    print('Saturated moisture am powder as fraction:'.ljust(tabs), f'{m_particle_am_sat:.5f}')
    print('Saturated moisture cryst powder as fraction:'.ljust(tabs), f'{m_particle_cryst_sat:.5f} \n')

    print('Total powder water content at beginning:'.ljust(tabs), f'{m_particle_tot_conc_initial:.5f} kg/m3')
    print('Total powder water content at saturation:'.ljust(tabs), f'{m_particle_tot_conc_sat:.5f} kg/m3\n')
    print('Total water content at saturation:'.ljust(tabs), f'{system_conc_sat:.5f}\n')
    print('Max time:'.ljust(tabs), hours, 'hours, = ', seconds, 'seconds.')
    print('Time step:'.ljust(tabs), f'{time_step:5f} s')
    print('')
print_info()

print('############################################# START COMPUTATION #############################################')
counter = 0
t_prev = 0
computation_start_time = time.time()
computed_system = odeint(compute_system, initials, discrete_time, args=(n_space_steps, step_length, True), mxstep=15000)
computation_tot_time = time.time() - computation_start_time
while_time_fraction = while_time/computation_tot_time

np.save(f'data_and_plots/wa {relative_humidity_gas_inlet}, T {temp_initial_celsius}, am {amorphous_material_initial}, '
        f'time {hours}, n {n_space_steps}', computed_system)

print('############################################# COMPUTATION COMPLETE #############################################')
total_water_all         = computed_system[:, 0:n_space_steps]
amorphicity_all         = computed_system[:, n_space_steps:n_space_steps*2]
energy_all              = computed_system[:, n_space_steps*2:n_space_steps*3]
water_activity_all      = computed_system[:, n_space_steps*3:n_space_steps*4]
temp_all_celsius        = computed_system[:, n_space_steps*5:n_space_steps*6]
temp_all_kelvin         = temp_all_celsius + kelvin

m_void_all              = np.zeros([resolution, n_space_steps])
m_void_all[0, :]        = m_void_initial

# temp_all_celsius        = np.zeros([resolution+1, n_space_steps])
# temp_all_celsius[0, :]  = temp_initial_celsius
# temp_all_kelvin         = np.zeros([resolution+1, n_space_steps])
# temp_all_kelvin[0, :]   = temp_initial

extra_water_weight_all  = np.zeros([resolution, n_space_steps])
m_particle_cryst_all    = np.zeros([resolution, n_space_steps])
m_particle_cryst_all[0, :] = m_particle_cryst_initial

m_particle_am_all       = np.zeros([resolution, n_space_steps])
m_particle_am_all[0, :] = m_particle_am_initial

m_powder_total_all      = np.zeros([resolution, n_space_steps])
m_powder_total_all[0, :]= m_particle_tot_initial

p_water_all      = np.zeros([resolution, n_space_steps])

density_gas_all         = np.zeros([resolution, n_space_steps]) + gas_density_initial

m_powder_conc_all       = np.zeros([resolution, n_space_steps])
water_W_all             = np.zeros([resolution, n_space_steps])
excess_water_curr             = np.zeros([resolution, n_space_steps])

# for t in range(1, 1+resolution):
for t in range(resolution):
    # m_void_all[t, :] = compute_H_from_water_activity_temp(water_activity_all[t - 1, :], temp_all_celsius[t - 1, :])
    m_void_all[t, :] = compute_H_from_water_activity_temp(water_activity_all[t, :])

    m_particle_cryst_all[t, :]  = compute_GAB_equilibrium_moisture_cryst(water_activity_all[t, :])
    m_particle_am_all[t, :]     = compute_GAB_equilibrium_moisture_am(water_activity_all[t, :], verbose=False)
    m_powder_total_all[t, :]    = m_particle_am_all[t, :] * amorphicity_all[t, :] + m_particle_cryst_all[t, :] * (1 - amorphicity_all[t, :])
    m_powder_conc_all[t, :]     = m_powder_total_all[t, :] * density_particle * (1 - porosity_powder)
    water_W_all[t, :]           = porosity_powder * m_void_all[t, :] * density_gas_all[t, :] + (1 - porosity_powder) * m_powder_total_all[t, :] * density_particle
    excess_water_curr[t, :]     = total_water_all[t, :] - water_W_all[t, :]

    # temp_all_kelvin[t, :], temp_all_celsius[t, :] = compute_temp_from_energy(m_void_all[t-1, :], m_powder_total_all[t-1, :], energy_all[t-1, :], excess_water_curr[t-1, :])
    # density_gas_all[t-1, :]           = compute_density_air(temp_all_celsius[t-1, :], m_void_all[t-1, :])
    density_gas_all[t, :]           = compute_density_air(m_void_all[t, :])

# m_void_all = m_void_all[1:, :]
# temp_all_kelvin, temp_all_celsius = temp_all_kelvin[1:, :], temp_all_celsius[1:, :]

m_void_conc_all         = m_void_all * porosity_powder * density_gas_all

tg_all                  = compute_glass_temp_mix(1 - m_particle_am_all, glass_temp_lactose, glass_temp_water_1)
t_tg_diff_all           = temp_all_kelvin - tg_all

hours_print = int(computation_tot_time/3600)
minutes_print = int( (computation_tot_time - hours_print*3600) / 60)
seconds_print = computation_tot_time - (minutes_print * 60 + hours_print * 3600)

print('Result:')
print('Total time:'.ljust(tabs), hours_print, 'h', minutes_print, 'min and', int(seconds_print), 's')
print('Time fraction used in while loop:'.ljust(tabs), f'{while_time_fraction:.2f}')

############################ AVERAGES OVER WHOLE CYLINDER ##############################################################
m_void_conc_avg = 0
m_void_avg      = 0
amorphicity_avg = 0
total_water_avg = 0
energy_avg      = 0
density_gas_avg = 0
temp_avg_kelvin = 0
water_activity_avg = 0
for n in range(n_space_steps):
    m_void_conc_avg += m_void_conc_all[:, n]
    m_void_avg      += m_void_all[:, n]
    amorphicity_avg += amorphicity_all[:, n]
    total_water_avg += total_water_all[:, n]
    energy_avg      += energy_all[:, n]
    density_gas_avg += density_gas_all[:, n]
    temp_avg_kelvin += temp_all_kelvin[:, n]
    water_activity_avg += water_activity_all[:, n]

m_void_conc_avg         /= n_space_steps
m_void_avg              /= n_space_steps
amorphicity_avg         /= n_space_steps
total_water_avg         /= n_space_steps
energy_avg              /= n_space_steps        # J/m3
density_gas_avg         /= n_space_steps
temp_avg_kelvin         /= n_space_steps
temp_avg_celsius        = temp_avg_kelvin - kelvin
water_activity_avg      /= n_space_steps

m_particle_cryst_avg    = compute_GAB_equilibrium_moisture_cryst(water_activity_avg)
m_particle_am_avg       = compute_GAB_equilibrium_moisture_am(water_activity_avg)
m_particle_tot_avg      = m_particle_am_avg * amorphicity_avg + m_particle_cryst_avg * (1-amorphicity_avg)
m_powder_avg            = m_particle_am_avg * amorphicity_avg + m_particle_cryst_avg * (1 - amorphicity_avg)
m_powder_conc_avg       = m_powder_avg * density_particle * (1 - porosity_powder)

tg_avg                  = compute_glass_temp_mix(1 - m_particle_am_avg, glass_temp_lactose, glass_temp_water_1)
t_tg_diff_avg           = temp_avg_kelvin - tg_avg

m_powder_diffs_avg      = m_powder_avg[1:] - m_powder_avg[:-1]                   # fraction, kg water/kg dry
m_powder_diffs_avg      = np.insert(m_powder_diffs_avg, 0, m_powder_diffs_avg[0], axis=0)

am_am_diffs = amorphicity_avg[1:] - amorphicity_avg[:-1]
am_am_diffs = np.insert(am_am_diffs, 0, 0, axis=0)
am_cryst_diffs = - am_am_diffs

heat_of_sorption = compute_heat_of_sorption(water_activity_avg, temp_avg_kelvin)   # J/kg
heat_of_sorption_vector = m_powder_diffs_avg * (heat_of_sorption) / time_step               # J/(kg s)
heat_of_sorption_vector /= 1000                                                         # J/(kg s) to J/(g s)

density_bed = density_powder + density_water * m_particle_tot_avg
heat_of_sorption_vector = energy_avg/(density_bed * 1000 * time_step)          # J/m3 to J/kg to J/g
heat_of_sorption_vector -= heat_of_sorption_vector[0]

heat_of_sorption_vector = (energy_avg[1:] - energy_avg[:-1])/(density_powder * 1000 * time_step)
heat_of_sorption_vector = np.insert(heat_of_sorption_vector, 0, 0)

heat_of_cryst_vector = am_cryst_diffs * (heat_of_crystallization/1000) / time_step      # J/kg to J/(g s)

total_energy_vector = heat_of_cryst_vector + heat_of_sorption_vector

########################################## PLOT ########################################################################
# amorphicity_avg = normalize_data(amorphicity_avg, False)
# amorphicity_all = normalize_data(amorphicity_all, False)
plot_tam()

