from initial_conditions import *
import time


######################################## FUNCTIONS #####################################################################
def compute_system(initials, t, n_space_steps, step_length, verbose=False):
    global counter
    global while_time

    global temp_celsius_prev
    global temp_kelvin_prev
    global m_void_prev
    global m_particle_tot_prev
    global water_activity_prev
    global water_activity_vector

    divider = 15000

    total_water     = initials[0:n_space_steps]
    amount_am       = initials[n_space_steps:n_space_steps * 2]
    total_energy    = initials[n_space_steps * 2:n_space_steps * 3]
    water_activity_vector = initials[n_space_steps * 3:n_space_steps * 4]

    amount_am = np.where(amount_am > amorphous_material_initial, amorphous_material_initial, amount_am)

    counter += 1
    ########################################### MOISTURE ###############################################################
    while_time_start = time.time()

    ##### ITERATION #####
    aw_sur_vector = np.zeros(n_space_steps) + water_activity_prev
    cons = [{"type": "ineq", 'fun': lambda x: x}, {"type": "ineq", 'fun': lambda x: 1 - x}]
    optimized_aw = scipy.optimize.minimize(compute_H_and_M_iteratively, aw_sur_vector,
                                           args=(total_water, temp_celsius_prev, amount_am),
                                           constraints=cons, method='SLSQP', options={'ftol': 1e-8})

    water_activity          = optimized_aw.x
    water_activity_change   = water_activity - water_activity_vector
    m_void = compute_H_from_aw_temp(water_activity, temp_celsius_prev)

    while_time += time.time() - while_time_start
    ### END ITERATION ###

    #### TABLES #####
    # m_void, extra_water_weight = compute_H_and_M_gas_tables(amount_am, total_water)
    # water_activity = compute_water_activity_from_m_void(m_void, p_saturated)

    # # water_activity, extra_water_weight = compute_H_and_M_gas_tables_new(amount_am, total_water, temp_celsius_prev)
    # # m_void = compute_H_from_aw_temp(water_activity, temp_celsius_prev)
    #
    # # if verbose and counter % divider == 0 and np.any(extra_water_weight > 0):
    # if np.any(extra_water_weight > 0):
    #     np.set_printoptions()
    #     # print('\nTotal water:'.ljust(tabs), total_water)
    #     print('Extra weight:'.ljust(tabs), extra_water_weight)
    #     print('Amount am:'.ljust(tabs), amount_am)
    #     print('m_void:'.ljust(tabs), m_void)
    #     print('Temp:'.ljust(tabs), temp_kelvin)
    #
    #     # print(abab)
    #
    # while_time_this = time.time() - while_time_start
    # while_time += while_time_this
    #
    # # if counter % divider2 == 0:
    # #     np.set_printoptions(formatter={'float': '{: .1e}'.format})
    # #     print('Table method')
    # #     print('m_void:'.ljust(tabs), m_void)
    # #     print('water act:'.ljust(tabs), water_activity)

    ## END TABLES ###

    # Update system for m_void

    m_particle_cryst    = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_particle_am       = compute_GAB_equilibrium_moisture_am(water_activity)
    m_particle_tot      = m_particle_cryst * (1-amount_am) + m_particle_am * amount_am

    # excess_water = np.where(water_activity >= 1, 1, 0)
    # total_water_new = m_void * porosity_powder * density_gas + m_particle_tot * (1 - porosity_powder) * density_particle
    ########################################### ENERGY #################################################################
    temp_kelvin, temp_celsius = compute_temp_from_energy(m_void_prev, m_particle_tot_prev, total_energy, verbose=False)
    # temp_kelvin, temp_celsius = compute_temp_from_energy(m_void, m_particle_tot, total_energy, verbose=False)

    enthalpy_vapor = heat_capacity_vapor * temp_celsius + heat_of_evaporation_water

    ################################## RESULTING CHANGES IN Q AND W ####################################################
    gas_density_curr = compute_density_air(temp_celsius, m_void)

    laplacian_conc = compute_laplacian(m_void * gas_density_curr, m_gas_sur * gas_density_initial, step_length)
    m_diffusion = laplacian_conc * diffusivity_eff                          # kg/m5 * m2/s = kg/(m3 s)

    # Energy, Higher temp = lower density -> surroundings will have higher density
    # Energy, Higher temp = higher enthalpy vapor -> surroundings will have lower enthalpy
    energy_diff_laplacian = compute_laplacian(enthalpy_vapor * m_void_prev * gas_density_curr, enthalpy_vapor_initial * m_gas_sur * gas_density_sur, step_length, False)
    change_energy_diffusion = diffusivity_eff * energy_diff_laplacian

    temp_laplacian = compute_laplacian(temp_celsius, temp_initial_celsius, step_length)
    conductivity_tot = porosity_powder * conductivity_gas + (1-porosity_powder) * conductivity_particle
    change_energy_conduction = conductivity_tot * temp_laplacian

    ###################################### AMORPHICITY #################################################################
    glass_temp = compute_glass_temp_mix(1 - m_particle_am, glass_temp_lactose, glass_temp_water_1)
    temp_diff = temp_kelvin - glass_temp
    change_amorph = compute_kinetics_avrami(amount_am, temp_diff)
    # change_amorph = np.zeros(n_space_steps)

    m_freed = (m_particle_am - m_particle_cryst) * change_amorph                    # kg/kg solid lactose
    m_freed *= (1-porosity_powder) * density_particle                               # kg/m3 water
    # m_crystallization = 0
    change_energy_evaporation       = m_freed * heat_of_evaporation_water              # Negative, kg/(m3 s) * J/kg
    change_energy_crystallization   = change_amorph * (1-porosity_powder) * density_particle * heat_of_crystallization         # Negative

    # change_energy = change_energy_diffusion + change_energy_conduction + change_energy_evaporation
    change_energy = change_energy_diffusion + change_energy_conduction + change_energy_evaporation + change_energy_crystallization

    if verbose and counter % divider == 0:
        hours_print = int(t / 3600)
        minutes_print = int((t - hours_print * 3600) / 60)
        seconds_print = int(t - (minutes_print * 60 + hours_print * 3600))

        print('\n')
        print('         TIME', hours_print, 'h;', minutes_print, 'min;', seconds_print, 's')
        np.set_printoptions(formatter={'float': '{: .2f}'.format})
        print('Temp C:'.ljust(tabs), temp_celsius)
        print('T-Tg:'.ljust(tabs), temp_diff)
        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        print('Total energy:'.ljust(tabs), total_energy)
        print('')

        print('m_void:'.ljust(tabs), m_void)
        print('m_gas_sur:'.ljust(tabs), f'{m_gas_sur:.1e}')
        print('M am:'.ljust(tabs), m_particle_am)

        np.set_printoptions(formatter={'float': '{: .3f}'.format})
        print('Water act:'.ljust(tabs), water_activity)
        np.set_printoptions(formatter={'float': '{: .5f}'.format})
        print('Am am:'.ljust(tabs), amount_am)

        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        print('Change am:'.ljust(tabs), change_amorph)
        print('')

        print('Change energy diff:'.ljust(tabs), change_energy_diffusion)
        print('Change energy cond:'.ljust(tabs), change_energy_conduction)
        print('Change energy cryst:'.ljust(tabs), change_energy_crystallization)
        print('Change energy evap:'.ljust(tabs), change_energy_evaporation)
        print('Change energy tot:'.ljust(tabs), change_energy)
        print('\n')

    change_moisture = m_diffusion
    system_change = np.concatenate([change_moisture, change_amorph, change_energy, water_activity_change]).astype('float64')
    m_void_prev = m_void
    m_particle_tot_prev = m_particle_tot
    temp_celsius_prev = temp_celsius
    temp_kelvin_prev = temp_kelvin
    water_activity_prev = water_activity
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
    x_min = -0.1
    x_min = 0
    x_max = hours + 0.1
    # x_max = 0.02
    fig.suptitle(f'Moisture sorption TAM at RH {relative_humidity_gas_inlet} and T {temp_initial-kelvin} C')
    # axs[0, 0].plot(time, m_void_conc_cryst, label='Crystalline', color='lightsteelblue', linestyle='--')
    # axs[0, 0].plot(time, m_void_conc_am, label='Amorhpous', color=m_color, linestyle='--')

    ########################################### MOISTURE ###############################################################
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

    ########################################### MISC. ##################################################################
    axs[0, 1].hlines(0, -1, hours+1, color=sat_color, linestyles='--')
    axs[0, 1].plot(time, t_tg_diff_avg,         label='Avg', color=temp_color)
    axs[0, 1].plot(time, t_tg_diff_all[:, 0],   label='1', color='rosybrown', linestyle='--')
    axs[0, 1].plot(time, t_tg_diff_all[:, 3],   label='4', color='darkgoldenrod', linestyle='--')
    axs[0, 1].legend()
    axs[0, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 1].set_xlim(x_min, x_max)
    axs[0, 1].set_title('Glass transition temp')
    axs[0, 1].set(xlabel='Time, hours', ylabel='Opeating temp - glass temp')

    axs[1, 1].plot(time, amorphicity_avg, label='Avg', color=am_color)
    axs[1, 1].plot(time, amorphicity_all[:, 0],                         label='1', color='green', linestyle='--')
    axs[1, 1].plot(time, amorphicity_all[:, int(n_space_steps/3)],      label=f'{1+int(n_space_steps/3)}', color='lightgreen', linestyle='--')
    axs[1, 1].plot(time, amorphicity_all[:, 2 * int(n_space_steps/3)],  label=f'{1+2*int(n_space_steps/3)}', color='limegreen', linestyle='--')
    axs[1, 1].plot(time, amorphicity_all[:, n_space_steps-1],           label=f'{n_space_steps}', color='darkgreen', linestyle='--')
    axs[1, 1].legend()
    axs[1, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 1].set_xlim(x_min, x_max)
    axs[1, 1].set_ylim(-0.1, 0.2)
    axs[1, 1].set_title('Amount amorhpous material')
    axs[1, 1].set(xlabel='Time, hours', ylabel='Am material, kg/kg lactose')

    # axs[2, 1].hlines(m_particle_cryst_sat, -1, hours+1, label='Saturated cryst', color=sat_color, linestyles='--')
    # axs[2, 1].plot(time, m_particle_cryst_avg, label='Crystalline', color=m_color)
    axs[2, 1].plot(time, energy_avg, label='Avg energy', color='darkorange')
    axs[2, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 1].set_xlim(x_min, x_max)
    axs[2, 1].set_title('Average energy')
    axs[2, 1].legend()
    axs[2, 1].set(xlabel='Time, hours', ylabel='Energy')

    ########################################### HEAT FLOW ##############################################################
    # axs[0, 2].plot(time, heat_of_sorption_vector, label='Heat of sorption', color='rosybrown', linestyle='--')
    axs[0, 2].plot(time, heat_of_cryst_vector, label='Crystallization', color='darkgoldenrod', linestyle='--')
    # axs[0, 2].plot(time, total_energy_vector, label='Total', color=temp_color)
    axs[0, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    # axs[1, 1].set_ylim(0, 0.9)
    axs[0, 2].set_xlim(x_min, x_max)
    axs[0, 2].set_title('Sorption')
    axs[0, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')
    axs[0, 2].legend()

    axs[1, 2].plot(time, water_activity_avg, label='Avg w_a', color=am_color)
    axs[1, 2].plot(time, water_activity_all[:, 0], label='1', color='green', linestyle='--')
    # axs[1, 2].plot(time, water_activity_all[:, int(n_space_steps / 3)], label=f'{1+int(n_space_steps/3)}',
    #                color='lightgreen', linestyle='--')
    # axs[1, 2].plot(time, water_activity_all[:, 2 * int(n_space_steps / 3)], label=f'{1+2*int(n_space_steps/3)}',
    #                color='limegreen', linestyle='--')
    axs[1, 2].plot(time, water_activity_all[:, n_space_steps - 1], label=f'{n_space_steps}', color='darkgreen',
                   linestyle='--')
    axs[1, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    # axs[1, 1].set_ylim(0, 0.9)
    axs[1, 2].set_xlim(x_min, x_max)
    axs[1, 2].set_title('Water activity over time')
    axs[1, 2].set(xlabel='Time, hours', ylabel='wa')
    axs[1, 2].legend()

    # max_cryst = np.max(heat_of_cryst_vector)
    # axs[2, 2].set_ylim(0, max_cryst)
    axs[2, 2].plot(time, temp_avg_celsius,                      label='Avg temp', color=temp_color)
    axs[2, 2].plot(time, temp_all_celsius[:, 0],                label='1', color='goldenrod', linestyle='--')
    axs[2, 2].plot(time, temp_all_celsius[:, n_space_steps-1],  label=f'{n_space_steps}', color='palegoldenrod', linestyle='--')
    axs[2, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 2].set_xlim(x_min, x_max)
    axs[2, 2].set_title('Temperature')
    axs[2, 2].set(xlabel='Time, hours', ylabel='T, deg C')
    axs[2, 2].legend()

    # plt.savefig('system_w_derivative4.pdf')
    plt.savefig('system_w_derivative5.pdf')
    plt.show()

################################ PARAMETERS & INITIAL CONDITIONS #######################################################

# Computing times
hours = 1.5
# hours = 0.02
seconds = hours * 60 * 60
resolution = 1000
# resolution = 15
time_step = seconds/resolution                      # each step in time is this many seconds, s
discrete_time = np.linspace(0, seconds, resolution)
while_time = 0

n_space_steps = 8
step_length = total_height_powder/n_space_steps

# m_gas_conc_initial_vector = np.zeros(n_space_steps) + m_void_conc_initial
# m_particle_tot_conc_initial_vector = np.zeros(n_space_steps) + m_particle_tot_conc_initial
m_tot_conc_vector = np.zeros(n_space_steps) + m_particle_tot_conc_initial + m_void_conc_initial
# m_tot_conc_vector = np.zeros(n_space_steps) + m_particle_cryst_conc_initial + m_void_conc_initial
initial_amorphicity_vector = np.zeros(n_space_steps) + amorphous_material_initial
energy_initial_vector = np.zeros(n_space_steps) + enthalpy_powder_initial
water_activity_initial_vector = np.zeros(n_space_steps) + water_activity_void_initial

initials = np.concatenate([m_tot_conc_vector, initial_amorphicity_vector, energy_initial_vector, water_activity_initial_vector])      # kg/m3, kg/kg

temp_kelvin_prev = np.zeros(n_space_steps) + temp_initial
temp_celsius_prev = temp_kelvin_prev - kelvin
m_void_prev = m_void_initial
m_particle_tot_prev = m_particle_tot_initial
water_activity_prev = water_activity_void_initial
tabs = 50
def print_info():
    print('\n')
    print('############################################# CONDITIONS #############################################')
    print('Total height of bed:', total_height_powder)

    print('Height of powder bed:'.ljust(tabs), f'{total_height_powder:.3f} m; {total_height_powder*1000:.2f} mm')
    print('Water activity surroundings:'.ljust(tabs), f'{water_activity_sur:.2f}')
    print('Moisture in the surrounding:'.ljust(tabs), f' {m_gas_sur:.5f}')
    print('Initial moisture gas as fraction:'.ljust(tabs), f'{m_void_initial:.5f}')
    print('Initial moisture gas as conc:'.ljust(tabs), f'{m_void_conc_initial:.5f} kg/m3')
    print('Saturated moisture gas as conc:'.ljust(tabs), f'{m_gas_conc_sat:.5f} kg/m3')
    print('Saturated moisture gas:'.ljust(tabs), f'{m_gas_sur:.5f} kg/m3\n')

    print('Initial moisture cryst powder as fraction:'.ljust(tabs), f'{m_particle_cryst_initial:.5f}')
    print('Initial moisture cryst powder as conc:'.ljust(tabs), f'{m_particle_cryst_conc_initial:.5f} kg/m3')

    print('Initial moisture am powder as fraction:'.ljust(tabs), f'{m_particle_am_initial:.5f}')
    print('Initial moisture am powder as conc:'.ljust(tabs), f'{m_particle_am_conc_initial:.5f} kg/m3\n')

    print('Total powder water content at beginning:'.ljust(tabs), f'{m_particle_tot_conc_initial:.5f}')
    print('Total powder water content at saturation:'.ljust(tabs), f'{m_particle_tot_conc_sat:.5f}\n')
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

print('############################################# COMPUTATION COMPLETE #############################################')
total_water_all         = computed_system[:, 0:n_space_steps]
# total_water_all         = np.insert(total_water_all, 0, 0)
amorphicity_all         = computed_system[:, n_space_steps:n_space_steps*2]
energy_all              = computed_system[:, n_space_steps*2:n_space_steps*3]
water_activity_all      = computed_system[:, n_space_steps*3:n_space_steps*4]

m_void_all              = np.zeros([resolution + 1, n_space_steps])
m_void_all[0, :]        = m_void_initial

temp_all_celsius        = np.zeros([resolution+1, n_space_steps])
temp_all_celsius[0, :]  = temp_initial_celsius
temp_all_kelvin         = np.zeros([resolution+1, n_space_steps])
temp_all_kelvin[0, :]   = temp_initial

# water_activity_all      = np.zeros([resolution+1, n_space_steps])
# water_activity_all[0, :]= water_activity_void_initial

extra_water_weight_all  = np.zeros([resolution+1, n_space_steps])
m_particle_cryst_all    = np.zeros([resolution, n_space_steps])
m_particle_cryst_all[0, :] = m_particle_cryst_initial

m_particle_am_all       = np.zeros([resolution, n_space_steps])
m_particle_am_all[0, :] = m_particle_am_initial

m_powder_total_all      = np.zeros([resolution, n_space_steps])
m_powder_total_all[0, :]= m_particle_tot_initial

p_water_all      = np.zeros([resolution, n_space_steps])

density_gas_all         = np.zeros([resolution, n_space_steps])

m_powder_conc_all       = np.zeros([resolution, n_space_steps])

cons = [{"type": "ineq", 'fun': lambda x: x}, {"type": "ineq", 'fun': lambda x: 1 - x}]
for t in range(1, 1+resolution):
    # water_activity_all[t, :], extra_water_weight_all[t, :] = compute_H_and_M_gas_tables_new(amorphicity_all[t-1, :], total_water_all[t-1, :], temp_all_celsius[t-1, :])

    # optimized_aw = scipy.optimize.minimize(compute_H_and_M_iteratively, water_activity_all[t-1, :],
    #                                        args=(total_water_all[t-1, :], temp_all_celsius[t-1, :], amorphicity_all[t-1, :]),
    #                                        constraints=cons, method='SLSQP', options={'ftol': 1e-8})    #TODO: changed from total_water_all[t-1, :]
    #
    # water_activity_all[t, :] = optimized_aw.x
    m_void_all[t, :] = compute_H_from_aw_temp(water_activity_all[t-1, :], temp_all_celsius[t-1, :])

    m_particle_cryst_all[t-1, :]  = compute_GAB_equilibrium_moisture_cryst(water_activity_all[t-1, :])
    m_particle_am_all[t-1, :]     = compute_GAB_equilibrium_moisture_am(water_activity_all[t-1, :], verbose=False)
    m_powder_total_all[t-1, :]    = m_particle_am_all[t-1, :] * amorphicity_all[t-1, :] + m_particle_cryst_all[t-1, :] * (1 - amorphicity_all[t-1, :])
    m_powder_conc_all[t-1, :]     = m_powder_total_all[t-1, :] * density_particle * (1 - porosity_powder)

    temp_all_kelvin[t, :], temp_all_celsius[t, :] = compute_temp_from_energy(m_void_all[t-1, :], m_powder_total_all[t-1, :], energy_all[t-1, :])
    density_gas_all[t-1, :]           = compute_density_air(temp_all_celsius[t-1, :], m_void_all[t-1, :])

m_void_all = m_void_all[1:, :]
# water_activity_all = water_activity_all[1:, :]
temp_all_kelvin, temp_all_celsius = temp_all_kelvin[1:, :], temp_all_celsius[1:, :]


m_void_conc_all         = m_void_all * porosity_powder * density_gas_all
# amorphicity_all         = normalize_data(amorphicity_all, zero_one=False)

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
m_void_avg = 0
amorphicity_avg = 0
total_water_avg = 0
energy_avg = 0
density_gas_avg = 0
temp_avg_kelvin = 0
# m_particle_tot_conc_total = 0
for n in range(n_space_steps):
    m_void_conc_avg += m_void_conc_all[:, n]
    m_void_avg      += m_void_all[:, n]
    amorphicity_avg += amorphicity_all[:, n]
    total_water_avg += total_water_all[:, n]
    energy_avg      += energy_all[:, n]
    density_gas_avg += density_gas_all[:, n]
    temp_avg_kelvin += temp_all_kelvin[:, n]

m_void_conc_avg         /= n_space_steps
m_void_avg              /= n_space_steps
amorphicity_avg         /= n_space_steps
total_water_avg         /= n_space_steps
energy_avg              /= n_space_steps
density_gas_avg         /= n_space_steps
temp_avg_kelvin         /= n_space_steps
temp_avg_celsius        = temp_avg_kelvin - kelvin

water_activity_avg      = compute_water_activity_from_m_void(m_void_avg, temp_avg_celsius)

m_particle_cryst_avg    = compute_GAB_equilibrium_moisture_cryst(water_activity_avg)
m_particle_am_avg       = compute_GAB_equilibrium_moisture_am(water_activity_avg)
m_powder_avg            = m_particle_am_avg * amorphicity_avg + m_particle_cryst_avg * (1 - amorphicity_avg)
m_powder_conc_avg       = m_powder_avg * density_particle * (1 - porosity_powder)

# temp_avg_kelvin, temp_avg_celsius = compute_temp_from_energy(m_void_avg, m_powder_avg, energy_avg)

tg_avg                  = compute_glass_temp_mix(1 - m_particle_am_avg, glass_temp_lactose, glass_temp_water_1)
t_tg_diff_avg           = temp_avg_kelvin - tg_avg

m_powder_diffs_avg      = m_powder_avg[1:] - m_powder_avg[:-1]                   # fraction, kg water/kg dry
m_powder_diffs_avg      = np.insert(m_powder_diffs_avg, 0, m_powder_diffs_avg[0], axis=0)

am_am_diffs = amorphicity_avg[1:] - amorphicity_avg[:-1]
am_am_diffs = np.insert(am_am_diffs, 0, 0, axis=0)
am_cryst_diffs = - am_am_diffs

# heat_of_sorption = compute_heat_of_sorption(water_activity_avg, temp_initial)   # J/kg
heat_of_sorption = compute_heat_of_sorption(water_activity_avg, temp_avg_kelvin)   # J/kg   TODO: added
# heat_of_sorption = 2.5
heat_of_sorption_vector = m_powder_diffs_avg * (heat_of_sorption) / time_step               # J/(kg s)
heat_of_sorption_vector /= 1000                                                         # J/(kg s) to J/(g s)

heat_of_cryst_vector = am_cryst_diffs * (heat_of_crystallization/1000) / time_step      # J/kg to J/(g s)

total_energy_vector = heat_of_cryst_vector + heat_of_sorption_vector

########################################## PLOT ########################################################################
plot_tam()

