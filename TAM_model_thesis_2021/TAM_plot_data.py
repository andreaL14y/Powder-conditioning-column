################################### FILE TO PLOT SIMULATION RESULTS ####################################################
from TAM_functions import*

# Choose what to plot, one of the following must be true
humidity_c    = False
humidity_c = True

weight_c      = False
# weight_c      = True

amorphicity_c = False
# amorphicity_c = True

low_hum_c     = False
# low_hum_c     = True

unhindered = False
# unhindered = True

if humidity_c:
    exp_list = np.array([16, 17, 18])
    exp_list = np.array([13, 14, 15])
    c_ind = 4, 6, -1
    labels = np.array(['RH 53', 'RH 57', 'RH 75'])
elif weight_c:
    exp_list = np.array([13, 20])
    c_ind = 2, 7
    labels = np.array(['150 mg', '300 mg'])
elif amorphicity_c:
    exp_list = np.array([14, 17])
    c_ind = 2, 6
    labels = np.array(['21%', '14%'])

elif low_hum_c:
    exp_list = np.array([24])
    exp_list = np.array([19])
    exp_list = np.array([22])
    c_ind = 0, 1
    if exp_list[0] == 24:
        labels = np.array(['RH 33%'])
    else:
        labels = np.array(['RH 38%'])
        c_ind = 2, 3

c_chosen = c_palette[c_ind, :]
grid_color = c_palette[1]

def plot_tam():
    fig, axs = plt.subplots(3, 3, figsize=(25, 15))
    time = discrete_time / 60

    grid_color = c_palette[1]
    x_min = 0
    x_max = hours * 60 + 5
    fig.suptitle(f'Exp: {exp}. RH {relative_humidity_gas_inlet}, T {temp_initial-kelvin} C and {amorphous_material_initial * 100} % am initially')

    ############################################# m_void ###############################################################
    # axs[0, 0].hlines(m_gas_sur, -1, hours+1, label='Sat sur', linestyles='--', color=c_palette[5])
    # axs[0, 0].plot(time, m_void_avg[:], label='Avg', color=c_palette[1])
    # axs[0, 0].plot(time, m_void_avg[:], label='Avg', color=c_palette[1])
    axs[0, 0].plot(time, -m_evap, label='Evap', color=c_palette[2])
    axs[0, 0].plot(time, m_powder_diffs_avg, label='Powder diff', color=c_palette[6])
    axs[0, 0].plot(time, 2*m_powder_diffs_avg - m_evap, label='Diff', color=c_palette[3])
    # axs[0, 0].plot(time, m_void_all[:, n_space_steps-1], label=f'{n_space_steps}', color=c_palette[4], linestyle='--')
    axs[0, 0].legend()
    axs[0, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 0].set_xlim(x_min, x_max)
    axs[0, 0].set_title('Void moisture concentration')
    axs[0, 0].set(xlabel='Time, hours', ylabel='Moisture conc void, kg/m3')

    ############################################# m_am #################################################################
    # axs[1, 0].hlines(m_particle_am_sat, -1, hours+1,                    label='Sat am', color=c_palette[5], linestyles='--')
    axs[1, 0].plot(time, m_particle_am_avg,                             label='Am', color=c_palette[1])
    axs[1, 0].plot(time, m_particle_am_all[:, 0],                       label='1', color=c_palette[2], linestyle='--')
    # axs[1, 0].plot(time, m_particle_am_all[:, int(n_space_steps / 3)],  label=f'{int(n_space_steps/3)}', color=c_palette[3], linestyle='--')
    axs[1, 0].plot(time, m_particle_am_all[:, n_space_steps-1],  label=f'{n_space_steps}', color=c_palette[4], linestyle='--')
    # axs[1, 0].plot(time, m_particle_am_all[:, int(2 * n_space_steps/3)], label=f'{int(2 * n_space_steps/3)}', color=c_palette[0], linestyle='--')
    # axs[1, 0].plot(time, m_particle_am_all[:, -1], label=f'{n_space_steps-1}', color=c_palette[1], linestyle='--')
    axs[1, 0].legend()
    axs[1, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 0].set_xlim(x_min, x_max)
    axs[1, 0].set_title('Moisture content amorphous')
    axs[1, 0].set(xlabel='Time, hours', ylabel='Moisture, kg/kg lactose')

    ########################################### Total Water W ##########################################################
    axs[2, 0].plot(time, total_water_avg,                           label='Avg', color=c_palette[1])
    axs[2, 0].plot(time, total_water_all[:, 0],                     label=f'{1}', color=c_palette[2], linestyle='--')
    # axs[2, 0].plot(time, total_water_all[:, int(n_space_steps/3)],  label=f'{int(n_space_steps/3)}', color=c_palette[3], linestyle='--')
    axs[2, 0].plot(time, total_water_all[:, n_space_steps-1],  label=f'{n_space_steps}', color=c_palette[4], linestyle='--')
    # axs[2, 0].plot(time, total_water_all[:, int(2 * n_space_steps/3)], label=f'{int(2 * n_space_steps/3)}', color=c_palette[0], linestyle='--')
    # axs[2, 0].plot(time, total_water_all[:, -1], label=f'{n_space_steps-1}', color=c_palette[0], linestyle='--')
    # axs[2, 0].plot(time, m_particle_cryst_avg, label='Cr', color=c_palette[2])
    axs[2, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 0].set_xlim(x_min, x_max)
    axs[2, 0].set_title('Total moisture system')
    axs[2, 0].legend()
    axs[2, 0].set(xlabel='Time, hours', ylabel='Moisture content, kg/m3')

    ########################################### T-Tg ###################################################################
    axs[0, 1].hlines(0, -1, hours+1, color=c_palette[-1], linestyles='--')
    axs[0, 1].plot(time, t_tg_diff_avg,         label='Avg', color=c_palette[-2])
    axs[0, 1].plot(time, t_tg_diff_all[:, 0],   label='1', color=c_palette[-3], linestyle='--')
    axs[0, 1].plot(time, t_tg_diff_all[:, n_space_steps-1],   label=f'{n_space_steps}', color=c_palette[-4], linestyle='--')
    axs[0, 1].legend()
    axs[0, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 1].set_xlim(x_min, x_max)
    axs[0, 1].set_title('Glass transition temp')
    axs[0, 1].set(xlabel='Time, hours', ylabel='Opeating temp - glass temp')

    ###################################### Amorphicity X ###############################################################
    axs[1, 1].plot(time, amorphicity_avg, label='Avg', color=c_palette[2])
    axs[1, 1].plot(time, amorphicity_all[:, 0],                         label='1', color=c_palette[3], linestyle='--')
    # axs[1, 1].plot(time, amorphicity_all[:, int(n_space_steps/3)],      label=f'{1+int(n_space_steps/3)}', color='green', linestyle='--')
    # axs[1, 1].plot(time, amorphicity_all[:, 2 * int(n_space_steps/3)],  label=f'{1+2*int(n_space_steps/3)}', color='limegreen', linestyle='--')
    axs[1, 1].plot(time, amorphicity_all[:, n_space_steps-1],           label=f'{n_space_steps}', color=c_palette[4], linestyle='--')
    axs[1, 1].legend()
    axs[1, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 1].set_xlim(x_min, x_max)
    axs[1, 1].set_ylim(-0.001, 0.004)
    # axs[1, 1].set_ylim(-0.1, np.max(amorphicity_avg) + 0.1)
    axs[1, 1].set_title('Amount amorhpous material')
    axs[1, 1].set(xlabel='Time, hours', ylabel='Am material, kg/kg lactose')

    ################################### Water Activity RH ##############################################################
    axs[2, 1].plot(time, water_activity_avg, label='Avg w_a', color=c_palette[3])
    axs[2, 1].plot(time, water_activity_all[:, 0], label='1', color=c_palette[4], linestyle='--')
    axs[2, 1].plot(time, water_activity_all[:, 1], label=f'{2}',
                   color=c_palette[5], linestyle='--')
    # axs[2, 1].plot(time, water_activity_all[:, 2 * int(n_space_steps / 3)], label=f'{1+2*int(n_space_steps/3)}',
    #                color='limegreen', linestyle='--')
    axs[2, 1].plot(time, water_activity_all[:, 2], label=f'{3}', color=c_palette[6],
                   linestyle='--')
    axs[2, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 1].set_xlim(0, x_max)
    # axs[2, 1].set_ylim(0.5, 1.1)
    axs[2, 1].set_title('Water activity over time')
    axs[2, 1].set(xlabel='Time, hours', ylabel='wa')
    axs[2, 1].legend()

    ################################### Heat of sorption ###############################################################
    axs[0, 2].plot(time, heat_of_sorption_vector, label='Heat of sorption', color=c_palette[8], linestyle='--')
    axs[0, 2].plot(time, heat_of_cryst_vector, label='Heat of cryst', color=c_palette[6], linestyle='--')
    axs[0, 2].plot(time, total_energy_vector, label='Total', color=c_palette[4])

    integral = simps(total_energy_vector, time * 3600)
    axs[0, 2].text(x_max / 2, 0.01, f'Integral {integral:0.2f}')

    axs[0, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 2].set_xlim(x_min, x_max)
    # axs[0, 2].set_ylim(-0.05, 0.03)
    axs[0, 2].set_title('Heat generated')
    axs[0, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')
    axs[0, 2].legend()

    ################################# Heat of crystallization ##########################################################
    integral = simps(heat_of_cryst_vector, time*3600)
    axs[1, 2].text(x_max/2, 0.01, f'Integral {integral:0.2f}')
    # axs[1, 2].plot(time, exc_water_avg, label='Excess water', color=c_palette[2]) #, linestyle='--')
    # axs[1, 2].plot(time, energy_avg, label='Heat flow', color=c_palette[5]) #, linestyle='--')
    # axs[1, 2].plot(time, heat_flow_vector, label='Heat flow', color=c_palette[5]) #, linestyle='--')
    axs[1, 2].plot(time * 60, heat_of_cryst_vector, label='Heat of cryst', color=c_palette[6], linestyle='--')

    # axs[1, 2].plot(time, total_energy_vector, label='Total', color=c_palette[6])
    axs[1, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 2].set_xlim(x_min, x_max * 60)
    # axs[1, 2].set_ylim(-0.004, 0.03)
    axs[1, 2].set_title('Heat generated')
    axs[1, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')
    axs[1, 2].legend()

    ####################################### Temperature ################################################################
    # max_cryst = np.max(heat_of_cryst_vector)
    # axs[2, 2].set_ylim(0, max_cryst)
    axs[2, 2].plot(time * 60, temp_avg_celsius,                      label='Avg temp', color=c_palette[-1])
    axs[2, 2].plot(time*60, temp_all_celsius[:, 0],                label='1', color=c_palette[4], linestyle='--')
    axs[2, 2].plot(time*60, temp_all_celsius[:, n_space_steps-1],  label=f'{n_space_steps}', color=c_palette[4], linestyle='--')
    axs[2, 2].hlines(temp_initial_celsius, x_min, x_max, label=f'{temp_initial_celsius}', color=c_palette[6], linestyle='--')
    axs[2, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 2].set_xlim(x_min, x_max*60)
    axs[2, 2].set_title('Temperature')
    axs[2, 2].set(xlabel='Time, hours', ylabel='T, deg C')
    axs[2, 2].legend()

    # plt.savefig(f'data_and_plots/exp {exp}, time {hours}, n {n_space_steps}.pdf')     # 1% am
    plt.show()

def plot_tam_4():
    # fig, axs = plt.subplots(2, 2, figsize=(25, 15))
    fig, axs = plt.subplots(2, 2, figsize=fig_size, dpi=100)
    # fig, axs = plt.subplots(2, 2, figsize=(10, 7), dpi=50, facecolor='w', edgecolor='k')
    time = discrete_time / 60
    minutes = hours * 60
    grid_color = c_palette[1]
    x_min = 0
    x_max = minutes + 5
    x_max = 200
    # fig.suptitle(f'Exp: {exp}. RH {relative_humidity_gas_inlet}, T {temp_initial-kelvin} C and {int(amorphous_material_initial * 100)} % am initially')

    ########################################### Total Water W ##########################################################
    axs[0, 0].plot(time, total_water_avg,                           label='Total water, W', color=c_palette[1])
    axs[0, 0].hlines(m_particle_tot_conc_initial, -1, minutes + 10, color=c_palette[3], linestyles='--', label='Initial')
    axs[0, 0].hlines(m_particle_cryst_conc_initial, -1, minutes + 10, color=c_palette[2], linestyles='--', label='Crystalline')
    # axs[2, 0].plot(time, total_water_all[:, 0],                     label=f'{1}', color=c_palette[2], linestyle='--')
    # axs[2, 0].plot(time, total_water_all[:, n_space_steps-1],  label=f'{n_space_steps}', color=c_palette[4], linestyle='--')
    axs[0, 0].grid(alpha = grid_alpha_min, color=grid_color)
    axs[0, 0].set_xlim(x_min, x_max)
    # axs[0, 0].set_title('Total moisture system')
    axs[0, 0].legend()
    axs[0, 0].set(xlabel='Time, min', ylabel='Moisture content, kg/m3')

    ########################################### T-Tg ###################################################################
    axs[0, 1].hlines(0, -1, minutes+10, color=c_palette[-4], linestyles='--')
    t_tg_diff_avg_pl = savgol_filter(t_tg_diff_avg, 11, 1)
    axs[0, 1].plot(time, t_tg_diff_avg_pl, color=c_palette[-2]) # ,         label=r'$T - T_G$'
    # axs[0, 1].plot(time, t_tg_diff_all[:, 0],   label='1', color=c_palette[-3], linestyle='--')
    # axs[0, 1].plot(time, t_tg_diff_all[:, n_space_steps-1],   label=f'{n_space_steps}', color=c_palette[-4], linestyle='--')
    # axs[0, 1].legend()
    axs[0, 1].grid(alpha = grid_alpha_min, color=grid_color)
    axs[0, 1].set_xlim(x_min, x_max)
    # axs[0, 1].set_title('Glass transition temp')
    axs[0, 1].set(xlabel='Time, min', ylabel=r'$T - T_G$, $^\circ$C')

    ###################################### Amorphicity X ###############################################################
    axs[1, 1].hlines(0, -1, minutes+10, color=c_palette[7], linestyles='--', label='0')
    axs[1, 1].hlines(amorphous_material_initial, -1, minutes+10, color=c_palette[4], linestyles='--', label='Initial')
    axs[1, 1].plot(time, amorphicity_avg, color=c_palette[6], label='Amorphicity')   #
    axs[1, 1].legend()
    axs[1, 1].grid(alpha=grid_alpha_min, color=grid_color)
    axs[1, 1].set_xlim(x_min, x_max)
    axs[1, 1].set_ylim(-0.01, np.max(amorphicity_avg) + 0.01)
    # axs[1, 1].set_title('Amorphicity, %')
    axs[1, 1].set(xlabel='Time, min', ylabel='Amorphicity, %')

    ############################################# Water activity #################################################################
    axs[1, 0].plot(time, water_activity_avg, color=c_palette[3])    # , label=r'Water activity $a_w$'
    # axs[1, 0].plot(time, water_activity_all[:, 0], label='1', color=c_palette[4], linestyle='--')
    # axs[1, 0].plot(time, water_activity_all[:, 1], label=f'{2}',
    #                color=c_palette[5], linestyle='--')
    # axs[2, 1].plot(time, water_activity_all[:, 2 * int(n_space_steps / 3)], label=f'{1+2*int(n_space_steps/3)}',
    #                color='limegreen', linestyle='--')
    # axs[1, 0].plot(time, water_activity_all[:, 2], label=f'{3}', color=c_palette[6],
    #                linestyle='--')
    # axs[1, 0].set_title('Water activity over time')
    axs[1, 0].set(xlabel='Time, hours', ylabel=r'Water activity, $a_w$')
    axs[1, 0].grid(alpha=grid_alpha_min, color=grid_color)
    axs[1, 0].set_xlim(x_min, x_max)
    # axs[1, 0].set_title('Water activity')
    axs[1, 0].set(xlabel='Time, min')   # , ylabel=r'Water activity, $a_w$'
    # axs[1, 0].legend(loc='best')

    plt.savefig(f'data_and_plots/211014-unhindered/exp {exp}, time {hours}, n {n_space_steps}.pdf')     # 1% am
    plt.show()


fig=plt.figure(figsize=fig_size, dpi= 100, facecolor='w', edgecolor='k')
for n_exp, exp in enumerate(exp_list):
    # one extra for the hindered version, to get smoother graphs
    n_space_steps = 4
    if unhindered:
        n_space_steps = 3
    batch, amorphous_material_initial, relative_humidity_gas_inlet, particle_diameter, weight, hours, resolution  = create_conditions(exp)

    # create conditions
    m_particle_tot_conc_initial, m_void_conc_initial, enthalpy_powder_initial, water_activity_void_initial,\
    m_void_initial, temp_initial_celsius, m_gas_sur, m_particle_cryst_initial, m_particle_am_initial, m_particle_tot_initial, \
    gas_density_initial, water_activity_sur, m_particle_cryst_conc_initial = \
        create_initial_conditions(relative_humidity_gas_inlet, amorphous_material_initial)
    seconds = hours * 60 * 60
    time_step = seconds/resolution
    discrete_time = np.linspace(0, seconds, resolution)
    if unhindered:
        computed_system = np.load(f'data_and_plots/211014-unhindered/exp {exp}, time {hours}, n {n_space_steps}.npy')
        save_string = f'data_and_plots/211014-unhindered/{exp_list}.pdf'
        computed_system = np.load(
            f'data_and_plots/211021-final/unhindered/exp {exp}, time {hours}, n {n_space_steps}.npy')
        save_string = \
            f'data_and_plots/211021-final/unhindered/{exp_list}.pdf'
    else:
        computed_system = np.load(
            f'data_and_plots/211021-final/hindered/exp {exp}, time {hours}, n {n_space_steps}.npy')
        save_string = \
            f'data_and_plots/211021-final/hindered/{exp_list}.pdf'

    # compute values from data; using for loop just to be able to hide it
    for do in range(1):
        total_water_all         = computed_system[:, 0:n_space_steps]
        amorphicity_all         = computed_system[:, n_space_steps:n_space_steps * 2]
        energy_all              = computed_system[:, n_space_steps * 2:n_space_steps * 3]
        water_activity_all      = computed_system[:, n_space_steps * 3:n_space_steps * 4]
        m_sur                   = computed_system[:, n_space_steps * 4:n_space_steps * 5]
        temp_all_celsius        = computed_system[:, n_space_steps * 5:n_space_steps * 6]
        m_particle_am_all       = computed_system[:, n_space_steps * 6:n_space_steps * 7]

        temp_all_kelvin         = temp_all_celsius + kelvin
        total_water_all_diffs   = total_water_all[1:, :] - total_water_all[:-1, :]

        m_void_all              = np.zeros([resolution, n_space_steps])
        m_void_all[0, :]        = m_void_initial

        extra_water_weight_all  = np.zeros([resolution, n_space_steps])
        m_particle_cryst_all    = np.zeros([resolution, n_space_steps])
        m_particle_cryst_all[0, :] = m_particle_cryst_initial

        m_powder_total_all      = np.zeros([resolution, n_space_steps])
        m_powder_total_all[0, :]= m_particle_tot_initial

        p_water_all             = np.zeros([resolution, n_space_steps])

        density_gas_all         = np.zeros([resolution, n_space_steps]) + gas_density_initial

        m_powder_conc_all       = np.zeros([resolution, n_space_steps])
        water_W_all             = np.zeros([resolution, n_space_steps])
        excess_water_curr             = np.zeros([resolution, n_space_steps])

        for t in range(resolution):
            m_void_all[t, :]            = compute_H_from_water_activity_temp(water_activity_all[t, :])

            m_particle_cryst_all[t, :]  = compute_GAB_equilibrium_moisture_cryst(water_activity_all[t, :])
            # m_particle_am_all[t, :]     = compute_GAB_equilibrium_moisture_am(water_activity_all[t, :], verbose=False)
            m_powder_total_all[t, :]    = m_particle_am_all[t, :] * amorphicity_all[t, :] + m_particle_cryst_all[t, :] * (1 - amorphicity_all[t, :])
            m_powder_conc_all[t, :]     = m_powder_total_all[t, :] * density_particle * (1 - porosity_powder)

            water_W_all[t, :]           = porosity_powder * m_void_all[t, :] * density_gas_all[t, :] + (1 - porosity_powder) * m_powder_total_all[t, :] * density_particle
            excess_water_curr[t, :]   = total_water_all[t, :] - water_W_all[t, :]

            # temp_all_kelvin[t, :], temp_all_celsius[t, :] = compute_temp_from_energy(m_void_all[t, :], m_powder_total_all[t, :], energy_all[t, :], excess_water_curr[t, :])
            density_gas_all[t, :]           = compute_density_air(m_void_all[t, :])

        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        m_void_conc_all         = m_void_all * porosity_powder * density_gas_all

        m_sur_tot = np.sum(m_sur, axis=1)/n_space_steps
        density_sur = compute_density_air(m_sur_tot)
        w_sur = m_sur_tot * density_sur

        total_water_all_tot = np.sum(total_water_all, axis = 1)/n_space_steps

        w_sur_diffs        = w_sur[1:] - w_sur[:-1]

        total_water_tot_diffs   = total_water_all_tot[1:] - total_water_all_tot[:-1]

        m_evap = w_sur_diffs + total_water_tot_diffs        # kg/m3
        m_evap = np.append(m_evap, 0)

        tg_all                  = compute_glass_temp_mix(1 - m_particle_am_all, glass_temp_lactose, glass_temp_water)
        t_tg_diff_all           = temp_all_kelvin - tg_all

        ############################ AVERAGES OVER WHOLE CYLINDER ##############################################################
        m_void_conc_avg = 0
        m_void_avg      = 0
        amorphicity_avg = 0
        total_water_avg = 0
        energy_avg      = 0
        density_gas_avg = 0
        temp_avg_kelvin = 0
        exc_water_avg   = 0
        m_particle_am_avg = 0
        water_activity_avg = 0
        total_water_all_diffs_avg = 0
        for n in range(n_space_steps):
            m_void_conc_avg += m_void_conc_all[:, n]
            m_void_avg      += m_void_all[:, n]
            amorphicity_avg += amorphicity_all[:, n]
            total_water_avg += total_water_all[:, n]
            energy_avg      += energy_all[:, n]
            density_gas_avg += density_gas_all[:, n]
            temp_avg_kelvin += temp_all_kelvin[:, n]
            exc_water_avg   += excess_water_curr[:, n]
            water_activity_avg += water_activity_all[:, n]
            # temp_celsius_sur_avg += temp_celsius_sur[:, n]
            total_water_all_diffs_avg   += total_water_all_diffs[:, n]
            m_particle_am_avg   += m_particle_am_all[:, n]

        m_void_conc_avg         /= n_space_steps
        m_void_avg              /= n_space_steps
        amorphicity_avg         /= n_space_steps
        total_water_avg         /= n_space_steps
        energy_avg              /= n_space_steps
        density_gas_avg         /= n_space_steps
        temp_avg_kelvin         /= n_space_steps
        exc_water_avg           /= n_space_steps
        water_activity_avg      /= n_space_steps
        m_particle_am_avg       /= n_space_steps
        total_water_all_diffs_avg/= n_space_steps
        total_water_all_diffs_avg = np.insert(total_water_all_diffs_avg, 0, 0)
        temp_avg_celsius        = temp_avg_kelvin - kelvin

        m_particle_cryst_avg    = compute_GAB_equilibrium_moisture_cryst(water_activity_avg)
        m_particle_tot_avg      = m_particle_am_avg * amorphicity_avg + m_particle_cryst_avg * (1-amorphicity_avg)
        m_powder_avg            = m_particle_am_avg * amorphicity_avg + m_particle_cryst_avg * (1 - amorphicity_avg)
        m_powder_conc_avg       = m_powder_avg * density_particle * (1 - porosity_powder)       # kg water/m3

        tg_avg                  = compute_glass_temp_mix(1 - m_particle_am_avg, glass_temp_lactose, glass_temp_water)
        t_tg_diff_avg           = temp_avg_kelvin - tg_avg

        m_powder_diffs_avg      = m_powder_avg[1:] - m_powder_avg[:-1]                   # fraction, kg water/kg dry
        m_powder_diffs_avg      = m_powder_conc_avg[1:] - m_powder_conc_avg[:-1]         # kg water/m3
        m_powder_diffs_avg      = np.insert(m_powder_diffs_avg, 0, m_powder_diffs_avg[0])

        am_am_diffs = amorphicity_avg[1:] - amorphicity_avg[:-1]    # fraction
        am_am_diffs = np.insert(am_am_diffs, 0, 0, axis=0)
        am_cryst_diffs = - am_am_diffs                              # fraction

    ################################################# HEAT #############################################################
    # to compute heat of sorption is only interesting up until crystallization starts; where sorption is positive and
    # water activity is not greater than salt solution RH
    m_powder_diffs_avg = np.where(m_powder_diffs_avg < 0, 0, m_powder_diffs_avg)
    m_powder_diffs_avg = np.where(water_activity_avg > water_activity_sur, 0, m_powder_diffs_avg)

    jump_step = 1
    heat_of_sorption_vector = m_powder_diffs_avg * heat_binding_diffs/(1000 * time_step)        # J/kg to J/g
    heat_of_sorption_vector = heat_of_sorption_vector[::jump_step]
    heat_of_evap = -m_evap * heat_of_evaporation_water/(1000 * time_step)                       # kg/m3 * J/kg

    heat_of_cryst_vector = am_cryst_diffs * heat_of_crystallization / (1000 * time_step)        # J/kg to J/(g s)
    heat_of_cryst_vector = heat_of_cryst_vector[::jump_step]

    # Use filters to smooth uneven graphs. These settings seem to work without affecting the curves too much.
    if unhindered:
        heat_of_sorption_vector = savgol_filter(heat_of_sorption_vector, 151, 1)
        heat_of_cryst_vector = savgol_filter(heat_of_cryst_vector, 5, 2)
    else:
        heat_of_sorption_vector = savgol_filter(heat_of_sorption_vector, 151, 1)
        heat_of_cryst_vector = savgol_filter(heat_of_cryst_vector, 5, 2)

    total_energy_vector = heat_of_cryst_vector + heat_of_sorption_vector
    if unhindered and not weight_c:
        total_energy_vector = savgol_filter(total_energy_vector, 41, 2)
    if unhindered == False:
        total_energy_vector = savgol_filter(total_energy_vector, 41, 2)

    discrete_time = discrete_time[::jump_step]
    minutes = discrete_time / 60
    seconds = minutes * 60
    x_start = np.where(minutes >= 15)[0][0]

    # choose which entity to plot
    # plt.plot(minutes, heat_of_cryst_vector, color=c_chosen[n_exp], label=f'Sim {exp}')
    # plt.plot(minutes, heat_of_sorption_vector, color=c_chosen[n_exp], label=f'Sim {exp}')
    plt.plot(minutes, total_energy_vector, color=c_chosen[n_exp], label=f'Sim {exp}: {labels[n_exp]}')
    integral = simps(total_energy_vector[x_start:], seconds[x_start:])
    plt.fill_between(minutes[x_start:], 0, total_energy_vector[x_start:], color=c_chosen[n_exp], alpha=0.2, label=f'Integral {integral:0.1f} J/g')
    x_index = np.where(heat_of_cryst_vector == np.max(heat_of_cryst_vector))[0]

    # print stuff where it looks good
    if humidity_c:
        if unhindered:
            plt.xlim(0, 160)
            plt.text(120 - 43 * (n_exp + 1 / 2), np.max(total_energy_vector[100:]) + 0.0002,
                     f'Peak {np.max(total_energy_vector[100:]):.2e} after {int(minutes[x_index])} minutes',
                     color=c_chosen[n_exp])
        else:
            plt.xlim(0, 300)

    elif weight_c:
        if unhindered:
            plt.text(120 * (1+ 1.5 * n_exp), np.max(total_energy_vector) + 0.0002, f'Peak {np.max(total_energy_vector):.2e} after {int(minutes[x_index])} minutes', color=c_chosen[n_exp])
            plt.xlim(0, 500)
        else:
            plt.text(120 * (1 + 1.5 * n_exp), np.max(total_energy_vector) + 0.0001,
                     f'Peak {np.max(total_energy_vector):.2e} after {int(minutes[x_index])} minutes',
                     color=c_chosen[n_exp])
            plt.xlim(0, 500)

    elif amorphicity_c:
        if unhindered:
            plt.xlim(0, 120)
            plt.text(20, np.max(total_energy_vector[100:]) + 0.0002,
                     f'Peak {np.max(total_energy_vector[100:]):.2e} after {int(minutes[x_index])} minutes',
                     color=c_chosen[n_exp])
        else:
            plt.xlim(0, 150)
            plt.text(40, np.max(total_energy_vector[100:]) + 0.0001,
                     f'Peak {np.max(total_energy_vector[100:]):.2e} after {int(minutes[x_index])} minutes',
                     color=c_chosen[n_exp])

    if humidity_c and unhindered == False:
        plt.hlines(np.max(total_energy_vector[300:]), 0, minutes[x_index], color=c_chosen[n_exp], linestyles='--')

    elif low_hum_c == False:
        plt.hlines(np.max(total_energy_vector[100:]), 0, minutes[x_index], color=c_chosen[n_exp], linestyles='--')

########################################## PLOT ########################################################################
plt.xlabel('Time, min')
plt.ylabel('Heat flow, W/g')
plt.legend()
plt.grid(color=grid_color, alpha=grid_alpha_min)

plt.savefig(save_string)
plt.show()
plt.close()

# function that plots 9 different entities, and several discretized sections: not very clear graphically, but good to
# investigate results further
# plot_tam()

# function that plots 4 interesting entities, clear and visual
# plot_tam_4()


