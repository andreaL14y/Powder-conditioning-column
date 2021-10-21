from TAM_functions import *
glass_temp_lactose2, glass_temp_lactose3 = 101 + kelvin, 97 + kelvin
############################################ FITTING ###################################################################
def plot_GAB():
    water_activities = np.linspace(0, 1, 1000)
    m_cryst = compute_GAB_equilibrium_moisture_cryst(water_activities)
    m_cryst_cheat = compute_GAB_equilibrium_moisture_cryst(water_activities, H_H=False)

    m_am = compute_GAB_equilibrium_moisture_am(water_activities, limit_aw = False)
    max_index = np.where(m_am >= 1)[0][0]
    wa_max = water_activities[max_index]
    fig, axs = plt.subplots(2, 1, figsize=(10, 12), dpi=100)

    change_index = np.where(m_cryst_cheat + m_cryst > 0.0015)[0][0]

    major_ticks = np.arange(0, 1.1, 0.1)
    axs[0].set_xticks(major_ticks)
    axs[1].set_xticks(major_ticks)
    axs[1].set_yticks(major_ticks)

    # fig.suptitle('Equilibrium moisture contents for crystalline and amorphous lactose by GAB')
    axs[0].plot(water_activities, m_cryst, label='Crystalline lactose with modification', color=c_palette[3])
    axs[0].plot(water_activities, m_cryst_cheat, label='Crystalline lactose without modification; H\' = H\'\' = 1',
                color=c_palette[2], linestyle='--')
    axs[0].vlines(water_activities[change_index], 0, m_cryst[change_index], color=c_palette[6], linestyle='--')
    axs[0].plot([water_activities[change_index]], m_cryst[change_index], 'x', ms=8, color=c_palette[6])

    ins = axs[0].inset_axes([0.7, 0.3, 0.2, 0.6])
    ins.plot(water_activities, m_cryst, color=c_palette[3])
    ins.plot(water_activities, m_cryst_cheat, color=c_palette[2], linestyle='--')
    ins.plot([water_activities[change_index]], m_cryst[change_index], 'x', ms=8, color=c_palette[6])
    ins.set_xlim(0.8, 1.01)
    ins.set_ylim(0.0004, 0.002)
    ins.grid(alpha=grid_alpha_min, color=grid_color)

    axs[0].legend()
    axs[0].grid(which='both', alpha=grid_alpha_min, color=grid_color)
    # axs[0].set_title('Crystalline lactose')

    axs[1].plot(water_activities, m_am, label='Amorhpous lactose', color=c_palette[2])

    axs[1].plot([wa_max], [m_am[max_index]], 'x', ms=8, color=c_palette[6])
    axs[1].vlines(wa_max, 0, m_am[max_index], linestyles='--', color=c_palette[6])
    axs[1].text(0.75, 1.01, f'{wa_max:.2f}', color=c_palette[0])  # , size='large'
    axs[1].legend()
    axs[1].grid(which='both', alpha=grid_alpha_min, color=grid_color)
    axs[1].set_ylim(0, 1.1)
    # axs[1].set_title('Amorhpous lactose')

    axs[0].set_xlabel('Water activity')
    axs[1].set_xlabel('Water activity')

    axs[0].set_ylabel('Moisture content, g water/g lactose')
    axs[1].set_ylabel('Moisture content, g water/g lactose')

    # plt.savefig('GAB_cryst_vs_am.pdf')
    plt.show()


def plot_glass_transition_curve(plot_experiments):
    weight_fractions_lactose = np.linspace(0, 1, 100)
    glass_temps = compute_glass_temp_mix(weight_fractions_lactose, glass_temp_lactose, glass_temp_water)

    ################################### EXP 1 ######################################
    rh_25 = np.array([0.11, 0.33, 0.3817, 0.5289, 0.5757, 0.753])
    n_RH = len(rh_25)

    moisture_contents_25 = compute_GAB_equilibrium_moisture_am(rh_25)  # water content g/100 g
    temp_25 = np.zeros(n_RH) + 25  # + kelvin

    glass_temps_C = glass_temps - kelvin
    glass_temp_lactose_C = glass_temp_lactose - kelvin
    glass_temp_water_C = glass_temp_water - kelvin

    fig = plt.figure(figsize=fig_size, dpi=100, facecolor='w', edgecolor='k')

    plt.plot(weight_fractions_lactose, glass_temps_C, color=c_palette[-1], label=r'$T_G$ solution')

    if plot_experiments:
        plt.plot(weight_fractions_lactose, glass_temps_C + 20, color=c_palette[4], label=r'$T_G$ + 20')
        plt.fill_between(weight_fractions_lactose, glass_temps_C, glass_temps_C + 20, color=c_palette[5])

        plt.plot(1 - moisture_contents_25, temp_25, 'o', markerfacecolor=c_palette[1], markeredgecolor='goldenrod',
                 fillstyle='full')
        plt.hlines(temp_25[0], 0, 1, linestyles='--', linewidth=1, color=c_palette[1],
                   label=f'$T$ experiment: {temp_25[0]:.0f}$^\circ$C')
        x = 1 - moisture_contents_25  # weight fraction lactose
        y_start = compute_glass_temp_mix(x, glass_temp_lactose, glass_temp_water) - kelvin
        plt.vlines(x, y_start, temp_25, linestyles='-', linewidth=0.5, color='darkgoldenrod')
        np.set_printoptions(formatter={'float': '{: .0f}'.format})
        print('At T', temp_25[0], r'distances $T-T_G$ are:')
        print(temp_25 - y_start, 'for RHs:')
        np.set_printoptions(formatter={'float': '{: .0f}'.format})
        print(rh_25 * 100)

    ax = fig.add_subplot(1, 1, 1)

    # Major ticks every 20, minor ticks every 5
    major_ticks = np.linspace(-1, 1.5, 11)
    minor_ticks = np.linspace(-1, 1.5, 101)
    major_ticks_y = np.arange(-140, 120, 20)
    minor_ticks_y = np.arange(-140, 120, 10)

    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks_y)
    ax.set_yticks(minor_ticks_y, minor=True)

    # Or if you want different settings for the grids:
    ax.grid(which='minor', alpha=grid_alpha_min, color=grid_color)
    ax.grid(which='major', alpha=grid_alpha_min, color=grid_color)

    plt.ylabel('Temperature, $^\circ$C')
    plt.xlabel('Weight fraction lactose')

    plt.hlines(glass_temp_lactose_C, -1, 2, linestyles='--', linewidth=1, color=c_palette[-3])
    plt.hlines(glass_temp_water_C, -1, 2, linestyles='--', linewidth=1, color=c_palette[-3])
    if plot_experiments:
        # plt.title('Experiment plan on $T-T_G$ map')
        plt.xlim(0.6, 1)
        plt.ylim(-140, 120)

    else:
        # plt.title('Glass transition temperature for lactose as a function of moisture content')
        plt.xlim(-0.02, 1.02)
        plt.text(0.2, glass_temp_lactose_C - 9, r'$T_G$ lactose', color=c_palette[-3])
        plt.text(0.9, glass_temp_water_C + 3, r'$T_G$ water', color=c_palette[-3])
        plt.ylim(-140, 120)
    plt.legend(loc='upper left')
    plt.show()

    # save file
    # if plot_experiments:
    #     fig.savefig('experiment_TTG_plot.pdf')
    # elif plot_all_tg:
    #     fig.savefig('different_tg_lactose.pdf')
    # else:
    #     fig.savefig('T_G-plot.pdf')


def plot_kinetics(am_amorph_initial, multiple_tg = False):
    tg_lactose_list = np.array([glass_temp_lactose, glass_temp_lactose2, glass_temp_lactose3])

    # aw_values = np.arange(0.46, 0.62, 0.02)
    aw_values = np.arange(0.54, 0.62, 0.02)
    n_TTG = len(aw_values)
    times_TTG = np.zeros([3, n_TTG])
    TTG_values = np.zeros([3, n_TTG])
    for i, tg_lactose in enumerate(tg_lactose_list):
        m_am = compute_GAB_equilibrium_moisture_am(aw_values)
        TGs = compute_glass_temp_mix(1 - m_am, tg_lactose, glass_temp_water)

        TTG_values[i, :] = (25 + kelvin) - TGs
        print('\nTTG values are:', TTG_values[i, :])

        for n, ttg in enumerate(TTG_values[i, :]):
            am_amorph = am_amorph_initial
            t = 0
            while am_amorph > 0.01 * am_amorph_initial:
                Y = am_amorph / am_amorph_initial
                Y = am_amorph

                am_diff = compute_kinetics_avrami(Y, ttg, single = True)[0]
                am_amorph += am_diff
                t += 1
            times_TTG[i, n] = t
            print(t, 's')

    fig = plt.figure(figsize=(10, 5), dpi=100)
    col_list = [0, 2, 3]

    if multiple_tg == True:
        number_sims = 3
        save_string = 'lactose_crystallization_times_tg_comparison.pdf'
    else:
        number_sims = 1
        save_string = 'lactose_crystallization_times_mine.pdf'

    for i in range(number_sims):
        TTG = TTG_values[i, :]
        times = times_TTG[i, :]
        plt.plot(aw_values, TTG, color=c_palette[4])
        plt.plot(aw_values, TTG, 'o', color=c_palette[col_list[i]],
                 label=f'$T_G$ lactose: {tg_lactose_list[i] - kelvin:.0f}$^\circ$C')

        dist = 0.005
        if tg_lactose_list[i] == np.max(tg_lactose_list):
            plt.vlines(aw_values, TTG, TTG + 20, linestyles=':', color=c_palette[7])

        dist2 = 4 * i + 12
        for n in range(n_TTG):
            if int(times[n] / 3600) > 24:
                plt.text(aw_values[n] - dist, TTG_values[-1, n] + dist2, f'{times[n]/(24 * 3600):0.0f} days')
            elif int(times[n] / 3600) > 0:
                plt.text(aw_values[n] - dist, TTG_values[-1, n] + dist2, f'{times[n]/3600:0.0f} h')
            elif int(times[n] / 60) > 0:
                plt.text(aw_values[n] - dist, TTG_values[-1, n] + dist2, f'{times[n]/60:0.0f} min')
            else:
                plt.text(aw_values[n] - dist, TTG_values[-1, n] + dist2, f'{times[n]:0.0f} s')

    ax = fig.add_subplot(1, 1, 1)

    # different settings for the grids:
    ax.grid(which='minor', alpha=grid_alpha_min, color=grid_color)
    ax.grid(which='major', alpha=grid_alpha_min, color=grid_color)

    plt.xlim(aw_values[0] - 0.01, aw_values[-1] + 0.01)
    plt.ylim(TTG_values[0, 0] - 5, TTG_values[-1, -1] + 25)
    plt.xlabel(r'Water activity, $a_w$')
    plt.ylabel(r'$T - T_G$')
    # plt.title('Time until 99% of fully amorphous lactose has crystallized \n at 25$^\circ$C and different values of $T - T_G$')
    plt.legend(loc='upper left')
    plt.show()
    # fig.savefig(save_string)


def plot_kinetics_K(am_amorph):
    min = -10
    max = 40
    n = 200
    temp_diffs = np.linspace(min, max, max - min + 1)
    temp_diffs = np.linspace(min, max, n)
    am_change = np.zeros(n)
    K = np.zeros(n)
    for t in range(n):
        am_change[t], K[t] = compute_kinetics_avrami(am_amorph, temp_diffs[t], single=True)

    fig = plt.figure(figsize=(10, 5), dpi=100)
    plt.plot(temp_diffs, np.log(K), color=c_palette[-3])

    plt.grid(alpha=grid_alpha_min, color=grid_color)
    # plt.title(f'Crystallization kinetics parameter K as a function of $T-T_G$, for X = {am_amorph_initial}')
    plt.ylabel(r'$log(K)$')
    plt.xlabel(r'$T - T_G$')
    # plt.savefig('kinetics_K_of_temperature_diffs.pdf')
    plt.show()

# plot_GAB()
# plot_glass_transition_curve(plot_experiments=True)
# plot_kinetics(multiple_tg=True, am_amorph_initial = 0.21)
# plot_kinetics_K(am_amorph = 0.21)