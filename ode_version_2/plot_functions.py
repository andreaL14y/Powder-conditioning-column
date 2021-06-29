from functions import *

############################################ FITTING ###################################################################
def plot_GAB():
    eq_air_am = equilibrium_curve_fit(moistures_am, alpha_parameter_am, N_am)
    fig, axs = plt.subplots(2, 2, figsize=(15, 15))

    fig.suptitle('Equilibrium GAB vs fitted curves')
    axs[0, 0].plot(relative_humidities, moistures_cryst, label='Moisture content')
    axs[0, 0].plot(eq_air_cryst, moistures_cryst, label='Fitted curve')
    axs[0, 0].legend()
    axs[0, 0].grid()
    axs[0, 0].set_title('Crystalline lactose')

    axs[1, 0].plot(relative_humidities, moistures_am, label='Moisture content')
    axs[1, 0].plot(eq_air_am, moistures_am, label='Fitted curve')
    axs[1, 0].legend()
    axs[1, 0].grid()
    axs[1, 0].set_ylim(0, 0.9)
    axs[1, 0].set_title('Amorhpous lactose')

    axs[1, 1].plot(relative_humidities, moistures_am_der, label='Moisture content')
    axs[1, 1].legend()
    axs[1, 1].grid()
    axs[1, 1].set_ylim(0, 0.9)
    axs[1, 1].set_title('Amorhpous lactose derivative')

    plt.xlabel('RH')
    plt.ylabel('Moisture, g/g lactose')

    plt.show()


def plot_glass_transition_curve():
    ############################################### DATA ###################################################################
    weight_fraction_lactose = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
    temperature = np.array([-135, -127, -112, -90, -50, 101]) + kelvin
    room_temp = 20 + kelvin

    ######################################## GORDON & TAYLOR ###############################################################
    weight_fractions_lactose = np.linspace(0, 1, 100)
    glass_temps_1 = compute_glass_temp_mix(weight_fractions_lactose, glass_temp_lactose, glass_temp_water_1)

    ############################################# DATA #####################################################################
    relative_humidities = np.array([0.37, 0.538, 0.662, 0.764, 0.858])
    n_RH = len(relative_humidities)

    moisture_contents = compute_GAB_equilibrium_moisture_am(relative_humidities)  # water content g/100 g
    exp_temp = np.zeros(n_RH) + 24 + kelvin

    moisture_content_RH30 = compute_GAB_equilibrium_moisture_am(0.3)
    moisture_content_RH80 = compute_GAB_equilibrium_moisture_am(0.8)

    RH30_temp = np.array([20, 60]) + kelvin
    RH80_temp = RH30_temp

    ############################################# PLOT #####################################################################
    fig = plt.figure(figsize=(14, 10), dpi=100, facecolor='w', edgecolor='k')

    # Plot glass temp as function of moisture content
    plt.plot(weight_fractions_lactose, glass_temps_1, color='brown', label=r'$T_G$ solution')

    # Add horizontal lines for some temperatures
    plt.hlines(glass_temp_lactose, 0, 1, linestyles='--', linewidth=1,
               color='brown')  # , label=f'T_G lactose: {glass_temp_lactose} K')
    plt.hlines(glass_temp_water_1, 0, 1, linestyles='--', linewidth=1,
               color='brown')  # , label=f'T_G water: {glass_temp_water_1} K')
    plt.hlines(room_temp, 0, 1, linestyles='--', linewidth=1,
               color='steelblue')  # , label=f'Room temp.: {room_temp} K = {room_temp - kelvin} C')
    plt.hlines(0 + kelvin, 0, 1, linestyles='--', linewidth=1, color='steelblue')  # , label='Freezing point water')

    # Add some data from article at constant temp but different RH
    plt.plot(1 - moisture_contents, exp_temp, 'o', markerfacecolor='black', markeredgecolor='goldenrod',
             fillstyle='full',
             label='RH: 85.8, 76.4, 66.2, 53.8, 37')

    # Add data from figure 4
    x = 1 - moisture_content_RH30
    y_start = RH30_temp[0]
    plt.vlines(x, y_start, RH30_temp[1], linestyles='-', linewidth=0.5, color='darkgoldenrod')

    x = 1 - moisture_content_RH80
    y_start = compute_glass_temp_mix(x, glass_temp_lactose, glass_temp_water_1)
    plt.vlines(x, y_start, RH30_temp[1], linestyles='-', linewidth=0.5, color='darkgoldenrod')

    size_m = 7
    plt.plot(1 - moisture_content_RH30, RH30_temp[0], 's', markerfacecolor='yellowgreen',
             markeredgecolor='darkolivegreen', fillstyle='full', markersize=size_m, label='Fig 4: RH 30, T 20 C')
    plt.plot(1 - moisture_content_RH30, RH30_temp[1], '^', markerfacecolor='red', markeredgecolor='darkred',
             fillstyle='full', markersize=size_m, label='Fig 4: RH 30, T 60 C')
    plt.plot(1 - moisture_content_RH80, RH80_temp[0], 'o', markerfacecolor='grey', markeredgecolor='black',
             fillstyle='full', markersize=size_m, label='Fig 4: RH 80, T 20 C')
    plt.plot(1 - moisture_content_RH80, RH80_temp[1], 'D', markerfacecolor='fuchsia', markeredgecolor='m',
             fillstyle='full', markersize=size_m, label='Fig 4: RH 80, T 60 C')

    x = 1 - moisture_contents
    y_start = compute_glass_temp_mix(x, glass_temp_lactose, glass_temp_water_1)
    plt.vlines(x, y_start, 24 + kelvin, linestyles='-', linewidth=0.5, color='darkgoldenrod')

    # Grid
    ax = fig.add_subplot(1, 1, 1)

    # Major ticks every 20, minor ticks every 5
    major_ticks = np.linspace(0, 1, 11)
    minor_ticks = np.linspace(0, 1, 101)
    major_ticks_y = np.arange(130, 400, 50)
    minor_ticks_y = np.arange(130, 400, 10)

    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks_y)
    ax.set_yticks(minor_ticks_y, minor=True)

    # Or if you want different settings for the grids:
    ax.grid(which='minor', alpha=0.1)
    ax.grid(which='major', alpha=0.4)

    # plt.grid()
    plt.xlim(0, 1)
    plt.ylim(130, 390)
    plt.ylabel('Temperature, K')
    plt.xlabel('Weight fraction lactose')
    plt.title('Glass transition temperature as a function of fraction of lactose in solution \n '
              '+ some equilibrium states for operating conditions from literature')
    plt.text(0.9, glass_temp_lactose - 7, r'$T_G$ lactose', color='brown')
    plt.text(0.9, glass_temp_water_1 + 3, r'$T_G$ water', color='brown')
    plt.text(0.02, 296, f'Room temperature, {room_temp} K', color='steelblue')
    plt.text(0.02, 266, f'Freezing point water, {0-kelvin} K', color='steelblue')
    plt.legend(loc=2)
    plt.show()
    # fig.savefig('T_G-plot.pdf')


# plot_GAB()
plot_glass_transition_curve()