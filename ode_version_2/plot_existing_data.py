from functions import*
from initial_conditions import*
n_space_steps = 3
hours = 2.1
hours = 8
seconds = hours * 60 * 60
resolution = 4500
# resolution = 50
time_step = seconds/resolution
discrete_time = np.linspace(0, seconds, resolution)

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
    axs[1, 1].set_ylim(-0.1, np.max(amorphicity_avg) + 0.1)
    axs[1, 1].set_title('Amount amorhpous material')
    axs[1, 1].set(xlabel='Time, hours', ylabel='Am material, kg/kg lactose')

    ################################### Water Activity RH ##############################################################
    axs[2, 1].plot(time, water_activity_avg, label='Avg w_a', color=am_color)
    axs[2, 1].plot(time, water_activity_all[:, 0], label='1', color='limegreen', linestyle='--')
    axs[2, 1].plot(time, water_activity_all[:, 1], label=f'{2}',
                   color='lightgreen', linestyle='--')
    # axs[2, 1].plot(time, water_activity_all[:, 2 * int(n_space_steps / 3)], label=f'{1+2*int(n_space_steps/3)}',
    #                color='limegreen', linestyle='--')
    axs[2, 1].plot(time, water_activity_all[:, 2], label=f'{3}', color='darkgreen',
                   linestyle='--')
    axs[2, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 1].set_xlim(1.5, x_max)
    axs[2, 1].set_ylim(0.5, 1.1)
    axs[2, 1].set_title('Water activity over time')
    axs[2, 1].set(xlabel='Time, hours', ylabel='wa')
    axs[2, 1].legend()

    ################################### Heat of sorption ###############################################################
    axs[0, 2].plot(time, heat_of_sorption_vector, label='Heat of sorption', color='rosybrown', linestyle='--')
    # axs[0, 2].plot(time, energy_avg, label='Energy', color='rosybrown')
    # axs[0, 2].plot(time, energy_all[:, 0], label='Energy 0', color='goldenrod', linestyle='--')
    # axs[0, 2].plot(time, energy_all[:, n_space_steps-1], label='Energy 3', color='fuchsia', linestyle='--')
    # axs[0, 2].plot(time, exc_water_avg, label='Excess water', color='rosybrown')
    # axs[0, 2].plot(time, excess_water_curr[:, 0], label='Excess water 1', color='goldenrod', linestyle='--')
    # axs[0, 2].plot(time, excess_water_curr[:, n_space_steps - 1], label='Excess water 3', color='fuchsia', linestyle='--')
    # axs[0, 2].plot(time, -evaporation_energy, label='Evaporation', color=temp_color)
    axs[0, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 2].set_xlim(x_min, x_max)
    axs[0, 2].set_ylim(-0.05, 0.03)
    axs[0, 2].set_title('Heat generated')
    axs[0, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')
    axs[0, 2].legend()

    ################################# Heat of crystallization ##########################################################
    # axs[1, 2].plot(time, exc_water_avg, label='Excess water', color='lightsteelblue') #, linestyle='--')
    # axs[1, 2].plot(time, energy_avg, label='Heat flow', color='navy') #, linestyle='--')
    # axs[1, 2].plot(time, heat_flow_vector, label='Heat flow', color='navy') #, linestyle='--')
    axs[1, 2].plot(time, heat_of_cryst_vector, label='Heat of cryst', color='navy', linestyle='--')
    # axs[1, 2].plot(time, energy_gas, label='Gas', color='lightblue') #, linestyle='--')
    # axs[1, 2].plot(time, energy_p, label='Particle', color='rosybrown') #, linestyle='--')
    # axs[1, 2].plot(time, total_energy_vector, label='Total', color=temp_color)
    axs[1, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 2].set_xlim(x_min, x_max)
    # axs[1, 2].set_ylim(-0.004, 0.03)
    axs[1, 2].set_title('Heat generated')
    axs[1, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')
    axs[1, 2].legend()

    ####################################### Temperature ################################################################
    # max_cryst = np.max(heat_of_cryst_vector)
    # axs[2, 2].set_ylim(0, max_cryst)
    # axs[2, 2].plot(time, temp_avg_celsius,                      label='Avg temp', color=temp_color)
    axs[2, 2].plot(time, energy_avg,                      label='Avg energy', color=temp_color)
    axs[2, 2].plot(time, energy_all[:, 0],                label='1', color='goldenrod', linestyle='--')
    # axs[2, 2].plot(time, energy_all[:, 1],                label='1', color='goldenrod', linestyle='--')
    # axs[2, 2].plot(time, energy_all[:, 2],                label='1', color='red', linestyle='--')
    # axs[2, 2].plot(time, temp_all_celsius[:, 0],                label='1', color='goldenrod', linestyle='--')
    # axs[2, 2].plot(time, temp_all_celsius[:, n_space_steps-1],  label=f'{n_space_steps}', color='palegoldenrod', linestyle='--')
    # axs[2, 2].hlines(temp_initial_celsius, x_min, x_max, label=f'{temp_initial_celsius}', color='firebrick', linestyle='--')
    axs[2, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 2].set_xlim(x_min, x_max)
    axs[2, 2].set_title('Temperature')
    axs[2, 2].set(xlabel='Time, hours', ylabel='T, deg C')
    axs[2, 2].legend()

    plt.savefig(f'data_and_plots/wa {relative_humidity_gas_inlet}, T {temp_initial_celsius}, am {amorphous_material_initial}, time {hours}, n {n_space_steps}.pdf')     # 1% am
    plt.show()


computed_system = np.load(f'data_and_plots/wa 0.574, T 10, am 0.16, time {hours}, n 3.npy')

total_water_all         = computed_system[:, 0:n_space_steps]
amorphicity_all         = computed_system[:, n_space_steps:n_space_steps*2]
energy_all              = computed_system[:, n_space_steps*2:n_space_steps*3]
water_activity_all      = computed_system[:, n_space_steps*3:n_space_steps*4]
temp_all_celsius        = computed_system[:, n_space_steps * 5:n_space_steps * 6]
temp_all_kelvin         = temp_all_celsius + kelvin
total_water_all_diffs = total_water_all[1:, :] - total_water_all[:-1, :]

m_void_all              = np.zeros([resolution, n_space_steps])
m_void_all[0, :]        = m_void_initial

extra_water_weight_all  = np.zeros([resolution, n_space_steps])
m_particle_cryst_all    = np.zeros([resolution, n_space_steps])
m_particle_cryst_all[0, :] = m_particle_cryst_initial

m_particle_am_all       = np.zeros([resolution, n_space_steps])
m_particle_am_all[0, :] = m_particle_am_initial

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
    m_particle_am_all[t, :]     = compute_GAB_equilibrium_moisture_am(water_activity_all[t, :], verbose=False)
    m_powder_total_all[t, :]    = m_particle_am_all[t, :] * amorphicity_all[t, :] + m_particle_cryst_all[t, :] * (1 - amorphicity_all[t, :])
    m_powder_conc_all[t, :]     = m_powder_total_all[t, :] * density_particle * (1 - porosity_powder)

    water_W_all[t, :]           = porosity_powder * m_void_all[t, :] * density_gas_all[t, :] + (1 - porosity_powder) * m_powder_total_all[t, :] * density_particle
    excess_water_curr[t, :]   = total_water_all[t, :] - water_W_all[t, :]

    temp_all_kelvin[t, :], temp_all_celsius[t, :] = compute_temp_from_energy(m_void_all[t, :], m_powder_total_all[t, :], energy_all[t, :], excess_water_curr[t, :])
    density_gas_all[t, :]           = compute_density_air(m_void_all[t, :])


m_void_conc_all         = m_void_all * porosity_powder * density_gas_all

tg_all                  = compute_glass_temp_mix(1 - m_particle_am_all, glass_temp_lactose, glass_temp_water_1)
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
temp_celsius_sur_avg = 0
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


print(temp_celsius_sur_avg)
m_void_conc_avg         /= n_space_steps
m_void_avg              /= n_space_steps
amorphicity_avg         /= n_space_steps
total_water_avg         /= n_space_steps
energy_avg              /= n_space_steps        # J/m3
# np.set_printoptions(formatter={'float': '{: .2f}'.format})
# print(energy_avg[discrete_time > 6600]/100000 - 20)
density_gas_avg         /= n_space_steps
temp_avg_kelvin         /= n_space_steps
exc_water_avg           /= n_space_steps
water_activity_avg      /= n_space_steps
# temp_celsius_sur_avg    /= n_space_steps
total_water_all_diffs_avg           /= n_space_steps
total_water_all_diffs_avg = np.insert(total_water_all_diffs_avg, 0, 0)
temp_avg_celsius        = temp_avg_kelvin - kelvin


m_particle_cryst_avg    = compute_GAB_equilibrium_moisture_cryst(water_activity_avg)
m_particle_am_avg       = compute_GAB_equilibrium_moisture_am(water_activity_avg)
m_particle_tot_avg      = m_particle_am_avg * amorphicity_avg + m_particle_cryst_avg * (1-amorphicity_avg)
m_powder_avg            = m_particle_am_avg * amorphicity_avg + m_particle_cryst_avg * (1 - amorphicity_avg)
m_powder_conc_avg       = m_powder_avg * density_particle * (1 - porosity_powder)

temp_avg_kelvin, temp_avg_celsius = compute_temp_from_energy(m_void_avg, m_particle_tot_avg, energy_avg, exc_water_avg)

tg_avg                  = compute_glass_temp_mix(1 - m_particle_am_avg, glass_temp_lactose, glass_temp_water_1)
t_tg_diff_avg           = temp_avg_kelvin - tg_avg

m_powder_diffs_avg      = m_powder_avg[1:] - m_powder_avg[:-1]                   # fraction, kg water/kg dry
m_powder_diffs_avg      = np.insert(m_powder_diffs_avg, 0, 0)

am_am_diffs = amorphicity_avg[1:] - amorphicity_avg[:-1]
am_am_diffs = np.insert(am_am_diffs, 0, 0, axis=0)
am_cryst_diffs = - am_am_diffs

############################################# TEMP AVG #################################################################
heat_capacity_bed = porosity_powder * heat_capacity_air + (1-porosity_powder) * heat_capacity_particle

# temp_diffs = temp_avg_celsius[1:] - temp_avg_celsius[:-1]
# temp_diffs = np.insert(temp_diffs, 0, 0)

heat_capacity_medium = 7

# temp_diffs = temp_avg_celsius - temp_avg_celsius[0]
# heat_flow_vector = temp_diffs * heat_capacity_medium/(1000 * time_step)

# heat_of_sorption = compute_heat_of_sorption(water_activity_avg, temp_avg_kelvin)   # J/kg
# heat_of_sorption_vector = m_powder_diffs_avg * (heat_of_sorption) / time_step               # J/(kg s)
# heat_of_sorption_vector /= 1000                                                         # J/(kg s) to J/(g s)

density_bed = density_powder + density_water * m_particle_tot_avg


# heat_flow_vector = conductivity_tot * contact_area_powder * temp_diffs / (weight * total_height_powder/2)   # J/kg
# heat_flow_vector /= 1000                                                                                    # J/g

heat_of_sorption_vector = (energy_avg - energy_avg[0])/(density_powder * 1000 * time_step)          # J/m3 to J/kg to J/g
heat_of_sorption_vector = (energy_avg[1:] - energy_avg[:-1])/(density_powder * 1000 * time_step)
heat_of_sorption_vector = np.insert(heat_of_sorption_vector, 0, 0)

heat_of_sorption_vector = m_powder_diffs_avg * heat_binding_diffs/(1000 * time_step)          # J/m3 to J/kg to J/g
# heat_of_sorption_vector = np.where(water_activity_avg > 0.52, 0, heat_of_sorption_vector)

evaporation_energy = total_water_all_diffs_avg * heat_of_evaporation_water/1000
# heat_of_sorption_vector -= evaporation_energy

heat_of_cryst_vector = am_cryst_diffs * (1-porosity_powder) * (heat_of_crystallization/1000) / time_step      # J/kg to J/(g s)
heat_of_cryst_vector = am_cryst_diffs * (heat_of_crystallization/1000) / time_step      # J/kg to J/(g s)
# heat_of_sorption_vector -= heat_of_cryst_vector


energy_gas = heat_capacity_air * temp_avg_celsius + m_void_avg * (heat_capacity_vapor * temp_avg_celsius + heat_of_evaporation_water)
energy_gas = (energy_gas - energy_gas[0])/(1000* time_step)

energy_p = heat_capacity_particle * temp_avg_celsius + m_particle_tot_avg * heat_capacity_water * temp_avg_celsius
energy_p = (energy_p - energy_p[0])/(1000 * time_step)

total_energy_vector = heat_of_cryst_vector + heat_of_sorption_vector

########################################## PLOT ########################################################################
# amorphicity_avg = normalize_data(amorphicity_avg, False)
# amorphicity_all = normalize_data(amorphicity_all, False)
plot_tam()