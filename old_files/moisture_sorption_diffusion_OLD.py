from initial_conditions import *
import time

######################################## FUNCTIONS #####################################################################
def compute_system(initials, t, n_space_steps, step_length, verbose=False):
    global counter
    global while_time
    # global m_particle_tot_conc
    counter += 1
    divider = 2000

    total_water = initials[0:n_space_steps]
    amount_am = initials[n_space_steps:n_space_steps * 2]
    amount_am = np.where(amount_am > amorphous_material_initial, amorphous_material_initial, amount_am)
    if verbose and counter % divider == 0:
        print('time', t)
        # print('Total water:'.ljust(tabs), total_water_avg)
        # print('amount_am:'.ljust(tabs), amount_am)

    # m_void_conc = initials[0:n_space_steps]
    # m_particle_tot_conc = initials[n_space_steps:n_space_steps*2]
    # amount_am = initials[n_space_steps*2:n_space_steps*3]
    # m_void = m_void_conc / (porosity_powder * density_gas)  # kg/kg water in void, max 0.0114
    # total_water_avg = m_void_conc + m_particle_tot_conc  # kg/m3 water in system

    ###################################### FRACTION ITERATION ##########################################################
    # m_void = np.where(m_void < 0, 0, m_void)                # Cannot be negative, or more than 1
    # m_void = np.where(m_void > 1, 1, m_void)
    #
    # max_diff = 0.00011                                        # Update m_void as powder sorbs. Total water remains.
    # gas_fraction_max = compute_gas_fraction_max(amount_am, m_particle_am_conc_sat, m_particle_cryst_conc_sat, m_gas_conc_sat)
    # while_time_start = time.time()
    # m_void = compute_H_and_M_gas_fractions(total_water_avg, max_diff, m_void, p_saturated, amount_am, m_gas_sur,
    #                                        n_space_steps, gas_fraction_min, gas_fraction_max, gas_fraction_initial)
    # while_time_this = time.time() - while_time_start
    # while_time += while_time_this
    #################################### END FRACTION ITERATION ########################################################

    while_time_start = time.time()
    ############################################# TABLES ###############################################################
    m_void = compute_H_and_M_gas_tables(amount_am, total_water)
    while_time_this = time.time() - while_time_start
    while_time += while_time_this
    ########################################## END TABLES ##############################################################

    if verbose and counter % divider == 0:
        print('m_void by table:'.ljust(tabs), m_void)
    # Update system for m_void
    water_activity = compute_water_activity_from_m_void(m_void, p_saturated)
    m_particle_cryst = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_particle_am = compute_GAB_equilibrium_moisture_am(water_activity)

    m_particle_tot = m_particle_cryst * (1-amount_am) + m_particle_am * amount_am
    m_particle_tot_conc_new = (1-porosity_powder) * density_particle * m_particle_tot
    m_void_conc_new = m_void * porosity_powder * density_gas

    # m_particle_tot_conc_change  = m_particle_tot_conc_new - m_particle_tot_conc
    # m_void_change_sorption     = m_void_conc_new - m_void_conc
    if verbose and counter % divider == 0:
        print('Total water before:'.ljust(tabs), total_water[0])
        print('Total water after:'.ljust(tabs), m_particle_tot_conc_new[0] + m_void_conc_new[0])
        # print('m_void_conc after[0]', m_void_conc_new[0])
        # print('m_gas_conc_sur', m_gas_conc_sat)
        # print('m_void[0]'.ljust(tabs), m_void[0])
        # print('m_gas_sur', m_gas_sur)
        # print('water_act[0]'.ljust(tabs), water_activity[0])
        # print('m_particle_am[0]'.ljust(tabs), m_particle_am[0])
        # print('m_particle_sat', m_particle_am_sat)
        # print('m_particle_tot_conc_change', m_particle_tot_conc_change[0])
        # print('m_void_change_sorption', m_void_change_sorption[0])
        print('')

    # OLD WAY NO FRACTION
    # m_void              = m_void_conc/(density_gas * porosity_powder)                               # fraction
    # water_activity      = compute_water_activity_from_m_void(m_void, p_saturated)
    #
    # m_particle_cryst    = compute_GAB_equilibrium_moisture_cryst(water_activity)
    # m_particle_am       = compute_GAB_equilibrium_moisture_am(water_activity)
    #     # END OLD WAY
    glass_temp          = compute_glass_temp_mix(1 - m_particle_am, glass_temp_lactose, glass_temp_water_1)
    temp_diff           = temp_initial - glass_temp
    change_amorph       = compute_kinetics_avrami(amount_am, temp_diff)

    laplacian_conc = compute_diffusion_laplacian(m_void, m_gas_sur, step_length)
    diffusion = laplacian_conc * moisture_diffusivity * porosity_powder * 0.5                       # kg/m5 * m2/s = kg/(m3 s)

    ##################################### DERIVATIVE ###################################################################
    # derivative_cryst = compute_derivative_m_GAB(water_activity)
    # derivative_am = compute_derivative_m_GAB(water_activity, cryst=False)
    # derivative_total = (amount_am - change_amorph) * derivative_am + \
    #                    (1 - amount_am) * derivative_cryst #+ \
    #                    # change_amorph * (m_particle_am - m_particle_cryst_sat) * 100
    #
    # derivative_total *= density_particle * (1-porosity_powder)                                      # kg/m3
    # m_void_change_cr = change_amorph * (m_particle_am - m_particle_cryst_sat)                       # fraction
    # m_void_change_cr *= density_particle * (1 - porosity_powder)                                    # kg/m3
    # m_void_change_cr /= 3000
    # top = diffusion #- m_void_change_cr #* (t - t_prev)
    # bottom = porosity_powder + (1-porosity_powder) * derivative_total
    # m_void_conc_change = diffusion + m_void_change_sorption
    # system_change = np.concatenate([m_void_conc_change, m_particle_tot_conc_change, change_amorph]).astype('float64')
    ##################################### DERIVATIVE END ###############################################################
    system_change = np.concatenate([diffusion, change_amorph]).astype('float64')
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
    axs[0, 0].hlines(m_gas_sat, -1, hours+1, label='Sat sur', linestyles='--', color=sat_color)
    # axs[0, 0].plot(time, m_void_conc_avg, label='Tot', color=m_color)
    axs[0, 0].plot(time, m_void_total[:], label='Tot', color=m_color)
    # axs[0, 0].plot(time, m_void_conc_all[:, 0], label='0', color=m_color_light, linestyle='--')
    axs[0, 0].plot(time, m_void_total_all[:, 0], label='0', color=m_color_light, linestyle='--')
    # axs[0, 0].plot(time, m_void_conc_all[:, int(n_space_steps/3)], label=f'{int(n_space_steps/3)}', color='cyan', linestyle='--')
    axs[0, 0].plot(time, m_void_total_all[:, int(n_space_steps/3)], label=f'{int(n_space_steps/3)}', color='cyan', linestyle='--')
    # axs[0, 0].plot(time, m_void_conc_all[:, int(2 * n_space_steps/3)], label=f'{int(2 * n_space_steps/3)}', color='darkcyan', linestyle='--')
    # axs[0, 0].plot(time, m_void_conc_all[:, n_space_steps-1], label=f'{n_space_steps-1}', color=m_color, linestyle='--')
    axs[0, 0].legend()
    axs[0, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 0].set_xlim(x_min, x_max)
    axs[0, 0].set_title('Void moisture concentration')
    axs[0, 0].set(xlabel='Time, hours', ylabel='Moisture conc void, kg/m3')

    axs[1, 0].hlines(m_particle_am_sat, -1, hours+1, label='Sat am', color=sat_color, linestyles='--')
    axs[1, 0].plot(time, m_particle_am_vector, label='Am', color=m_color)
    axs[1, 0].plot(time, m_particle_am_all_vector[:, 0], label='0', color=m_color_light, linestyle='--')
    axs[1, 0].plot(time, m_particle_am_all_vector[:, int(n_space_steps/3)], label=f'{int(n_space_steps/3)}', color='cyan', linestyle='--')
    # axs[1, 0].plot(time, m_particle_am_all[:, int(2 * n_space_steps/3)], label=f'{int(2 * n_space_steps/3)}', color='darkcyan', linestyle='--')
    # axs[1, 0].plot(time, m_particle_am_all[:, -1], label=f'{n_space_steps-1}', color=m_color, linestyle='--')
    axs[1, 0].legend()
    axs[1, 0].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 0].set_xlim(x_min, x_max)
    axs[1, 0].set_title('Moisture content amorphous')
    axs[1, 0].set(xlabel='Time, hours', ylabel='Moisture, kg/kg lactose')

    # axs[2, 0].hlines(m_particle_am_conc_sat, -1, hours+1, color=sat_color, linestyles='--') #label='Sat am, kg/m3',
    axs[2, 0].hlines(m_particle_cryst_conc_sat, -1, hours+1, color=sat_color, linestyles='--') #label='Sat cryst, kg/m3',
    axs[2, 0].plot(time, total_water, label='Tot', color=m_color)
    axs[2, 0].plot(time, total_water_all[:, 0], label=f'{0}', color=m_color_light, linestyle='--')
    axs[2, 0].plot(time, total_water_all[:, int(n_space_steps/3)], label=f'{int(n_space_steps/3)}', color='cyan', linestyle='--')
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
    axs[0, 1].plot(time, t_tg_diff, label='T-Tg', color=temp_color)
    axs[0, 1].plot(time, t_tg_diff_all[:, 0], label='0', color=temp_color, linestyle='--')
    axs[0, 1].legend()
    axs[0, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[0, 1].set_xlim(x_min, x_max)
    axs[0, 1].set_title('Glass transition temp')
    axs[0, 1].set(xlabel='Time, hours', ylabel='Opeating temp - glass temp')

    axs[1, 1].plot(time, amorphicity_norm, label='Amorphicity', color=am_color)
    # amorphicity_all[:, 29] = normalize_data(amorphicity_all[:, 29], zero_one=False)
    axs[1, 1].plot(time, amorphicity_all[:, 0], label='0', color='green', linestyle='--')
    axs[1, 1].plot(time, amorphicity_all[:, int(n_space_steps/3)], label=f'{int(n_space_steps/3)}', color='lightgreen', linestyle='--')
    axs[1, 1].plot(time, amorphicity_all[:, 2 * int(n_space_steps/3)], label='20', color=f'{2*int(n_space_steps/3)}', linestyle='--')
    axs[1, 1].plot(time, amorphicity_all[:, n_space_steps-1], label=f'{n_space_steps}', color='darkgreen', linestyle='--')
    axs[1, 1].legend()
    axs[1, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[1, 1].set_xlim(x_min, x_max)
    axs[1, 1].set_ylim(-0.1, 1.1)
    axs[1, 1].set_title('Amount amorhpous material')
    axs[1, 1].set(xlabel='Time, hours', ylabel='Am material, kg/kg lactose')

    # axs[2, 1].hlines(m_particle_cryst_sat, -1, hours+1, label='Saturated cryst', color=sat_color, linestyles='--')
    # axs[2, 1].plot(time, m_particle_cryst_avg, label='Crystalline', color=m_color)
    axs[2, 1].plot(time, m_powder_diffs, label='Diff moisture', color=m_color)
    axs[2, 1].grid(color=grid_color, linestyle='-', linewidth=0.2)
    axs[2, 1].set_xlim(x_min, x_max)
    axs[2, 1].set_title('Total moisture content powder')
    axs[2, 1].set_xlim(x_min, x_max)
    axs[2, 1].legend()
    axs[2, 1].set(xlabel='Time, hours', ylabel='Moisture content, kg/kg')

    ########################################### HEAT FLOW ##############################################################
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

    max_cryst = np.max(heat_of_cryst_vector)
    axs[2, 2].plot(time, total_energy_vector, label='Total energy', color=temp_color)
    axs[2, 2].legend()
    axs[2, 2].grid(color=grid_color, linestyle='-', linewidth=0.2)
    # axs[2, 2].set_ylim(0, max_cryst)
    axs[2, 2].set_xlim(x_min, x_max)
    axs[2, 2].set_title('Total energy')
    axs[2, 2].set(xlabel='Time, hours', ylabel='Heat, J/(g s)')

    # plt.savefig('system_w_W_fraction.pdf')
    plt.savefig('system_w_derivative.pdf')
    plt.show()

################################ PARAMETERS & INITIAL CONDITIONS #######################################################
n_space_steps = 4
density_powder = 300
while_time = 0
# porosity_powder = 1 - density_powder/density_particle
# porosity_powder = density_powder/density_particle

weight = 0.00015                            # kg, 150 mg
weight = 0.00015                            # kg, 150 mg
volume = weight/density_powder              # m3
diameter_vial = 0.005                       # m, 5 mm
diameter_vial = 0.011                       # m, measured as 11 mm
area_vial = np.pi * (diameter_vial/2)**2    # m2
total_height = volume/area_vial
print('Total height of bed:', total_height)
step_length = total_height/n_space_steps

# Computing times
hours = 2
seconds = hours * 60 * 60
resolution = 1000
time_step = seconds/resolution                      # each step in time is this many seconds, s
discrete_time = np.linspace(0, seconds, resolution)

m_gas_conc_initial_vector = np.zeros(n_space_steps) + m_void_conc_initial
initial_amorphicity_vector = np.zeros(n_space_steps) + amorphous_material_initial
m_particle_tot_conc_initial_vector = np.zeros(n_space_steps) + m_particle_tot_conc_initial
# initials = np.concatenate([m_gas_conc_initial_vector, m_particle_tot_conc_initial_vector, initial_amorphicity_vector])      # kg/m3, kg/kg
initials = np.concatenate([m_gas_conc_initial_vector + m_particle_tot_conc_initial_vector, initial_amorphicity_vector])      # kg/m3, kg/kg
# print(initials)
tabs = 50
def print_info():
    print('\n')
    print('############################################# CONDITIONS #############################################')
    print('Water activity surroundings:'.ljust(tabs), f'{water_activity_sur:.2f}')
    print('Moisture in the surrounding as fraction:'.ljust(tabs), f' {m_gas_sur:.5f}')
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
computed_system = odeint(compute_system, initials, discrete_time, args=(n_space_steps, step_length, True), mxstep=2000)
computation_tot_time = time.time() - computation_start_time
while_time_fraction = while_time/computation_tot_time

print('############################################# COMPUTATION COMPLETE #############################################')
# m_void_conc_all       = computed_system[:, 0:n_space_steps]
# m_particle_tot_conc_all     = computed_system[:, n_space_steps:n_space_steps*2]
# amorphicity_all             = computed_system[:, n_space_steps*2:n_space_steps*3]
# amorphicity_all             = normalize_data(amorphicity_all, zero_one=False)

total_water_all             = computed_system[:, 0:n_space_steps]
amorphicity_all             = computed_system[:, n_space_steps:n_space_steps*2]


m_void_total_all = np.zeros([resolution, n_space_steps])
for t in range(resolution):
    if t % 500 == 0:
        print('')
        print('Am:', amorphicity_all[t, :])
        print('W:', total_water_all[t, :])
    # if t % 100 == 0:
    #     print('Calculating step ', t, 'of', resolution)
    m_void_total_all[t, :] = compute_H_and_M_gas_tables(amorphicity_all[t, :], total_water_all[t, :])

m_void_conc_total_all = m_void_total_all * porosity_powder * density_gas
amorphicity_all             = normalize_data(amorphicity_all, zero_one=False)

print('Result:')
print('Total time:'.ljust(tabs), computation_tot_time)
print('Time fraction used in while loop:'.ljust(tabs), while_time_fraction)
# print('m_void_conc_all', computed_system[0, 0:n_space_steps])
# print('m_void_conc_all', computed_system[1, 0:n_space_steps], '\n')
#
# print('m_particle_tot_conc_all', computed_system[0, n_space_steps:n_space_steps*2])
# print('m_particle_tot_conc_all', computed_system[1, n_space_steps:n_space_steps*2])
# print('m_particle_tot_conc_all', computed_system[2, n_space_steps:n_space_steps*2])
# print('m_particle_tot_conc_all', computed_system[3, n_space_steps:n_space_steps*2], '\n')
#
# print('amorphicity_all', computed_system[0, n_space_steps*2:n_space_steps*3])
# print('amorphicity_all', computed_system[1, n_space_steps*2:n_space_steps*3], '\n')

m_void_conc_total = 0
amorphicity_total = 0
total_water = 0
# m_particle_tot_conc_total = 0
for n in range(n_space_steps):
    if n == 0:
        m_void_conc_total += m_void_conc_total_all[:, n]
        amorphicity_total += amorphicity_all[:, n]
        # m_particle_tot_conc_total += m_particle_tot_conc_all[:, n]
        total_water += total_water_all[:, n]

    else:
        m_void_conc_total += m_void_conc_total_all[:, n]
        amorphicity_total += amorphicity_all[:, n]
        # m_particle_tot_conc_total += m_particle_tot_conc_all[:, n]
        total_water += total_water_all[:, n]


m_void_conc_total /= (n_space_steps)
m_void_total = m_void_conc_total/(porosity_powder * density_gas)
amorphicity_total /= n_space_steps
total_water /= n_space_steps
# m_particle_tot_conc_total /= n_space_steps

m_void_tot_vector               = m_void_conc_total/(density_gas * porosity_powder)
m_void_all_vector               = m_void_conc_total_all/(density_gas * porosity_powder)

water_activity_void_vector      = compute_water_activity_from_m_void(m_void_tot_vector, p_saturated)
water_activity_all_vector       = compute_water_activity_from_m_void(m_void_all_vector, p_saturated)

m_particle_cryst_vector         = compute_GAB_equilibrium_moisture_cryst(water_activity_void_vector)
m_particle_cryst_all_vector     = compute_GAB_equilibrium_moisture_cryst(water_activity_all_vector)

m_particle_am_vector            = compute_GAB_equilibrium_moisture_am(water_activity_void_vector)
m_particle_am_all_vector        = compute_GAB_equilibrium_moisture_am(water_activity_all_vector)

m_powder_total                  = m_particle_am_vector * amorphicity_total + m_particle_cryst_vector * (1-amorphicity_total)
m_powder_conc_total             = m_powder_total * density_particle * (1 - porosity_powder)
m_powder_total_all_vector       = m_particle_am_all_vector * amorphicity_all + m_particle_cryst_all_vector * (1-amorphicity_all)
m_powder_conc_all_vector        = m_powder_total_all_vector * density_particle * (1 - porosity_powder)

tg_s            = compute_glass_temp_mix(1 - m_particle_am_vector, glass_temp_lactose, glass_temp_water_1)
tg_s_all        = compute_glass_temp_mix(1 - m_particle_am_all_vector, glass_temp_lactose, glass_temp_water_1)
t_tg_diff       = temp_initial - tg_s
t_tg_diff_all   = temp_initial - tg_s_all

m_powder_diffs = m_powder_total[1:] - m_powder_total[:-1]                   # fraction, kg water/kg dry
m_powder_diffs = np.insert(m_powder_diffs, 0, m_powder_diffs[0], axis=0)
# m_powder_diffs_avg[m_powder_diffs_avg < 0] = 0
# print(m_powder_diffs_avg[0:10])

am_am_diffs = amorphicity_total[1:] - amorphicity_total[:-1]
am_am_diffs = np.insert(am_am_diffs, 0, 0, axis=0)
am_cryst_diffs = - am_am_diffs

heat_of_sorption = compute_heat_of_sorption(water_activity_void_vector, temp_initial)   # J/kg
# heat_of_sorption = 2.5
heat_of_sorption_vector = m_powder_diffs * (heat_of_sorption) / time_step               # J/(kg s)
heat_of_sorption_vector /= 1000                                                         # J/(kg s) to J/(g s)

heat_of_cryst_vector = am_cryst_diffs * (heat_of_crystallization/1000) / time_step      # J/kg to J/(g s)

total_energy_vector = heat_of_cryst_vector + heat_of_sorption_vector

amorphicity_norm = amorphicity_total
# amorphicity_norm = normalize_data(amorphicity, zero_one=False)

########################################## PLOT ########################################################################
plot_tam()

