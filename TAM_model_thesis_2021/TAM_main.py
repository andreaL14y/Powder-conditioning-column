#################################### MAIN FILE TO RUN SIMULATIONS ######################################################
from TAM_functions import *

# main function; computes the full system via the differential equations
def compute_system(initials, t, n_space_steps, step_length, verbose=False):
    global counter
    global while_time
    global water_activity_just_used

    divider = 25000                     # interval at which to print status
    counter += 1                        # iteration no

    total_water         = initials[0:n_space_steps]
    amount_am           = initials[n_space_steps:n_space_steps * 2]
    total_energy        = initials[n_space_steps * 2:n_space_steps * 3]
    water_activity_prev = initials[n_space_steps * 3:n_space_steps * 4]
    m_sur               = initials[n_space_steps * 4:n_space_steps * 5]
    temp_celsius_prev   = initials[n_space_steps * 5:n_space_steps * 6]

    if hindered_am == False:
        m_particle_am_prev  = initials[n_space_steps * 6:n_space_steps * 7]
    else:
        m_particle_am       = initials[n_space_steps * 6:n_space_steps * 7]

    water_activity_sur_vector = compute_water_activity_from_m_void(m_sur, temp_initial_celsius)
    amount_am                 = np.where(amount_am > amorphous_material_initial, amorphous_material_initial, amount_am)

    ########################################### MOISTURE ###############################################################
    while_time_start = time.time()                      # compute time fraction taken up by optimization

    # iterate to find optimum water activity. ftol can be decreased to get more precise results, but will be slower
    optimized_aw = scipy.optimize.minimize(compute_H_and_M_iteratively, water_activity_just_used, method='SLSQP',
                                           args=(total_water, amount_am),
                                           bounds=((0, 1),) * n_space_steps, options={'ftol': 1e-07})
    water_activity          = optimized_aw.x

    water_activity_change   = water_activity - water_activity_prev
    m_void                  = compute_H_from_water_activity_temp(water_activity)
    while_time += time.time() - while_time_start        # compute time fraction taken up by optimization

    m_particle_cryst            = compute_GAB_equilibrium_moisture_cryst(water_activity)
    m_particle_am_max           = compute_GAB_equilibrium_moisture_am(water_activity)

    # if no hindered diffusion in amorpohus -> everything in equilibrium
    if hindered_am == False:
        m_particle_am = m_particle_am_max
        m_particle_am_change = m_particle_am - m_particle_am_prev

    # otherwise compute how moisture in amorphous particles change
    else:
        m_particle_am_laplacian = np.zeros(n_space_steps)
        for n in range(n_space_steps):
            m_particle_am_laplacian[n]     = compute_laplacian(m_particle_am[n], m_particle_am_max[n], am_layer_thickness, double_bc=False)
        m_particle_am_change_diff   = m_particle_am_laplacian * am_diffusivity

    m_particle_tot      = m_particle_cryst * (1-amount_am) + m_particle_am * amount_am
    m_particle_cryst_prev = compute_GAB_equilibrium_moisture_cryst(water_activity_prev)

    if hindered_am == False:
        m_particle_tot_prev = m_particle_cryst_prev * (1 - amount_am) + m_particle_am_prev * amount_am
    else:
        m_particle_tot_prev = m_particle_cryst_prev * (1 - amount_am) + m_particle_am * amount_am

    ########################################### ENERGY #################################################################
    gas_density     = compute_density_air(m_void)
    water_W         = porosity_powder * m_void * gas_density + (1 - porosity_powder) * m_particle_tot * density_particle
    excess_water    = total_water - water_W
    excess_water    = np.where((excess_water > 0) & (water_activity == 1), excess_water, 0)

    temp_kelvin, temp_celsius   = compute_temp_from_energy(m_void, m_particle_tot, total_energy, excess_water, verbose=False)
    temp_change                 = temp_celsius - temp_celsius_prev
    energy_vapor                =  m_void * compute_enthalpy_vapor(temp_celsius) * gas_density

    ######################################### AIR IN AMPOULE ###########################################################
    # compute how the humidity of the air above the powder bed changes
    gas_density_sur = compute_density_air(m_sur)
    m_sur_bounds    = np.array([m_gas_sur, m_void[0]])
    m_sur_change    = compute_laplacian(m_sur, m_sur_bounds, step_length_air, double_bc=True, inflow=True)

    ################################## RESULTING CHANGES IN Q AND W ####################################################
    # compute moisture diffusion into powder bed
    laplacian_conc  = compute_laplacian(m_void * gas_density, m_sur[-1] * gas_density_sur[-1], step_length)
    m_diffusion     = laplacian_conc * diffusivity_eff                  # kg/m5 * m2/s = kg/(m3 s)
    change_moisture = m_diffusion                                       # total water in bed only changes via diffusion

    # compute energy diffusion into powder bed
    energy_vapor_sur = m_sur[-1] * compute_enthalpy_vapor(temp_initial_celsius) * gas_density_sur[-1]
    energy_diff_laplacian = compute_laplacian(energy_vapor, energy_vapor_sur, step_length)
    change_energy_diffusion = diffusivity_eff * energy_diff_laplacian   # m2/s * (J/kg * kg/m3)/m2 = J/(m3 s)

    # compute how heat conducts through powder bed
    temp_bounds = np.array([temp_initial_celsius, temp_initial_celsius])
    temp_laplacian = compute_laplacian(temp_celsius, temp_bounds, step_length, double_bc=True)
    change_energy_conduction = conductivity_tot * temp_laplacian

    ###################################### AMORPHICITY #################################################################
    glass_temp      = compute_glass_temp_mix(1 - m_particle_am, glass_temp_lactose, glass_temp_water)
    temp_diff       = temp_kelvin - glass_temp
    change_amorph   = compute_kinetics_avrami(amount_am, temp_diff)                                   # negative

    change_energy_crystallization = change_amorph * density_particle * (1 - porosity_powder) * heat_of_crystallization  # Negative, 1/s * kg/m3 * J/kg
    change_energy_crystallization = - change_energy_crystallization         # Positive

    change_energy_evap      = change_moisture * heat_of_evaporation_water           # positive during sorption kg/m3 * J/kg = J/m3

    m_diffs_p               = m_particle_tot - m_particle_tot_prev                  # positive during sorption, then neg
    m_diffs_p_conc          = m_diffs_p * (1 - porosity_powder) * density_particle
    change_energy_binding   = m_diffs_p_conc * heat_binding_diffs / density_water

    change_energy = change_energy_conduction + change_energy_crystallization + change_energy_diffusion - change_energy_evap + change_energy_binding

    if hindered_am == True:
        # moisture content of amorphous particles changes by diffusion AND by the water suddenly available from
        # crystallization. The last one is already in the particles; does not need to diffuse in
        m_particle_am_change = m_particle_am_change_diff - change_amorph * (m_particle_am - m_particle_cryst)

    # print status
    if verbose and counter % divider == 0:
        hours_print = int(t / 3600)
        minutes_print = int((t - hours_print * 3600) / 60)
        seconds_print = int(t - (minutes_print * 60 + hours_print * 3600))
        time_percentage = t/seconds
        print('')
        print(f'        Exp {exp} {time_percentage*100:.1f} % completed: TIME', hours_print, 'h;', minutes_print, 'min;', seconds_print, 's')
        np.set_printoptions(formatter={'float': '{: .2f}'.format})
        print('T-Tg:'.ljust(tabs), temp_diff)
        print('Temp C:'.ljust(tabs), temp_celsius)
        # print('Temp C sur:'.ljust(tabs), temp_celsius_sur)
        np.set_printoptions(formatter={'float': '{: .1e}'.format})
        # print('Total energy:'.ljust(tabs), total_energy)
        print('')

        # print('m_void:'.ljust(tabs), m_void)
        print('M am:'.ljust(tabs), m_particle_am)
        print('M am sat:'.ljust(tabs), m_particle_am_max)
        # print('Excess water:'.ljust(tabs), excess_water)

        np.set_printoptions(formatter={'float': '{: .3f}'.format})
        print('Water act:'.ljust(tabs), water_activity)
        print('Water act sur:'.ljust(tabs), water_activity_sur_vector)
        np.set_printoptions(formatter={'float': '{: .5f}'.format})
        print('Am am:'.ljust(tabs), amount_am)

        # print('Change m_am diffusion:'.ljust(tabs), m_particle_am_change_diff)
        print('Change m_am crystallization:'.ljust(tabs), - change_amorph * (m_particle_am - m_particle_cryst))

        # np.set_printoptions(formatter={'float': '{: .1e}'.format})
        # print('Change am:'.ljust(tabs), change_amorph)
        # print('')
        # #
        # print('Change energy diff:'.ljust(tabs), change_energy_diffusion)
        # print('Change energy cond:'.ljust(tabs), change_energy_conduction)
        # print('Change energy cryst:'.ljust(tabs), change_energy_crystallization)
        # print('Change energy tot:'.ljust(tabs), change_energy)

        print('N iterations in optimize:'.ljust(tabs), optimized_aw.nfev)
        print(f'        {time_percentage*100:.1f} % completed: TIME', hours_print, 'h;', minutes_print, 'min;', seconds_print, 's')
        print('')

    # save this just to help optimizer reach solution faster
    water_activity_just_used = water_activity

    system_change       = np.concatenate([change_moisture, change_amorph, change_energy, water_activity_change,
                                          m_sur_change, temp_change, m_particle_am_change]).astype('float64')
    return system_change

# loop through the experiments
for exp in exp_list:
    print(f'############################################ SETUP EXP {exp} #############################################')
    counter = 0
    batch, amorphous_material_initial, relative_humidity_gas_inlet, particle_diameter, weight, hours, resolution = \
        create_conditions(exp)

    volume              = weight/density_powder                                 # volume of powder bed, m3
    volume_ampoule      = in_area_TAM_cyl * height_air_TAM_cyl                  # total volume ampoule
    total_height_powder = volume / in_area_TAM_cyl                              # how high is powder filled
    height_air_TAM_cyl  = height_air_TAM_cyl - total_height_powder

    seconds = hours * 60 * 60
    time_step = seconds / resolution                        # each step in time is this many seconds, s
    discrete_time = np.linspace(0, seconds, resolution)     # create output times
    while_time = 0

    step_length = total_height_powder / n_space_steps       # length of each step in powder bed
    step_length_air = height_air_TAM_cyl / n_space_steps    # length of each step in air above powder bed

    if hindered_am == True:
        save_string = f'data_and_plots/211021-final/hindered/exp {exp}, time {hours}, n {n_space_steps}'
    else:
        save_string = f'data_and_plots/211021-final/unhindered/exp {exp}, time {hours}, n {n_space_steps}, dens 170'

    m_particle_tot_conc_initial, m_void_conc_initial, enthalpy_powder_initial, water_activity_void_initial, \
    m_void_initial, temp_initial_celsius, m_gas_sur, m_particle_cryst_initial, m_particle_am_initial, \
    m_particle_tot_initial, gas_density_initial, water_activity_sur, m_particle_cryst_conc_initial = \
        create_initial_conditions(relative_humidity_gas_inlet, amorphous_material_initial)

    m_tot_conc_vector           = np.zeros(n_space_steps) + m_particle_tot_conc_initial + m_void_conc_initial
    initial_amorphicity_vector  = np.zeros(n_space_steps) + amorphous_material_initial
    energy_initial_vector       = np.zeros(n_space_steps) + enthalpy_powder_initial
    water_activity_initial_vector = np.zeros(n_space_steps) + water_activity_void_initial
    m_sur_initial_vector        = np.zeros(n_space_steps) + m_void_initial
    temp_initial_vector         = np.zeros(n_space_steps) + temp_initial_celsius
    m_am_initial_vector         = np.zeros(n_space_steps) + m_particle_am_initial

    initials = np.concatenate([m_tot_conc_vector, initial_amorphicity_vector, energy_initial_vector,
                               water_activity_initial_vector, m_sur_initial_vector, temp_initial_vector, m_am_initial_vector])

    water_activity_just_used = np.zeros(n_space_steps) + water_activity_void_initial

    print(f'############################################ START EXP {exp} #############################################')
    print('Max time:'.ljust(tabs), hours, 'hours, = ', seconds, 'seconds.')
    print('Time step:'.ljust(tabs), f'{time_step:5f} s')
    computation_start_time = time.time()

    computed_system = odeint(compute_system, initials, discrete_time, args=(n_space_steps, step_length, True),
                             mxstep=25000)

    # check computation time and print it, along with the fraction used in optimization
    computation_tot_time = time.time() - computation_start_time
    while_time_fraction = while_time / computation_tot_time
    hours_print = int(computation_tot_time/3600)
    minutes_print = int( (computation_tot_time - hours_print*3600) / 60)
    seconds_print = computation_tot_time - (minutes_print * 60 + hours_print * 3600)

    print('Result:')
    print('Total time:'.ljust(tabs), hours_print, 'h', minutes_print, 'min and', int(seconds_print), 's')
    print('Time fraction used in while loop:'.ljust(tabs), f'{while_time_fraction:.2f}')
    np.save(save_string, computed_system)
    print('###################################### COMPUTATION EXP {exp} COMPLETE #####################################')
