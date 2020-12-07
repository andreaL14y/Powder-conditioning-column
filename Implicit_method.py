import numpy as np
import math
from define_functions import*
from test_main import*

#from control_volume_computation_andrea import *
time_steps = 100
number_measure_points = 100
space_step = bed_length/number_measure_points
time_step = space_step / superficial_velocity

############### INITIAL CONDITIONS ############################################################
moisture_particle = np.zeros(number_measure_points) + moisture_particle_initial
temperature_particle = np.zeros(number_measure_points) + temp_initial
moisture_gas = np.zeros(number_measure_points) + moisture_gas_initial_in
moisture_gas_2 = np.zeros(number_measure_points) + moisture_gas_initial_in                   
temperature_gas = np.zeros(number_measure_points) + temp_initial

############### INITIAL PARAMETERS THAT CHANGE DIRECTLY WITH THE TEMPERATURE ###################
#moisture_density = np.zeros(number_measure_points)+moisture_density                                    # TODO: for simplifications this will be considered as constant
pressure_saturated = np.zeros(number_measure_points)+pressure_saturated_initial

############### INITIAL PARAMETERS THAT CHANGE INDIRECTLY WITH THE TEMPERATURE #################
molar_concentration_moisture = np.zeros(number_measure_points)+molar_concentration_moisture_initial #dep. on moisture_density
relative_humidity = np.zeros(number_measure_points)+relative_humidity_gas_initial #dep. on pressure saturated
partial_pressure_moisture = np.zeros(number_measure_points)+partial_pressure_moisture_initial #dep. on pressure saturated and RH_gas

################ INITIAL PARAMETERS THAT CHANGE WITH MOISTURE & TEMPERATURE (PARTICLE/GAS)#####

constant = np.zeros(number_measure_points)+constant_initial #dep. on pressure saturated, k_GP
k_GP = np.zeros(number_measure_points)+k_GP_initial #dep. on molar_concentration_moisture

f_X = np.zeros(number_measure_points)


################ SYSTEM OF LIN EQ FOR du/dt=a*Laplacian*u-b*u+r ################################
################ Moisture Particle                              ################################
M_mp = np.zeros((number_measure_points, number_measure_points))
rhs_mp = np.zeros(number_measure_points)
b_mp = np.zeros(number_measure_points)
r_mp = np.zeros(number_measure_points)

################ Moisture Gas                                   ################################
M_mg = np.zeros((number_measure_points, number_measure_points))
rhs_mg = np.zeros(number_measure_points)
a_mg = moisture_diffusivity
c_mg = gas_velocity #TODO: is this the right velocity?
r_mg_old = np.zeros(number_measure_points) #TODO really necessary?
r_mg_new = np.zeros(number_measure_points)

################ Temperature Particle                           ################################
K_tp = np.zeros((number_measure_points, number_measure_points))
rhs_tp = np.zeros(number_measure_points)
b_tp = np.zeros(number_measure_points)
r_tp = np.zeros(number_measure_points)

################ Temperature Gas                                ################################
K_tg = np.zeros((number_measure_points, number_measure_points))
rhs_tg = np.zeros(number_measure_points)
b_tg = np.zeros(number_measure_points)
r_tg = np.zeros(number_measure_points)

for t in range(time_steps):
    at_most_bed_length = min(number_measure_points, t+1)
    ##### Saving current values
    moisture_particle_current = moisture_particle
    moisture_gas_current = moisture_gas
    temperature_particle_current = temperature_particle
    temperature_gas_current = temperature_gas


    for i in range(at_most_bed_length):
        ##### Moisture particle 
        f_X [i] = compute_equilibrium_moisture(alpha_parameter, moisture_particle_current[i], N)

        b_mp[i] = -1*k_GP[i]*specific_surface_area*pressure_saturated[i]/pressure_ambient
        r_mp[i] = k_GP[i]*specific_surface_area*pressure_saturated[i]*relative_humidity[i]/pressure_ambient
        if i==0:
            M_mp[i,i]=1 #TODO: initial cond
            rhs_mp[i] = time_step*(b_mp[i]/2*moisture_particle_current[i])+time_step*(r_mp[i+1]+r_mp[i])/2 #TODO: BC? WHAT IS r_MP[i], r_mp[i+1], r_mp[i-1]?
        elif i==(number_measure_points-1):
            M_mp[i,i]=1 #TODO BC
            rhs_mp[i] = time_step*b_mp[i]/2*moisture_particle_current[i]+time_step*(r_mp[i]) #TODO
        else: 
            M_mp[i,i]=1-time_step*b_mp[i]/2
            rhs_mp[i] = time_step*b_mp[i]/2*moisture_particle_current[i]+time_step*(r_mp[i]+r_mp[i+1])/2

        moisture_particle[:i]=np.linalg.solve(M_mp[:i,:i], rhs_mp[:i])

        ##### Moisture Gas
        r_mg_old[i] = (-k_GP[i]*specific_surface_area*pressure_saturated[i]*particle_density* \
            porosity_powder/(pressure_ambient*gas_density*(1-porosity_powder)))*(relative_humidity[i]- \
                compute_equilibrium_moisture(alpha_parameter, moisture_particle_current[i], N))

        r_mg_new[i] = (-k_GP[i]*specific_surface_area*pressure_saturated[i]*particle_density* \
            porosity_powder/(pressure_ambient*gas_density*(1-porosity_powder)))*(relative_humidity[i]- \
                compute_equilibrium_moisture(alpha_parameter, moisture_particle[i], N))

        if i==0:
            M_mg[i,i] = time_step*a_mg/space_step**2 +1 #TODO: BC
            M_mg[i,i+1]= -time_step/2*(a_mg/space_step**2-c_mg/(2*space_step)) #TODO BC
        elif i==(number_measure_points-1):
            M_mg[i,i-1] = -time_step/2*(a_mg/space_step**2+c_mg/(2*space_step)) #TODO BC
            M_mg[i,i]= time_step*a_mg/space_step**2 +1 #TODO BC
        else: 
            M_mg[i,i-1] = -time_step/2*(a_mg/space_step**2+c_mg/(2*space_step))
            M_mg[i,i]= time_step*a_mg/space_step**2 +1
            M_mg[i,i+1]= -time_step/2*(a_mg/space_step**2-c_mg/(2*space_step))
            rhs_mg[i-1] = time_step/2*(a_mg/space_step**2+c_mg/(2*space_step))*moisture_gas_current[i-1]
            rhs_mg[i] = -time_step*a_mg/space_step**2*moisture_gas_current[i]
            rhs_mg[i+1] = time_step/2*(a_mg/space_step**2-c_mg/(2*space_step))*moisture_gas_current[i+1]

        r_mg = time_step/2*(r_mg_old+r_mg_new)
        moisture_gas[:i] = np.linalg.solve(M_mg[:i,:i], rhs_mg[:i]+r_mg[:i])

        
        ############ UPDATE PARAMETERS 
        pressure_saturated[i]=compute_p_saturated(A, B, temperature_gas[i], C)
        
        relative_humidity[i]=compute_relative_humidity_from_Y(molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[i], pressure_saturated[i])
        
        molar_concentration_moisture[i]=compute_molar_concentration(relative_humidity[i], pressure_saturated[i], R_gas_constant, temperature_gas[i])
        
        partial_pressure_moisture[i]=compute_partial_pressure_moisture(molar_concentration_moisture[i], R_gas_constant, temperature_gas[i])  

        k_GP[i] = compute_mass_transfer_coefficient(
            moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, 
            flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i])[3]
        
        constant[i]=k_GP[i] * specific_surface_area * pressure_saturated[i] / pressure_ambient
    # UPDATE PARAMETERS

