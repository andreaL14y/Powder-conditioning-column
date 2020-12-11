import numpy as np
import math
from define_functions import*
from input_parameters import*
from copy import deepcopy

#from control_volume_computation_andrea import *
time_steps = 24
number_measure_points = 101
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


################ SYSTEM OF LIN EQ FOR du/dt=a*u+r               ################################
################ Moisture Particle                              ################################
M_mp = np.zeros((number_measure_points, number_measure_points))
rhs_mp = np.zeros(number_measure_points)
b_mp = np.zeros(number_measure_points)
r_mp = np.zeros(number_measure_points)

################ SYSTEM OF LIN EQ FOR du/dt=a*Laplacian*u-c*Gradient*u+r #######################
################ Moisture Gas                                   ################################
M_mg = np.zeros((number_measure_points, number_measure_points))
rhs_mg = np.zeros(number_measure_points)
a_mg = moisture_diffusivity
c_mg = gas_velocity #TODO: is this the right velocity?
r_mg_old = np.zeros(number_measure_points) #TODO really necessary?
r_mg_new = np.zeros(number_measure_points)

################ SYSTEM OF LIN EQ FOR du/dt=a*Laplacian*u-b*u+r ################################
################ Temperature Particle                           ################################
K_tp = np.zeros((number_measure_points, number_measure_points))
rhs_tp = np.zeros(number_measure_points)
a_tp = conductivity_particle/(particle_density*particle_heat_capacity)
b_tp = -heat_transfer_coefficient_initial*specific_surface_area*1/particle_heat_capacity
r_tp = np.zeros(number_measure_points)
r_tp_helper = np.zeros(number_measure_points)

################ Temperature Gas                                ################################
K_tg = np.zeros((number_measure_points, number_measure_points))
rhs_tg = np.zeros(number_measure_points)
a_tg = conductivity_gas/(gas_density*gas_heat_capacity)
b_tg = np.zeros(number_measure_points)
r_tg = np.zeros(number_measure_points)
r_tg_helper = np.zeros(number_measure_points)

for t in range(time_steps):
    print(t)
    at_most_bed_length = min(number_measure_points, t+1)
    ##### Saving current values
    moisture_particle_current = deepcopy(moisture_particle)
    moisture_gas_current = deepcopy(moisture_gas)
    temperature_particle_current = deepcopy(temperature_particle)
    temperature_gas_current = deepcopy(temperature_gas)

    for i in range(at_most_bed_length):
        ##### Moisture particle 
        b_mp[i] = -1*k_GP[i]*specific_surface_area*pressure_saturated[i]/pressure_ambient
        r_mp[i] = k_GP[i]*specific_surface_area*pressure_saturated[i]*relative_humidity[i]/pressure_ambient
        if i==0:
            M_mp[i,i]= 1 #TODO: initial cond
            rhs_mp[i] = moisture_particle_initial #time_step*(b_mp[i]/2*moisture_particle_current[i])+time_step*(r_mp[i+1]+r_mp[i])/2 #TODO: BC? WHAT IS r_MP[i], r_mp[i+1], r_mp[i-1]?
        elif i==(number_measure_points-1):
            M_mp[i,i] = 1 #TODO BC
            rhs_mp[i] = 0 #time_step*b_mp[i]/2*moisture_particle_current[i]+time_step*(r_mp[i]) #TODO
        else: 
            M_mp[i,i]= 1-time_step*b_mp[i]/2
            rhs_mp[i] = time_step*b_mp[i]/2*moisture_particle_current[i]+time_step*(r_mp[i]+r_mp[i+1])/2

        moisture_particle[:i]=np.linalg.solve(M_mp[:i,:i], rhs_mp[:i])

        # ##### Moisture Gas
        # r_mg_old[i] = (-k_GP[i]*specific_surface_area*pressure_saturated[i]*particle_density* \
        #     porosity_powder/(pressure_ambient*gas_density*(1-porosity_powder)))*(relative_humidity[i]- \
        #         compute_equilibrium_moisture(alpha_parameter, moisture_particle_current[i], N))

        # r_mg_new[i] = (-k_GP[i]*specific_surface_area*pressure_saturated[i]*particle_density* \
        #     porosity_powder/(pressure_ambient*gas_density*(1-porosity_powder)))*(relative_humidity[i]- \
        #         compute_equilibrium_moisture(alpha_parameter, moisture_particle[i], N))

        # if i==0:
        #     M_mg[i,i] = 1 #TODO: BC
        #     rhs_mg[i]=moisture_gas_initial_in
        #     #M_mg[i,i+1]= -time_step/2*(a_mg/space_step**2-c_mg/(2*space_step)) #TODO BC
        # elif i==(number_measure_points-1):
        #     # M_mg[i,i-1] = -time_step/2*(a_mg/space_step**2+c_mg/(2*space_step)) #TODO BC
        #     M_mg[i,i] = 1 #TODO BC
        #     rhs_mg = time_step/2*(r_mg_old[i]+r_mg_new[i])+moisture_gas_current[i]
        # else: 
        #     M_mg[i,i]= time_step*a_mg/space_step**2 +1
        #     rhs_mg[i] = -time_step*a_mg/space_step**2*moisture_gas_current[i]
        #     if i != 1:
        #         M_mg[i,i-1] = -time_step/2*(a_mg/space_step**2+c_mg/(2*space_step))
        #         rhs_mg[i-1] = time_step/2*(a_mg/space_step**2+c_mg/(2*space_step))*moisture_gas_current[i-1]
        #     if i != number_measure_points-2:
        #         M_mg[i,i+1]= -time_step/2*(a_mg/space_step**2-c_mg/(2*space_step))
        #         rhs_mg[i+1] = time_step/2*(a_mg/space_step**2-c_mg/(2*space_step))*moisture_gas_current[i+1]

        # r_mg = time_step/2*(r_mg_old+r_mg_new)
        # moisture_gas[:i] = np.linalg.solve(M_mg[:i,:i], rhs_mg[:i]+r_mg[:i])

        # ##### Temperature particle
        # r_tp_helper[i]=constant[i]*(relative_humidity[i]-compute_equilibrium_moisture(alpha_parameter, moisture_particle[i], N))*\
        #     heat_of_vaporization*1/particle_heat_capacity+heat_transfer_coefficient_initial/particle_heat_capacity*temperature_gas[i] #TODO is it moisture particle current here? or the new X_P? same question for T_G. IF it is the "old" values I could write all in one if loop
        
        # if i==0:
        #     K_tp[i,i]= 1 #TODO BC
        #     #K_tp[i,i+1]= - time_step*a_tp/(2*space_step**2)  #TODO BC
        #     rhs_tp[i]=temp_initial
        #     r_tp[i] = r_tp_helper[i] #TODO
        # elif i==(number_measure_points-1):
        #     #K_tp[i,i-1] = - time_step*a_tp/(2*space_step**2) #TODO BC
        #     K_tp[i,i]= 1-time_step*b_tp/2  #TODO BC
        #     rhs_tp[i]= time_step/2*((b_tp+1)*temperature_particle_current[i]+r_tp_helper[i]+r_tp_helper[i+1])
        #     r_tp[i]= time_step/2*(r_tp_helper[i-1]+r_tp_helper[i])
        # else: 
        #     K_tp[i,i]= 1-time_step-b_tp
        #     r_tp[i]= time_step/2*(r_tp_helper[i-1]+r_tp_helper[i])
        #     if i != 1: 
        #         K_tp[i,i-1] = - time_step*a_tp/(2*space_step**2)
        #         rhs_tp[i-1] = time_step*a_tp/(2*space_step**2)*temperature_particle_current[i-1]
        #     rhs_tp[i] = (time_step+b_tp-1)*temperature_particle_current[i]
        #     if i != (number_measure_points-2):
        #         rhs_tp[i+1] = time_step*a_tp/(2*space_step**2)*temperature_particle_current[i+1]
        #         K_tp[i,i+1]= - time_step*a_tp/(2*space_step**2)

        # temperature_particle[:i] = np.linalg.solve(K_tp[:i,:i], rhs_tp[:i]+r_tp[:i])
        
        # ##### Temperature gas
        # r_tg_helper[i] = (particle_density*porosity_powder*temperature_particle[i]/(gas_density*(1-porosity_powder*gas_heat_capacity)))*\
        #     (-constant[i]*(relative_humidity[i]-compute_equilibrium_moisture(alpha_parameter, moisture_particle[i], N))*\
        #         moisture_vapor_heat_capacity+heat_transfer_coefficient_initial*specific_surface_area) #TODO is it moisture particle current here? or the new X_P? same question for T_G. IF it is the "old" values I could write all in one if loop
        
        # b_tg[i] = particle_density*porosity_powder/(gas_density*(1-porosity_powder)*gas_heat_capacity)*(constant[i]*\
        #     (relative_humidity[i]-compute_equilibrium_moisture(alpha_parameter, moisture_particle[i], N))*\
        #         moisture_vapor_heat_capacity-heat_transfer_coefficient_initial*specific_surface_area)
        
        # if i==0:
        #     K_tg[i,i]= 1 #TODO BC
        #     #K_tg[i,i+1]= - time_step*a_tg/(2*space_step**2)  #TODO BC
        #     r_tg[i] = r_tg_helper[i] #TODO
        #     rhs_tg[i] = temp_initial
        # elif i==(number_measure_points-1):
        #     # K_tg[i,i-1] = - time_step*a_tg/(2*space_step**2) #TODO BC
        #     K_tg[i,i]= 1-time_step*b_tg/2  #TODO BC
        #     rhs_tg[i] = time_step/2*((b_tg+1)*temperature_gas_current[i]+r_tg_helper[i]+r_tg_helper[i+1])
        #     r_tg[i]= time_step/2*(r_tg_helper[i-1]+r_tg_helper[i])
        # else: 
        #     K_tg[i,i]= 1-time_step-b_tg[i]
        #     rhs_tg[i] = (time_step+b_tg[i]-1)*temperature_particle_current[i]
        #     r_tg[i]= time_step/2*(r_tg_helper[i-1]+r_tg_helper[i])
        #     if i != 1:
        #         K_tg[i,i-1] = - time_step*a_tg/(2*space_step**2)
        #         rhs_tg[i-1] = time_step*a_tg/(2*space_step**2)*temperature_particle_current[i-1]
        #     if i != (number_measure_points-2):
        #         rhs_tg[i+1] = time_step*a_tg/(2*space_step**2)*temperature_particle_current[i+1]
        #         K_tg[i,i+1]= - time_step*a_tg/(2*space_step**2)
            

        # temperature_gas[:i] = np.linalg.solve(K_tg[:i,:i], rhs_tg[:i]+r_tg[:i])


        ############ UPDATE PARAMETERS 
        pressure_saturated[i]=compute_p_saturated(A, B, temperature_gas[i], C)
        
        relative_humidity[i]=compute_relative_humidity_from_Y(molar_mass_dry_air, molar_mass_moisture, pressure_ambient, moisture_gas[i], pressure_saturated[i])
        
        molar_concentration_moisture[i]=compute_molar_concentration(relative_humidity[i], pressure_saturated[i], R_gas_constant, temperature_gas[i])
        
        partial_pressure_moisture[i]=compute_partial_pressure_moisture(molar_concentration_moisture[i], R_gas_constant, temperature_gas[i])  

        k_GP[i] = compute_mass_transfer_coefficient(
            moisture_diffusivity, gas_viscosity, column_diameter, porosity_powder, gas_density, particle_density, 
            flow_rate, particle_diameter, molar_mass_moisture, superficial_velocity, molar_concentration_moisture[i])[3]
        
        constant[i]=k_GP[i] * specific_surface_area * pressure_saturated[i] / pressure_ambient

print('moisture gas', moisture_gas)
print('moisture particle', moisture_particle)
print('temperature particle', temperature_particle)
print('temperature gas', temperature_gas)

