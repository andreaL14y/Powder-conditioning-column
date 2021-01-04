from PyQt5 import QtCore, QtGui, QtWidgets
from pcc import*
from input_class import*
from pcc_define_functions import*
from pcc_ode_solver import*

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    pcc = QtWidgets.QDialog()
    ui = Ui_pcc()
    ui.setupUi(pcc)
    pcc.setWindowTitle("Powder conditioning column")
    pcc.setWindowIcon(QtGui.QIcon('astra-logo'))

    ####### VALIDATOR (ONLY NUMBER INPUT) #########################################################################
    dvalidator=QtGui.QDoubleValidator()
    dvalidator.setBottom(0.001)

    ivalidator=QtGui.QIntValidator()


    ui.R_gas_constant.setValidator(dvalidator)
    ui.antoine_constants_A.setValidator(dvalidator)
    ui.antoine_constant_B.setValidator(dvalidator)
    ui.antoine_constant_C.setValidator(dvalidator)
    ui.bed_length.setValidator(dvalidator)
    ui.column_diameter.setValidator(dvalidator)
    ui.volumetric_flow_rate.setValidator(dvalidator)
    ui.temp_init.setValidator(dvalidator)
    ui.temp_walls.setValidator(dvalidator)
    ui.relative_humidity_bed_init.setValidator(dvalidator)
    ui.relative_humidity_gas_initial.setValidator(dvalidator)
    ui.relative_humidity_gas_end.setValidator(dvalidator)
    ui.pressure_ambient.setValidator(dvalidator)
    ui.porosity_powder.setValidator(dvalidator)
    ui.N_parameter.setValidator(dvalidator)
    ui.alpha_parameter.setValidator(dvalidator)
    ui.particle_density.setValidator(dvalidator)
    ui.gas_density.setValidator(dvalidator)
    ui.particle_diameter.setValidator(dvalidator)
    ui.heat_of_vaporization.setValidator(dvalidator)
    ui.gas_viscosity.setValidator(dvalidator)
    ui.moisture_diffusivity.setValidator(dvalidator)
    ui.molar_mass_moisture.setValidator(dvalidator)
    ui.molar_mass_dry_air.setValidator(dvalidator)
    ui.moisture_vapor_heat_capacity.setValidator(dvalidator)
    ui.moisture_liquid_heat_capacity.setValidator(dvalidator)
    ui.particle_heat_capacity.setValidator(dvalidator)
    ui.gas_heat_capacity.setValidator(dvalidator)
    ui.conductivity_particle.setValidator(dvalidator)
    ui.conductivity_gas.setValidator(dvalidator)
    ui.boiling_temperature.setValidator(dvalidator)

    ui.max_time.setValidator(ivalidator)
    ui.n_space_steps.setValidator(ivalidator)
    ui.n_height_steps.setValidator(ivalidator)
    ui.resolution.setValidator(ivalidator)

    in_param = input_parameter()
    dis_param = discretization_parameter()
    init_param = initial_conditions()

    def compute():
        ####### SET INPUT PARAMETERS ######################################################################
        in_param.R_gas_constant =float(ui.R_gas_constant.text())
        in_param.antoine_constant_A=float(ui.antoine_constants_A.text())
        in_param.antoine_constant_B=float(ui.antoine_constant_B.text())
        in_param.antoine_constant_C=float(ui.antoine_constant_C.text())
        in_param.bed_length=float(ui.bed_length.text())
        in_param.column_diameter=float(ui.column_diameter.text())
        in_param.volumetric_flow_rate=float(ui.volumetric_flow_rate.text())
        in_param.temp_initial=float(ui.temp_init.text())
        in_param.temp_walls=float(ui.temp_walls.text())
        in_param.relative_humidity_bed_init=float(ui.relative_humidity_bed_init.text())
        in_param.relative_humidity_gas_initial=float(ui.relative_humidity_gas_initial.text())
        in_param.relative_humidity_gas_end=float(ui.relative_humidity_gas_end.text())
        in_param.pressure_ambient=float(ui.pressure_ambient.text())
        in_param.porosity_powder=float(ui.porosity_powder.text())
        in_param.N_parameter=float(ui.N_parameter.text())
        in_param.alpha_parameter=float(ui.alpha_parameter.text())
        in_param.particle_density=float(ui.particle_density.text())
        in_param.gas_density=float(ui.gas_density.text())
        in_param.particle_diameter=float(ui.particle_diameter.text())
        in_param.heat_of_vaporization=float(ui.heat_of_vaporization.text())
        in_param.gas_viscosity=float(ui.gas_viscosity.text())
        in_param.moisture_diffusivity =float(ui.moisture_diffusivity.text())
        in_param.molar_mass_moisture=float(ui.molar_mass_moisture.text())
        in_param.molar_mass_dry_air=float(ui.molar_mass_dry_air.text())
        in_param.moisture_vapor_heat_capacity=float(ui.moisture_vapor_heat_capacity.text())
        in_param.moisture_liquid_heat_capacity=float(ui.moisture_liquid_heat_capacity.text())
        in_param.particle_heat_capacity=float(ui.particle_heat_capacity.text())
        in_param.gas_heat_capacity=float(ui.gas_heat_capacity.text())
        in_param.conductivity_particle=float(ui.conductivity_particle.text())
        in_param.conductivity_gas=float(ui.conductivity_gas.text())
        in_param.boiling_temperature=float(ui.boiling_temperature.text())

        ###### SET DISCRETIZATION PARAMETERS #####################################################
        dis_param.max_time=int(ui.max_time.text())
        dis_param.n_space_steps=int(ui.n_space_steps.text())
        dis_param.n_height_steps=int(ui.n_height_steps.text())
        dis_param.resolution=int(ui.resolution.text())
        dis_param.height_of_interest=int(dis_param.n_height_steps/2)+1

        ###### SET INITIAL CONDITIONS #############################################################
        init_param.kelvin = 273.15
        init_param.cross_sectional_area = np.pi * (in_param.column_diameter/2)**2
        init_param.gas_velocity=compute_velocity(in_param)
        init_param.specific_surface_area = compute_specific_surface_area(in_param)                      # m2/kg, SSA
        init_param.flow_rate = volumetric_flow_rate_m3_per_second(in_param)                             # m3/s
        init_param.superficial_velocity = init_param.flow_rate / (np.pi*(in_param.column_diameter/2)**2)         # superficial velocity U in m/s
        init_param.pressure_saturated_initial = compute_p_saturated_vector(in_param.temp_initial, in_param, init_param)
        init_param.partial_pressure_moisture_initial = init_param.pressure_saturated_initial * in_param.relative_humidity_gas_initial
        init_param.molar_concentration_moisture_initial = compute_molar_concentration_vector(in_param.relative_humidity_gas_initial, \
            init_param.pressure_saturated_initial, in_param.temp_initial, in_param)
        init_param.moisture_gas_initial_bed = compute_Y_from_RH(in_param.relative_humidity_bed_init, init_param.pressure_saturated_initial, \
            in_param)
        init_param.moisture_gas_initial_in = compute_Y_from_RH(in_param.relative_humidity_gas_initial, init_param.pressure_saturated_initial, \
            in_param)
        init_param.moisture_gas_end = compute_Y_from_RH(in_param.relative_humidity_gas_end, init_param.pressure_saturated_initial, \
            in_param)
        init_param.moisture_particle_initial = compute_moisture_particle_from_RH(in_param.relative_humidity_bed_init, in_param)
        init_param.moisture_particle_saturated = compute_moisture_particle_from_RH(in_param.relative_humidity_gas_initial, in_param)
        init_param.k_GP_initial = compute_mass_transfer_coefficient_vector(init_param.molar_concentration_moisture_initial, in_param, \
            init_param)[3]
        init_param.constant_initial = init_param.k_GP_initial * init_param.specific_surface_area * init_param.pressure_saturated_initial/in_param.pressure_ambient
        init_param.heat_transfer_coefficient = compute_heat_transfer_coefficient(init_param.molar_concentration_moisture_initial, in_param, init_param)
        init_param.temp_min = min(in_param.temp_initial, in_param.temp_walls)

        #pcc.close()
        computation_of_system(in_param, dis_param, init_param)
        
    ui.pushButton.clicked.connect(compute)
    pcc.show()
    sys.exit(app.exec_())