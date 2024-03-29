# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\Andi L14y\Desktop\UNI GOT\P6\MatMod\GUI\PCCt.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_pcc(object):
    def setupUi(self, pcc):
        pcc.setObjectName("pcc")
        pcc.resize(585, 711)
        self.verticalLayout = QtWidgets.QVBoxLayout(pcc)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.line = QtWidgets.QFrame(pcc)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout.addWidget(self.line, 1, 2, 1, 1)
        self.formLayout_4 = QtWidgets.QFormLayout()
        self.formLayout_4.setObjectName("formLayout_4")
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_4.setItem(0, QtWidgets.QFormLayout.LabelRole, spacerItem)
        self.label_29 = QtWidgets.QLabel(pcc)
        self.label_29.setObjectName("label_29")
        self.formLayout_4.setWidget(2, QtWidgets.QFormLayout.SpanningRole, self.label_29)
        self.label_27 = QtWidgets.QLabel(pcc)
        self.label_27.setObjectName("label_27")
        self.formLayout_4.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_27)
        self.bed_length = QtWidgets.QLineEdit(pcc)
        self.bed_length.setObjectName("bed_length")
        self.formLayout_4.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.bed_length)
        self.label_35 = QtWidgets.QLabel(pcc)
        self.label_35.setObjectName("label_35")
        self.formLayout_4.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.label_35)
        self.column_diameter = QtWidgets.QLineEdit(pcc)
        self.column_diameter.setObjectName("column_diameter")
        self.formLayout_4.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.column_diameter)
        self.label_36 = QtWidgets.QLabel(pcc)
        self.label_36.setObjectName("label_36")
        self.formLayout_4.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label_36)
        self.volumetric_flow_rate = QtWidgets.QLineEdit(pcc)
        self.volumetric_flow_rate.setObjectName("volumetric_flow_rate")
        self.formLayout_4.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.volumetric_flow_rate)
        self.label_33 = QtWidgets.QLabel(pcc)
        self.label_33.setObjectName("label_33")
        self.formLayout_4.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.label_33)
        self.temp_init = QtWidgets.QLineEdit(pcc)
        self.temp_init.setObjectName("temp_init")
        self.formLayout_4.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.temp_init)
        self.label_34 = QtWidgets.QLabel(pcc)
        self.label_34.setObjectName("label_34")
        self.formLayout_4.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.label_34)
        self.temp_walls = QtWidgets.QLineEdit(pcc)
        self.temp_walls.setObjectName("temp_walls")
        self.formLayout_4.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.temp_walls)
        self.label_31 = QtWidgets.QLabel(pcc)
        self.label_31.setObjectName("label_31")
        self.formLayout_4.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.label_31)
        self.relative_humidity_bed_init = QtWidgets.QLineEdit(pcc)
        self.relative_humidity_bed_init.setObjectName("relative_humidity_bed_init")
        self.formLayout_4.setWidget(9, QtWidgets.QFormLayout.FieldRole, self.relative_humidity_bed_init)
        self.label_32 = QtWidgets.QLabel(pcc)
        self.label_32.setObjectName("label_32")
        self.formLayout_4.setWidget(10, QtWidgets.QFormLayout.LabelRole, self.label_32)
        self.relative_humidity_gas_initial = QtWidgets.QLineEdit(pcc)
        self.relative_humidity_gas_initial.setObjectName("relative_humidity_gas_inlet")
        self.formLayout_4.setWidget(10, QtWidgets.QFormLayout.FieldRole, self.relative_humidity_gas_initial)
        self.label_30 = QtWidgets.QLabel(pcc)
        self.label_30.setObjectName("label_30")
        self.formLayout_4.setWidget(11, QtWidgets.QFormLayout.LabelRole, self.label_30)
        self.relative_humidity_gas_end = QtWidgets.QLineEdit(pcc)
        self.relative_humidity_gas_end.setObjectName("relative_humidity_gas_end")
        self.formLayout_4.setWidget(11, QtWidgets.QFormLayout.FieldRole, self.relative_humidity_gas_end)
        self.label_28 = QtWidgets.QLabel(pcc)
        self.label_28.setObjectName("label_28")
        self.formLayout_4.setWidget(12, QtWidgets.QFormLayout.LabelRole, self.label_28)
        self.pressure_ambient = QtWidgets.QLineEdit(pcc)
        self.pressure_ambient.setObjectName("pressure_ambient")
        self.formLayout_4.setWidget(12, QtWidgets.QFormLayout.FieldRole, self.pressure_ambient)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_4.setItem(3, QtWidgets.QFormLayout.LabelRole, spacerItem1)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_4.setItem(1, QtWidgets.QFormLayout.LabelRole, spacerItem2)
        self.gridLayout.addLayout(self.formLayout_4, 0, 2, 1, 1)
        self.formLayout_5 = QtWidgets.QFormLayout()
        self.formLayout_5.setObjectName("formLayout_5")
        self.label_18 = QtWidgets.QLabel(pcc)
        self.label_18.setObjectName("label_18")
        self.formLayout_5.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_18)
        self.molar_mass_moisture = QtWidgets.QLineEdit(pcc)
        self.molar_mass_moisture.setObjectName("molar_mass_moisture")
        self.formLayout_5.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.molar_mass_moisture)
        self.label_19 = QtWidgets.QLabel(pcc)
        self.label_19.setObjectName("label_19")
        self.formLayout_5.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_19)
        self.molar_mass_dry_air = QtWidgets.QLineEdit(pcc)
        self.molar_mass_dry_air.setObjectName("molar_mass_dry_air")
        self.formLayout_5.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.molar_mass_dry_air)
        self.label_20 = QtWidgets.QLabel(pcc)
        self.label_20.setObjectName("label_20")
        self.formLayout_5.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.label_20)
        self.moisture_vapor_heat_capacity = QtWidgets.QLineEdit(pcc)
        self.moisture_vapor_heat_capacity.setObjectName("heat_capacity_vapor")
        self.formLayout_5.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.moisture_vapor_heat_capacity)
        self.label_21 = QtWidgets.QLabel(pcc)
        self.label_21.setObjectName("label_21")
        self.formLayout_5.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_21)
        self.moisture_liquid_heat_capacity = QtWidgets.QLineEdit(pcc)
        self.moisture_liquid_heat_capacity.setObjectName("heat_capacity_water")
        self.formLayout_5.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.moisture_liquid_heat_capacity)
        self.label_22 = QtWidgets.QLabel(pcc)
        self.label_22.setObjectName("label_22")
        self.formLayout_5.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.label_22)
        self.particle_heat_capacity = QtWidgets.QLineEdit(pcc)
        self.particle_heat_capacity.setObjectName("heat_capacity_particle")
        self.formLayout_5.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.particle_heat_capacity)
        self.label_23 = QtWidgets.QLabel(pcc)
        self.label_23.setObjectName("label_23")
        self.formLayout_5.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label_23)
        self.gas_heat_capacity = QtWidgets.QLineEdit(pcc)
        self.gas_heat_capacity.setObjectName("heat_capacity_air")
        self.formLayout_5.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.gas_heat_capacity)
        self.label_24 = QtWidgets.QLabel(pcc)
        self.label_24.setObjectName("label_24")
        self.formLayout_5.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.label_24)
        self.conductivity_particle = QtWidgets.QLineEdit(pcc)
        self.conductivity_particle.setObjectName("conductivity_particle")
        self.formLayout_5.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.conductivity_particle)
        self.label_25 = QtWidgets.QLabel(pcc)
        self.label_25.setObjectName("label_25")
        self.formLayout_5.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.label_25)
        self.conductivity_gas = QtWidgets.QLineEdit(pcc)
        self.conductivity_gas.setObjectName("conductivity_gas")
        self.formLayout_5.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.conductivity_gas)
        self.label_26 = QtWidgets.QLabel(pcc)
        self.label_26.setObjectName("label_26")
        self.formLayout_5.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.label_26)
        self.boiling_temperature = QtWidgets.QLineEdit(pcc)
        self.boiling_temperature.setObjectName("boiling_temperature")
        self.formLayout_5.setWidget(9, QtWidgets.QFormLayout.FieldRole, self.boiling_temperature)
        spacerItem3 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_5.setItem(0, QtWidgets.QFormLayout.LabelRole, spacerItem3)
        self.gridLayout.addLayout(self.formLayout_5, 2, 2, 1, 1)
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.label = QtWidgets.QLabel(pcc)
        self.label.setObjectName("label")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.label_37 = QtWidgets.QLabel(pcc)
        self.label_37.setObjectName("label_37")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_37)
        self.max_time = QtWidgets.QLineEdit(pcc)
        self.max_time.setObjectName("max_time")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.max_time)
        self.label_38 = QtWidgets.QLabel(pcc)
        self.label_38.setObjectName("label_38")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_38)
        self.n_space_steps = QtWidgets.QLineEdit(pcc)
        self.n_space_steps.setObjectName("n_space_steps")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.n_space_steps)
        self.label_39 = QtWidgets.QLabel(pcc)
        self.label_39.setObjectName("label_39")
        self.formLayout_2.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.label_39)
        self.label_40 = QtWidgets.QLabel(pcc)
        self.label_40.setObjectName("label_40")
        self.formLayout_2.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_40)
        self.resolution = QtWidgets.QLineEdit(pcc)
        self.resolution.setObjectName("resolution")
        self.formLayout_2.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.resolution)
        spacerItem4 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_2.setItem(5, QtWidgets.QFormLayout.LabelRole, spacerItem4)
        self.label_41 = QtWidgets.QLabel(pcc)
        self.label_41.setObjectName("label_41")
        self.formLayout_2.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label_41)
        spacerItem5 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_2.setItem(7, QtWidgets.QFormLayout.LabelRole, spacerItem5)
        self.label_2 = QtWidgets.QLabel(pcc)
        self.label_2.setObjectName("label_2")
        self.formLayout_2.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.R_gas_constant = QtWidgets.QLineEdit(pcc)
        self.R_gas_constant.setObjectName("r_gas_constant")
        self.formLayout_2.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.R_gas_constant)
        self.label_3 = QtWidgets.QLabel(pcc)
        self.label_3.setObjectName("label_3")
        self.formLayout_2.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.label_5 = QtWidgets.QLabel(pcc)
        self.label_5.setObjectName("label_5")
        self.formLayout_2.setWidget(10, QtWidgets.QFormLayout.LabelRole, self.label_5)
        self.antoine_constants_A = QtWidgets.QLineEdit(pcc)
        self.antoine_constants_A.setObjectName("antoine_constants_A")
        self.formLayout_2.setWidget(10, QtWidgets.QFormLayout.FieldRole, self.antoine_constants_A)
        self.label_6 = QtWidgets.QLabel(pcc)
        self.label_6.setObjectName("label_6")
        self.formLayout_2.setWidget(11, QtWidgets.QFormLayout.LabelRole, self.label_6)
        self.antoine_constant_B = QtWidgets.QLineEdit(pcc)
        self.antoine_constant_B.setObjectName("antoine_constant_B")
        self.formLayout_2.setWidget(11, QtWidgets.QFormLayout.FieldRole, self.antoine_constant_B)
        self.label_7 = QtWidgets.QLabel(pcc)
        self.label_7.setObjectName("label_7")
        self.formLayout_2.setWidget(12, QtWidgets.QFormLayout.LabelRole, self.label_7)
        self.antoine_constant_C = QtWidgets.QLineEdit(pcc)
        self.antoine_constant_C.setObjectName("antoine_constant_C")
        self.formLayout_2.setWidget(12, QtWidgets.QFormLayout.FieldRole, self.antoine_constant_C)
        self.height_of_interest = QtWidgets.QComboBox(pcc)
        self.height_of_interest.setObjectName("height_of_interest")
        self.height_of_interest.addItem("")
        self.height_of_interest.addItem("")
        self.formLayout_2.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.height_of_interest)
        self.gridLayout.addLayout(self.formLayout_2, 0, 0, 1, 1)
        self.line_2 = QtWidgets.QFrame(pcc)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayout.addWidget(self.line_2, 1, 0, 1, 1)
        self.formLayout_3 = QtWidgets.QFormLayout()
        self.formLayout_3.setObjectName("formLayout_3")
        self.label_8 = QtWidgets.QLabel(pcc)
        self.label_8.setObjectName("label_8")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.SpanningRole, self.label_8)
        self.label_10 = QtWidgets.QLabel(pcc)
        self.label_10.setObjectName("label_10")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_10)
        self.porosity_powder = QtWidgets.QLineEdit(pcc)
        self.porosity_powder.setObjectName("porosity_powder")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.porosity_powder)
        self.label_11 = QtWidgets.QLabel(pcc)
        self.label_11.setObjectName("label_11")
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.label_11)
        self.N_parameter = QtWidgets.QLineEdit(pcc)
        self.N_parameter.setObjectName("N_parameter")
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.N_parameter)
        self.label_9 = QtWidgets.QLabel(pcc)
        self.label_9.setObjectName("label_9")
        self.formLayout_3.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_9)
        self.label_13 = QtWidgets.QLabel(pcc)
        self.label_13.setObjectName("label_13")
        self.formLayout_3.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.label_13)
        self.particle_density = QtWidgets.QLineEdit(pcc)
        self.particle_density.setObjectName("density_particle")
        self.formLayout_3.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.particle_density)
        self.label_14 = QtWidgets.QLabel(pcc)
        self.label_14.setObjectName("label_14")
        self.formLayout_3.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label_14)
        self.gas_density = QtWidgets.QLineEdit(pcc)
        self.gas_density.setObjectName("density_gas")
        self.formLayout_3.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.gas_density)
        self.label_12 = QtWidgets.QLabel(pcc)
        self.label_12.setObjectName("label_12")
        self.formLayout_3.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.label_12)
        self.particle_diameter = QtWidgets.QLineEdit(pcc)
        self.particle_diameter.setObjectName("particle_diameter")
        self.formLayout_3.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.particle_diameter)
        self.label_15 = QtWidgets.QLabel(pcc)
        self.label_15.setObjectName("label_15")
        self.formLayout_3.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.label_15)
        self.heat_of_vaporization = QtWidgets.QLineEdit(pcc)
        self.heat_of_vaporization.setObjectName("heat_of_sorption")
        self.formLayout_3.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.heat_of_vaporization)
        self.label_16 = QtWidgets.QLabel(pcc)
        self.label_16.setObjectName("label_16")
        self.formLayout_3.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.label_16)
        self.gas_viscosity = QtWidgets.QLineEdit(pcc)
        self.gas_viscosity.setObjectName("gas_viscosity")
        self.formLayout_3.setWidget(9, QtWidgets.QFormLayout.FieldRole, self.gas_viscosity)
        self.label_17 = QtWidgets.QLabel(pcc)
        self.label_17.setObjectName("label_17")
        self.formLayout_3.setWidget(10, QtWidgets.QFormLayout.LabelRole, self.label_17)
        self.moisture_diffusivity = QtWidgets.QLineEdit(pcc)
        self.moisture_diffusivity.setObjectName("moisture_diffusivity")
        self.formLayout_3.setWidget(10, QtWidgets.QFormLayout.FieldRole, self.moisture_diffusivity)
        self.alpha_parameter = QtWidgets.QLineEdit(pcc)
        self.alpha_parameter.setObjectName("alpha_parameter")
        self.formLayout_3.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.alpha_parameter)
        spacerItem6 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_3.setItem(1, QtWidgets.QFormLayout.FieldRole, spacerItem6)
        self.gridLayout.addLayout(self.formLayout_3, 2, 0, 1, 1)
        self.line_4 = QtWidgets.QFrame(pcc)
        self.line_4.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.gridLayout.addWidget(self.line_4, 0, 1, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.line_3 = QtWidgets.QFrame(pcc)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.verticalLayout.addWidget(self.line_3)
        self.pushButton = QtWidgets.QPushButton(pcc)
        self.pushButton.setObjectName("pushButton")
        self.verticalLayout.addWidget(self.pushButton)

        self.retranslateUi(pcc)
        QtCore.QMetaObject.connectSlotsByName(pcc)
        pcc.setTabOrder(self.max_time, self.n_space_steps)
        pcc.setTabOrder(self.n_space_steps, self.height_of_interest)
        pcc.setTabOrder(self.height_of_interest, self.resolution)
        pcc.setTabOrder(self.resolution, self.R_gas_constant)
        pcc.setTabOrder(self.R_gas_constant, self.antoine_constants_A)
        pcc.setTabOrder(self.antoine_constants_A, self.antoine_constant_B)
        pcc.setTabOrder(self.antoine_constant_B, self.antoine_constant_C)
        pcc.setTabOrder(self.antoine_constant_C, self.bed_length)
        pcc.setTabOrder(self.bed_length, self.column_diameter)
        pcc.setTabOrder(self.column_diameter, self.volumetric_flow_rate)
        pcc.setTabOrder(self.volumetric_flow_rate, self.temp_init)
        pcc.setTabOrder(self.temp_init, self.temp_walls)
        pcc.setTabOrder(self.temp_walls, self.relative_humidity_bed_init)
        pcc.setTabOrder(self.relative_humidity_bed_init, self.relative_humidity_gas_initial)
        pcc.setTabOrder(self.relative_humidity_gas_initial, self.relative_humidity_gas_end)
        pcc.setTabOrder(self.relative_humidity_gas_end, self.pressure_ambient)
        pcc.setTabOrder(self.pressure_ambient, self.porosity_powder)
        pcc.setTabOrder(self.porosity_powder, self.N_parameter)
        pcc.setTabOrder(self.N_parameter, self.alpha_parameter)
        pcc.setTabOrder(self.alpha_parameter, self.particle_density)
        pcc.setTabOrder(self.particle_density, self.gas_density)
        pcc.setTabOrder(self.gas_density, self.particle_diameter)
        pcc.setTabOrder(self.particle_diameter, self.heat_of_vaporization)
        pcc.setTabOrder(self.heat_of_vaporization, self.gas_viscosity)
        pcc.setTabOrder(self.gas_viscosity, self.moisture_diffusivity)
        pcc.setTabOrder(self.moisture_diffusivity, self.molar_mass_moisture)
        pcc.setTabOrder(self.molar_mass_moisture, self.molar_mass_dry_air)
        pcc.setTabOrder(self.molar_mass_dry_air, self.moisture_vapor_heat_capacity)
        pcc.setTabOrder(self.moisture_vapor_heat_capacity, self.moisture_liquid_heat_capacity)
        pcc.setTabOrder(self.moisture_liquid_heat_capacity, self.particle_heat_capacity)
        pcc.setTabOrder(self.particle_heat_capacity, self.gas_heat_capacity)
        pcc.setTabOrder(self.gas_heat_capacity, self.conductivity_particle)
        pcc.setTabOrder(self.conductivity_particle, self.conductivity_gas)
        pcc.setTabOrder(self.conductivity_gas, self.boiling_temperature)
        pcc.setTabOrder(self.boiling_temperature, self.pushButton)

    def retranslateUi(self, pcc):
        _translate = QtCore.QCoreApplication.translate
        pcc.setWindowTitle(_translate("pcc", "Dialog"))
        self.label_29.setText(_translate("pcc", "Column and flow specific parameters:"))
        self.label_27.setText(_translate("pcc", "Bed length (m): "))
        self.bed_length.setText(_translate("pcc", "0.2"))
        self.label_35.setText(_translate("pcc", "Column diameter (m):"))
        self.column_diameter.setText(_translate("pcc", "0.1"))
        self.label_36.setText(_translate("pcc", "Volumetric flow rate (l/min):"))
        self.volumetric_flow_rate.setText(_translate("pcc", "1"))
        self.label_33.setText(_translate("pcc", "Initial temperature (K):"))
        self.temp_init.setText(_translate("pcc", "293.15"))
        self.label_34.setText(_translate("pcc", "Temperature walls (K):"))
        self.temp_walls.setText(_translate("pcc", "293.15"))
        self.label_31.setText(_translate("pcc", "Initial relative humidity (bed):"))
        self.relative_humidity_bed_init.setText(_translate("pcc", "0.2"))
        self.label_32.setText(_translate("pcc", "Initial relative humidity (gas):"))
        self.relative_humidity_gas_initial.setText(_translate("pcc", "0.9"))
        self.label_30.setText(_translate("pcc", "End relative humidity (gas):"))
        self.relative_humidity_gas_end.setText(_translate("pcc", "0.2"))
        self.label_28.setText(_translate("pcc", "Ambient pressure (Pa):"))
        self.pressure_ambient.setText(_translate("pcc", "101325"))
        self.label_18.setText(_translate("pcc", "Molar mass of moisture (kg/mol):"))
        self.molar_mass_moisture.setText(_translate("pcc", "0.018"))
        self.label_19.setText(_translate("pcc", "Molar mass of dry air (kg/mol):"))
        self.molar_mass_dry_air.setText(_translate("pcc", "0.02897"))
        self.label_20.setText(_translate("pcc", "Moisture vapor heat capacity (J/(kg*K)):"))
        self.moisture_vapor_heat_capacity.setText(_translate("pcc", "2000"))
        self.label_21.setText(_translate("pcc", "Moisture liquid heat capacity (J/(kg*K)):"))
        self.moisture_liquid_heat_capacity.setText(_translate("pcc", "4000"))
        self.label_22.setText(_translate("pcc", "Particle heat capacity (J/(kg*K)):"))
        self.particle_heat_capacity.setText(_translate("pcc", "1000"))
        self.label_23.setText(_translate("pcc", "Gas heat capacity (J/(kg*K)):"))
        self.gas_heat_capacity.setText(_translate("pcc", "1000"))
        self.label_24.setText(_translate("pcc", "Conductivity particle (W/(m*K)):"))
        self.conductivity_particle.setText(_translate("pcc", "0.1"))
        self.label_25.setText(_translate("pcc", "Conductivity gas (W/(m*K)):"))
        self.conductivity_gas.setText(_translate("pcc", "0.01"))
        self.label_26.setText(_translate("pcc", "Boiling temperature (K):"))
        self.boiling_temperature.setText(_translate("pcc", "373.15"))
        self.label.setText(_translate("pcc", "Discretization parameters:"))
        self.label_37.setText(_translate("pcc", "Maximal time (s):"))
        self.max_time.setText(_translate("pcc", "550000"))
        self.label_38.setText(_translate("pcc", "Steps in lengths:"))
        self.n_space_steps.setText(_translate("pcc", "10"))
        self.label_39.setText(_translate("pcc", "Height of interest:"))
        self.label_40.setText(_translate("pcc", "Resolution"))
        self.resolution.setText(_translate("pcc", "1000"))
        self.label_41.setText(_translate("pcc", "Air specific parameters:"))
        self.label_2.setText(_translate("pcc", "Gas constant:"))
        self.R_gas_constant.setText(_translate("pcc", "8.314"))
        self.label_3.setText(_translate("pcc", "Antoine constants:"))
        self.label_5.setText(_translate("pcc", "A:"))
        self.antoine_constants_A.setText(_translate("pcc", "8.07131 "))
        self.label_6.setText(_translate("pcc", "B:"))
        self.antoine_constant_B.setText(_translate("pcc", "1730.630"))
        self.label_7.setText(_translate("pcc", "C:"))
        self.antoine_constant_C.setText(_translate("pcc", "233.426"))
        self.height_of_interest.setItemText(0, _translate("pcc", "in the middle"))
        self.height_of_interest.setItemText(1, _translate("pcc", "at the wall"))
        self.label_8.setText(_translate("pcc", "Material specific parameters:"))
        self.label_10.setText(_translate("pcc", "Porosity powder:"))
        self.porosity_powder.setText(_translate("pcc", "0.6"))
        self.label_11.setText(_translate("pcc", "N:"))
        self.N_parameter.setText(_translate("pcc", "1"))
        self.label_9.setText(_translate("pcc", "Alpha (10<Alpha<100):"))
        self.label_13.setText(_translate("pcc", "Particle density (kg/m^3): "))
        self.particle_density.setText(_translate("pcc", "1500"))
        self.label_14.setText(_translate("pcc", "Gas density (kg/m^3):"))
        self.gas_density.setText(_translate("pcc", "1"))
        self.label_12.setText(_translate("pcc", "Particle diameter (m):"))
        self.particle_diameter.setText(_translate("pcc", "0.00001"))
        self.label_15.setText(_translate("pcc", "Heat of vaporization (J/kg):"))
        self.heat_of_vaporization.setText(_translate("pcc", "1000000"))
        self.label_16.setText(_translate("pcc", "Gas viscosity (kg/(m*s)):"))
        self.gas_viscosity.setText(_translate("pcc", "0.00001"))
        self.label_17.setText(_translate("pcc", "Moisture diffusivity (m^2/s): "))
        self.moisture_diffusivity.setText(_translate("pcc", "0.00001"))
        self.alpha_parameter.setText(_translate("pcc", "25"))
        self.pushButton.setText(_translate("pcc", "Compute"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    pcc = QtWidgets.QDialog()
    ui = Ui_pcc()
    ui.setupUi(pcc)
    pcc.show()
    sys.exit(app.exec_())
