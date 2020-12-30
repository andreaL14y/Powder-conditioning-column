class input_parameter():
    def __init__(self):
        self.R_gas_constant : float
        self.temp_kelvin : float
        self.antoine_constant_A: float
        self.antoine_constant_B: float
        self.antoine_constant_C: float
        self.bed_length: float
        self.column_diameter: float
        self.volumetric_flow_rate: float
        self.temp_initial: float
        self.temp_walls: float
        self.relative_humidity_bed_init: float
        self.relative_humidity_gas_initial: float
        self.relative_humidity_gas_end: float
        self.pressure_ambient: float
        self.porosity_powder: float
        self.N_parameter: float
        self.alpha_parameter: float
        self.particle_density: float
        self.gas_density: float
        self.particle_diameter: float
        self.heat_of_vaporization: float
        self.gas_viscosity: float
        self.moisture_diffusivity : float
        self.molar_mass_moisture: float
        self.molar_mass_dry_air: float
        self.moisture_vapor_heat_capacity: float
        self.moisture_liquid_heat_capacity: float
        self.particle_heat_capacity: float
        self.gas_heat_capacity: float
        self.conductivity_particle: float
        self.conductivity_gas: float
        self.boiling_temperature: float

class discretization_parameter():
    def __init__(self):
        self.max_time: int
        self.n_space_steps: int
        self.n_height_steps: int
        self.resolution: int
        self.height_of_interest: int

class initial_conditions():
    def __init__(self):
        self.gas_velocity: float
        self.specific_surface_area: float
        self.flow_rate: float
        self.superficial_velocity: float
        self.pressure_saturated_initial: float
        self.partial_pressure_moisture_initial: float
        self.molar_concentration_moisture_initial: float
        self.moisture_gas_initial_bed: float
        self.moisture_gas_initial_in: float
        self.moisture_gas_end: float
        self.moisture_particle_initial: float
        self.moisture_particle_saturated: float
        self.k_GP_initial: float
        self.constant_initial: float
        self.heat_transfer_coefficient: float
        self.temp_min: float
        self.kelvin: float
        self.cross_sectional_area: float