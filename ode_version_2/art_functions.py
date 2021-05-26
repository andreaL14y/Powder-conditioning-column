import numpy as np

# Universal constants
h_planck = 6.62607004 * 10 **(-34)      # m2 * kg / s
k_boltzmann = 1.38064852 * 10 **(-23)      # m2 * kg / (s2 * K)
r_gas_constant = 8.314462

# lactose
m1 = 21
c1 = -16.7
k1 = 4.67

m2 = -95
c2 = c1
k2 = k1

def compute_entropy_of_activation_S(moisture):
    entropy_of_activation_S = m2 * c2 * k2 * moisture/((1 - k2 * moisture) * (1 - k2 * moisture + c2 * k2 * moisture))
    return entropy_of_activation_S


def compute_enthalpy_of_activation_H(moisture):
    enthalpy_of_activation_H = m1 * c1 * k1 * moisture/((1 - k1 * moisture) * (1 - k1 * moisture + c1 * k1 * moisture))
    return enthalpy_of_activation_H

def compute_crystal_growth_rate():
    pass