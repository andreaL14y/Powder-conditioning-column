
def derivative_m_isotherm(water_activity):
    M0 = 0.0488
    c = 3.23
    f = 1.16
    top = c * f * M0 * (c * f**2 * water_activity**2 - f**2 * water_activity**2 +1)
    bottom = (f * water_activity - 1)**2 * (c * f * water_activity - f * water_activity + 1)**2

    derivative = top/bottom
    return derivative


def compute_moisture_change (water_activity, porosity, diffusion, laplacian):
    derivative = derivative_m_isotherm(water_activity)

    top = diffusion * laplacian
    bottom = (1 - porosity) * derivative + porosity
    moisture_change = top/bottom

    return moisture_change

