# Relative humidity


def relative_humidity(vapor_pressure, temp_kelvin):
    p_saturated=compute_p_saturated(A, B, temp_kelvin, C)
    RH = vapor_pressure/p_saturated
    return RH

def vapor_pressure_moisture(RH, temp_kelvin):
    p_saturated=compute_p_saturated(A, B, temp_kelvin, C)
    vapor_pressure_moisture=RH*p_saturated
    return vapor_pressure_moisture

