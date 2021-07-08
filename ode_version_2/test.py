from input_parameters import *

weight = 0.00015                            # kg, 150 mg
density_powder = 500                        # kg/m3
volume = weight/density_powder              # m3
diameter_vial = 0.005                       # m, 5 mm
area_vial = np.pi * (diameter_vial/2)**2    # m2

height = volume/area_vial

print('Height would be:\n', height, 'm\n', height*1000, 'mm')


# Want to finish at 0.99998, 0.0000
M = 0
H = 1

diff = H-M
while abs(diff) > 0.001:
    diff = 0.99998 - M
    if diff > 0:
        H = H * 0.5
    else:
        H = H * 1.01
    M = 1-H


gas_fraction = np.array([0, 0.5, 1, -1, 0])
diff = np.array([1, 0.5, 1, 0, 0])
max_diff = 0.5
gas_to_min = 0.2

print(gas_fraction)
gas_fraction = np.where(diff > max_diff, gas_fraction - gas_to_min/2, gas_fraction)
print(gas_fraction)
