# Running Simplified Simulation
The `simplified_simulation_LUXE.py` script takes a  .yaml file as input.
This .yaml file is divided into the sections [detector](#detector), [dipole magnet](#dipole&20magnet) and 
[settings](#settings). All values, but the dipole field are converted to SI-conform units.

## detector
* `radiation length`: [kg / m²]
* `rho`: [kg / m³]
* `thickness`: [m]
* `num pixel x:` number pixels in x
* `num pixel y:` number pixels in y

## dipole magnet
* `dipole_field:` list of values describing the dipole field in z [T]
* `dipole start:` [m]
* `dipole end:` [m]

## settings
* `scattering:` True/False
