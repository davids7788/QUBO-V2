# Running Simplified Simulation
The Simplified Simulations take a  .yaml file as input.
This .yaml file is divided into the sections [detector](#detector), [dipole magnet](#dipole&20magnet) and [settings](#settings). All numbers are converted to SI-conform units.

## detector
* `radiation length`: [kg / m²]
* `rho`: [kg /m³]
* `thickness`: [m]
* `num pixel x:` number pixels in x
* `num pixel y:` number pixels in y
* `layer positions x_start:` list of values [m]
* `layer positions x_end:` list of values [m]
* `layer positions y_start:` list of values [m]
* `layer positions y_end:` list of values [m]
* `layer positions z:` only needed for simplified setup, list of values [m] 
* `layers right z positions:` only needed for full setup, list of values [m] 
* `layers left z positions:` only needed for full setup, list of values [m] 

## dipole magnet
* `dipole_field:` list of values describing the dipole field in z [T]
* `dipole start:` [m]
* `dipole end:` [m]

## settings
* `scattering:` True/False
