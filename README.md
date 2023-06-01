#  Workflow 

This framework provides a chain of scripts for track reconstruction at [LUXE](https://arxiv.org/abs/2102.02032).

## Simplified simulation
A simplified simulation is employed to compute the path of positrons through a dipole magnet and the interaction with 
a tracking detector. To run the simulation the following arguments are needed:
   * `config:` configuration of the simulation, see [here](docs/simplified_simulation.input.md)
   * `particles:` particles (positrons) from output file of ptarmigan in .h5 format, see [here](https://github.com/tgblackburn)
   * `geometry:` file of the detector configuration in .csv format
   * `target folder:` folder to which results are saved

```bash
python simplified_simulation.py [-h] [--config_file CONFIG_FILE] [--ptarmigan_file PTARMIGAN_FILE] [--geometry_file GEOMETRY_FILE] [--target_folder TARGET_FOLDER]
```

## Preselection and QUBO-preparation
In the preselection, doublet and triplets are created. To reduce the computational costs of this 
combinatorial task, constraints on the creation of multiplets are set. Additionally, coefficients for the
pattern recognition task in the form of a Quadratic Unconstrained Binary Optimisation (QUBO) are applied.
Generator level multiplets from truth information are computed and saved into a single .npy file as a reference
to enable the possibility of calculating efficiency and fake rate.
To run the preselection and QUBO preparation the following arguments are needed:
   * `configuration:` configuration of the preselection, see [here](docs/qubo_preparation.md)
   * `tracking data:` detector hit information in .csv similar [TrackML challenge](https://www.kaggle.com/c/trackml-particle-identification)
   * `geometry:` file of the detector configuration in .csv format
   * `target folder:` folder to which results are saved

```bash
python qubo_preparation.py [-h] [--config_file CONFIG_FILE] [--tracking_data TRACKING_DATA] [--geometry_file GEOMETRY_FILE] [--target_folder TARGET_FOLDER]
```

At the end of the preselection, plots can be created with truth information to check if parameters are set well and results 
make sense.

## QUBO solve
The qubo is solved and results are to a on-the-fly created subdirectory of the folder containing a .npy file with 
triplet objects. The result is saved as .npy file  with various information about the solving process. From these results, xplets are reconstructed and checked if there are ambiguities, this is solved and 
an integrated track reconstruction efficiency and fake rate is calculated and displayed in the terminal. The following arguments are needed:
   * `configuration:` configuration of the qubo, see [here](docs/qubo_solve_input.md) 
   * `qubo folder:` folder with a prepared triplet list

```bash
python qubo_solve.py [-h] [--config_file CONFIG_FILE] [--qubo_folder QUBO_FOLDER]

```






