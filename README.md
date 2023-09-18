#  Workflow 

This framework provides a chain of scripts for track reconstruction at [LUXE](https://arxiv.org/abs/2102.02032).

## Simplified Simulation
A simplified simulation is employed to compute the path of positrons through a dipole magnet and the interaction with 
a tracking detector. To run the simulation the following arguments are needed:
   * `config_file:` configuration of the simulation, see [here](docs/simplified_simulation_input.md)
   * `ptarmigan_file:` particles (positrons) form output file of ptarmigan in .h5 format, see [here](https://github.com/tgblackburn)
   * `geometry_file:` file of the detector configuration in .csv format
   * `target_folder:` folder to which results are saved

```bash
python simplified_simulation.py --config_file <configuration> --ptarmigan_file  <ptarmigan_file> --geometry_file <geometry_file> --target_folder <target_folder>
```

## Pattern Building (QUBO preparation)
In the preselection, possible parts of track candidates, doublets and triplets, are created. To reduce the computational
costs of this combinatorial task, constraints on the creation of these multiplets-plets are set. Additionally, coefficients for the
pattern recognition task in the form of a Quadratic Unconstrained Binary Optimisation (QUBO) are set.
To run the preselection the following arguments are needed:
   * `config_file:` configuration of the preselection, see [here](docs/qubo_preparation_input.md)
   * `tracking_data:` detector hit information in .csv similar [TrackML challenge](https://www.kaggle.com/c/trackml-particle-identification)
   * `geometry_file:` file of the detector configuration in .csv format
   * `target_folder:` folder to which results are saved
   * `sample_composition:` signal, signal+background, blinded
   * `simulation_tool:` simplified_simulation, key4hep_csv 

```bash
python qubo_preparation --config_file <config_file> --ptarmigan_file <ptarmigan_file> --geometry_file <geometry_file> --target_folder <target_folder> --smaple_composition <sample_composition> --simulation_tool <simulation_tool>
```

At the end of the preselection, plots are created with truth information to check if parameters are set well and results 
make sense.

## QUBO solve
The QUBO cost function is minimised. Various backends and quantum circuit compositions can be built with the
specified variables in the configuration_files. Returned are efficiency and a .npy file with tracks build out of the
remaining triplet candidates.
   * `config_file:` configuration of qubo_solve, see [here](docs/qubo_solve_input.md)
   * `qubo_folder:` folder with a .npy triplet list file

```bash
python qubo_solve.py --config_file <configuration> --qubo_folder <qubo_folder>
```






