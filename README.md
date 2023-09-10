#  Workflow 

This framework provides a chain of scripts for track reconstruction at [LUXE](https://arxiv.org/abs/2102.02032).

## Simplified simulation
A simplified simulation is employed to compute the path of positrons through a dipole magnet and the interaction with 
a tracking detector. To run the simulation the following arguments are needed:
   * `configuration:` configuration of the simulation, see [here](docs/simplified_simulation_input.md)
   * `particles:` particles (positrons) form output file of ptarmigan in .h5 format, see [here](https://github.com/tgblackburn)
   * `geometry:` file of the detector configuration in .csv format
   * `save_to:` folder to which results are saved

```bash
python simplified_simulation.py configuration particles geometry save_to
```

## Preselection
In the preselection, possible parts of track candidates, doublets and triplets, are created. To reduce the computational
costs of this combinatorial task, constraints on the creation of these x-plets are set. Additionally, coefficients for the
pattern recognition task in the form of a Quadratic Unconstrained Binary Optimisation (QUBO) are set.
To run the preselection the following arguments are needed:
   * `configuration:` configuration of the preselection, see [here](docs/make_triplets_input.md)
   * `tracking_data:` detector hit information in .csv similar [TrackML challenge](https://www.kaggle.com/c/trackml-particle-identification)
   * `geometry:` file of the detector configuration in .csv format
   * `save_to:` folder to which results are saved

```bash
python create_triplets.py configuration tracking_data geometry save_to
```

At the end of the preselection, plots are created with truth information to check if parameters are set well and results 
make sense.

## Pattern reconstruction
In the pattern recognition, triplets are marked as kept (1) and discarded (0). At the moment, a bit_flip approach, and
a VQE approach are implemented to solve the QUBO. To run the pattern recognition, the following arguments are needed:
   * `configuration:` configuration of the preselection, see [here](docs/qubo_solve_input.md)
   * `qubo_folder:` folder with a .npy triplet list file

```bash
python solve_qubo.py configuration qubo_folder
```
From the results of the qubo, xplets can be created on reconstruction level. The following arguments are needed:
   * `qubo_folder:` folder with a solved qubo (has a 9-digit number in the front by default)
   * `geometry:` file of the detector configuration in .csv format

```bash
python create_reco_xplets.py qubo_folder geometry
```

To create xplets on generator level (truth xplets), the following arguments are needed:
   * `tracking_data:` detector hit information in .csv similar [TrackML challenge](https://www.kaggle.com/c/trackml-particle-identification)
   * `geometry:` file of the detector configuration in .csv format
   * `save_to:` folder to which results are saved
```bash
python gen_multiplets.py tracking_data geometry save_to
```






