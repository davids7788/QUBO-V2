#  Workflow 

This framework provides a chain of scripts for track reconstruction at [LUXE](https://arxiv.org/abs/2102.02032).

## Simplified simulation
A simplified simulation is employed to compute the path of positrons through a dipole magnet and the interaction with 
a tracking detector. To run the simulation the following arguments are needed:
   * `configuration:` configuration of the simulation, see [here](docs/simplified_simulation_input.md)
   * `particles:` particles (positrons) from output file of ptarmigan in .h5 format, see [here](https://github.com/tgblackburn)
   * `geometry:` file of the detector configuration in .csv format
   * `save_to:` folder to which results are saved

```bash
python simplified_simulation_LUXE.py configuration particles geometry save_to
```

## Preselection and QUBO-building
In the preselection, possible parts of track candidates, doublets and triplets, are created. To reduce the computational
costs of this combinatorial task, constraints on the creation of these Xplets are set. Additionally, coefficients for the
pattern recognition task in the form of a Quadratic Unconstrained Binary Optimisation (QUBO) are set.
Additionally, generator level Xplets from truth information are computed and saved into a single .npy file.
To run the preselection the following arguments are needed:
   * `configuration:` configuration of the preselection, see [here](docs/qubo_preparation_LUXE.md)
   * `tracking_data:` detector hit information in .csv similar [TrackML challenge](https://www.kaggle.com/c/trackml-particle-identification)
   * `geometry:` file of the detector configuration in .csv format
   * `save_to:` folder to which results are saved

```bash
python qubo_preparation_LUXE.py configuration tracking_data geometry save_to
```

At the end of the preselection, plots can be created with truth information to check if parameters are set well and results 
make sense.

## QUBO-solving
The qubo is solved and results are saved into a .npy file with a dictionary with various information about the 
solving process. From these results, xplets are reconstructed and checked if there are ambiguities, this is solved and 
an integrated track reconstruction efficiency and fake rate is calculated on the fly The following arguments are needed:
   * `configuration:` 
   * `folder_with_triplet_list:` folder with a solved qubo (has a 9-digit number in the front by default)

```bash
python qubo_solve.py configuration folder_with_triplet_list
```






