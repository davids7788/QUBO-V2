
##  Simplified Simulation 

The Simplified Simulation is a simulation based on statistics for propagating particles from the IP to the tracking system.
As input it takes files from [ptarmigan](https://github.com/tgblackburn/ptarmigan).

There are two modes: 
* simple setup (sl): four equidistant detector layers resembling 18 chips without gaps between the chips
* full setup (fl): four rows of two partly overlapping detector layers, consisting of 9 chips with a small gap between the chips

Additionally some parameters can be tuned, e.g:
* Scattering ON/OFF
* (In)homogenous magnetic field strength
* Detector parameters (thickness, effective radiation length, density)

The setup is made via an yaml file, an example is given in the examples folder.
The result is in the form of a nested dictionary inside the ouput .npy file.
With make_csv.py one can convert the results into the csv format similar to the [TrackML challenge](https://www.kaggle.com/c/trackml-particle-identification).




