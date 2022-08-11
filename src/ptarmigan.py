import h5py


class PtargmiganSimData:
    """Stores the .hdf5/.h5 file and provides a method to access information about
    momenta and positions of the stored particles.
    :param
        file_name: file name as string, has to be an .hdf5 or h5 file
    """
    def __init__(self,
                 file_name):
        self.data_from_file = h5py.File(file_name, 'r')

    def get_particle_info(self, species, info, particle_number=None):
        """Returns the position/momentum of the particles as a 4-vector.
        :param
            species          : 'photon','electron' or 'positron'
            info             : 'position', 'momentum'
            particle_number  : number of particle in order of the dataset

        :return 4-vector is returned:
            position: ct, x, y, z || original [mm] --> converted to [m]
            momentum: E/c, px, py, pz || original [GeV] -> [MeV]
        """
        if info == 'momentum':
            if particle_number is not None:
                return 1000 * self.data_from_file['final-state/' + species + '/' + info][particle_number]
            else:
                return 1000 * self.data_from_file['final-state/' + species + '/' + info][:]
            
        elif info == 'position':
            if particle_number is not None:
                return 1 / 1000 * self.data_from_file['final-state/' + species + '/' + info][particle_number]
            else:
                return 1 / 1000 * self.data_from_file['final-state/' + species + '/' + info][:]
        
        else:
            if particle_number is not None:
                return self.data_from_file['final-state/' + species + '/' + info][particle_number]
            else:
                return self.data_from_file['final-state/' + species + '/' + info][:]
