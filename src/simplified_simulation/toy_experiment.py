from simplified_simulation.particle import Particle


class MCToyExperiment:
    """Class for accumulating infos about the experiment and executing it. Results are stored inside the result object.
    :param
        name                       : Toy Experiment name, e.g. "test series 1: gauss distributed"
        particle_source            : particle source object : hdf5 file in a Ptarmigan class wrapper
        species                    : 'photon, 'electron' or 'positron
        detector_plane_list        : list of detectors, list has to be ordered by z-coordinate of the detectors
        dipole_magnet              : dipole magnet object
        result                     : object for storing results
        scattering                 : True or False
    """
    def __init__(self,
                 name,
                 particle_source,
                 species,
                 detector_plane_list,
                 dipole_magnet,
                 result,
                 scattering):

        self.name = name
        self.particle_source = particle_source
        self.species = species
        self.detector_plane_list = detector_plane_list
        self.dipole_magnet = dipole_magnet
        self.result = result
        self.scattering = scattering

    def start_experiment(self):
        """Starts the experiment. Th particle moves first to the dipole, the trajectory gets bend and it moves
        to the detector planes."""

        particles = len(self.particle_source.get_particle_info(self.species, 'position'))

        for i in range(particles):

            particle_momentum_log = {}
            particle_position_log = {}

            particle = Particle(position=self.particle_source.get_particle_info(self.species, 'position', i)[1:],
                                momentum=self.particle_source.get_particle_info(self.species, 'momentum', i)[1:],
                                particle_id=str(i))
            
            # particle related dictionaries
            particle_momentum_log.update({'IP': particle.momentum})    
            particle_position_log.update({'IP': particle.position}) 

            # particle moves to dipole location
            z_dist_dipole = abs(particle.position[2] - self.dipole_magnet.dipole_start)
            particle.move_along_beam(z_dist_dipole)
            particle_momentum_log.update({'Magnetic Dipole Start': particle.momentum})
            particle_position_log.update({'Magnetic Dipole Start': particle.position})

            self.dipole_magnet.traversing_through_dipole_magnet(particle)
            particle_momentum_log.update({'Magnetic Dipole End': particle.momentum})
            particle_position_log.update({'Magnetic Dipole End': particle.position})

            for z_detector in sorted(set([plane.z_position for plane in self.detector_plane_list])):
                z_dist = z_detector - particle.position[2]
                particle.move_along_beam(z_dist)

                # Particle interaction with the detector planes
                for k, plane in enumerate(self.detector_plane_list):
                    if plane.z_position != z_detector:
                        continue
                        
                    if any([particle.position[0] < plane.limits_x[0],
                            particle.position[0] > plane.limits_x[1],
                            particle.position[1] < plane.limits_y[0],
                            particle.position[1] > plane.limits_y[1]]):
                        continue

                    # Gaussian scattering part
                    if self.scattering:
                        particle.momentum = plane.gaussian_scattering(particle.momentum)
                    else:
                        pass                        
                    particle_momentum_log.update({f"Plane {k}": particle.momentum})
                    particle_position_log.update({f"Plane {k}": particle.position})

                    # Adds entries with true, smeared and pixelated hits to dictionaries, stored in layer objects
                    plane.true_hits_dictionary.update({particle.particle_ID: particle.position})

            self.result.particle_momentum_history.update({particle.particle_ID: particle_momentum_log})
            self.result.particle_position_history.update({particle.particle_ID: particle_position_log})

        for k, plane in enumerate(self.detector_plane_list):
            self.result.list_of_planes.update({f"Plane {k}": plane})
            self.result.true_detector_hits_dictionary.update({f"Plane {k}": plane.true_hits_dictionary})
