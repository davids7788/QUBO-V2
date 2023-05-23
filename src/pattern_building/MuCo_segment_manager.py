import os
import csv
import numpy as np

from pattern_building.MuCo_segment import MuCoDetectorSegment


class MuCoSegmentManager:
    def __init__(self,
                 configuration,
                 geometry_folder):
        """Class for handling segments for doublet and triplet creation for Muon Collider data.
        :param configuration: pattern building configuration file
        :param geometry_folder: folder address with files providing information about the MuonCollider detector geometry
        """
        # num segments in z-direction for barrel layers, in r-direction for endcap layers and in phi-direction for
        # both barrel and endcap layers
        self.num_segments_barrel_z = configuration['segments']['barrel z']
        self.num_segments_endcap_r = configuration['segments']['endcap r']
        self.num_segments_phi = configuration['segments']['phi']

        # max allowed difference in phi for doublets [rad]
        self.doublet_phi_max_vxd = configuration['doublet']['phi']['VXDTrackerBarrel']
        self.doublet_phi_max_inner = configuration['doublet']['phi']['InnerTrackerBarrel']
        self.doublet_phi_max_outer = configuration['doublet']['phi']['OuterTrackerBarrel']

        # angle in phi [rad]
        self.segment_size_phi = 2 * np.pi / self.num_segments_phi

        # folder containing geometry .csv files for the muon collider
        self.geometry_folder = geometry_folder

        # dictionary for storing vxd tracker segments, key: layer_name, value: segment_list
        self.vxd_tracker_barrel_segments = {}
        self.vxd_tracker_endcap_segments = {}

        # dictionary for storing inner tracker segments, key: layer_name, value: segment_list
        self.inner_tracker_barrel_segments = {}
        self.inner_tracker_endcap_segments = {}

        # dictionary for storing outer tracker segments, key: layer_name, value: segment_list
        self.outer_tracker_barrel_segments = {}
        self.outer_tracker_endcap_segments = {}

        # predefining lookup table for layer mapping
        self.possible_layer_mapping = {}

        # save information about layers in a dictionary
        self.detector_layer_information = {}

        self.segment_mapping = {}
        self.segment_mapping_key = {}

    def create_MuCo_segments(self) -> None:
        """Segments are created according to their phi, z and r distance from IP coordinates.
        """
        tracker_files = os.listdir(self.geometry_folder)
        for tracker_file in tracker_files:
            print(f'Processing geometry file: {tracker_file}')
            detector_layers = {}
            with open(f'{self.geometry_folder}/{tracker_file}', 'r') as file:
                csv_reader = csv.reader(file)
                next(csv_reader)  # header structure: ["name", "r_start", "r_end", "z_start", "z_end"]
                for row in csv_reader:
                    detector_layers.update({row[0]: [float(r) for r in (row[1:])]})
            self.detector_layer_information.update({tracker_file: detector_layers})
            # save information in class attribute
            if "Endcap" in tracker_file:
                self.process_endcap_geometry(detector_layers, tracker_file)
            else:
                self.process_barrel_geometry(detector_layers, tracker_file)

    def process_endcap_geometry(self,
                                detector_layers,
                                tracker_file_name) -> None:
        """Creates segments from an endcap geometry file
        :param detector_layers: dictionary with boundary information about the detector layers
        :param tracker_file_name: name of the processed tracker file
        """
        for layer_name, det_limits in detector_layers.items():
            segment_list = []
            segment_size_r = (det_limits[1] - det_limits[0]) / self.num_segments_endcap_r
            for phi in range(self.num_segments_phi):
                for r in range(self.num_segments_endcap_r):
                    new_endcap_segment = MuCoDetectorSegment(name=f'{layer_name}_{r}_{phi}',
                                                             phi_start=- np.pi + phi * self.segment_size_phi,
                                                             phi_end=- np.pi + (phi + 1) * self.segment_size_phi,
                                                             z_start=det_limits[2],
                                                             z_end=det_limits[3],
                                                             r_start=r * segment_size_r + det_limits[0],
                                                             r_end=(r + 1) * segment_size_r + det_limits[0])
                    segment_list.append(new_endcap_segment)
            if tracker_file_name == 'VXDTrackerEndcap.csv':
                self.vxd_tracker_endcap_segments.update({f'{layer_name}': segment_list})
            if tracker_file_name == 'ITrackerEndcap.csv':
                self.inner_tracker_endcap_segments.update({f'{layer_name}': segment_list})
            if tracker_file_name == 'OTrackerEndcap.csv':
                self.outer_tracker_endcap_segments.update({f'{layer_name}': segment_list})

    def process_barrel_geometry(self,
                                detector_layers,
                                tracker_file_name) -> None:
        """Creates segments from an endcap geometry file
        :param detector_layers: dictionary with boundary information about the detector layers
        :param tracker_file_name: name of the processed tracker file
        """
        for layer_name, det_limits in detector_layers.items():
            segment_list = []
            segment_size_z = (det_limits[3] - det_limits[2]) / self.num_segments_barrel_z
            for phi in range(self.num_segments_phi):
                for z in range(self.num_segments_barrel_z):
                    new_barrel_segment = MuCoDetectorSegment(name=f'{layer_name}_{z}_{phi}',
                                                             phi_start=- np.pi + phi * self.segment_size_phi,
                                                             phi_end=- np.pi + (phi + 1) * self.segment_size_phi,
                                                             z_start=z * segment_size_z + det_limits[2],
                                                             z_end=(z + 1) * segment_size_z + det_limits[2],
                                                             r_start=det_limits[0],
                                                             r_end=det_limits[1])
                    segment_list.append(new_barrel_segment)
            if tracker_file_name == 'VXDTracker.csv':
                self.vxd_tracker_barrel_segments.update({f'{layer_name}': segment_list})
            if tracker_file_name == 'ITracker.csv':
                self.inner_tracker_barrel_segments.update({f'{layer_name}': segment_list})
            if tracker_file_name == 'OTracker.csv':
                self.outer_tracker_barrel_segments.update({f'{layer_name}': segment_list})

    def get_segment_at_known_xyz_value(self,
                                       x: float,
                                       y: float,
                                       z: float,
                                       file_name: str,
                                       layer_number: str) -> MuCoDetectorSegment:
        """Returns segment of layer at given coordinates.
        :param x: x value of hit
        :param y: y value of hit
        :param z: z value of hit
        :param file_name: name of file
        :param layer_number: detector layer
        """
        if "_" in layer_number:
            layer_number = int(layer_number.split("_")[0])
        else:
            layer_number = int(float(layer_number))

        # Get detector layer name
        detector_layer_name = MuCoSegmentManager.get_detector_layer_name(z,
                                                                         file_name,
                                                                         layer_number)
        # Get information about the dimensions in r , z and phi of the detector
        detector_file = '.'.join([detector_layer_name.split('_')[0], 'csv'])
        detector_layer_information = self.detector_layer_information[detector_file][detector_layer_name]

        # Compute
        phi = np.arctan2(y, x)

        if "Endcap" in file_name:
            r_start = detector_layer_information[0]
            r_end = detector_layer_information[1]
            segment_size_r = (r_end - r_start) / self.num_segments_endcap_r
            r = np.sqrt(x**2 + y**2)
            segment_position = self.get_endcap_segment_index(segment_size_r,
                                                             r,
                                                             phi,
                                                             r_start)
            if "VXD" in file_name:
                return self.vxd_tracker_endcap_segments[detector_layer_name][segment_position]
            if "ITracker" in file_name:
                return self.inner_tracker_endcap_segments[detector_layer_name][segment_position]
            if "OTracker" in file_name:
                return self.outer_tracker_endcap_segments[detector_layer_name][segment_position]
        else:
            z_start = detector_layer_information[2]
            z_end = detector_layer_information[3]
            segment_size_z = (z_end - z_start) / self.num_segments_barrel_z
            segment_position = self.get_barrel_segment_index(segment_size_z,
                                                             z,
                                                             phi,
                                                             z_start)
            if "VXD" in file_name:
                return self.vxd_tracker_barrel_segments[detector_layer_name][segment_position]
            if "ITracker" in file_name:
                return self.inner_tracker_barrel_segments[detector_layer_name][segment_position]
            if "OTracker" in file_name:
                return self.outer_tracker_barrel_segments[detector_layer_name][segment_position]

    def get_barrel_segment_index(self,
                                 segment_size_z: float,
                                 z: float,
                                 phi: float,
                                 z_start: float) -> int:
        """Returns the index of the barrel segment in the corresponding segment list.
        :param segment_size_z: size of segment in [mm]
        :param z: distance in z from IP [mm]
        :param phi: phi value of hit, computed from xy coordinates [rad]
        :param z_start: start of detector layer in z [mm]
        """
        phi_index = int((phi + np.pi) / self.segment_size_phi)
        z_index = int((z - z_start) / segment_size_z)
        return self.num_segments_barrel_z * phi_index + z_index

    def get_endcap_segment_index(self,
                                 segment_size_r: float,
                                 r: float,
                                 phi: float,
                                 r_start: float) -> int:
        """Returns the index of the endcap segment in the corresponding segment list.
        :param segment_size_r: size of segment in [mm]
        :param r: radial distance from IP [mm]
        :param phi: phi value of hit, computed from xy coordinates [rad]
        :param r_start: start of detector layer in r [mm]
        """
        phi_index = int((phi + np.pi) / self.segment_size_phi)
        r_index = int((r - r_start) / segment_size_r)
        return self.num_segments_phi * r_index + phi_index

    @staticmethod
    def get_detector_layer_name(z: float,
                                file_name: str,
                                layer_number: int) -> str:
        """Returns detector layer for finding it in the segment storage dictionary.
        :param z: z value of hit
        :param file_name: name of file
        :param layer_number: detector layer
        """
        if 'Endcap' in file_name:
            if z < 0:
                add_string = 'Endcap_l'
            else:
                add_string = 'Endcap_r'
        else:
            add_string = '_'
        if 'VXD' in file_name:
            return f'{"VXDTracker"}{add_string}{layer_number}'
        if 'ITracker' in file_name:
            return f'{"ITracker"}{add_string}{layer_number}'
        if 'OTracker' in file_name:
            return f'{"OTracker"}{add_string}{layer_number}'

    def layer_mapping(self,
                      no_endcaps) -> None:
        """Chose the layer mapping based on if endcaps are considered or not.
        """
        if no_endcaps:
            # for the simple mapping also just consecutive layers are considered
            self.simple_layer_mapping()
        else:
            # for the advanced mapping consecutive layers of consecutive layers are also considered
            self.layer_mapping_vxd_tracker_barrel()
            self.layer_mapping_vxd_tracker_endcap()
            self.layer_mapping_inner_tracker_barrel()
            self.layer_mapping_inner_tracker_endcap()
            self.layer_mapping_outer_tracker_barrel()
            self.layer_mapping_outer_tracker_endcap()

    def simple_layer_mapping(self) -> None:
        # Vertex Detector consists of 8 layers, last layer is mapped to first layer the Inner Tracker Barrel
        for i in range(7):
            self.possible_layer_mapping.update({f'VXDTracker_{i}': [f'VXDTracker_{i + 1}']})
        self.possible_layer_mapping.update({f'VXDTracker_7': ['ITracker_0']})

        # Inner Tracker Barrel consists of 3 layers, last layer is mapped to the first layer of the Outer Tracker Barrel
        for i in range(2):
            self.possible_layer_mapping.update({f'ITracker_{i}': [f'ITracker_{i + 1}']})
        self.possible_layer_mapping.update({f'ITracker_2': ['OTracker_0']})

        # Outer Tracker Barrel consists of 3 layers, last layer is mapped to an empty list
        for i in range(2):
            self.possible_layer_mapping.update({f'OTracker_{i}': [f'OTracker_{i + 1}']})
        self.possible_layer_mapping.update({f'OTracker_2': []})

    def layer_mapping_vxd_tracker_barrel(self) -> None:
        # Vertex Detector Barrel consists of 8 layers, last layer is mapped to first layer the Inner Tracker Barrel
        for i in range(6):
            self.possible_layer_mapping.update({f'VXDTracker_{i}': [f'VXDTracker_{i + 1}',
                                                                    f'VXDTracker_{i + 2}',
                                                                    'VXDTrackerEndcap_r0',
                                                                    'VXDTrackerEndcap_r1',
                                                                    'VXDTrackerEndcap_l0',
                                                                    'VXDTrackerEndcap_l1']})
            self.possible_layer_mapping.update({'VXDTracker_6': ['VXDTracker_7',
                                                                 'ITracker_0',
                                                                 'VXDTrackerEndcap_r0',
                                                                 'VXDTrackerEndcap_r1',
                                                                 'VXDTrackerEndcap_l0',
                                                                 'VXDTrackerEndcap_l1']})
            self.possible_layer_mapping.update({'VXDTracker_7': ['ITracker_0',
                                                                 'ITracker_1',
                                                                 'VXDTrackerEndcap_r0',
                                                                 'VXDTrackerEndcap_r1',
                                                                 'VXDTrackerEndcap_l0',
                                                                 'VXDTrackerEndcap_l1']})

    def layer_mapping_vxd_tracker_endcap(self) -> None:
        # Vertex Detector Barrel consists of 8 layers on each side
        for side in ['r', 'l']:
            for i in range(6):
                self.possible_layer_mapping.update({f'VXDTrackerEndcap_{side}{i}': [f'VXDTrackerEndcap_{side}{i + 1}',
                                                                                    f'VXDTrackerEndcap_{side}{i + 2}',
                                                                                    'ITracker_0',
                                                                                    'ITracker_1']})
            self.possible_layer_mapping.update({f'VXDTrackerEndcap_{side}6': [f'VXDTrackerEndcap_{side}7',
                                                                              'ITracker_0',
                                                                              f'ITrackerEndcap_{side}0']})
            self.possible_layer_mapping.update({f'VXDTrackerEndcap_{side}7': ['ITracker_0',
                                                                              f'ITrackerEndcap_{side}0',
                                                                              f'ITrackerEndcap_{side}1']})

    def layer_mapping_inner_tracker_barrel(self):
        self.possible_layer_mapping.update({f'ITracker_0': ['ITracker_1',
                                                            'ITracker_2',
                                                            'ITrackerEndcap_l0',
                                                            'ITrackerEndcap_l1',
                                                            'ITrackerEndcap_r0',
                                                            'ITrackerEndcap_r1']})
        self.possible_layer_mapping.update({f'ITracker_1': ['ITracker_2',
                                                            'OTracker_1',
                                                            'ITrackerEndcap_l0',
                                                            'ITrackerEndcap_r0']})
        self.possible_layer_mapping.update({f'ITracker_2': ['OTracker_0',
                                                            'OTracker_1',
                                                            'ITrackerEndcap_l0',
                                                            'ITrackerEndcap_l1',
                                                            'ITrackerEndcap_r0',
                                                            'ITrackerEndcap_r1']})

    def layer_mapping_inner_tracker_endcap(self):
        for side in ['r', 'l']:
            self.possible_layer_mapping.update({f'ITrackerEndcap_{side}0': ['ITracker_2',
                                                                            'OTracker_0',
                                                                            f'ITrackerEndcap_{side}1',
                                                                            f'ITrackerEndcap_{side}2']})
            self.possible_layer_mapping.update({f'ITrackerEndcap_{side}1': ['OTracker_0',
                                                                            f'ITrackerEndcap_{side}2',
                                                                            f'ITrackerEndcap_{side}3',
                                                                            f'OTrackerEndcap_{side}0',
                                                                            f'OTrackerEndcap_{side}1']})
            self.possible_layer_mapping.update({f'ITrackerEndcap_{side}2': [f'ITrackerEndcap_{side}3',
                                                                            f'ITrackerEndcap_{side}4',
                                                                            f'OTrackerEndcap_{side}0',
                                                                            f'OTrackerEndcap_{side}1']})
            self.possible_layer_mapping.update({f'ITrackerEndcap_{side}3': [f'ITrackerEndcap_{side}4',
                                                                            f'ITrackerEndcap_{side}5',
                                                                            f'OTrackerEndcap_{side}1',
                                                                            f'OTrackerEndcap_{side}2']})
            self.possible_layer_mapping.update({f'ITrackerEndcap_{side}4': [f'ITrackerEndcap_{side}5',
                                                                            f'ITrackerEndcap_{side}6',
                                                                            f'OTrackerEndcap_{side}2',
                                                                            f'OTrackerEndcap_{side}3']})
            self.possible_layer_mapping.update({f'ITrackerEndcap_{side}5': [f'ITrackerEndcap_{side}6',
                                                                            f'OTrackerEndcap_{side}3']})
            self.possible_layer_mapping.update({f'ITrackerEndcap_{side}6': []})

    def layer_mapping_outer_tracker_barrel(self):
        self.possible_layer_mapping.update({f'OTracker_0': ['OTracker_1',
                                                            'OTracker_2',
                                                            'OTrackerEndcap_l0',
                                                            'OTrackerEndcap_l1',
                                                            'OTrackerEndcap_r0',
                                                            'OTrackerEndcap_r1']})
        self.possible_layer_mapping.update({f'OTracker_1': ['OTracker_2',
                                                            'OTrackerEndcap_l0',
                                                            'OTrackerEndcap_l1',
                                                            'OTrackerEndcap_r0',
                                                            'OTrackerEndcap_r1']})
        self.possible_layer_mapping.update({f'OTracker_2': []})

    def layer_mapping_outer_tracker_endcap(self):
        for side in ['r', 'l']:
            for i in range(2):
                self.possible_layer_mapping.update({f'OTrackerEndcap_{side}{i}': [f'OTrackerEndcap_{side}{i + 1}',
                                                                                  f'OTrackerEndcap_{side}{i + 2}']})
            self.possible_layer_mapping.update({f'OTrackerEndcap_{side}2': [f'OTrackerEndcap_{side}3']})
            self.possible_layer_mapping.update({f'OTrackerEndcap_{side}3': []})

    def create_segment_mapping(self,
                               no_endcaps=True) -> None:
        if no_endcaps:
            segment_dictionaries = [self.vxd_tracker_barrel_segments,
                                    self.inner_tracker_barrel_segments,
                                    self.outer_tracker_barrel_segments]
        else:
            segment_dictionaries = [self.vxd_tracker_barrel_segments,
                                    self.vxd_tracker_endcap_segments,
                                    self.inner_tracker_barrel_segments,
                                    self.inner_tracker_endcap_segments,
                                    self.outer_tracker_barrel_segments,
                                    self.outer_tracker_endcap_segments]

        for segment_dictionary in segment_dictionaries:

            # segment dictionaries have following structure --> key: layer_name, value: segment_list
            for detector_layer in segment_dictionary.keys():
                for segment in segment_dictionary[detector_layer]:

                    # Adding key to segment_mapping dictionary
                    self.segment_mapping.update({segment.name: []})
                    self.segment_mapping_key.update({segment.name: segment})

                    # Calculating angle from IP to segment boundaries
                    min_zr, max_zr = MuCoSegmentManager.get_zr_cone_for_segment(segment)

                    # look up target detector from predefined look up table
                    for target_detector in self.possible_layer_mapping[detector_layer]:

                        # Get target detector
                        dictionary_storage = self.get_detector_dictionary(target_detector)

                        # loop over target segments in target detector
                        for target_segment in dictionary_storage[target_detector]:

                            if not self.check_phi_compatibility(segment, target_segment):
                                continue

                            # Calculating target segment angle from IP to segment boundaries
                            target_min_zr, target_max_zr = MuCoSegmentManager.get_zr_cone_for_segment(target_segment)
                            if target_max_zr >= min_zr and target_min_zr <= max_zr:
                                self.segment_mapping[segment.name].append(target_segment)

    def check_phi_compatibility(self,
                                start_segment: MuCoDetectorSegment,
                                target_segment: MuCoDetectorSegment) -> bool:
        """Checks the compatibility of segments of consecutive detector layers with respect to phi
        :param start_segment: segment object
        :param target_segment: segment object
        """
        start_segment_phi_index = int(start_segment.name.split('_')[-1])
        target_segment_phi_index = int(target_segment.name.split('_')[-1])

        # number of allowed difference in
        if 'VXD' in target_segment.name:
            doublet_phi_max = self.doublet_phi_max_vxd
        elif 'ITracker' in target_segment.name:
            doublet_phi_max = self.doublet_phi_max_inner
        elif 'OTracker' in target_segment.name:
            doublet_phi_max = self.doublet_phi_max_outer
        else:
            print('Problem occurred while setting max_phi difference for segments!')
            exit()
        max_distance_to_neighbour_segment = int(doublet_phi_max / self.segment_size_phi) + 1

        # think of minimum distance of segments in a circle
        if abs(start_segment_phi_index - target_segment_phi_index) > self.num_segments_phi / 2:
            if abs(abs(start_segment_phi_index - target_segment_phi_index) - self.num_segments_phi) > \
                    max_distance_to_neighbour_segment:
                return False
            return True
        else:
            if abs(start_segment_phi_index - target_segment_phi_index) > max_distance_to_neighbour_segment:
                return False
            return True

    @staticmethod
    def get_zr_cone_for_segment(segment) -> tuple[float, float]:
        """Calculates min and max angle with respect to the IP for the given segment
        """
        max_zr = np.arctan2(segment.z_end, segment.r_start)
        min_zr = np.arctan2(segment.z_start, segment.r_end)
        return min_zr, max_zr

    def get_detector_dictionary(self,
                                string) -> dict:
        """Returns the corresponding dictionary for the given string
        """
        if "VXDTracker" in string:
            if "Endcap" in string:
                return self.vxd_tracker_endcap_segments
            else:
                return self.vxd_tracker_barrel_segments

        if "ITracker" in string:
            if "Endcap" in string:
                return self.inner_tracker_endcap_segments
            else:
                return self.inner_tracker_barrel_segments

        if "OTracker" in string:
            if "Endcap" in string:
                return self.outer_tracker_endcap_segments
            else:
                return self.outer_tracker_barrel_segments
