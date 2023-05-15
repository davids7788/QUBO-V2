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
        self.segments_barrel_z = configuration['segments']['barrel']['z']
        self.segments_barrel_phi = configuration['segments']['barrel']['phi']
        self.segments_endcap_r = configuration['segments']['endcap']['r']
        self.segments_endcap_phi = configuration['segments']['endcap']['phi']

        self.segment_size_phi = 2 * np.pi / self.segments_barrel_phi

        self.geometry_folder = geometry_folder

        self.vxd_tracker_barrel_segments = {}
        self.vxd_tracker_endcap_segments = {}

        self.inner_tracker_barrel_segments = {}
        self.inner_tracker_endcap_segments = {}

        self.outer_tracker_barrel_segments = {}
        self.outer_tracker_endcap_segments = {}

        self.possible_layer_mapping = {}

        self.detector_layer_information = {}

        self.segment_mapping = {}
        self.segment_mapping_key = {}

    def create_MuCo_segments(self) -> None:
        """Segments are created according to their phi, theta and distance from IP coordinates.
        The name of the segments gives information about their position and layer.
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

    def process_endcap_geometry(self, detector_layers, tracker_file_name) -> None:
        """Creates segments from an endcap geometry file
        :param detector_layers: dictionary with boundary information about the detector layers
        :param tracker_file_name: name of the processed tracker file
        """
        for layer_name, det_limits in detector_layers.items():
            segment_list = []
            segment_size_r = (det_limits[1] - det_limits[0]) / self.segments_endcap_r
            for phi in range(self.segments_endcap_phi):
                for r in range(self.segments_endcap_r):
                    new_barrel_segment = MuCoDetectorSegment(name=f'{layer_name}_{r}_{phi}',
                                                             phi_start=- np.pi + phi * self.segment_size_phi,
                                                             phi_end=- np.pi + (phi + 1) * self.segment_size_phi,
                                                             z_start=det_limits[2],
                                                             z_end=det_limits[3],
                                                             r_start=r * segment_size_r + det_limits[0],
                                                             r_end=(r + 1) * segment_size_r + det_limits[0])
                    segment_list.append(new_barrel_segment)
            if tracker_file_name == 'VXDTrackerEndcap.csv':
                self.vxd_tracker_endcap_segments.update({f'{layer_name}': segment_list})
            if tracker_file_name == 'ITrackerEndcap.csv':
                self.inner_tracker_endcap_segments.update({f'{layer_name}': segment_list})
            if tracker_file_name == 'OTrackerEndcap.csv':
                self.outer_tracker_endcap_segments.update({f'{layer_name}': segment_list})

    def process_barrel_geometry(self, detector_layers, tracker_file_name) -> None:
        """Creates segments from an endcap geometry file
        :param detector_layers: dictionary with boundary information about the detector layers
        :param tracker_file_name: name of the processed tracker file
        """
        for layer_name, det_limits in detector_layers.items():
            segment_list = []
            segment_size_z = (det_limits[3] - det_limits[2]) / self.segments_barrel_z
            for phi in range(self.segments_endcap_phi):
                for z in range(self.segments_barrel_z):
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
        """Gets segment of layer at given coordinates.
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

        detector_layer_name = MuCoSegmentManager.get_detector_layer_name_at_known_xyz_value(z,
                                                                                            file_name,
                                                                                            layer_number)
        detector_file = '.'.join([detector_layer_name.split('_')[0], 'csv'])
        detector_layer_information = self.detector_layer_information[detector_file][detector_layer_name]

        phi = np.arctan2(y, x)

        if "Endcap" in file_name:
            segment_size_phi = 2 * np.pi / self.segments_endcap_phi
            r_start = detector_layer_information[0]
            r_end = detector_layer_information[1]
            segment_size_r = (r_end - r_start) / self.segments_endcap_r
            r = np.sqrt(x**2 + y**2)
            segment_position = self.get_endcap_segment(segment_size_phi,
                                                       segment_size_r,
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
            segment_size_phi = 2 * np.pi / self.segments_barrel_phi

            z_start = detector_layer_information[2]
            z_end = detector_layer_information[3]
            segment_size_z = (z_end - z_start) / self.segments_barrel_z
            segment_position = self.get_barrel_segment(segment_size_phi,
                                                       segment_size_z,
                                                       z,
                                                       phi,
                                                       z_start)
            if "VXD" in file_name:
                return self.vxd_tracker_barrel_segments[detector_layer_name][segment_position]
            if "ITracker" in file_name:
                return self.inner_tracker_barrel_segments[detector_layer_name][segment_position]
            if "OTracker" in file_name:
                return self.outer_tracker_barrel_segments[detector_layer_name][segment_position]

    def get_barrel_segment(self,
                           segment_size_phi: float,
                           segment_size_z: float,
                           z: float,
                           phi: float,
                           z_start: float):
        """Returns the index of the barrel segment in the corresponding segment list"""
        phi_index = int((phi + np.pi) / segment_size_phi)
        z_index = int((z - z_start) / segment_size_z)
        return self.segments_barrel_z * phi_index + z_index

    def get_endcap_segment(self,
                           segment_size_phi: float,
                           segment_size_r: float,
                           r: float,
                           phi: float,
                           r_start: float):
        """Returns the index of the endcap segment in the corresponding segment list"""
        phi_index = int((phi + np.pi) / segment_size_phi)
        r_index = int((r - r_start) / segment_size_r)
        return self.segments_endcap_phi * r_index + phi_index

    @staticmethod
    def get_detector_layer_name_at_known_xyz_value(z: float,
                                                   file_name: str,
                                                   layer_number: int) -> str:
        """Gets detector layer for finding it in the segment storage dictionary.
        :param z: z value of hit
        :param file_name: name of file
        :param layer_number: detector layer
        """
        if "Endcap" in file_name:
            if z < 0:
                add_string = 'Endcap_l'
            else:
                add_string = 'Endcap_r'
        else:
            add_string = '_'
        if "VXD" in file_name:
            return f'{"VXDTracker"}{add_string}{layer_number}'
        if "ITracker" in file_name:
            return f'{"ITracker"}{add_string}{layer_number}'
        if "OTracker" in file_name:
            return f'{"OTracker"}{add_string}{layer_number}'

    def layer_mapping(self, no_endcaps):
        if no_endcaps:
            self.simple_layer_mapping()
        else:
            self.layer_mapping_vxd_tracker_barrel()
            self.layer_mapping_vxd_tracker_endcap()
            self.layer_mapping_inner_tracker_barrel()
            self.layer_mapping_inner_tracker_endcap()
            self.layer_mapping_outer_tracker_barrel()
            self.layer_mapping_outer_tracker_endcap()

    def simple_layer_mapping(self):
        for i in range(7):
            self.possible_layer_mapping.update({f'VXDTracker_{i}': [f'VXDTracker_{i + 1}']})
        self.possible_layer_mapping.update({f'VXDTracker_7': ['ITracker_0']})
        for i in range(2):
            self.possible_layer_mapping.update({f'ITracker_{i}': [f'ITracker_{i + 1}']})
        self.possible_layer_mapping.update({f'ITracker_2': ['OTracker_0']})
        for i in range(2):
            self.possible_layer_mapping.update({f'OTracker_{i}': [f'OTracker_{i + 1}']})
        self.possible_layer_mapping.update({f'OTracker_2': []})

    def layer_mapping_vxd_tracker_barrel(self):
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

    def layer_mapping_vxd_tracker_endcap(self):
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

    def create_segment_mapping(self, no_endcaps=True):
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
            for detector_layer in segment_dictionary.keys():
                for segment in segment_dictionary[detector_layer]:
                    self.segment_mapping.update({segment.name: set()})
                    self.segment_mapping_key.update({segment.name: segment})
                    max_zr = np.arctan2(segment.z_end, segment.r_start)
                    min_zr = np.arctan2(segment.z_start, segment.r_end)
                    for target_detector in self.possible_layer_mapping[detector_layer]:
                        dictionary_storage = self.get_detector_dictionary(target_detector)
                        for target_segment in dictionary_storage[target_detector]:
                            if abs(segment.phi_start) == abs(target_segment.phi_start) or \
                                    abs(segment.phi_start) == abs(target_segment.phi_end) or \
                                    abs(segment.phi_end) == abs(target_segment.phi_start):
                                target_max_zr = np.arctan2(target_segment.z_end, target_segment.r_start)
                                target_min_zr = np.arctan2(target_segment.z_start, target_segment.r_end)
                                if target_max_zr >= min_zr and target_min_zr <= max_zr or \
                                        max_zr >= target_min_zr and min_zr <= target_max_zr or \
                                    target_max_zr <= max_zr and target_min_zr >= min_zr or \
                                    max_zr <= target_max_zr and min_zr >= target_min_zr:
                                    self.segment_mapping[segment.name].add(target_segment)

    def get_detector_dictionary(self, string):
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
