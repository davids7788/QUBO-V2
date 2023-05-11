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

        self.geometry_folder = geometry_folder

        self.vxd_tracker_barrel_segments = {}
        self.vxd_tracker_endcap_segments = {}

        self.inner_tracker_barrel_segments = {}
        self.inner_tracker_endcap_segments = {}

        self.outer_tracker_barrel_segments = {}
        self.outer_tracker_endcap_segments = {}

        self.segment_mapping = {}

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
            if "Endcap" in tracker_file:
                self.process_endcap_geometry(detector_layers, tracker_file)
            else:
                self.process_barrel_geometry(detector_layers, tracker_file)

    def process_endcap_geometry(self, detector_layers, tracker_file_name) -> None:
        """Creates segments from an endcap geometry file
        :param detector_layers: dictionary with boundary information about the detector layers
        :param tracker_file_name: name of the processed tracker file
        """
        segment_size_phi = 2 * np.pi / self.segments_endcap_phi
        for name, det_limits in detector_layers.items():
            segment_list = []
            segment_size_r = (det_limits[1] - det_limits[0]) / self.segments_endcap_r
            for phi in range(self.segments_endcap_phi):
                for r in range(self.segments_endcap_r):
                    new_barrel_segment = MuCoDetectorSegment(name=f'{tracker_file_name}_{r}_{phi}',
                                                             phi_start=phi * segment_size_phi,
                                                             phi_end=(phi + 1) * segment_size_phi,
                                                             z_start=det_limits[2],
                                                             z_end=det_limits[3],
                                                             r_start=r * segment_size_r + det_limits[0],
                                                             r_end=(r + 1) * segment_size_r + det_limits[0])
                    segment_list.append(new_barrel_segment)
                    if tracker_file_name == 'VXDTrackerEndcaps.csv':
                        self.vxd_tracker_endcap_segments.update({f'{det_limits[2]}_{det_limits[3]}': segment_list})
                    if tracker_file_name == 'ITrackerEndcaps.csv':
                        self.inner_tracker_endcap_segments.update({f'{det_limits[2]}_{det_limits[3]}': segment_list})
                    if tracker_file_name == 'OTrackerEndcaps.csv':
                        self.outer_tracker_endcap_segments.update({f'{det_limits[2]}_{det_limits[3]}': segment_list})

    def process_barrel_geometry(self, detector_layers, tracker_file_name) -> None:
        """Creates segments from an endcap geometry file
        :param detector_layers: dictionary with boundary information about the detector layers
        :param tracker_file_name: name of the processed tracker file
        """
        segment_size_phi = 2 * np.pi / self.segments_barrel_phi
        for name, det_limits in detector_layers.items():
            segment_list = []
            segment_size_z = (det_limits[3] - det_limits[2]) / self.segments_barrel_z
            for phi in range(self.segments_endcap_phi):
                for z in range(self.segments_endcap_r):
                    new_barrel_segment = MuCoDetectorSegment(name=f'{tracker_file_name}_{z}_{phi}',
                                                             phi_start=phi * segment_size_phi,
                                                             phi_end=(phi + 1) * segment_size_phi,
                                                             z_start=z * segment_size_z + det_limits[0],
                                                             z_end=(z + 1) * segment_size_z + det_limits[0],
                                                             r_start=det_limits[0],
                                                             r_end=det_limits[1])
                    segment_list.append(new_barrel_segment)
                    if tracker_file_name == 'VXDTracker.csv':
                        self.vxd_tracker_barrel_segments.update({f'{det_limits[0]}_{det_limits[1]}': segment_list})
                    if tracker_file_name == 'ITracker.csv':
                        self.inner_tracker_barrel_segments.update({f'{det_limits[0]}_{det_limits[1]}': segment_list})
                    if tracker_file_name == 'OTracker.csv':
                        self.outer_tracker_barrel_segments.update({f'{det_limits[0]}_{det_limits[1]}': segment_list})
