import pandas as pd
import csv

from src.segment import Segment

from typing import List


class SegmentManager:
    def __init__(self,
                 configuration: dict,
                 detector_geometry: str):
        """Class for handling segments for doublet and triplet creation
        :param configuration: information needed for the segments:
                {
                doublet: {dx/x0: <value>,
                          eps:   <value>},
                triplet: {angle diff x: <value>,
                          angle diff y: <value>}
                binning: {num bins x: <value>}
                qubo parameters: {b_ij conflict: <value>,
                                  b_ij match: value,
                                  a_i: <value>}
                scale range parameters: {z_scores: <value>,
                                         quality: <value>,
                                         interaction: <value>}
                }
        :param detector_geometry: .csv detector layer file geometry file
        """

        self.configuration = configuration   # Dictionary of selection criteria
        self.detector_layers = []
        with open(detector_geometry, 'r') as file:
            csv_reader = csv.reader(file)
            next(csv_reader)    # skip header, csv files should consist of one line
            for row in csv_reader:
                self.detector_layers.append([float(r) for r in (row[1:])])
        self.detector_layers.sort(key=lambda entry: entry[-1])   # sort by z-value

        self.num_layers = len(set([layer[-1] for layer in self.detector_layers]))
        self.segment_list = []   # List of segment objects
        self.segment_mapping = {}   # <segment> (key): [<segment_name_0>, <segment_name_1>, ...] (value)
        self.reference_layer_z = None   # reference layer for dx/x0 criterion
        self.layer_z_values = [entry[-1] for entry in self.detector_layers]
        self.layer_x_min_values = [entry[0] for entry in self.detector_layers]
        self.layer_x_max_values = [entry[1] for entry in self.detector_layers]

        # create segments and mapping
        if self.num_layers == 4:
            self.create_segments_simplified_model()
            self.segment_mapping_simplified_model()
            self.setup = "Simplified LUXE"
        elif self.num_layers == 8:
            self.setup = "Full LUXE"
            pass  # to be implemented --> FullLUXE

    def target_segments(self,
                        name: str):
        """Takes the name of a segment and the target segments are returned
        :param name: name of segment
        :return:
            target segments
        """
        segment_names = [segment_name for segment_name in self.segment_mapping[name]]
        return [segment for segment in self.segment_list if segment.name in segment_names]

    def get_segment_at_known_xz_value(self,
                                      x: float,
                                      z: float):
        """Finds the correct segment to store a hit. The segment list is a flattened 2d array with
        <number layers> * <num segments per layer> entries
        :param x: x value of
        :param z: z value of hit
        :return:
            index of corresponding segment in segment list
        """
        z_index = self.layer_z_values.index(z)
        x_range = self.detector_layers[z_index][1] - self.detector_layers[z_index][0]
        x_bin_size = x_range / self.configuration["binning"]["num bins x"]
        x_index = int((x - self.detector_layers[z_index][0]) / x_bin_size)
        return z_index * self.configuration["binning"]["num bins x"] + x_index

    def segment_mapping_simplified_model(self):
        """Maps the segments according to the doublet preselection criteria.
        Stores the result inside the segment mapping attribute.
        """
        # set reference layer as first layer counting from the IP
        self.reference_layer_z = min(self.layer_z_values)

        # get maximum detector dimension in x
        min_x_detector_dimension = min(self.layer_x_min_values)
        max_x_detector_dimension = max(self.layer_x_max_values)

        for segment in self.segment_list:
            target_list = []
            for target_segment in self.segment_list:
                if target_segment.layer != segment.layer + 1:   # only consecutive layers
                    continue
                if target_segment.x_end < segment.x_start:   # heavy scattering excluded
                    continue

                # max and minimum dx values
                max_dx = target_segment.x_end - segment.x_start
                min_dx = max([target_segment.x_start - segment.x_end, 0])

                # max x_0 range on reference screen
                x0_max = self.x0_at_z_ref(target_segment.x_start, segment.x_end, target_segment.z_position,
                                          segment.z_position)

                x0_min = self.x0_at_z_ref(target_segment.x_end, segment.x_start, target_segment.z_position,
                                          segment.z_position)

                # correct for detector dimensions
                if x0_min < min_x_detector_dimension:
                    x0_min = min_x_detector_dimension
                    if x0_max < min_x_detector_dimension:
                        continue
                if x0_max > max_x_detector_dimension:
                    x0_max = max_x_detector_dimension

                dx_x0_interval = pd.Interval(self.configuration["doublet"]["dx/x0"] -
                                             self.configuration["doublet"]["eps"],
                                             self.configuration["doublet"]["dx/x0"] +
                                             self.configuration["doublet"]["eps"])

                max_dx_interval = pd.Interval(min_dx / x0_max, max_dx / x0_min)

                if dx_x0_interval.overlaps(max_dx_interval):
                    target_list.append(target_segment.name)

            self.segment_mapping.update({segment.name: target_list})

    def x0_at_z_ref(self,
                    x_end: float,
                    x_start: float,
                    z_end: float,
                    z_start: float):
        """Help function for calculation x position of doublet at a z-reference value, which was set before
        as a class attribute
        :param x_end: x-position of target segment
        :param x_start: x-position of segment
        :param z_end: z-position of target segment
        :param z_start: z-position of segment
        :return:
            x_position at the reference layer
        """
        dx = x_end - x_start
        dz = z_end - z_start
        return x_end - dx * (z_end - self.reference_layer_z) / dz

    def create_segments_simplified_model(self):
        """Splitting data into num bins segment to reduce the combinatorial computational costs.
        For the simplified model only.
        """
        for i, layer in enumerate(self.detector_layers):   # ordered in z value from lowest to highest
            x_max = layer[1]
            x_min = layer[0]
            segment_size = (x_max - x_min) / int(self.configuration["binning"]["num bins x"])
            for j in range(int(self.configuration["binning"]["num bins x"])):
                self.segment_list.append(Segment(f"L{i}_S{j}",  # Layer i segment j
                                                 i,
                                                 x_min + j * segment_size,
                                                 x_min + (j + 1) * segment_size,
                                                 self.detector_layers[i][-1]))
