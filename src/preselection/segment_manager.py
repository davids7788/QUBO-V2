import pandas as pd
import csv

from preselection.segment import Segment


class SegmentManager:
    def __init__(self,
                 configuration: dict,
                 detector_geometry: str):
        """Class for handling segments for doublet and triplet creation
        :param configuration: information needed for the segments, delivered by loading yaml file as a nested
               python dictionary:
                {
                doublet: {dx/x0: float,
                          eps:   float>,
                          dy/y0: float},
                triplet: {max scattering: float}
                binning: {num bins x: int,
                          num bins y: int}
                qubo parameters: {b_ij conflict: float
                                  b_ij match: <name of an implemented function> or float,
                                  a_i: <name of an implemented function> or float}
                scale range parameters: {z_scores: True or False,
                                         quality: null (= None or [float, float],
                                         interaction: null (= None) or [float, float]}
                }
        :param detector_geometry: string with address of a.csv geometry file
        """

        self.configuration = configuration   # dictionary of selection criteria
        self.detector_layers = []   # list containing coordinate information about detector layers

        # load detector layers / chips information and add them to the detector layers list
        with open(detector_geometry, 'r') as file:
            csv_reader = csv.reader(file)
            next(csv_reader)    # skip header, csv files should consist of one line
            for row in csv_reader:
                self.detector_layers.append([float(r) for r in (row[1:])])

        # information string
        print(f"Loaded file {detector_geometry.split('/')[-1]} for segmentation algorithm\n")

        self.segment_mapping = {}   # <segment> (key): [<segment_name_0>, <segment_name_1>, ...] (value)
        self.segment_list = []   # List of segment objects

    def create_segments_simplified_LUXE(self):
        """Splitting data into num bins segment to reduce the combinatorial computational costs.
        For the simplified model only.
        """
        for i, layer in enumerate(self.detector_layers):   # ordered in z value from lowest to highest
            x_min = layer[0]
            x_max = layer[1]
            y_min = layer[2]
            y_max = layer[3]

            segment_size_x = (x_max - x_min) / int(self.configuration["binning"]["num bins x"])
            segment_size_y = (y_max - y_min) / int(self.configuration["binning"]["num bins x"])
            for j in range(int(self.configuration["binning"]["num bins x"])):
                for k in range(int(self.configuration["binning"]["num bins y"])):
                    self.segment_list.append(Segment(f"L{i}_S{j}_{k}",  # Layer i segment-x j segment-y k
                                                     i,
                                                     x_min + j * segment_size_x,
                                                     x_min + (j + 1) * segment_size_x,
                                                     y_min + k * segment_size_y,
                                                     y_min + (k + 1) * segment_size_y,
                                                     self.detector_layers[i][-1]))

    def segment_mapping_simplified_model(self):
        """Maps the segments according to the doublet preselection criteria.
        Stores the result inside the segment mapping attribute.
        """

        max_x_detector, max_y_detector, min_x_detector, min_y_detector = self.get_detector_dimensions_LUXE()

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

                max_dy = target_segment.y_end - segment.y_start
                min_dy = min([target_segment.y_end - segment.y_start, 0])

                # max x_0 range on reference screen
                x0_max = self.x0_at_z_ref(target_segment.x_start, segment.x_end, target_segment.z_position,
                                          segment.z_position)

                x0_min = self.x0_at_z_ref(target_segment.x_end, segment.x_start, target_segment.z_position,
                                          segment.z_position)

                # correct for detector dimensions
                if x0_min < min_x_detector:
                    x0_min = min_x_detector
                    if x0_max < min_x_detector:
                        continue
                if x0_max > max_x_detector:
                    x0_max = max_x_detector

                dx_x0_interval = pd.Interval(self.configuration["doublet"]["dx/x0"] -
                                             self.configuration["doublet"]["eps"],
                                             self.configuration["doublet"]["dx/x0"] +
                                             self.configuration["doublet"]["eps"])

                dy_x0_interval = pd.Interval(self.configuration["doublet"]["dy/x0"] -
                                             self.configuration["doublet"]["eps"],
                                             self.configuration["doublet"]["dy/x0"] +
                                             self.configuration["doublet"]["eps"])

                max_dx_interval = pd.Interval(min_dx / x0_max, max_dx / x0_min)
                max_dy_interval = pd.Interval(min_dy / x0_max, max_dy / x0_min)

                if dx_x0_interval.overlaps(max_dx_interval):
                    if dy_x0_interval.overlaps(max_dy_interval):
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
        return x_end - dx * (z_end - self.get_z_reference_layer_LUXE()) / dz

    def target_segments(self,
                        name: str):
        """Takes the name of a segment and the target segments are returned
        :param name: name of segment
        :return:
            target segments
        """
        segment_names = [segment_name for segment_name in self.segment_mapping[name]]
        return [segment for segment in self.segment_list if segment.name in segment_names]

    def get_segment_at_known_xyz_value(self,
                                       x: float,
                                       y: float,
                                       z: float):
        """Finds the correct segment to store a hit. The segment list is a flattened 3d array with
        <number layers> * <num segments in x per layer> * <num segments in y per layer> entries.
        Only possible if all the layers have the exact same length!
        :param x: x value of hit
        :param y: y value of hit
        :param z: z value of hit
        :return:
            index of corresponding segment in segment list
        """
        for i, segment in enumerate(self.segment_list):
            if segment.z_position != z:
                continue
            if segment.x_start > x:
                continue
            if segment.x_end < x:
                continue
            if segment.y_start > y:
                continue
            if segment.y_end < y:
                continue
            return i

    def get_z_reference_layer_LUXE(self):
        """For LUXE the reference layer is always the layer with the lowest distance to the IP,
        thus with the lowest z-value.
        :return
            z-value of reference layer
        """
        return min([layer[-1] for layer in self.detector_layers])

    def get_detector_dimensions_LUXE(self):
        """Extracts detector dimensions from the LUXE setup files.
        :return
            max_x_detector, max_y_detector, min_x_detector, min_y_detector
        """
        min_x = min([layer[1] for layer in self.detector_layers])
        max_x = max([layer[2] for layer in self.detector_layers])
        min_y = min([layer[3] for layer in self.detector_layers])
        max_y = max([layer[4] for layer in self.detector_layers])

        return max_x, max_y, min_x, min_y
