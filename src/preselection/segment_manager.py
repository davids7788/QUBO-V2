import pandas as pd
import csv

from preselection.segment import LUXEDetectorSegment


class SegmentManager:
    def __init__(self,
                 binning: tuple[int, int],
                 detector_geometry: str):
        """Class for handling segments for doublet and triplet creation
        :param binning: [x_bins, y_bins] number of bins in x and y direction
        :param detector_geometry: string with address of a.csv geometry file
        """

        self.binning = binning   # dictionary of selection criteria
        self.detector_chips = []   # list containing coordinate information about detector layers
        self.setup = None

        # load detector layers / chips information and add them to the detector layers list
        with open(detector_geometry, 'r') as file:
            csv_reader = csv.reader(file)
            next(csv_reader)    # skip header, csv files should consist of one line
            for row in csv_reader:
                self.detector_chips.append([float(r) for r in (row[1:])])

        print(f"Using geometry file {detector_geometry.split('/')[-1]} for segmentation algorithm\n")

        # z-position -> layer number, sets are automatically ordered in python -> values ordered from lowest to highest
        self.z_position_to_layer = list(set([chip[4] for chip in self.detector_chips]))

        # check setup, corresponds to LUXE_sl.csv = simplified, LUXE_fl.csv = full
        if len(self.z_position_to_layer) == 4:
            self.setup = "simplified"
        if len(self.z_position_to_layer) == 8:
            self.setup = "full"
        else:
            print("No valid LUXE setup was chosen!")

        self.segment_mapping = {}   # <segment> (key): [<segment_name_0>, <segment_name_1>, ...] (value)
        self.segment_storage = {}   # dictionary for storing and organising segment objects

    def create_LUXE_segments(self):
        """Segments are created according to their x, y and z coordinates. The name of the segments gives
        information about their position and layer.
        """
        for chip in self.detector_chips:   # ordered in z value from lowest to highest
            x_min = chip[0]
            x_max = chip[1]
            y_min = chip[2]
            y_max = chip[3]

            layer_number = self.z_position_to_layer.index(chip[4])

            # creating dictionary key for each layer, value is a list
            self.segment_storage.update({layer_number: []})

            segment_size_x = (x_max - x_min) / int(self.binning[0])
            segment_size_y = (y_max - y_min) / int(self.binning[1])
            for j in range(int(self.binning[0])):
                for k in range(int(self.binning[1])):
                    # name of the segment consists of layer number, segment numbering in x and in y
                    new_segment = LUXEDetectorSegment(f"L{layer_number}_SX{j}_SY{k}",
                                                      layer_number,
                                                      x_min + j * segment_size_x,
                                                      x_min + (j + 1) * segment_size_x,
                                                      y_min + k * segment_size_y,
                                                      y_min + (k + 1) * segment_size_y,
                                                      self.detector_chips[layer_number][-1])
                    self.segment_storage[layer_number].insert(len(self.segment_storage[layer_number]), new_segment)

    def segment_mapping_LUXE(self):
        """Maps the segments according to the doublet preselection criteria. That means, that if there are hits inside
        the area, defined by the segment, that should be considered for creating doublets, a connection to the target
        segment is stored inside the segment mapping attribute.
        """
        max_x_detector, max_y_detector, min_x_detector, min_y_detector = self.get_detector_dimensions_LUXE()

        for index, z_position in enumerate(self.z_position_to_layer):
            for segment in self.segment_list:
                target_list = []
                for target_segment in self.segment_list:
                    if self.setup == "final"
                    if target_segment.layer != segment.layer + 1:   # only consecutive layers
                        continue
                    if target_segment.x_end < segment.x_start:   # heavy scattering excluded
                        continue

                    min_dy = min([abs(target_segment.y_end - segment.y_start),
                                  abs(target_segment.y_start - segment.y_end),
                                  abs(target_segment.y_end - segment.y_end),
                                  abs(target_segment.y_start - segment.y_start)])

                    # max and minimum dx values
                    max_dx = target_segment.x_end - segment.x_start
                    min_dx = max([target_segment.x_start - segment.x_end, 0])

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

                    max_dx_interval = pd.Interval(min_dx / x0_max, max_dx / x0_min)

                    if dx_x0_interval.overlaps(max_dx_interval):
                        if min_dy / x0_max < self.configuration["doublet"]["dy/x0"]:
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
        min_x = min([layer[0] for layer in self.detector_layers])
        max_x = max([layer[1] for layer in self.detector_layers])
        min_y = min([layer[2] for layer in self.detector_layers])
        max_y = max([layer[3] for layer in self.detector_layers])

        return max_x, max_y, min_x, min_y
