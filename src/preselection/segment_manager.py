import pandas as pd
import csv

from preselection.segment import LUXEDetectorSegment


class SegmentManager:
    def __init__(self,
                 binning: list[int, int],
                 mapping_criteria: dict,
                 detector_geometry: str):
        """Class for handling segments for doublet and triplet creation
        :param configuration: information needed for the segments:
                {
                doublet: {dx/x0: <value>,
                          eps:   <value>,
                          dy/x0: <value>},
                triplet: {angle diff x: <value>,
                          angle diff y: <value>}
                binning: {num bins x: <value>,
                          num bins y: <value>}
                qubo parameters: {b_ij conflict: <value>,
                                  b_ij match: value,
                                  a_i: <value>}
                scale range parameters: {z_scores: <value>,
                                         quality: <value>,
                                         interaction: <value>}
                }
        :param detector_geometry: .csv detector layer file geometry file
        """

        self.mapping_criteria = mapping_criteria
        self.binning = binning
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
        self.detector_range_x = [min([chip[1] for chip in self.detector_chips]),
                                 max([chip[2] for chip in self.detector_chips])]

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
        for index, z_position in enumerate(self.z_position_to_layer):
            for segment in self.segment_storage[index]:
                target_list = []
                if self.setup == "full":
                    target_list = self.segment_storage[index + 1] + self.segment_storage[index + 2]
                if self.setup == "simplified":
                    target_list = self.segment_storage[index + 1]

                for target_segment in self.segment_storage[index + 1] + self.segment_storage[index + 2]:
                    check_compatibility = self.is_compatible_with_target_LUXE_segment(segment, target_segment)
                    if check_compatibility:
                        target_list.append(target_segment.name)

                self.segment_mapping.update({segment.name: target_list})

    @staticmethod
    def get_min_dy_of_two_segments(source_y: list[float, float],
                                   target_y: list[float, float]) -> float:
        """Calculating the minimum difference in y of two segments depending on their location on the detector.
        :param source_y: y edge coordinates of source segment
        :param target_y: y edge coordinates of source segment
        :return
            minimum distance of segments in y
        """
        if source_y[0] == target_y[0] and source_y[1] == target_y[1]:
            return 0.0
        elif source_y[0] > target_y[1]:
            return source_y[0] - target_y[1]
        else:
            return target_y[0] - source_y[1]


    def is_compatible_with_target_LUXE_segment(self,
                                               source_segment: LUXEDetectorSegment,
                                               target_segment: LUXEDetectorSegment) -> bool:
        """Takes two segments and checks if they are compatible with the preselection criteria.
        :param source_segment: segment from which the mapping starts
        :param target_segment: segment considered as a target
        :return:
            True if compatible, else false
        """

        # exclude heavy scattering in x-direction
        if target_segment.x_end < source_segment.x_start:
            return False

        min_dy = SegmentManager.get_min_dy_of_two_segments([source_segment.y_start, source_segment.y_end],
                                                           [target_segment.y_start, target_segment.y_end])


        # max x_0 range on reference segment
        x0_max = self.x0_at_z_ref(target_segment.x_start,
                                  source_segment.x_end,
                                  target_segment.z_position,
                                  source_segment.z_position)

        #  exclude heavy scattering in y-direction
        if min_dy / x0_max > self.mapping_criteria["dy/x0"]:
            return False

        x0_min = self.x0_at_z_ref(target_segment.x_end,
                                  source_segment.x_start,
                                  target_segment.z_position,
                                  source_segment.z_position)





        # max and minimum dx values
        max_dx = target_segment.x_end - segment.x_start
        min_dx = max([target_segment.x_start - segment.x_end, 0])



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
        return None

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
        return x_end - dx * (z_end - self.z_position_to_layer[0]) / dz

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

