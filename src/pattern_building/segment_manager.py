import csv

from math_functions.geometry import x0_at_z_ref
from pattern_building.segment import DetectorSegment


class SegmentManager:
    def __init__(self,
                 configuration,
                 detector_geometry: str):
        """Class for handling segments for doublet and triplet creation for LUXE data.
        :param configuration: pattern building configuration file
        :param detector_geometry: .csv detector layer file geometry file
        """
        self.mapping_criteria = configuration['doublet']
        self.binning = [configuration['binning']['num bins x'], configuration['binning']['num bins y']]
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
        self.z_position_to_layer.sort()
        # check setup, corresponds to LUXE_sl.csv = simplified, LUXE_fl.csv = full

        if len(self.z_position_to_layer) == 4:
            self.setup = "simplified"
        elif len(self.z_position_to_layer) == 8:
            self.setup = "full"
        else:
            print("No valid LUXE setup was chosen!")

        self.segment_mapping = {}   # <segment> (key): [<target_segment_0>, <target_segment_1>, ...] (value)
        self.segment_storage = {}   # organising segment objects, subdicts ordered by construction in the following way:
        # <layer> (key):
        # [<segment_x0_y0>, <segment_x0_y1>, ..., <segment_x0_yn>
        #  <segment_x1_y0>, <segment_x1_y1>, ...,
        #  ...,                                 , <segment_xm_yn>]
        self.layer_ranges = {}

    def create_LUXE_segments(self) -> None:
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
                    new_segment = DetectorSegment(f"L{layer_number}_SX{j}_SY{k}",
                                                  layer_number,
                                                  x_min + j * segment_size_x,
                                                  x_min + (j + 1) * segment_size_x,
                                                  y_min + k * segment_size_y,
                                                  y_min + (k + 1) * segment_size_y,
                                                  self.z_position_to_layer[layer_number])
                    self.segment_storage[layer_number].append(new_segment)
            # fortunately x and y are arranged in a way that the min and max x a nd y values can be accessed easily

            if layer_number not in self.layer_ranges.keys():
                self.layer_ranges.update({layer_number: [[x_min,
                                                          x_max,
                                                          y_min,
                                                          y_max]]})
            else:
                self.layer_ranges[layer_number].append([x_min,
                                                        x_max,
                                                        y_min,
                                                        y_max])

    def check_if_full_overlap(self,
                              source_segment: DetectorSegment,
                              target_segment: DetectorSegment) -> bool:
        """Checks if segments come from not consecutive layers e.g. 1 and 3, 4, and 6 if there is a segment in the
        in between layer, which covers the path completely.
        :param source_segment: source segment
        :param target_segment: target segment
        :return
            True if there is an in-between segment covering path of two segments, else False
        """
        dz_segments = target_segment.z_position - source_segment.z_position
        # since segments are the same size, this applies for the end as well
        dx_segments = target_segment.x_start - source_segment.x_start
        source_layer_index = self.z_position_to_layer.index(source_segment.z_position)
        target_layer_index = self.z_position_to_layer.index(target_segment.z_position)
        if target_layer_index == source_layer_index + 2:
            in_between_index = int((target_layer_index + source_layer_index) / 2)
            dz_ratio = (target_segment.z_position - self.z_position_to_layer[in_between_index]) / dz_segments
            x_start_in_between = target_segment.x_start - (dx_segments * dz_ratio)
            x_end_in_between = target_segment.x_end - (dx_segments * dz_ratio)
            for chip in self.layer_ranges[in_between_index]:
                if chip.x_start < x_start_in_between and chip.x_end > x_end_in_between:
                    return True
        return False

    def segment_mapping_LUXE(self) -> None:
        """Maps the segments according to the doublet pattern_building criteria. That means, that if there are hits inside
        the area, defined by the segment, that should be considered for creating doublets, a connection to the target
        segment is stored inside the segment mapping attribute.
        """

        for index, z_position in enumerate(self.z_position_to_layer):
            if self.setup == "full":
                if index > len(self.z_position_to_layer) - 3:
                    continue
            if self.setup == "simplified":
                if index > len(self.z_position_to_layer) - 2:
                    continue
            for segment in self.segment_storage[index]:
                target_list = []
                if self.setup == "full":
                    target_list = self.segment_storage[index + 1] + self.segment_storage[index + 2]
                if self.setup == "simplified":
                    target_list = self.segment_storage[index + 1]
                for target_segment in target_list:
                    # if self.z_position_to_layer.index(target_segment.z_position) - index == 2:
                    #     check_if_overlapping_in_between_layer = self.check_if_full_overlap(segment, target_segment)
                    #     if check_if_overlapping_in_between_layer:
                    #         continue
                    check_compatibility = self.is_compatible_with_target_LUXE_segment(segment, target_segment)
                    if check_compatibility:
                        if segment.name in self.segment_mapping.keys():
                            self.segment_mapping[segment.name].append(target_segment)
                        else:
                            self.segment_mapping.update({segment.name: [target_segment]})


    @staticmethod
    def get_min_dy_of_two_segments(source_y: list[float],
                                   target_y: list[float]) -> float:
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
                                               source_segment: DetectorSegment,
                                               target_segment: DetectorSegment) -> bool:
        """Takes two segments and checks if they are compatible with the pattern_building criteria.
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

        detector_range_x_at_z_ref = [self.segment_storage[0][0].x_start,
                                     self.segment_storage[0][-1].x_end]

        # max x_0 on reference segment, only points on existing detector parts make sense, strictly positive value
        x0_max = x0_at_z_ref(target_segment.x_start,
                             source_segment.x_end,
                             target_segment.z_position,
                             source_segment.z_position,
                             self.z_position_to_layer[0])

        if x0_max < detector_range_x_at_z_ref[0]:
            x0_max = detector_range_x_at_z_ref[0]
        if x0_max > detector_range_x_at_z_ref[1]:
            x0_max = detector_range_x_at_z_ref[1]

        #  exclude heavy scattering in y-direction
        if min_dy / x0_max > self.mapping_criteria["dy/x0 eps"]:
            return False

        # min x_0 on reference segment, only points on existing detector parts make sense, strictly positive value
        x0_min = x0_at_z_ref(target_segment.x_end,
                             source_segment.x_start,
                             target_segment.z_position,
                             source_segment.z_position,
                             self.z_position_to_layer[0])

        if x0_min < detector_range_x_at_z_ref[0]:
            x0_min = detector_range_x_at_z_ref[0]
        if x0_min > detector_range_x_at_z_ref[1]:
            x0_min = detector_range_x_at_z_ref[1]

        max_dx = target_segment.x_end - source_segment.x_start
        min_dx = max([target_segment.x_start - source_segment.x_end, 0])

        if max_dx / x0_min < self.mapping_criteria["dx/x0"] - self.mapping_criteria["dx/x0 eps"]:
            return False

        elif min_dx / x0_max > self.mapping_criteria["dx/x0"] + self.mapping_criteria["dx/x0 eps"]:
            return False
        else:
            return True

    def target_segments(self,
                        name: str) -> list[DetectorSegment]:
        """Takes the name of a segment and the target segments are returned
        :param name: name of segment
        :return:
            target segments
        """
        return self.segment_mapping[name]

    def get_segment_at_known_xyz_value(self,
                                       x: float,
                                       y: float,
                                       z: float) -> DetectorSegment:
        """Finds the correct segment of a hit, given via x,y,z value.
        :param x: x value of hit
        :param y: y value of hit
        :param z: z value of hit
        :return:
            index of corresponding segment in segment list
        """

        print(x, y, z)
        subdict_key = self.z_position_to_layer.index(z)
        print(subdict_key)

        if self.setup == 'full':
            subdict = self.layer_ranges[subdict_key]
            count = 0
            for s_dict in subdict:
                print(s_dict)
                if s_dict[0] <= x <= s_dict[1] and s_dict[2] <= y <= s_dict[3]:
                    print(count)

        exit()

        if self.setup == 'simplified':
            x_index = int((x - self.layer_ranges[subdict][0]) /
                          ((self.layer_ranges[subdict][1] - self.layer_ranges[subdict][0]) / self.binning[0]))

            y_index = int((y - self.layer_ranges[subdict][2]) /
                          ((self.layer_ranges[subdict][3] - self.layer_ranges[subdict][2]) / self.binning[1]))
            print(x_index, y_index)

        return self.segment_storage[subdict][self.binning[1] * x_index + y_index]
