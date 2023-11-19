import csv

from math_functions.geometry import x0_at_z_ref
from pattern_building.segment import DetectorSegment
from pattern.detector_hit import DetectorHit


class LUXESegmentManager:
    """Managing segmentation algorithm for the given detector geometry.
    """
    def __init__(self,
                 configuration,
                 detector_geometry: str):
        """Setting fields according to the configuration and detector geometry.
        :param configuration: pattern building configuration file
        :param detector_geometry: .csv detector layer file geometry file
        """
        self.mapping_criteria = configuration['doublet']
        print(f'Segment mapping based on doublet preselection criteria:\n'
              f'dx/x0: {self.mapping_criteria["dx/x0"]}\n'
              f'dx/x0 eps: {self.mapping_criteria["dx/x0 eps"]}\n'
              f'dy/x0: {self.mapping_criteria["dy/x0"]}\n'
              f'dy/x0 eps : {self.mapping_criteria["dy/x0 eps"]}\n')

        self.binning = [configuration['binning']['num bins x'],
                        configuration['binning']['num bins y']]
        print(f'Segmentation binning:\n'
              f'num bins x: {self.binning[0]}\n'
              f'num bins y: {self.binning[1]}\n')

        # list containing coordinate information about detector layers
        self.detector_chips = []
        self.setup = None

        # load detector layers / chips information and add them to the detector layers list
        with open(detector_geometry, 'r') as file:
            csv_reader = csv.reader(file)
            next(csv_reader)    # skip header, csv files should consist of one line
            for row in csv_reader:
                self.detector_chips.append([float(r) for r in (row[1:])])

        print(f"Using geometry file {detector_geometry.split('/')[-1]} for segmentation algorithm")

        # here the first layer is the one with the lowest z position, in the usual geometry file, this is different!
        self.z_position_to_layer = list(set([chip[4] for chip in self.detector_chips]))
        self.z_position_to_layer.sort()

        # check setup
        if len(self.z_position_to_layer) == 4:
            self.setup = "simplified"
            self.chips_per_layer = 1
            print(f'--> simplified geometry setup...\n')
        elif len(self.z_position_to_layer) == 8:
            self.setup = "full"
            self.chips_per_layer = 9
            print(f'--> full geometry setup...\n')
        else:
            print("--> on valid LUXE setup set!")
            exit()

        # <segment> (key): [<target_segment_0>, <target_segment_1>, ...] (value)
        self.segment_mapping = {}
        self.segment_storage = {}

        # organising segment objects, subdicts ordered by construction in the following way:
        # <layer> (key):
        # [<segment_x0_y0>, <segment_x0_y1>, ..., <segment_x0_yn>
        #  <segment_x1_y0>, <segment_x1_y1>, ...,
        #  ...,                                 , <segment_xm_yn>]
        self.layer_ranges = {}

    def create_LUXE_segments(self) -> None:
        """Create segments according to their x, y and z coordinates. The name of the segments gives
        information about their position and layer.
        """
        for layer_number, stave in enumerate(self.z_position_to_layer):
            x_vals = [chip[0] for chip in self.detector_chips if chip[4] == stave] + \
                     [chip[1] for chip in self.detector_chips if chip[4] == stave]
            y_vals = [chip[2] for chip in self.detector_chips if chip[4] == stave] + \
                     [chip[3] for chip in self.detector_chips if chip[4] == stave]

            # extract spatial covering in x,y of the chips in the detector
            min_x, max_x, min_y, max_y = min(x_vals), max(x_vals), min(y_vals), max(y_vals)
            self.layer_ranges.update({layer_number: [min_x, max_x, min_y, max_y]})

            # creating dictionary key for each layer, value is a list
            if layer_number not in self.segment_storage.keys():
                self.segment_storage.update({layer_number: []})
            segment_size_x = (max_x - min_x) / int(self.binning[0])
            segment_size_y = (max_y - min_y) / int(self.binning[1])

            # segment creation
            for j in range(int(self.binning[0])):
                for k in range(int(self.binning[1])):

                    # name of the segment consists of layer number, segment numbering in x and in y
                    new_segment = DetectorSegment(f"L{layer_number}_SX{j}_SY{k}",
                                                  layer_number,
                                                  min_x + j * segment_size_x,
                                                  min_x + (j + 1) * segment_size_x,
                                                  min_y + k * segment_size_y,
                                                  min_y + (k + 1) * segment_size_y,
                                                  self.z_position_to_layer[layer_number] - 0.1e-4,  # 100 mu
                                                  self.z_position_to_layer[layer_number] + 0.1e-4)  # 100 mu

                    self.segment_storage[layer_number].append(new_segment)

    def segment_mapping_LUXE(self) -> None:
        """Maps the segments according to the doublet pattern_building criteria. That means, that if there are hits
        inside the area, defined by the segment, that should be considered for creating doublets, a connection to the
        target  segment is stored inside the segment mapping attribute.
        """
        for index, z_position in enumerate(self.z_position_to_layer):
            print(f'Finding target segments of segments from layer {index}')
            if index > len(self.z_position_to_layer) - 2:
                continue
            for segment in self.segment_storage[index]:
                target_list = []
                if self.setup == "full":
                    if index == len(self.z_position_to_layer) - 2:
                        target_list = self.segment_storage[index + 1]
                    else:
                        target_list = self.segment_storage[index + 1] + self.segment_storage[index + 2]

                if self.setup == "simplified":
                    target_list = self.segment_storage[index + 1]

                for target_segment in target_list:
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
        :param target_y: y edge coordinates of target segment
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
            True if compatible, else False
        """
        # exclude heavy scattering in x-direction
        if target_segment.x_end < source_segment.x_start:
            return False

        min_dy = LUXESegmentManager.get_min_dy_of_two_segments([source_segment.y_start, source_segment.y_end],
                                                               [target_segment.y_start, target_segment.y_end])

        # max x_0 on reference layer
        x0_max = x0_at_z_ref(target_segment.x_start,
                             source_segment.x_end,
                             target_segment.z_start,
                             source_segment.z_start,
                             self.z_position_to_layer[0])

        #  exclude heavy scattering in y-direction
        if min_dy / x0_max / (target_segment.z_start - source_segment.z_start) > self.mapping_criteria["dy/x0 eps"]:
            return False

        # min x_0 on reference segment, only points on existing detector parts make sense, strictly positive value
        x0_min = x0_at_z_ref(target_segment.x_end,
                             source_segment.x_start,
                             target_segment.z_start,
                             source_segment.z_start,
                             self.z_position_to_layer[0])

        max_dx = target_segment.x_end - source_segment.x_start
        min_dx = max([target_segment.x_start - source_segment.x_end, 0])

        max_dx_x0 = abs(max_dx / x0_min / (target_segment.z_start - source_segment.z_start))
        min_dx_x0 = abs(min_dx / x0_max / (target_segment.z_start - source_segment.z_start))

        if min_dx_x0 > self.mapping_criteria["dx/x0"] + self.mapping_criteria["dx/x0 eps"]:
            return False
        if max_dx_x0 < self.mapping_criteria["dx/x0"] - self.mapping_criteria["dx/x0 eps"]:
            return False

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
                                       hit: DetectorHit) -> DetectorSegment | None:
        """Finds the correct segment of a hit, given via x,y,z value.
        :param hit: DetectorHit object

        :return:
            index of corresponding segment in segment list
        """

        for segment in self.segment_storage[self.z_position_to_layer.index(hit.z)]:
            if segment.x_start <= hit.x <= segment.x_end and segment.y_start <= hit.y <= segment.y_end:
                return segment
        return None
