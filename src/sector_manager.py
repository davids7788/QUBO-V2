import pandas as pd
from typing import List


class SectorManager:
    def __init__(self,
                 selection_criteria: dict,
                 detector_layers):
        """Class for handling sectors for doublet and triplet creation
        :param selection_criteria: selection criteria for doublet and triplet selection:
                {doublet: {dx/x0: <value>,
                           eps:   <value>},
                 triplet: {angle diff x: <value>,
                           angle diff y: <value>}}
        :param detector_geometry: list of detector layer dimensions:
               [[l0_x0, l0_y0, l0_x1, l0_y1, l0_z], [l1_x0, l1_y0, l1_x1, l1_y1, l1_z], ...]
        """
        self.selection_criteria = selection_criteria   # Dictionary of selection criteria
        self.detector_layers = detector_layers
        self.sector_list = []   # List of sector objects
        self.sector_mapping = {}   # Map for sectors as <sector>: [<sector_name_0>, <sector_name_1>, ...]
        self.reference_layer_z = None

    def target_sector(self, name):
        """Takes the name of a sector and the target sectors are returned
        :param name: name of sector
        :return: target sectors
        """
        sector_names = [sector_name for sector_name in self.sector_mapping[name]]
        return [sector for sector in self.sector_list if sector.name in sector_names]

    def sector_mapping_simplified_model(self):
        """Maps the sectors according to the doublet preselection criteria.
        Stores the result inside the sector mapping attribute.
        """
        # set reference layer as first layer counting from the IP
        self.reference_layer_z = min([entry[-1] for entry in self.detector_layers])

        # get maximum detector dimension in x
        min_x_detector_dimension = min([entry[0] for entry in self.detector_layers])
        max_x_detector_dimension = max([entry[0] for entry in self.detector_layers])

        for sector in self.sector_list:
            target_list = []
            for target_sector in self.sector_list:
                if target_sector.layer != sector.layer + 1:   # only consecutive layers
                    continue

                # max and minimum dx values
                max_dx = target_sector.x_end - sector.x_start
                min_dx = max([target_sector.x_start - sector.x_end, 0])

                # max x_0 range on reference screen
                x0_max = self.z_at_x0(target_sector.x_start,
                                      sector.x_end,
                                      target_sector.z_position,
                                      sector.z_position)

                x0_min = self.z_at_x0(target_sector.x_end,
                                      sector.x_start,
                                      target_sector.z_position,
                                      sector.z_position)

                # correct for detector dimensions
                if x0_min < min_x_detector_dimension:
                    x0_min = min_x_detector_dimension
                    if x0_max < min_x_detector_dimension:
                        continue
                if x0_max > max_x_detector_dimension:
                    x0_max = max_x_detector_dimension

                dx_x0_interval = pd.Interval(self.selection_criteria["doublet"]["dx/x0"] -
                                             self.selection_criteria["doublet"]["eps"],
                                             self.selection_criteria["doublet"]["dx/x0"] +
                                             self.selection_criteria["doublet"]["eps"])

                max_dx_interval = pd.Interval(min_dx / x0_max, max_dx / x0_min)

                if dx_x0_interval.overlaps(max_dx_interval):
                    target_list.append(target_sector.name)

            self.sector_mapping.update({sector.name: target_list})

    def z_at_x0(self, x_end, x_start, z_end, z_start):
        """
        Help function for calculation x position of doublet at a z-reference value, which was set before
        as a class attribute
        :param x_end: x-position of target sector
        :param x_start: x-position of sector
        :param z_end: z-position of target sector
        :param z_start: z-position of sector
        :return: x_position at the reference layer
        """
        dx = x_end - x_start
        dz = z_end - z_start
        return x_end - dx * abs(z_end - self.reference_layer_z) / dz
