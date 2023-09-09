import unittest
import numpy as np
import sys
sys.path.insert(0, "../src")

from pattern.doublet import Doublet


class TestDoublet(unittest.TestCase):

    def test_is_correct_match(self):
        test_doublet_0 = Doublet(hit_1_particle_key=0,
                                 hit_2_particle_key=0,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[1, 1, 1],
                                 hit_1_id=0,
                                 hit_2_id=1,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_1 = Doublet(hit_1_particle_key=1,
                                 hit_2_particle_key=2,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[1, 1, 1],
                                 hit_1_id=2,
                                 hit_2_id=3,
                                 energy_1=0,
                                 energy_2=0)

        self.assertTrue(test_doublet_0.from_same_particle)
        self.assertFalse(test_doublet_1.from_same_particle)

    def test_xz_angle(self):
        test_doublet_2 = Doublet(hit_1_particle_key=0,
                                 hit_2_particle_key=0,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 1],
                                 hit_1_id=0,
                                 hit_2_id=1,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_3 = Doublet(hit_1_particle_key=1,
                                 hit_2_particle_key=2,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[1, 1, 1],
                                 hit_1_id=2,
                                 hit_2_id=3,
                                 energy_1=0,
                                 energy_2=0)

        self.assertEqual(test_doublet_2.xz_angle(), 0)
        self.assertEqual(test_doublet_3.xz_angle(), np.pi/4)

    def test_yz_angle(self):
        test_doublet_4 = Doublet(hit_1_particle_key=0,
                                 hit_2_particle_key=0,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 1],
                                 hit_1_id=0,
                                 hit_2_id=1,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_5 = Doublet(hit_1_particle_key=1,
                                 hit_2_particle_key=2,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[1, 1, 1],
                                 hit_1_id=2,
                                 hit_2_id=3,
                                 energy_1=0,
                                 energy_2=0)

        self.assertEqual(test_doublet_4.yz_angle(), 0)
        self.assertEqual(test_doublet_5.yz_angle(), np.pi/4)


if __name__ == '__main__':
    unittest.main()
