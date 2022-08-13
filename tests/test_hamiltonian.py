import unittest
import sys
sys.path.insert(0, "../src")

from pattern.doublet import Doublet
from pattern.triplet import Triplet
from qubo.hamiltonian import Hamiltonian


class TestHamiltonian(unittest.TestCase):

    def test_linear_interactions(self):
        test_doublet_0 = Doublet(hit_1_particle_key=0,
                                 hit_2_particle_key=0,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=0,
                                 hit_2_id=1,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_1 = Doublet(hit_1_particle_key=1,
                                 hit_2_particle_key=1,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=2,
                                 hit_2_id=3,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_2 = Doublet(hit_1_particle_key=2,
                                 hit_2_particle_key=2,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=4,
                                 hit_2_id=5,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_3 = Doublet(hit_1_particle_key=3,
                                 hit_2_particle_key=3,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=6,
                                 hit_2_id=7,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_4 = Doublet(hit_1_particle_key=4,
                                 hit_2_particle_key=4,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=8,
                                 hit_2_id=9,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_5 = Doublet(hit_1_particle_key=5,
                                 hit_2_particle_key=5,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=10,
                                 hit_2_id=11,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_6 = Doublet(hit_1_particle_key=6,
                                 hit_2_particle_key=6,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=12,
                                 hit_2_id=13,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_7 = Doublet(hit_1_particle_key=7,
                                 hit_2_particle_key=7,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=14,
                                 hit_2_id=15,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_8 = Doublet(hit_1_particle_key=8,
                                 hit_2_particle_key=8,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=16,
                                 hit_2_id=17,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_9 = Doublet(hit_1_particle_key=9,
                                 hit_2_particle_key=9,
                                 hit_1_position=[0, 0, 0],
                                 hit_2_position=[0, 0, 0],
                                 hit_1_id=18,
                                 hit_2_id=19,
                                 energy_1=0,
                                 energy_2=0)
        test_doublet_10 = Doublet(hit_1_particle_key=10,
                                  hit_2_particle_key=10,
                                  hit_1_position=[0, 0, 0],
                                  hit_2_position=[0, 0, 0],
                                  hit_1_id=20,
                                  hit_2_id=21,
                                  energy_1=0,
                                  energy_2=0)
        test_doublet_11 = Doublet(hit_1_particle_key=11,
                                  hit_2_particle_key=11,
                                  hit_1_position=[0, 0, 0],
                                  hit_2_position=[0, 0, 0],
                                  hit_1_id=22,
                                  hit_2_id=23,
                                  energy_1=0,
                                  energy_2=0)
        test_doublet_12 = Doublet(hit_1_particle_key=12,
                                  hit_2_particle_key=12,
                                  hit_1_position=[0, 0, 0],
                                  hit_2_position=[0, 0, 0],
                                  hit_1_id=24,
                                  hit_2_id=25,
                                  energy_1=0,
                                  energy_2=0)
        test_doublet_13 = Doublet(hit_1_particle_key=13,
                                  hit_2_particle_key=13,
                                  hit_1_position=[0, 0, 0],
                                  hit_2_position=[0, 0, 0],
                                  hit_1_id=26,
                                  hit_2_id=27,
                                  energy_1=0,
                                  energy_2=0)
        test_doublet_14 = Doublet(hit_1_particle_key=14,
                                  hit_2_particle_key=14,
                                  hit_1_position=[0, 0, 0],
                                  hit_2_position=[0, 0, 0],
                                  hit_1_id=28,
                                  hit_2_id=29,
                                  energy_1=0,
                                  energy_2=0)
        test_doublet_15 = Doublet(hit_1_particle_key=15,
                                  hit_2_particle_key=15,
                                  hit_1_position=[0, 0, 0],
                                  hit_2_position=[0, 0, 0],
                                  hit_1_id=30,
                                  hit_2_id=31,
                                  energy_1=0,
                                  energy_2=0)

        test_triplet_0 = Triplet(test_doublet_0,
                                 test_doublet_1,
                                 0)
        test_triplet_1 = Triplet(test_doublet_2,
                                 test_doublet_3,
                                 1)
        test_triplet_2 = Triplet(test_doublet_4,
                                 test_doublet_5,
                                 2)
        test_triplet_3 = Triplet(test_doublet_6,
                                 test_doublet_7,
                                 3)
        test_triplet_4 = Triplet(test_doublet_8,
                                 test_doublet_9,
                                 4)
        test_triplet_5 = Triplet(test_doublet_10,
                                 test_doublet_11,
                                 5)
        test_triplet_6 = Triplet(test_doublet_12,
                                 test_doublet_13,
                                 6)
        test_triplet_7 = Triplet(test_doublet_14,
                                 test_doublet_15,
                                 7)

        test_triplet_0.quality = -0.0
        test_triplet_1.quality = -0.1
        test_triplet_2.quality = -0.2
        test_triplet_3.quality = -0.3
        test_triplet_4.quality = -0.4
        test_triplet_5.quality = -0.5
        test_triplet_6.quality = -0.6
        test_triplet_7.quality = -0.7

        test_triplet_0.interactions.update({1: -0.9, 2: -0.8, 3: -0.7, 4: 1, 5: 1})
        test_triplet_1.interactions.update({2: -0.8, 3: -0.7, 4: -0.6, 5: 1, 6: 1})
        test_triplet_2.interactions.update({3: -0.7, 4: -0.6, 5: -0.5, 6: 1, 7: 1})
        test_triplet_3.interactions.update({4: -0.6, 5: -0.5, 6: -0.4, 7: 1, 0: 1})
        test_triplet_4.interactions.update({5: -0.5, 6: -0.4, 7: -0.3, 0: 1, 1: 1})
        test_triplet_5.interactions.update({6: -0.4, 7: -0.3, 0: -0.2, 1: 1, 2: 1})
        test_triplet_6.interactions.update({7: -0.3, 8: -0.2, 1: -0.1, 2: 1, 3: 1})
        test_triplet_7.interactions.update({0: -0.3, 1: -0.2, 2: -0.1, 3: 1, 4: 1})

        test_hamiltonian_0 = Hamiltonian([test_triplet_0,
                                          test_triplet_1,
                                          test_triplet_2,
                                          test_triplet_3,
                                          test_triplet_4,
                                          test_triplet_5,
                                          test_triplet_6,
                                          test_triplet_7],
                                         [0, 1, 0, 1, 0, 1, 0, 1, 0],
                                         rescaling=None)

        self.assertSequenceEqual(test_hamiltonian_0.linear_term().tolist(),
                                 [0.0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7])

        test_hamiltonian_1 = Hamiltonian([test_triplet_0,
                                          test_triplet_1,
                                          test_triplet_3,
                                          test_triplet_7],
                                         [0, 1, 0, 1, 0, 1, 0, 1, 0],
                                         rescaling=None)

        self.assertSequenceEqual(test_hamiltonian_1.linear_term().tolist(), [1.0, 0.9, -0.8, -0.7])


if __name__ == '__main__':
    unittest.main()
