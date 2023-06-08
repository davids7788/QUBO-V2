import numpy as np

from pattern.multiplet import Multiplet
from math_functions.geometry import x0_at_z_ref


def plot_dx_x0(z_0, gen_xplet_file):
    dx_x0 = []
    gen_multiplet = np.load(gen_xplet_file, allow_pickle=True)
    for multiplet in gen_multiplet:
        for i in range(len(multiplet.hit_ids) - 1):
            if multiplet.z[i + 1] == multiplet.z[i]:
                continue
            dx = multiplet.x[i + 1] - multiplet.x[i]
            dx_x0.append(dx / x0_at_z_ref(multiplet.x[i + 1],
                                          multiplet.x[i],
                                          multiplet.z[i + 1],
                                          multiplet.z[i],
                                          z_0))
    print(np.mean(dx_x0), np.std(dx_x0))


