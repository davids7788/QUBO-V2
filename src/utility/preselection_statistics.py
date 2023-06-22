import numpy as np
import matplotlib.pyplot as plt

from pattern.multiplet import Multiplet
from math_functions.geometry import x0_at_z_ref


def plot_dx_x0(z_0, gen_xplet_file):
    dx_x0 = []
    dy_x0 = []
    gen_multiplet = np.load(gen_xplet_file, allow_pickle=True)
    for multiplet in gen_multiplet:
        for i in range(len(multiplet.hit_ids) - 1):
            if multiplet.z[i + 1] == multiplet.z[i]:
                continue
            dx = multiplet.x[i + 1] - multiplet.x[i]
            dy = multiplet.y[i + 1] - multiplet.y[i]
            x0 = x0_at_z_ref(multiplet.x[i + 1],
                             multiplet.x[i],
                             multiplet.z[i + 1],
                             multiplet.z[i],
                             z_0)
            dx_x0_v = dx / x0
            if 0.05 < dx_x0_v < 0.55:
                dx_x0.append(dx_x0_v)
            dy_x0.append(dy / x0)
    print(f'\ndx/x0 from truth: {np.around(np.mean(dx_x0),4)} +/- {np.around(np.std(dx_x0),4)}')
    print(f'dy/y0 from truth: {np.around(np.mean(dy_x0), 4)} +/- {np.around(np.std(dy_x0), 4)}\n')

    # plt.hist(dy_x0, bins=50, range=(-0.002, 0.002))
    # plt.show()



