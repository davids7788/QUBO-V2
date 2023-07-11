import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np
import matplotlib.pyplot as plt
import errno
import os

from scipy.optimize import curve_fit
from scipy.stats import norm
from matplotlib.patches import Rectangle


class Visualisation:
    """Class for visualizing results
    :param source: beam source object used in the experiment
    :param result: file with experimental results
    """
    def __init__(self,
                 source,
                 result):
        self.source = source
        self.result = result

    def beam_source_visualisation(self):
        """Visualisation of the position of the particles, stemming from the beam source. Also shows the energy level
        which is described via color. z direction is the beam axis."""

        if not os.path.isdir('beam_source_visualisation'):
            os.mkdir('beam_source_visualisation')

        position = self.source.get_particle_info('positron', 'position')
        momentum = self.source.get_particle_info('positron', 'momentum')
        positron_weights = self.source.get_particle_info('positron', 'weight')
        
        energy_level = np.array([value[3] for value in momentum])

        x_data = np.array([value[1] for value in position])
        y_data = np.array([value[2] for value in position])
        z_data = np.array([value[3] for value in position])

        fig = plt.figure(dpi=200)
        ax = fig.add_subplot(111, projection='3d')

        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')
        ax.set_title('Point of origin distribution of the beam.')
        scat = ax.scatter3D(x_data,
                            y_data,
                            z_data,
                            c=1e-3 * energy_level,
                            cmap='Greens',
                            label='Momentum [GeV]')
        ax.view_init(30, 205)
        fig.colorbar(scat, ax=ax, shrink=0.4, aspect=3)
        ax.legend(loc='best')
        plt.savefig('beam_source_visualisation/positron_source_position_distribution.pdf')
        plt.close()

        # momentum_distribution
        num_bins = (100, 100, 100)

        px_data = np.array([value[1] for value in momentum])
        py_data = np.array([value[2] for value in momentum])
        pt_data = np.array([np.sqrt(value[1]**2 + value[2]**2) for value in momentum])

        plt.figure(figsize=(10, 10))

        # best fit of data
        (mu_x, sigma_x) = norm.fit(px_data)
        (mu_y, sigma_y) = norm.fit(py_data)
        (mu_pt, sigma_pt) = norm.fit(pt_data)

        # x angle distribution - projected on y axis
        plt.subplot(3, 1, 1)
        plt.title(f'px, py and pt distribution of beam source', fontsize=14)

        plt.hist(px_data,
                 weights=positron_weights,
                 bins=num_bins[0],
                 color='red',
                 linewidth=2.5,
                 histtype='step',
                 alpha=0.75,
                 label=f'$\mu$ = {np.round(mu_x, 4)}'
                       f'\n$\sigma$ = {np.round(sigma_x, 4)}' + '\nnum_bins = ' +
                       str(num_bins[0]))

        plt.xlabel('px [MeV]')
        plt.ylabel('hit counts')
        plt.legend(loc='best', fontsize=12)

        # y angle distribution - projected on y axis
        plt.subplot(3, 1, 2)

        plt.hist(py_data,
                 weights=positron_weights,
                 bins=num_bins[1],
                 color='red',
                 linewidth=2.5,
                 histtype='step',
                 alpha=0.75,
                 label=f'$\mu$ = {np.round(mu_y, 4)}'
                       f'\n$\sigma$ = {np.round(sigma_y, 4)}' + '\nnum_bins = ' +
                       str(num_bins[1]))

        plt.xlabel('py [MeV]')
        plt.ylabel('hit counts')
        plt.legend(loc='best', fontsize=12)

        plt.subplot(3, 1, 3)

        plt.hist(pt_data,
                 weights=positron_weights,
                 bins=num_bins[1],
                 color='red',
                 linewidth=2.5,
                 histtype='step',
                 alpha=0.75,
                 label=f'$\mu$ = {np.round(mu_pt, 4)}'
                       f'\n$\sigma$ = {np.round(sigma_pt, 4)}' + '\nnum_bins = ' +
                       str(num_bins[2]))

        plt.xlabel('pt [MeV]')
        plt.ylabel('hit counts')
        plt.legend(loc='best', fontsize=12)
        plt.savefig('beam_source_visualisation/px_py_pt_beam_source.pdf')
        plt.close()

    def plot_experiment_with_results(self,
                                     x_lim,
                                     y_lim,
                                     z_lim,
                                     number_of_displayed_tracks=20):
        """Plotting a visualisation of the experiment.
        :param x_lim: plot limits x
        :param y_lim: plot limits y
        :param z_lim: plot limits z
        :param number_of_displayed_tracks: number of single tracks from one plane to another shown in the plot
        """
        if not os.path.isdir('tracks_example'):
            os.mkdir('tracks_example')
            
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(projection='3d')

        ax.set_xlim(z_lim[0], z_lim[1])  # for better visualisation switch x and z axis
        ax.set_ylim(y_lim[0], y_lim[1])
        ax.set_zlim(x_lim[0], x_lim[1])

        for plane in self.result[()]['Plane list'].values():
            # plot plane
            p_x = plane.limits_x
            p_y = plane.limits_y
            p = Rectangle((p_x[0], p_y[0]), (p_x[-1] - p_x[0]), (p_y[-1] - p_y[0]), color='grey', alpha=0.3)
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=plane.z_position, zdir="x")

            # plot entries
            for value in enumerate(plane.true_hits_dictionary.values()):
                x_hit = value[1][0]
                y_hit = value[1][1]

                if any([x_hit >= max(plane.limits_x),
                        y_hit >= max(plane.limits_y),
                        x_hit <= min(plane.limits_x),
                        y_hit <= min(plane.limits_y)]):
                    continue

                pixel = plane.bin_edges_for_hit([x_hit, y_hit, plane.z_position])

                dx = pixel[0][1] - pixel[0][0]
                dy = pixel[1][1] - pixel[1][0]

                hit = Rectangle((pixel[0][0], pixel[1][0]), dx, dy, color='red')
                ax.add_patch(hit)
                art3d.pathpatch_2d_to_3d(hit, z=plane.z_position, zdir="x")

        for i in range(number_of_displayed_tracks):
            list_of_hits = []
            for plane in self.result[()]['Plane list'].values():
                if str(i) in plane.true_hits_dictionary.keys():
                    list_of_hits.append(plane.true_hits_dictionary[str(i)])

            if len(list_of_hits) > 1:
                for start, end in zip(list_of_hits[0:-1], list_of_hits[1:]):
                    ax.plot([start[0],
                             end[0]],
                            [start[1],
                             end[1]],
                            [start[2],
                             end[2]],
                            zdir='x',
                            color='blue',
                            linewidth=1)
        plt.savefig('tracks_example/plotted_example_tracks.pdf')
        plt.close()
