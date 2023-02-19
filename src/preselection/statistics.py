def plot_and_save_statistics(self,
                             num_particles,
                             preselection_statistic_dx_x0,
                             preselection_statistic_scattering):
    """This functions plots and saves various statistics to the same folder where the triplet list is saved to.
    The results are saved as 'triplet_coefficients_statistics.pdf' and 'triplet_interactions.pdf' into the
    target folder.
    :param num_particles: number of particles in the current tracking file
    :param preselection_statistic_dx_x0: doublet preselection statistic
    :param preselection_statistic_scattering: scattering / triplet angle statistics
    """
    # Number of interactions with other triplets
    interactions_list = []
    for t in self.triplet_list:
        interactions_list.append(len(list(t.interactions.keys())))
    n, bins, _ = plt.hist(interactions_list, bins=max(interactions_list))
    plt.figure(figsize=(12, 9))
    plt.hist(np.array(interactions_list),
             weights=1 / sum(n) * np.ones(len(interactions_list)),
             bins=max(interactions_list),
             edgecolor="firebrick",
             linewidth=3,
             histtype='step',
             label=f"Number of particles: {num_particles}\n"
                   f"Number of triplets: {len(self.triplet_list)}\n")
    plt.yscale("log")
    plt.legend(loc="best", fontsize=20)

    plt.xlabel("Number of interactions with other triplets", fontsize=20)
    plt.ylabel("Fraction of counts", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(f"{self.save_to_folder}/triplet_interactions.pdf")
    plt.close()

    # distribution of coefficients
    fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 12))

    ax1.hist([self.quality_correct_match_list, self.quality_wrong_match_list],
             bins=50,
             label=[r"quality correct match", r"quality wrong match"],
             edgecolor='k',
             color=["goldenrod", "royalblue"],
             align="left")
    ax1.set_yscale('log')
    ax1.legend(loc='best', fontsize=20)
    ax1.tick_params(axis='both', labelsize=20)
    ax1.set_ylabel("counts", fontsize=20)
    ax1.set_xlabel("[a.u]", fontsize=20)
    n2, bins2, patches2 = ax2.hist([self.connectivity_correct_match_list, self.connectivity_wrong_match_list],
                                   bins=50,
                                   label=[r"interaction correct match",
                                          r"interaction wrong match"],
                                   edgecolor='k',
                                   color=["goldenrod", "royalblue"],
                                   align="left",
                                   rwidth=1)

    width = patches2[1][-1].get_width()
    height = patches2[1][-1].get_height()

    patches2[1][-1].remove()
    ax2.bar(1, align="center", height=height, width=width, color='red', edgecolor="k", label="conflicts")

    ax2.set_yscale('log')
    ax2.legend(loc='best', fontsize=20)
    ax2.set_ylabel("counts", fontsize=20)
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_xlabel("[a.u]", fontsize=20)
    plt.savefig(f"{self.save_to_folder}/qubo_coefficients_statistics.pdf")
    plt.close()

    # preselection statistics of truth doublets / triplets
    fig, (ax3, ax4) = plt.subplots(2, figsize=(12, 12))
    ax3.set_title(r"dx/$x_0$ truth doublets", fontsize=20, loc="left")
    ax3.hist(np.array(preselection_statistic_dx_x0) * 1e2,
             bins=50,
             label=f"$\mu$ = {np.around(np.mean(preselection_statistic_dx_x0), 5)}\n"
                   f"$\sigma$ = {np.around(np.std(preselection_statistic_dx_x0), 5)}",
             edgecolor="blue",
             linewidth=3,
             histtype='step',
             align="left")
    ax3.legend(loc='best', fontsize=20)
    ax3.tick_params(axis='both', labelsize=20)
    ax3.set_ylabel("counts", fontsize=20)
    ax3.set_xlabel("[$10^{-2}$ a.u]", fontsize=20, loc="right")

    ax4.set_title("Scattering angle truth triplets", fontsize=20, loc="left")
    ax4.hist(np.array(preselection_statistic_scattering) * 1e3,
             range=(min(preselection_statistic_scattering) * 1e3,
                    max(preselection_statistic_scattering) * 1e3),
             bins=50,
             label=r"$\sqrt{angle_{xz}^2 + angle_{yz}^2}$:"
                   f"\n$\mu$ = {np.around(np.mean(preselection_statistic_scattering), 5)}\n"
                   f"$\sigma$ = {np.around(np.std(preselection_statistic_scattering), 5)}\n",
             edgecolor="goldenrod",
             linewidth=3,
             histtype='step',
             align="left")

    ax4.legend(loc='best', fontsize=20)
    ax4.set_ylabel("counts", fontsize=20)
    ax4.tick_params(axis='both', labelsize=20)
    ax4.set_xlabel("[$10^{-3}$ rad]", fontsize=20, loc="right")
    plt.savefig(f"{self.save_to_folder}/preselection_truth_statistics.pdf")
