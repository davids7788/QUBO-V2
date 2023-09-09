import matplotlib.pyplot as plt
import numpy as np

from pattern_building.qubo_coefficients import QuboCoefficients


def plot_coefficients_statistics(num_particles: float,
                                 qubo_coefficients: QuboCoefficients) -> None:

    """This functions plots and saves various statistics to the same folder where the triplet list is saved to.
    The results are saved as 'triplet_coefficients_statistics.pdf' and 'triplet_interactions.pdf' into the
    target folder.
    :param num_particles: number of particles in the current tracking file
    :param qubo_coefficients: QuboCoefficients object
    """
    print("Plot statistics ...")
    quality_wrong_match_list = []
    connectivity_wrong_match_list = []
    quality_correct_match_list = []
    connectivity_correct_match_list = []
    num_conflicts = 0

    qubo_coefficients.triplet_list = list(qubo_coefficients.triplet_list)
    t_mapping = {key: value for key, value in zip([t_n.triplet_id for t_n in qubo_coefficients.triplet_list],
                                                  np.arange(len(qubo_coefficients.triplet_list)))}

    for i, t1 in enumerate(qubo_coefficients.triplet_list):
        if t1.from_same_particle():
            quality_correct_match_list.append(t1.quality)
        else:
            quality_wrong_match_list.append(t1.quality)

        for i_key, i_value in t1.interactions.items():
            if t_mapping[i_key] > i:
                continue
            t2 = qubo_coefficients.triplet_list[t_mapping[i_key]]
            if t1.from_same_particle() and t2.is_correct_match():  # need to have 2 overlap hits, so this is enough
                connectivity_correct_match_list.append(i_value)
            else:
                if i_value >= 0:
                    num_conflicts += 1
                else:
                    connectivity_wrong_match_list.append(i_value)

    interactions_list = []
    for t in qubo_coefficients.triplet_list:
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
                   f"Number of triplets: {len(qubo_coefficients.triplet_list)}\n"
                   f"Number of interactions: {'%.2e' %(sum(interactions_list) / 2)}")
    plt.yscale("log")
    plt.legend(loc="best", fontsize=20)

    plt.xlabel("Number of interactions with other triplets", fontsize=20)
    plt.ylabel("Fraction of counts", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(f"{qubo_coefficients.save_to_folder}/triplet_interactions.pdf")
    plt.close()

    # distribution of coefficients
    fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 12))

    ax1.hist(quality_correct_match_list,
             bins=50,
             range=(min(min(quality_wrong_match_list),
                        min(quality_correct_match_list)),
                    max(max(quality_wrong_match_list),
                        max(quality_correct_match_list))),
             label="quality correct match",
             histtype='step',
             linewidth=2,
             edgecolor="goldenrod",
             align="left")
    ax1.hist(quality_wrong_match_list,
             range=(min(min(quality_wrong_match_list),
                        min(quality_correct_match_list)),
                    max(max(quality_wrong_match_list),
                        max(quality_correct_match_list))),
             bins=50,
             label="quality wrong match",
             histtype='step',
             linewidth=2,
             edgecolor="royalblue",
             align="left")
    ax1.set_yscale('log')
    ax1.legend(loc='best', fontsize=20)
    ax1.tick_params(axis='both', labelsize=20)
    ax1.set_ylabel("counts", fontsize=20)
    ax1.set_xlabel("[a.u]", fontsize=20)
    ax2.hist(connectivity_correct_match_list,
             range=(min(min(connectivity_correct_match_list),
                        min(connectivity_wrong_match_list)),
                    max(max(connectivity_correct_match_list),
                        max(connectivity_wrong_match_list))),
             bins=50,
             label="interaction correct match",
             histtype="step",
             linewidth=2,
             edgecolor="goldenrod",
             align="left")
    ax2.hist(connectivity_wrong_match_list,
             range=(min(min(connectivity_correct_match_list),
                        min(connectivity_wrong_match_list)),
                    max(max(connectivity_correct_match_list),
                        max(connectivity_wrong_match_list))),
             bins=50,
             label="interaction wrong match",
             histtype="step",
             linewidth=2,
             edgecolor="royalblue",
             align="left")
    ax2.hist([],
             label=f"Number of conflicts: {num_conflicts}",
             range=(min(min(connectivity_correct_match_list),
                        min(connectivity_wrong_match_list)),
                    max(max(connectivity_correct_match_list),
                        max(connectivity_wrong_match_list))),
             color="white")
    ax2.set_yscale('log')
    ax2.legend(loc='best', fontsize=20)
    ax2.set_ylabel("counts", fontsize=20)
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_xlabel("[a.u]", fontsize=20)
    plt.savefig(f"{qubo_coefficients.save_to_folder}/qubo_coefficients_statistics.pdf")
    plt.close()
