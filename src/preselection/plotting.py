import matplotlib.pyplot as plt
import numpy as np

from preselection.qubo_coefficients import QuboCoefficients


def plot_and_save_statistics(num_particles: float,
                             qubo_coefficients: QuboCoefficients) -> None:

    """This functions plots and saves various statistics to the same folder where the triplet list is saved to.
    The results are saved as 'triplet_coefficients_statistics.pdf' and 'triplet_interactions.pdf' into the
    target folder.
    :param num_particles: number of particles in the current tracking file
    :param qubo_coefficients: QuboCoefficients object
    """
    # Number of interactions with other triplets

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
                   f"Number of triplets: {len(qubo_coefficients.triplet_list)}\n")
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

    ax1.hist(qubo_coefficients.quality_correct_match_list,
             bins=50,
             range=(min(min(qubo_coefficients.quality_wrong_match_list),
                        min(qubo_coefficients.quality_correct_match_list)),
                    max(max(qubo_coefficients.quality_wrong_match_list),
                        max(qubo_coefficients.quality_correct_match_list))),
             label="quality correct match",
             histtype='step',
             linewidth=2,
             edgecolor="goldenrod",
             align="left")
    ax1.hist(qubo_coefficients.quality_wrong_match_list,
             range=(min(min(qubo_coefficients.quality_wrong_match_list),
                        min(qubo_coefficients.quality_correct_match_list)),
                    max(max(qubo_coefficients.quality_wrong_match_list),
                        max(qubo_coefficients.quality_correct_match_list))),
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
    ax2.hist(qubo_coefficients.connectivity_correct_match_list,
             range=(min(min(qubo_coefficients.connectivity_correct_match_list),
                        min(qubo_coefficients.connectivity_wrong_match_list)),
                    max(max(qubo_coefficients.connectivity_correct_match_list),
                        max(qubo_coefficients.connectivity_wrong_match_list))),
             bins=50,
             label="interaction correct match",
             histtype="step",
             linewidth=2,
             edgecolor="goldenrod",
             align="left")
    ax2.hist(qubo_coefficients.connectivity_wrong_match_list,
             bins=50,
             label="interaction wrong match",
             histtype="step",
             linewidth=2,
             edgecolor="royalblue",
             align="left")
    ax2.set_yscale('log')
    ax2.legend(loc='best', fontsize=20)
    ax2.set_ylabel("counts", fontsize=20)
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_xlabel("[a.u]", fontsize=20)
    plt.savefig(f"{qubo_coefficients.save_to_folder}/qubo_coefficients_statistics.pdf")
    plt.close()
