import matplotlib.pyplot as plt
from matplotlib.legend import Legend

xi = [4, 5, 7]
xi_v2 = [7]

c1_random_guess_impact_eff = [0.701, 0.472, 0.094]
c1_random_guess_impact_frate = [0.027, 0.082, 0.198]
c1_random_guess_impact_eff_error = [0.039, 0.012, 0.07]
c1_random_guess_impact_frate_error = [0.005,0.006, 0.08]

c2_random_guess_impact_eff = [0.913, 0.743, 0.294]
c2_random_guess_impact_frate = [0.02, 0.082, 0.399]
c2_random_guess_impact_eff_error = [0.007, 0.008, 0.07]
c2_random_guess_impact_frate_error = [0.002, 0.003, 0.08]

c3_random_guess_impact_eff = [0.942, 0.845, 0.449]
c3_random_guess_impact_frate = [0.017, 0.064, 0.39]
c3_random_guess_impact_eff_error = [0.006, 0.005, 0.013]
c3_random_guess_impact_frate_error = [0.004, 0.005, 0.007]

c3_ones_impact_reversed_eff = [0.98, 0.927, 0.657]
c3_ones_impact_reversed_frate = [0.006, 0.033, 0.248]
c3_ones_impact_reversed_eff_error = [0.003, 0.003, 0.004]
c3_ones_impact_reversed_frate_error = [0.002, 0.002, 0.004]

c3_v2_eff = [0.654]
c3_v2_frate = [0.238]
c3_v2_eff_error = [0.004]
c3_v2_frate_error = [0.004]

c4_random_guess_impact_eff = [0.939, 0.832, 0.449]
c4_random_guess_impact_frate = [0.017, 0.072, 0.398]
c4_random_guess_impact_eff_error = [0.005, 0.005, 0.009]
c4_random_guess_impact_frate_error = [0.003, 0.003, 0.006]

c4_ones_impact_reversed_eff = [0.98, 0.93, 0.667]
c4_ones_impact_reversed_frate = [0.005, 0.03, 0.235]
c4_ones_impact_reversed_eff_error = [0.003, 0.003, 0.004]
c4_ones_impact_reversed_frate_error = [0.001, 0.002, 0.005]

c4_v2_eff = [0.66]
c4_v2_frate = [0.22]
c4_v2_eff_error = [0.004]
c4_v2_frate_error = [0.005]

c5_random_guess_impact_eff = [0.913, 0.859, 0.57]
c5_random_guess_impact_frate = [0.008, 0.033, 0.264]
c5_random_guess_impact_eff_error = [0.004, 0.004, 0.023]
c5_random_guess_impact_frate_error = [0.003, 0.003, 0.009]

c5_ones_impact_reversed_eff = [0.923, 0.883, 0.668]
c5_ones_impact_reversed_frate = [0.004, 0.019, 0.16]
c5_ones_impact_reversed_eff_error = [0.003, 0.004, 0.004]
c5_ones_impact_reversed_frate_error = [0.002, 0.002, 0.003]

c5_v2_eff = [0.652]
c5_v2_frate = [0.152]
c5_v2_eff_error = [0.003]
c5_v2_frate_error = [0.003]

plt.figure(figsize=(8, 6), dpi=200)
# plt.errorbar(xi, c1_random_guess_impact_eff, yerr=c1_random_guess_impact_eff_error,
#              label="baseline", color="grey", marker="o", markersize=8.5, linewidth=3.5)
# plt.errorbar(xi, c2_random_guess_impact_eff, yerr=c2_random_guess_impact_eff_error,
#              label="small", color="navy", alpha=0.25, marker="o", markersize=8.5, linewidth=3.5)
# plt.errorbar(xi, c3_random_guess_impact_eff, yerr=c3_random_guess_impact_eff_error,
#              label="medium", color="navy", alpha=0.33, marker="o", markersize=8.5, linewidth=3.5, linestyle="--")
# plt.errorbar(xi, c4_random_guess_impact_eff, yerr=c4_random_guess_impact_eff_error,
#              label="large", color="navy", alpha=0.66, marker="o", markersize=8.5, linewidth=3.5, linestyle="--")
# plt.errorbar(xi, c5_random_guess_impact_eff, yerr=c5_random_guess_impact_eff_error,
#              label="huge", color="navy", alpha=0.99, marker="o", markersize=8.5, linewidth=3.5, linestyle="--")
plt.errorbar(xi, c3_ones_impact_reversed_eff, yerr=c3_ones_impact_reversed_eff_error,
             label="medium", color="seagreen", alpha=0.33, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi, c4_ones_impact_reversed_eff, yerr=c4_ones_impact_reversed_eff_error,
             label="large", color="seagreen", alpha=0.66, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi, c5_ones_impact_reversed_eff, yerr=c5_ones_impact_reversed_eff_error,
             label="huge", color="seagreen", alpha=0.99, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi_v2, c3_v2_eff, yerr=c3_v2_eff_error,
             label="medium", color="crimson", alpha=0.33, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi_v2, c4_v2_eff, yerr=c4_v2_eff_error,
             label="large", color="crimson", alpha=0.66, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi_v2, c5_v2_eff, yerr=c5_v2_eff_error,
             label="huge", color="crimson", alpha=0.66, marker="o", markersize=8.5, linewidth=3.5)
# plt.text(3.985, 0.3, "[random bin. vector]\nlow $\longrightarrow$ high impact",
#          fontsize=12, bbox=dict(facecolor="navy", alpha=0.5))
# plt.text(5.025, 0.3, "[1 1 1 1 ...1 1 1 1]\nhigh $\longrightarrow$ low impact",
#          fontsize=12, bbox=dict(facecolor="seagreen", alpha=0.5))
plt.legend(loc="lower left", fontsize=16, ncol=2)
plt.xlabel(r"$\xi$", fontsize=16, loc="right")
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel("Efficiency", fontsize=16, loc="top")
plt.ylim(0, 1.05)
plt.grid()
plt.show()

plt.figure(figsize=(8, 6), dpi=200)
# plt.errorbar(xi, c1_random_guess_impact_frate, yerr=c1_random_guess_impact_frate_error,
#              label="baseline", color="grey", marker="o", markersize=8.5, linewidth=3.5)
# plt.errorbar(xi, c2_random_guess_impact_frate, yerr=c2_random_guess_impact_frate_error,
#              label="small", color="navy", alpha=0.25, marker="o", markersize=8.5, linewidth=3.5)
# plt.errorbar(xi, c3_random_guess_impact_frate, yerr=c3_random_guess_impact_frate_error,
#              label="medium", color="navy", alpha=0.33, marker="o", markersize=8.5, linewidth=3.5, linestyle="--")
# plt.errorbar(xi, c4_random_guess_impact_frate, yerr=c4_random_guess_impact_frate_error,
#              label="large", color="navy", alpha=0.66, marker="o", markersize=8.5, linewidth=3.5, linestyle="--")
# plt.errorbar(xi, c5_random_guess_impact_frate, yerr=c5_random_guess_impact_frate_error,
#              label="huge", color="navy", alpha=0.99, marker="o", markersize=8.5, linewidth=3.5, linestyle="--")
plt.errorbar(xi, c3_ones_impact_reversed_frate, yerr=c3_ones_impact_reversed_frate_error,
             label="medium", color="seagreen", alpha=0.33, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi, c4_ones_impact_reversed_frate, yerr=c4_ones_impact_reversed_frate_error,
             label="large", color="seagreen", alpha=0.66, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi, c5_ones_impact_reversed_frate, yerr=c5_ones_impact_reversed_frate_error,
             label="huge", color="seagreen", alpha=0.99, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi_v2, c3_v2_frate, yerr=c3_v2_frate_error,
             label="medium", color="crimson", alpha=0.33, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi_v2, c4_v2_frate, yerr=c4_v2_frate_error,
             label="large", color="crimson", alpha=0.66, marker="o", markersize=8.5, linewidth=3.5)
plt.errorbar(xi_v2, c5_v2_frate, yerr=c5_v2_frate_error,
             label="huge", color="crimson", alpha=0.99, marker="o", markersize=8.5, linewidth=3.5)
plt.legend(loc="upper left", fontsize=16, ncol=2)
# plt.text(3.985, 0.36, "[random bin. vector]\nlow $\longrightarrow$ high impact",
#          fontsize=12, bbox=dict(facecolor="navy", alpha=0.5))
# plt.text(5.025, 0.36, "[1 1 1 1 ...1 1 1 1]\nhigh $\longrightarrow$ low impact",
#          fontsize=12, bbox=dict(facecolor="seagreen", alpha=0.5))
plt.xlabel(r"$\xi$", fontsize=16, loc="right")
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel("Fake rate", fontsize=16, loc="top")
plt.ylim(0, 0.55)
plt.grid()
plt.show()