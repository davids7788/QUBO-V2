import matplotlib.pyplot as plt
from matplotlib.legend import Legend

xi = [4, 5, 7]


qubo_analytical_eff = [0.974, 0.906, 0.583]
qubo_analytical_eff_error = [0.003, 0.004, 0.005]
qubo_analytical_frate = [0.003, 0.017, 0.119]
qubo_analytical_frate_error = [0.001, 0.002, 0.002]

ACTS_eff = [0.980226, 0.90095, 0.633594]
ACTS_eff_error = [0.00350474, 0.00295252, 0.00187182]
ACTS_frate = [0.000959693, 0.009932, 0.0948692]
ACTS_frate_error = [0.000619849, 0.000997547, 0.00135573]

plt.figure(figsize=(8, 6), dpi=200)

plt.errorbar(xi, qubo_analytical_eff, yerr=qubo_analytical_eff_error,
             label="VQE (Eigensolver)", color="blue", alpha=0.85, marker="o", markersize=9.5, linewidth=3.)
plt.errorbar(xi, ACTS_eff, yerr=ACTS_eff_error,
             label="CKF", color="black", alpha=0.85, marker="v", markersize=9.5, linewidth=3.)

plt.legend(loc="lower left", fontsize=16, ncol=2)
plt.xlabel(r"$\xi$", fontsize=16, loc="right")
plt.xticks([4,5,6,7], fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel("Efficiency", fontsize=16, loc="top")
plt.ylim(0, 1.05)
plt.grid()
plt.show()

plt.figure(figsize=(8, 6), dpi=200)

plt.errorbar(xi, qubo_analytical_frate, yerr=qubo_analytical_frate_error,
             label="VQE (Eigensolver)", color="blue", alpha=0.85, marker="o", markersize=9.5, linewidth=3.)
plt.errorbar(xi, ACTS_frate, yerr=ACTS_frate_error,
             label="CKF", color="black", alpha=0.85, marker="v", markersize=9.5, linewidth=3.)

plt.legend(loc="upper left", fontsize=16, ncol=2)
plt.xlabel(r"$\xi$", fontsize=16, loc="right")
plt.xticks([4,5,6,7], fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel("Fake rate", fontsize=16, loc="top")
plt.ylim(0.0, 1.05)
plt.grid()
plt.show()
