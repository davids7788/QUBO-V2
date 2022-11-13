import matplotlib.pyplot as plt
from matplotlib.legend import Legend

xi = [4, 5, 7]


qubo_analytical_eff = [0.9698, 0.8930, 0.5341]
qubo_analytical_eff_error = [0.003, 0.004, 0.005]
qubo_analytical_frate = [0.0014, 0.0096, 0.0619]
qubo_analytical_frate_error = [0.001, 0.002, 0.002]

ACTS_eff = [0.980226, 0.90095, 0.633594]
ACTS_eff_error = [0.00350474, 0.00295252, 0.00187182]
ACTS_frate = [0.000959693, 0.009932, 0.0948692]
ACTS_frate_error = [0.000619849, 0.000997547, 0.00135573]

plt.figure(figsize=(6, 6), dpi=250)

plt.plot(xi, qubo_analytical_eff, label="VQE", color="blue", marker="o", markersize=9.5, linewidth=3.)
plt.plot(xi, ACTS_eff, label="CKF", color="black", marker="v", markersize=9.5, linewidth=3.)

plt.legend(loc="best", fontsize=16)
plt.xlabel(r"$\xi$", fontsize=16, loc="right")
plt.xticks([4, 5, 6, 7], fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel("Efficiency", fontsize=16)
plt.ylim(0.5, 1.0)
plt.grid()
plt.savefig("fig1", bbox_inches='tight')

plt.figure(figsize=(6, 6), dpi=250)

plt.plot(xi, qubo_analytical_frate, label="VQE", color="blue", marker="o", markersize=9.5,  linewidth=3.)
plt.plot(xi, ACTS_frate, label="CKF", color="black", marker="v", markersize=9.5, linewidth=3.)

plt.legend(loc="best", fontsize=16)
plt.xlabel(r"$\xi$", fontsize=16, loc="right")
plt.xticks([4, 5, 6, 7], fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel("Fake rate", fontsize=16)
plt.ylim(0.0, 0.1)
plt.grid()
plt.savefig("fig2", bbox_inches='tight')
