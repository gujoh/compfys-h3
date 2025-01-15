# Imports
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200

data = np.genfromtxt("../data/task1.csv", delimiter=",")
positions = np.genfromtxt("../data/positions_task1.csv", delimiter=",")
dt = 0.02
t = np.arange(data.shape[0]) * dt

plt.plot(t, data[:, 0])
mean_walkers = np.mean(data[5000:, 0])
plt.plot(t, np.zeros(data.shape[0]) + mean_walkers, "--", label=f"Mean number of walkers = {mean_walkers:.1f}")
plt.legend()
plt.ylabel("Number of walkers")
plt.xlabel(r"Time ($\tau$)")
plt.show()

plt.plot(t, data[:, 1], label="$E_T$")
mean_e = np.mean(data[5000:, 1])
plt.plot(t, np.zeros(data.shape[0]) + mean_e, "--", label=f"Mean $E_T$ = {mean_e:.3f}")
plt.ylabel("$E_H$")
plt.xlabel(r"Time ($\tau$)")
plt.legend()
plt.show()

x = np.linspace(-5, 20, 1000)
y = np.exp(- np.exp(- x) - x / 2) / np.sqrt(np.pi)
plt.hist(positions[5000:], bins=50, density=True, edgecolor="black", alpha=0.7, label="Experimental result")
plt.plot(x, y, linewidth=3, label="Analytical result")
plt.xlabel("x [$a_0$]")
plt.ylabel(r"$\Psi$(x)")
plt.legend()
plt.show()

print(positions[np.abs(positions) < 0.00000001].shape)
print(np.min(np.abs(positions)))


