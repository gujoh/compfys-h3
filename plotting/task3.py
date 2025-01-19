# Imports
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200

data = np.genfromtxt("data/task3.csv", delimiter=",")
dt = 0.1
t = np.arange(data.shape[0]) * dt

plt.plot(t, data[:, 0])
mean_walkers = np.mean(data[2000:, 0])
print(f'Standard deviation, number of walkers: {np.std(data[2000:,0])}')
plt.plot(t, np.zeros(data.shape[0]) + mean_walkers, "--", label=f"Mean number of walkers = {mean_walkers:.1f}")
plt.legend()
plt.ylabel("Number of walkers")
plt.xlabel(r"Time ($\tau$)")
plt.show()

plt.plot(t, data[:, 1], label="$E_T$")
mean_e = np.mean(data[2000:, 1])
print(f'Standard deviation, reference energy: {np.std(data[2000:,1])}')
plt.plot(t, np.zeros(data.shape[0]) + mean_e, "--", label=f"Mean $E_T$ = {mean_e:.3f}")
plt.ylabel("$E_H$")
plt.xlabel(r"Time ($\tau$)")
plt.legend()
plt.show()

