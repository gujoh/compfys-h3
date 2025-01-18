import numpy as np 
import matplotlib.pyplot as plt 
plt.rcParams['figure.dpi'] = 200

data = np.genfromtxt("data/task5.csv", delimiter=",")
x1 = np.linspace(0.01, 0.4, 20)
x = np.linspace(0, 0.4, 20)

slope, intercept = np.polyfit(x, data[:, 0], 1)
fit_line = slope * x + intercept
plt.scatter(x1, data[:, 0], label="Data")
plt.plot(x, fit_line, color='red', label=f"Fit: y = {slope:.3f}x {intercept:.3f}")
plt.scatter(x[0], fit_line[0], label=f"Ground state estimate: {fit_line[0]:.3f}")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()

coef = np.polyfit(x, data[:, 1], 2)
quadratic_fit = coef[0] * x ** 2 + coef[1] * x + coef[2]
plt.scatter(x1, data[:, 1], label="Data")
plt.scatter(x[0], quadratic_fit[0], label=f"Ground state estimate: {quadratic_fit[0]:.3f}")
plt.plot(x, quadratic_fit, color='red', label=f"Fit: $y = {coef[0]:.3f}x^2 {coef[1]:.3f}x {coef[2]:.3f}$")
plt.plot(x, np.zeros(20) - 2.903, "--", color="black", label="$E_0$ = 2.903")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()


