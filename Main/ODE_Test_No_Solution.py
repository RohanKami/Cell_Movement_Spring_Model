import matplotlib.pyplot as plt
import numpy as np

# Initialize variables for the first timestep
deltaT = 0.00001
steps = int(1 / deltaT)
t_values = np.arange(0, 1, deltaT)
y_numerical = np.zeros(steps)

# Calculate y_numerical values for the first timestep using vectorized operations
integrand = np.exp(-t_values**2)
y_numerical = np.cumsum(integrand) * deltaT

# Initialize variables for the second timestep
deltaT2 = 0.0001
steps2 = int(1 / deltaT2)
t_values2 = np.arange(0, 1, deltaT2)
y_numerical2 = np.zeros(steps2)

# Calculate y_numerical values for the second timestep using vectorized operations
integrand2 = np.exp(-t_values2**2)
y_numerical2 = np.cumsum(integrand2) * deltaT2

# Interpolate y_numerical2 to match the t_values of y_numerical
y_numerical2_interpolated = np.interp(t_values, t_values2, y_numerical2)

# Calculate the absolute difference
abs_difference = np.abs(y_numerical - y_numerical2_interpolated)

# Create the plot
plt.figure(figsize=(12, 8))

# Plot y_numerical
plt.plot(t_values, y_numerical, label=f'Numerical Integration (deltaT={deltaT})', linewidth=2)

# Plot y_numerical2
plt.plot(t_values, y_numerical2_interpolated, label=f'Numerical Integration (deltaT={deltaT2})', linestyle='--', linewidth=2)

# Plot the absolute difference
plt.plot(t_values, abs_difference, label='Absolute Difference', linestyle='-.', linewidth=2)

# Add labels and title
plt.xlabel('t')
plt.ylabel('Value')
plt.title('Comparison of Numerical Integrations and Their Difference')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
