import matplotlib.pyplot as plt
import numpy as np

# Initialize variables
deltaT = 0.001
steps = int(1 / deltaT)
t_values = np.arange(0, 1, deltaT)
y_numerical = np.zeros(steps)
y_analytical = np.zeros(steps)

# Calculate y_numerical and y_analytical values
y = 0
for i in range(steps):
    t = i * deltaT
    y += deltaT * -t**2
    y_numerical[i] = y
    y_analytical[i] = -1/3 * t**3
    
error = abs(y_analytical[steps-1] - y_numerical[steps-1])
print(error)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(t_values, y_numerical, label='Numerical Integration (y_numerical)')
plt.plot(t_values, y_analytical, label='Analytical Solution (y_analytical)', linestyle='--')
plt.xlabel('t')
plt.ylabel('Value')
plt.title('Comparison of Numerical Integration and Analytical Solution')
plt.legend()
plt.grid(True)
plt.show()