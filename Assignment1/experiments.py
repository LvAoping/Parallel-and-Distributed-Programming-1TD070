import numpy as np
import matplotlib.pyplot as plt

# Replace these sample execution times with your actual data
processor_counts = np.array([1, 2, 4, 8, 12, 16, 20, 24, 28, 32])
execution_times = np.array([1.037985, 0.587072, 1.161410, 1.345415, 1.514724, 0.303111, 0.313265, 0.308031, 0.268239, 0.294812])

# Calculate the ideal speedup based on the number of processors
ideal_speedup = processor_counts

# Calculate the actual speedup using the execution times
actual_speedup = execution_times[0] / execution_times
print(actual_speedup)

# Plot the actual speedup and ideal speedup on the same graph
plt.plot(processor_counts, actual_speedup, marker='o', label='Actual Speedup')
plt.plot(processor_counts, ideal_speedup, linestyle='--', label='Ideal Speedup')

# Add labels and a legend
plt.xlabel('Number of Processors')
plt.ylabel('Speedup')
plt.legend()
plt.title('Strong Scaling Speedup')

# Display the graph
plt.show()


