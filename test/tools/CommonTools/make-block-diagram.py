import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal for <class 'numpy.float64'> type is zero.")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

# Define the control, validation, and signal regions
control_region = {'MET': [0, 70], 'MT': [0, float('inf')]}
validation_region = {'MET': [0, float('inf')], 'MT': [0, 100]}
signal_region = {'MET': [120, float('inf')], 'MT': [100, float('inf')]}

# Create a meshgrid for plotting
MET = np.linspace(0, 200, 400)
MT = np.linspace(0, 200, 400)
MET, MT = np.meshgrid(MET, MT)

# Initialize the region array
region = np.zeros(MET.shape)

# Assign control region
region[(MET < 70)] = 1  # Control region

# Plot the control region
plt.figure(figsize=(10, 8))
plt.contourf(MET, MT, region, levels=[0.5, 1.5], colors='blue', alpha=0.3)

# Assign validation region
region[(MT < 100)] = 2  # Validation region

# Plot the validation region with hatched very light grey lines only where it overlaps with control region
plt.contourf(MET, MT, region, levels=[1.5, 2.5], colors='none', hatches=['//'], extend='neither')

# Assign signal region
region[(MET > 120) & (MT > 100)] = 3  # Signal region

# Plot the signal region with green filled color
plt.contourf(MET, MT, region, levels=[2.5, 3.5], colors='green', alpha=0.3)

# Add individual legends for the regions in the corresponding regions with white background and increased font size
plt.text(35, 150, 'Control Region: \nMET < 70 GeV', fontsize=12, ha='center', color='black', bbox=dict(facecolor='white', alpha=1))
plt.text(100, 50, 'Validation Region: \nMT < 100 GeV', fontsize=14, ha='center', color='black', bbox=dict(facecolor='white', alpha=1))
plt.text(160, 150, 'Signal Region: \nMET > 120 GeV\n MT > 100 GeV', fontsize=12, ha='center', color='black', bbox=dict(facecolor='white', alpha=1))

# Increase the x-axis and y-axis title font size and add arrows to indicate extension beyond 200 GeV
plt.xlabel('MET (GeV) ', fontsize=14)
plt.ylabel('MT (GeV) ', fontsize=14)

# Set x and y limits to show that the axes extend to infinity
plt.xlim(0, 220)
plt.ylim(0, 220)

# Add arrows to indicate extension beyond the plot limits
plt.annotate('', xy=(220, 0), xytext=(200, 0), arrowprops=dict(arrowstyle="->"))
plt.annotate('', xy=(0, 220), xytext=(0, 200), arrowprops=dict(arrowstyle="->"))

# Remove top and right spines to indicate axes extend to infinity
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.grid(False)

# Save the plot as a PDF file
plt.savefig("met_mt_phase_space.pdf")

# Show the plot
plt.show()
