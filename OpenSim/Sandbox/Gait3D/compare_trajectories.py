import opensim as osim
import os
import matplotlib.pyplot as plt

build_path = '/Users/nbianco/Developer/opensim/worktrees/gait3d/gait3d/build'   

# Load the states
states1 = osim.TimeSeriesTable(os.path.join(build_path, "states_pointpath_zero_control0.1.sto"))
states2 = osim.TimeSeriesTable(os.path.join(build_path, "states_hyfydy_zero_control0.1_pointbasedpath.sto"))

time1 = states1.getIndependentColumn()
time2 = states2.getIndependentColumn()

labels1 = states1.getColumnLabels()
labels2 = states2.getColumnLabels()

# Plot the first 16 columns of each table, each column on the same subplot
# Use a 4x4 grid for the plots and use the labels as the titles
side = 5
fig, axs = plt.subplots(side, side, figsize=(10, 10))

# Set the titles of the subplots
for i in range(side):
    for j in range(side):
        axs[i, j].set_title(labels1[i*side + j + 1]) 
        axs[i, j].plot(time1, states1.getDependentColumn(labels1[i*side + j + 1]), label='PointPathMuscle')
        axs[i, j].plot(time2, states2.getDependentColumn(labels2[i*side + j + 1]), label='HyfydyMuscle')
        axs[i, j].legend()
plt.show()
