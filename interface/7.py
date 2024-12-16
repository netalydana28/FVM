import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Set the figure size
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True

# Random data for the contour plot
data = np.random.randn(800).reshape(10, 10, 8)

# Create a figure and a set of subplots
fig, ax = plt.subplots()

# Method to change the contour data points
def animate(i):
    ax.clear()
    ax.contourf(data[:, :, i], cmap='plasma')

# Call animate method
ani = animation.FuncAnimation(fig, animate, 5, interval=50, blit=False)

# Display the plot
plt.show()