#Imports
import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt
import sys
from matplotlib.animation import FuncAnimation

#Set font size
import matplotlib
matplotlib.rcParams.update({'font.size': 15})

#Load the data cube
cubic = pa.cx_cube()
cubic.load("cube.bin")

#Save slices of data cube to appropriate arrays
cub_arr = np.sqrt(np.real(np.array(pa.conj(cubic)@cubic)))
re_cub_arr = np.real(np.array(cubic))
im_cub_arr = np.imag(np.array(cubic))


#Fill arrays with total probability values
p_arr = np.zeros(len(cub_arr))
i_arr = np.zeros(len(cub_arr))

dt = 2.5e-5
for i in range(len(cub_arr)):
    i_arr[i] = i*dt
    p_arr[i] = np.sum(cub_arr[i]**2)



"""
The following code for animation the probability function over times is directly copied from the project description,
and therefore not authored by me
"""

z_data_list = []
for i in range(len(cub_arr)):
    z_data_list.append(cub_arr[i])

h = 0.005
x_points = np.arange(0, 1, h)
y_points = np.arange(0, 1, h)
x, y = np.meshgrid(x_points, y_points, sparse=True)


# Array of time points
dt = 2.5e-5
t_points = np.arange(0, 0.008+dt, dt)

# Some settings
fontsize = 12
t_min = t_points[0]
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]

norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

# Create figure
fig = plt.figure()
ax = plt.gca()

# Create a colour scale normalization according to the max z value in the first frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

# Plot the first frame
img = ax.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

# Axis labels
plt.xlabel("x", fontsize=fontsize)
plt.ylabel("y", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label("$\sqrt{p(x,y,t)}$", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white",
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
    img.set_norm(norm)

    # Update z data
    img.set_data(z_data_list[i])

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img


#Plots
fig1, ax1 = plt.subplots()
ax1.plot(i_arr,p_arr)
ax1.hlines(1,i_arr[0],i_arr[-1],color='k',linestyles='--')
ax1.set_xlabel('Time [1]')
ax1.set_ylabel('Total probability [1]')
fig1.tight_layout()

from sklearn import preprocessing

fig2, ax2 = plt.subplots()
ax2.plot(y_points[1:-1],preprocessing.normalize([cub_arr[79,:,160]**2])[0])
ax2.set_xlabel("y [1]")
ax2.set_ylabel("Detection probability [1]")
fig2.tight_layout()

fig3,ax3 = plt.subplots()
img3 = ax3.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)
ax3.set_xlabel('x [1]')
ax3.set_ylabel('y [1]')
cbar3 = fig3.colorbar(img3, ax=ax3)
cbar3.set_label("$\sqrt{p(x,y,t=0)}$", fontsize=fontsize)
cbar3.ax.tick_params(labelsize=fontsize)
fig3.tight_layout()


norm4 = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[40]))

fig4,ax4 = plt.subplots()
img4 = ax4.imshow(z_data_list[40], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm4)
ax4.set_xlabel('x [1]')
ax4.set_ylabel('y [1]')
cbar4 = fig4.colorbar(img4, ax=ax4)
cbar4.set_label("$\sqrt{p(x,y,,t=0.001)}$", fontsize=fontsize)
cbar4.ax.tick_params(labelsize=fontsize)
fig4.tight_layout()

norm5 = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[79]))

fig5,ax5 = plt.subplots()
img5 = ax5.imshow(z_data_list[79], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm5)
ax5.set_xlabel('x [1]')
ax5.set_ylabel('y [1]')
cbar5 = fig5.colorbar(img5, ax=ax5)
cbar5.set_label("$\sqrt{p(x,y,t=0.002)}$", fontsize=fontsize)
cbar5.ax.tick_params(labelsize=fontsize)
fig5.tight_layout()

norm6 = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(im_cub_arr[0]))

fig6,ax6 = plt.subplots()
img6 = ax6.imshow(im_cub_arr[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm6)
ax6.set_xlabel('x [1]')
ax6.set_ylabel('y [1]')
cbar6 = fig6.colorbar(img6, ax=ax6)
cbar6.set_label("Im{$u(x,y,t=0)$}", fontsize=fontsize)
cbar6.ax.tick_params(labelsize=fontsize)
fig6.tight_layout()

norm7 = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(im_cub_arr[40]))

fig7,ax7 = plt.subplots()
img7 = ax7.imshow(im_cub_arr[40], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm7)
ax7.set_xlabel('x [1]')
ax7.set_ylabel('y [1]')
cbar7 = fig7.colorbar(img7, ax=ax7)
cbar7.set_label("Im{$u(x,y,t=0.001)$}", fontsize=fontsize)
cbar7.ax.tick_params(labelsize=fontsize)
fig7.tight_layout()

norm8 = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(im_cub_arr[79]))

fig8,ax8 = plt.subplots()
img8 = ax8.imshow(im_cub_arr[79], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm8)
ax8.set_xlabel('x [1]')
ax8.set_ylabel('y [1]')
cbar8 = fig8.colorbar(img7, ax=ax8)
cbar8.set_label("Im{$u(x,y,t=0.002)$}", fontsize=fontsize)
cbar8.ax.tick_params(labelsize=fontsize)
fig8.tight_layout()

norm9 = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(re_cub_arr[0]))

fig9,ax9 = plt.subplots()
img9 = ax9.imshow(re_cub_arr[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm9)
ax9.set_xlabel('x [1]')
ax9.set_ylabel('y [1]')
cbar9 = fig9.colorbar(img9, ax=ax9)
cbar9.set_label("Re{$u(x,y,t=0)$}", fontsize=fontsize)
cbar9.ax.tick_params(labelsize=fontsize)
fig9.tight_layout()

norm10 = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(re_cub_arr[40]))

fig10,ax10 = plt.subplots()
img10 = ax10.imshow(re_cub_arr[40], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm10)
ax10.set_xlabel('x [1]')
ax10.set_ylabel('y [1]')
cbar10 = fig10.colorbar(img10, ax=ax10)
cbar10.set_label("Re{$u(x,y,t=0.001)$}", fontsize=fontsize)
cbar10.ax.tick_params(labelsize=fontsize)
fig10.tight_layout()

norm11 = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(re_cub_arr[79]))

fig11,ax11 = plt.subplots()
img11 = ax11.imshow(re_cub_arr[79], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm11)
ax11.set_xlabel('x [1]')
ax11.set_ylabel('y [1]')
cbar11 = fig11.colorbar(img11, ax=ax11)
cbar11.set_label("Re{$u(x,y,t=0.002)$}", fontsize=fontsize)
cbar11.ax.tick_params(labelsize=fontsize)
fig11.tight_layout()

anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=False, blit=0)
# Use matplotlib.animation.FuncAnimation to put it all together

# Run the animation!
plt.show()

# # Save the animation
anim.save('./animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)
