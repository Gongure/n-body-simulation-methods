
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib.backend_tools import ToolBase
from matplotlib.animation import FuncAnimation
import matplotlib


def animate(data, time_steps, bboxes):
    for x in data:
        for i in range(len(x['position'])):
            x['position'][i] = x['position'][i].value

    PLOT_MAX = 5
    ARROW_LENGTH = 0.5

    plt.ion()
    fig = plt.figure()
    ax = Axes3D(fig)
    ax = fig.add_subplot(111, projection='3d')
    fig.suptitle("BruteForce", fontsize=12)

    colors = [body['color'] for body in data]

    for line in bboxes:
        for bbox in line:
            ooo = bbox[0]
            iii = bbox[1]
            ooi = [ooo[0], ooo[1], iii[2]]
            oio = [ooo[0], iii[1], ooo[2]]
            oii = [ooo[0], iii[1], iii[2]]
            ioo = [iii[0], ooo[1], ooo[2]]
            ioi = [iii[0], ooo[1], iii[2]]
            iio = [iii[0], iii[1], ooo[2]]
    
            
            

    for i in range(time_steps):
        plt.cla()

        ax.set_xlabel('X in AE')
        ax.set_ylabel('Y in AE')
        ax.set_zlabel('Z in AE')

        # Set the limits of the plot
        ax.set_xlim(-PLOT_MAX, PLOT_MAX)
        ax.set_ylim(-PLOT_MAX, PLOT_MAX)
        ax.set_zlim(-PLOT_MAX, PLOT_MAX)

        current_positions_x = []
        current_positions_y = []
        current_positions_z = []

        # umscrheiben: vorher schon unterteilenen und dann nur noch die richtigen nehmen   
        for body in data:
            current_positions_x += [body['position'][i][0]]
            current_positions_y += [body['position'][i][1]]
            current_positions_z += [body['position'][i][2]]

        ax.scatter(
            xs=current_positions_x, ys=current_positions_y, zs=current_positions_z, c=colors)

        """
        if showName:
            ax.text(current_positions_x, current_positions_y,
                    current_positions_z, names)
        """

        """
        for body in data:

            ax.scatter(
                xs=['position'][i][0], ys=body['position'][i][1], zs=body['position'][i][2])

            if showName:
                ax.text(body['position'][i][0], body['position'][i]
                        [1], body['position'][i][2], body['name'])

            if show_velocityVector:
                ax.quiver(body['position'][i][0], body['position'][i][1], body['position'][i][2], body['velocity'][i][0],
                          body['velocity'][i][1], body['velocity'][i][2], normalize=True, length=ARROW_LENGTH, color="red")

            ax.scatter(xs=PLOT_MAX, ys=PLOT_MAX, zs=PLOT_MAX, alpha=0.0)
            ax.scatter(xs=-PLOT_MAX, ys=-PLOT_MAX, zs=-PLOT_MAX, alpha=0.0)
        """
        plt.pause(0.001)


def visualize_planetary_motionEndPic(data, time_steps):
    for x in data:
        for i in range(len(x['position'])):
            x['position'][i] = x['position'][i].value

    fig = plt.figure()
    ax = plt.axes(projection="3d")

    PLOT_MAX = 5

    ax.set_xlim(-PLOT_MAX, PLOT_MAX)
    ax.set_ylim(-PLOT_MAX, PLOT_MAX)
    ax.set_zlim(-PLOT_MAX, PLOT_MAX)

    for planet in data:
        name = planet["name"]
        positions = planet["position"]

        x = [pos[0] for pos in positions]
        y = [pos[1] for pos in positions]
        z = [pos[2] for pos in positions]

        ax.plot3D(x, y, z, label=name)

    plt.legend()
    plt.show()


def animate_solar_system(data, interval=50):
    for x in data:
        for i in range(len(x['position'])):
            x['position'][i] = x['position'][i].value

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Keep track of the current time step
    time_step = 0
    max_time_step = max([len(planet['position']) for planet in data])

    # Create a scatter plot for each planet
    scatters = []
    for planet in data:
        scatters.append(ax.scatter([], [], [], label=planet['name']))

    # Set up axis labels and limits
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    def update(frame):
        nonlocal time_step
        for i, planet in enumerate(data):
            if time_step >= len(planet['position']):
                continue
            x, y, z = planet['position'][time_step]
            scatters[i].set_offsets(np.c_[x, y, z])
        time_step = (time_step + 1) % max_time_step
        return scatters

    anim = FuncAnimation(fig, update, frames=max_time_step, interval=interval)
    # anim.save('solar_system.gif', writer='imagemagick')
    plt.show()
