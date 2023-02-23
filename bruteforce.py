

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
matplotlib.rcParams["toolbar"] = "toolmanager"


# Julian date calculator
# Input is one string in the format DD.MM.YYYY
# Output is a string in the format '2453599.01542'


def julian_date(date):
    date = date.split('.')
    day = int(date[0])
    month = int(date[1])
    year = int(date[2])
    if month == 1 or month == 2:
        month += 12
        year -= 1
    A = int(year / 100)
    B = 2 - A + int(A / 4)
    if year < 0:
        C = int((365.25 * year) - 0.75)
    else:
        C = int(365.25 * year)
    D = int(30.6001 * (month + 1))
    JD = B + C + D + day + 1720994.5
    JD = str(JD)
    JD = JD.split('.')
    JD = JD[0] + '.' + JD[1][:5]
    return JD


def fetch_data(date):

    # Define a list of objects to retrieve data for
    objects = [{'name': 'Sun', 'mass': 1.989e+30 * u.kg, 'id': '10', 'color': 'yellow'},  {'name': 'Mercury', 'mass': 3.3022e+23 * u.kg, 'id': '199', 'color': 'gray'},  {'name': 'Venus', 'mass': 4.8685e24 * u.kg, 'id': '299', 'color': 'yellow'},  {'name': 'Earth', 'mass': 5.97237e24 *
                                                                                                                                                                                                                                                        u.kg, 'id': '399', 'color': 'blue'},  {'name': 'Mars', 'mass': 6.4185e23 * u.kg, 'id': '499', 'color': 'red'},               {'name': 'Jupiter', 'mass': 1.8986e+27 * u.kg, 'id': '599', 'color': 'orange'},  {'name': 'Saturn', 'mass': 5.6846e+26 * u.kg, 'id': '699', 'color': 'yellow'}]

    data = []

    # Iterate over the objects
    for obj in objects:
        # Query the JPL Horizons database using Astroquery

        # id_type='majorbody' ; get_raw_response=True

        #epochs = {'start': '2010-01-01', 'stop': '2010-03-01', 'step': '10d'}
        result = Horizons(
            id=obj['id'], location='500@10', epochs={'start': '2010-01-01', 'stop': '2010-03-01', 'step': '10d'}).vectors()

        # Extract the position and velocity data from the result and convert to astropy units
        # Turn the data into a numpy array
        velocity = np.array([result['vx'][0], result['vy']
                            [0], result['vz'][0]]) * u.au/u.day

        position = np.array(
            [result['x'][0], result['y'][0], result['z'][0]]) * u.au

        # Append the data to the list
        data.append({'name': obj['name'], 'mass': obj['mass'],
                    'position': [position], 'velocity': [velocity], 'color': obj['color']})

    """
    # Print the data in a more readable format
    for obj in data:
        print(obj['name'])
        print(obj['position'])
        print(obj['velocity'])
        print(obj['mass'])
        print('')
    """
    return data


def get_end_date(start_date, total_time):
    total_time = total_time.value
    start_date = start_date.split('.')
    start_date[2] = str(int(start_date[2]) + int(total_time))
    # reverse the split
    start_date = '.'.join(start_date)

    return julian_date(start_date)


def evaluate_results(simulation_results, end_conditions):
    for i in range(len(simulation_results)):
        print(simulation_results[i]['name'])
        # print(simulation_results[i]['position'][-1])
        # print(end_conditions[i]['position'])
        a = simulation_results[i]['position'][-1]
        b = end_conditions[i]['position'][0]
        ab = a - b

        distance = np.linalg.norm(ab)

        print(distance)
        print('')


def setup(start_date, total_time, time_step_size):
    end_date = get_end_date(start_date, total_time)
    end_conditions = fetch_data(end_date)

    start_date = julian_date(start_date)
    time_steps = int(total_time / time_step_size)

    inital_conditions = fetch_data(start_date)

    simulation_results = simulate(
        time_steps, time_step_size, inital_conditions)

    evaluate_results(simulation_results, end_conditions)

    # Animates the results
    #visualize_planetary_motionEndPic(simulation_results, time_steps)
    animate(simulation_results, time_steps)


def simulate(time_steps, time_step_size, initial_conditions):

    # Setup the initial conditions
    current_conditions = initial_conditions

    # Calculate the gravitational constant
    G = 6.674e-11 * u.m**3 * u.kg**-1 * u.s**-2

    # Iterate over the time steps
    for i in range(time_steps):

        for body in current_conditions:

            otherPart = np.array([0, 0, 0]) * u.kg / u.au**2

            # F = Sum of f  = Gm * sum of m * (position of other body - position of body) / distance^3

            for other_body in current_conditions:
                if body['name'] != other_body['name']:
                    # Calculate the connection vector between the bodies
                    connectionVector = other_body['position'][-1] - \
                        body['position'][-1]

                    # Calculate the length of the connection vector
                    distance = np.linalg.norm(connectionVector)

                    # Normalize the connection vector
                    direction = connectionVector / distance

                    '''# Calculate the gravitational force between the bodies
                    force = G * (body['mass'] *
                                 other_body['mass']) / (distance**2)'''

                    thingi = other_body['mass'] / (distance**2)
                    otherPart += thingi * direction

                    '''
                    # Calculate the resultant force
                    resultingForce += force * direction'''

            resultingForce = G * body['mass'] * otherPart

            # Calculate the acceleration of the body
            acceleration = resultingForce / body['mass']

            # Calculate the new velocity of the body
            body['velocity'].append(
                body['velocity'][-1] + (acceleration * time_step_size))

        for body in current_conditions:
            # Calculate the new position of the body
            body['position'].append(
                body['position'][-1] + body['velocity'][-1] * time_step_size)

        print(str(i) + ' / ' +
              str(time_steps))

    return current_conditions


def animate(data, time_steps):
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


# setup(start_date, total_time in years, time_step_size)
setup('16.08.2005', 3 * u.year, 6 * u.day)
