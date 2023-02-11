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


def bruteForceSimulation(time_steps, time_step_size, initial_conditions):

    # Setup the initial conditions
    current_conditions = initial_conditions

    # Calculate the gravitational constant
    G = 6.674e-11 * u.m**3 * u.kg**-1 * u.s**-2

    # Iterate over the time steps
    for i in range(time_steps):

        for body in current_conditions:

            resultingForce = np.array([0, 0, 0]) * u.N

            for other_body in current_conditions:
                if body['name'] != other_body['name']:
                    # Calculate the connection vector between the bodies
                    connectionVector = other_body['position'][-1] - \
                        body['position'][-1]

                    # Calculate the length of the connection vector
                    distance = np.linalg.norm(connectionVector)

                    # Normalize the connection vector
                    direction = connectionVector / distance

                    # Calculate the gravitational force between the bodies
                    force = G * (body['mass'] *
                                 other_body['mass']) / (distance**2)

                    # Calculate the resultant force
                    resultingForce += force * direction

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
