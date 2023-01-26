from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.jplhorizons import Horizons
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


PLOT_MAX = 1

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
    objects = [{'name': 'Sun', 'mass': 1.989e+30 * u.kg, 'id': '10'},  {'name': 'Mercury', 'mass': 3.3022e+23 * u.kg, 'id': '199'},  {'name': 'Venus', 'mass': 4.8685e24 * u.kg, 'id': '299'},  {'name': 'Earth', 'mass': 5.97237e24 * u.kg, 'id': '399'},  {'name': 'Mars', 'mass': 6.4185e23 * u.kg, 'id': '499'},
               {'name': 'Jupiter', 'mass': 1.8986e+27 * u.kg, 'id': '599'},  {'name': 'Saturn', 'mass': 5.6846e+26 * u.kg, 'id': '699'}]
    # Initialize an empty list to store the data
    data = []

    # Iterate over the objects
    for obj in objects:
        # Query the JPL Horizons database using Astroquery

        # id_type='majorbody' ; get_raw_response=True
        result = Horizons(
            id=obj['id'], location='500@10', epochs=date).vectors()

        print(result)

        # Extract the position and velocity data from the result and convert to astropy units
        # Turn the data into a numpy array
        velocity = np.array([result['vx'][0], result['vy']
                            [0], result['vz'][0]]) * u.au/u.day

        position = np.array(
            [result['x'][0], result['y'][0], result['z'][0]]) * u.au

        # Append the data to the list
        data.append({'name': obj['name'], 'mass': obj['mass'],
                    'position': [position], 'velocity': [velocity]})

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


def setup(start_date, total_time, time_step_size):
    start_date = julian_date(start_date)
    time_steps = int(total_time / time_step_size)

    inital_conditions = fetch_data(start_date)

    simulate(time_steps, time_step_size, inital_conditions)


def simulate(time_steps, time_step_size, initial_conditions):

    # Setup the initial conditions
    current_conditions = initial_conditions

    for body in current_conditions:
        distanceToSun = np.linalg.norm(
            body['position'][-1] - current_conditions[0]['position'][-1])

        print(body['name'] + " : " + str(distanceToSun))

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

    # For debugging purposes, print the distance between every planet to the sun
    for body in current_conditions:
        distanceToSun = np.linalg.norm(
            body['position'][-1] - current_conditions[0]['position'][-1])

        print(body['name'] + " : " + str(distanceToSun))

    # animate(results)


def animate(results):
    pass


# setup(start_date, total_time, time_step_size)
setup('16.08.2005', 1 * u.year, 1 * u.day)
