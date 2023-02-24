import numpy as np
from astropy import units as u


def brute_force_simulation(time_steps, time_step_size, initial_conditions):
    # tim_step_size defines the time between each iteration

    # Set up the initial conditions
    current_conditions = initial_conditions

    '''
    # Calculate the gravitational constant
    G = 6.674e-11 * u.m**3 * u.kg**-1 * u.s**-2
    '''

    # Iterate over the time steps
    for i in range(time_steps):

        for body in current_conditions:

            integration_result = np.array([0, 0, 0]) * u.kg / u.au**2

            # F = Sum of f  = Gm * sum of m * (position of other body - position of body) / distance^3

            for other_body in current_conditions:
                if body['name'] != other_body['name']:
                    # Calculate the connection vector between the bodies
                    connection_vector = other_body['position'][-1] - \
                        body['position'][-1]

                    # Calculate the length of the connection vector
                    distance = np.linalg.norm(connection_vector)

                    # Normalize the connection vector
                    direction = connection_vector / distance

                    '''# Calculate the gravitational force between the bodies
                    force = G * (body['mass'] *
                                 other_body['mass']) / (distance**2)'''

                    integration_force_step = other_body['mass'] / (distance**2)
                    integration_result += integration_force_step * direction

                    '''
                    # Calculate the resultant force
                    resulting_force += force * direction'''

            resulting_force = body['GM'] * integration_result

            # Calculate the acceleration of the body
            acceleration = resulting_force / body['mass']

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
