import numpy as np
from astropy import units as u


def brute_force_simulation(time_steps, time_step_size, initial_conditions):
    current_conditions = initial_conditions

    # Iteriere über die Zeitschritte
    for i in range(time_steps):

        # Iteriere über die Körper
        for body in current_conditions:

            integration_result = np.array([0, 0, 0]) * u.kg / u.au**2

            # Gehe über alle anderen Körper
            for other_body in current_conditions:
                if body['name'] != other_body['name']:
                    # Berechne den Vektor zwischen den beiden Körpern
                    connection_vector = other_body['position'][-1] - \
                        body['position'][-1]

                    # Berechne die Distanz zwischen den beiden Körpern
                    distance = np.linalg.norm(connection_vector)

                    # Normalisiere den Vektor
                    direction = connection_vector / distance

                    # Berechne die Kraft zwischen den beiden Körpern
                    integration_force_step = other_body['mass'] / (distance**2)
                    integration_result += integration_force_step * direction

            # F = Sum of f  = GM * sum of m * (position of other body - position of body) / distance^3
            resulting_force = body['GM'] * integration_result

            # Berechne die Beschleunigung des Körpers mit a = F / m
            acceleration = resulting_force / body['mass']

            # Berechne die neue Geschwindigkeit des Körpers mit v = v0 + a * dt
            body['velocity'].append(
                body['velocity'][-1] + (acceleration * time_step_size))

        for body in current_conditions:
            # Berechne die neue Position des Körpers mit p = p0 + v * dt
            body['position'].append(
                body['position'][-1] + body['velocity'][-1] * time_step_size)

        print(str(i) + ' / ' +
              str(time_steps))

    return current_conditions
