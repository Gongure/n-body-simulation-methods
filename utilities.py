
import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord
from astropy import units as u


def evaluate_results(simulation_results, end_conditions):

    # array of distances
    distances = []

    for i in range(len(simulation_results)):
        print(simulation_results[i]['name'])
        # print(simulation_results[i]['position'][-1])
        # print(end_conditions[i]['position'])
        a = simulation_results[i]['position'][-1]
        b = end_conditions[i]['position'][0]
        ab = a - b

        distance = np.linalg.norm(ab)

        distances.append(distance)

        print(distance)
        print('')
    return distances


def fetch_data(date):

    # Define a list of objects to retrieve data for
    objects = [{'name': 'Sun', 'mass': 1.989e+30 * u.kg, 'id': '10', 'color': 'yellow'},  {'name': 'Mercury', 'mass': 3.3022e+23 * u.kg, 'id': '199', 'color': 'gray'},  {'name': 'Venus', 'mass': 4.8685e24 * u.kg, 'id': '299', 'color': 'yellow'},  {'name': 'Earth', 'mass': 5.97237e24 *
                                                                                                                                                                                                                                                        u.kg, 'id': '399', 'color': 'blue'},  {'name': 'Mars', 'mass': 6.4185e23 * u.kg, 'id': '499', 'color': 'red'},               {'name': 'Jupiter', 'mass': 1.8986e+27 * u.kg, 'id': '599', 'color': 'orange'},  {'name': 'Saturn', 'mass': 5.6846e+26 * u.kg, 'id': '699', 'color': 'yellow'}]

    data = []

    # Iterate over the objects
    for obj in objects:
        # Query the JPL Horizons database using Astroquery

        # id_type='majorbody' ; get_raw_response=True
        result = Horizons(
            id=obj['id'], location='500@10', epochs=date).vectors()

        # Extract the position and velocity data from the result and convert to astropy units
        # Turn the data into a numpy array
        velocity = np.array([result['vx'][0], result['vy']
                            [0], result['vz'][0]]) * u.au/u.day

        position = np.array(
            [result['x'][0], result['y'][0], result['z'][0]]) * u.au

        # Append the data to the list
        data.append({'name': obj['name'], 'mass': obj['mass'],
                    'position': [position], 'velocity': [velocity], 'color': obj['color']})
    return data

    # summ up how the array "data" is structured. Leave out the color:
    #data = [{'name': 'Sun', 'mass': 1.989e+30 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Mercury', 'mass': 3.3022e+23 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Venus', 'mass': 4.8685e24 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Earth', 'mass': 5.97237e24 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Mars', 'mass': 6.4185e23 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Jupiter', 'mass': 1.8986e+27 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Saturn', 'mass': 5.6846e+26 * u.kg, 'position': [position], 'velocity': [velocity]}]


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
