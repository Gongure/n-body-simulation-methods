
import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time


def evaluate_results(simulation_results, end_conditions):

    # the average deviation of the simulation is the sum of the deviations of each object
    # the deviation of an object is the distance between the final position of the object and the final position of the object in the simulation divided by the last position minus the penultimate position

    total_deviation = 0

    for i in range(len(simulation_results)):

        # the deviation of the ith object
        deviation = 0

        # the final position of the ith object
        final_position = end_conditions[i]['position'][-1]

        # the final position of the ith object in the simulation
        final_simulation_position = simulation_results[i]['position'][-1]

        # the penultimate position of the ith object
        penultimate_position = end_conditions[i]['position'][-2]

        # is this the same as the following?
        deviation = np.linalg.norm((final_position - final_simulation_position).value) / \
            np.linalg.norm((final_position - penultimate_position).value + 1)

        # add the deviation of the ith object to the total deviation
        total_deviation += deviation

    # the average deviation is the total deviation divided by the number of objects
    average_deviation = total_deviation / len(simulation_results)

    return average_deviation


def fetch_data(start_date, end_date, time_step_size):
    objects = [{'name': 'Sun', 'mass': 1.989e+30 * u.kg, 'id': '10', 'color': 'yellow'},  {'name': 'Mercury', 'mass': 3.3022e+23 * u.kg, 'id': '199', 'color': 'gray'},  {'name': 'Venus', 'mass': 4.8685e24 * u.kg, 'id': '299', 'color': 'yellow'},  {'name': 'Earth', 'mass': 5.97237e24 *
                                                                                                                                                                                                                                                        u.kg, 'id': '399', 'color': 'blue'},  {'name': 'Mars', 'mass': 6.4185e23 * u.kg, 'id': '499', 'color': 'red'},               {'name': 'Jupiter', 'mass': 1.8986e+27 * u.kg, 'id': '599', 'color': 'orange'},  {'name': 'Saturn', 'mass': 5.6846e+26 * u.kg, 'id': '699', 'color': 'yellow'}]

    data = []

    if start_date == 'on':
        for obj in objects:
            # Query the JPL Horizons database using Astroquery

            if start_date == 'on':
                # id_type='majorbody' ; get_raw_response=True
                result = Horizons(
                    id=obj['id'], location='500@10', epochs=end_date).vectors()

            # Extract the position and velocity data from the result and convert to astropy units
            # Turn the data into a numpy array
            velocity = np.array([result['vx'][0], result['vy']
                                [0], result['vz'][0]]) * u.au/u.day

            position = np.array(
                [result['x'][0], result['y'][0], result['z'][0]]) * u.au

            # Append the data to the list
            data.append({'name': obj['name'], 'mass': obj['mass'],
                        'position': [position], 'velocity': [velocity], 'color': obj['color']})
    else:
        if time_step_size.unit == u.day:
            time_step_size = str(str(int(time_step_size.value)) + 'd')
        elif time_step_size.unit == u.hour:
            time_step_size = str(str(int(time_step_size.value)) + 'h')
        elif time_step_size.unit == u.minute:
            time_step_size = str(str(int(time_step_size.value)) + 'm')
        elif time_step_size.unit == u.second:
            time_step_size = str(str(int(time_step_size.value)) + 's')
        for obj in objects:
            # Query the JPL Horizons database using Astroquery

            isostart = Time(start_date, format='jd').iso
            isoend = Time(end_date, format='jd').iso
            result = Horizons(
                id=obj['id'], location='500@10', epochs={'start': isostart, 'stop': isoend, 'step': time_step_size}).vectors()

            # Extract the position and velocity data from the result and convert to astropy units
            # Turn the data into a numpy array
            velocity = np.array([result['vx'][0], result['vy']
                                [0], result['vz'][0]]) * u.au/u.day

            position = []
            for i in range(len(result['x'])):
                position.append(
                    np.array([result['x'][i], result['y'][i], result['z'][i]]) * u.au)

            # für v unwichitg
            # Append the data to the list
            data.append({'name': obj['name'], 'mass': obj['mass'],
                        'position': position, 'velocity': [velocity], 'color': obj['color']})
    return data

    # summ up how the array "data" is structured. Leave out the color:
    #data = [{'name': 'Sun', 'mass': 1.989e+30 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Mercury', 'mass': 3.3022e+23 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Venus', 'mass': 4.8685e24 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Earth', 'mass': 5.97237e24 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Mars', 'mass': 6.4185e23 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Jupiter', 'mass': 1.8986e+27 * u.kg, 'position': [position], 'velocity': [velocity]}, {'name': 'Saturn', 'mass': 5.6846e+26 * u.kg, 'position': [position], 'velocity': [velocity]}]
'''

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
    return JD'''
