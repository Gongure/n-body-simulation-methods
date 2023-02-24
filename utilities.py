
import numpy as np
from astropy import units as u
from astropy.time import Time
from astroquery.jplhorizons import Horizons


def evaluate_results(simulation_results, end_conditions):
    total_deviation = 0

    for i in range(len(simulation_results)):

        # the deviation of the ith object

        # the final position of the ith object
        final_position = end_conditions[i]['position'][-1]

        # the final position of the ith object in the simulation
        final_simulation_position = simulation_results[i]['position'][-1]

        # deviation
        deviation = np.linalg.norm(
            (final_position - final_simulation_position).value)

        '''
        # the penultimate position of the ith object
        penultimate_position = end_conditions[i]['position'][-2]

        # is this the same as the following?
        deviation = np.linalg.norm((final_position - final_simulation_position).value) / \
            np.linalg.norm((final_position - penultimate_position).value + 1)
        '''

        # add the deviation of the ith object to the total deviation
        total_deviation += deviation

    # the average deviation is the total deviation divided by the number of objects
    average_deviation = total_deviation / len(simulation_results)

    return average_deviation


def fetch_data(start_date, end_date, time_step_size):

    objects = [
        {'name': 'Sun', 'mass': 1.989e+30 * u.kg,
            'id': '10', 'color': 'yellow', 'GM': 1.3271244004193938E+11},
        {'name': 'Mercury', 'mass': 3.3022e+23 * u.kg, 'id': '199',
            'color': 'gray', 'GM': 2.2031780000000021E+04},
        {'name': 'Venus', 'mass': 4.8685e24 * u.kg, 'id': '299',
            'color': 'yellow', 'GM': 3.2485859200000006E+05},
        {'name': 'Earth', 'mass': 5.97237e24 * u.kg, 'id': '399',
            'color': 'blue', 'GM': 3.9860043543609598E+05},
        {'name': 'Mars', 'mass': 6.4185e23 * u.kg, 'id': '499',
            'color': 'red', 'GM': 4.282837362069909E+04},
        {'name': 'Jupiter', 'mass': 1.8986e+27 * u.kg, 'id': '599',
            'color': 'orange', 'GM': 1.266865349218008E+08},
        {'name': 'Saturn', 'mass': 5.6846e+26 * u.kg, 'id': '699',
            'color': 'yellow', 'GM': 3.793120749865224E+07},
        {'name': 'Uranus', 'mass': None,
            'id': '799', 'color': 'blue', 'GM': 5.793951322279009E+06},
        {'name': 'Neptune', 'mass': None,
            'id': '899', 'color': 'blue', 'GM': 6.835099502439672E+06},
        {'name': 'Pluto', 'mass': None,
            'id': '999', 'color': 'red', 'GM': 8.696138177608748E+02},
        {'name': 'Luna', 'mass': None, 'id': '301',
            'color': 'white', 'GM': 4.9028000661637961E+03},
        {'id': '401', 'GM': 0.0007087546066894452,
            'name': 'Phobos', 'mass': None, 'color': 'gray'},
        {'id': '402', 'GM': 9.615569648120313e-05,
            'name': 'Deimos', 'mass': None, 'color': 'gray'},
        {'id': '501', 'GM': 5959.916033410404, 'name': 'Io',
            'mass': None, 'color': 'yellow'},
        {'id': '502', 'GM': 3202.738774922892, 'name': 'Europa',
            'mass': None, 'color': 'white'},
        {'id': '503', 'GM': 9887.834453334144, 'name': 'Ganymede',
            'mass': None, 'color': 'white'},
        {'id': '504', 'GM': 7179.28936139727, 'name': 'Callisto',
            'mass': None, 'color': 'white'},
        {'id': '505', 'GM': 0.1378480571202615,
            'name': 'Amalthea', 'mass': None, 'color': 'red'},
        {'id': '601', 'GM': 2.503522884661795,
            'name': 'Mimas', 'mass': None, 'color': 'white'},
        {'id': '602', 'GM': 7.211292085479989, 'name': 'Enceladus',
            'mass': None, 'color': 'white'},
        {'id': '603', 'GM': 41.21117207701302,
            'name': 'Tethys', 'mass': None, 'color': 'white'},
        {'id': '604', 'GM': 73.11635322923193, 'name': 'Dione',
            'mass': None, 'color': 'white'},
        {'id': '605', 'GM': 153.9422045545342, 'name': 'Rhea',
            'mass': None, 'color': 'white'},
        {'id': '606', 'GM': 8978.138845307376, 'name': 'Titan',
            'mass': None, 'color': 'orange'},
        {'id': '607', 'GM': 0.3718791714191668,
            'name': 'Hyperion', 'mass': None, 'color': 'gray'},
        {'id': '608', 'GM': 120.5134781724041, 'name': 'Iapetus',
            'mass': None, 'color': 'black'},
        {'id': '609', 'GM': 0.5531110414633374,
            'name': 'Phoebe', 'mass': None, 'color': 'black'},
        {'id': '610', 'GM': 0.1266231296945636,
            'name': 'Janus', 'mass': None, 'color': 'white'},
        {'id': '611', 'GM': 0.03513977490568457,
            'name': 'Epimetheus', 'mass': None, 'color': 'gray'},
        {'name': 'Atlas', 'mass': None, 'id': '615',
            'color': 'gray', 'GM': 3.759718886965353e-04},
        {'name': 'Prometheus', 'mass': None, 'id': '616',
            'color': 'gray', 'GM': 1.066368426666134e-02},
        {'name': 'Pandora', 'mass': None, 'id': '617',
            'color': 'gray', 'GM': 9.1037683110543e-03},
        {'name': 'Ariel', 'mass': None, 'id': '701',
            'color': 'gray', 'GM': 83.46344431770477},
        {'name': 'Umbriel', 'mass': None, 'id': '702',
            'color': 'gray', 'GM': 85.09338094489388},
        {'name': 'Titania', 'mass': None, 'id': '703',
            'color': 'gray', 'GM': 226.9437003741248},
        {'name': 'Oberon', 'mass': None, 'id': '704',
            'color': 'gray', 'GM': 205.3234302535623},
        {'name': 'Miranda', 'mass': None, 'id': '705',
            'color': 'gray', 'GM': 4.3195168992321},
        {'name': 'Charon', 'mass': None, 'id': '901',
            'color': 'gray', 'GM': 105.8799888601881},
        {'name': 'Nix', 'mass': None, 'id': '902',
            'color': 'gray', 'GM': 0.00304817564816976},
        {'name': 'Hydra', 'mass': None, 'id': '903',
            'color': 'gray', 'GM': 0.003211039206155255}
    ]
    # Fehlenden Massen hinzu fügen
    # GM in richtige Einheit umrechnen
    G = 6.674e-20 * u.km**3 * u.kg**-1 * u.s**-2
    for obj in objects:
        obj['GM'] = obj['GM'] * u.km**3 / u.s**2
        obj['mass'] = obj['GM'] / G
        obj['GM'] = obj['GM'].to(u.au**3 / u.day**2)
        # G in km^3 / kg * s^2

    data = []

    # Die Funktion unterscheidet zwischen zwei Fällen: Wenn start_date == 'on' ist, wird die Position und Geschwindigkeit der Objekte zu einem bestimmten Zeitpunkt berechnet. Wenn start_date == 'off' ist, wird die Position und Geschwindigkeit der Objekte zu einem bestimmten Zeitraum berechnet.
    if start_date == 'on':
        for obj in objects:

            # Die Daten werden von JPL Horizons geholt
            result = Horizons(
                id=obj['id'], location='500@10', epochs=end_date).vectors()

            # Die velocity aus den Daten wird in ein Array umgewandelt
            velocity = np.array([result['vx'][0], result['vy']
                                [0], result['vz'][0]]) * u.au/u.day

            # Die position aus den Daten wird in ein Array umgewandelt
            position = np.array(
                [result['x'][0], result['y'][0], result['z'][0]]) * u.au

            # Die Neuen Daten werden in ein Array gespeichert
            data.append({'name': obj['name'], 'mass': obj['mass'], 'GM': obj['GM'],
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

            iso_start = Time(start_date, format='jd').iso
            iso_end = Time(end_date, format='jd').iso
            result = Horizons(
                id=obj['id'], location='500@10', epochs={'start': iso_start, 'stop': iso_end, 'step': time_step_size}
            ).vectors()

            # Extract the position and velocity data from the result and convert to astropy units
            # Turn the data into a numpy array

            position = []
            for i in range(len(result['x'])):
                position.append(
                    np.array([result['x'][i], result['y'][i], result['z'][i]]) * u.au)

            data.append({'name': obj['name'], 'mass': obj['mass'], 'GM': obj['GM'],
                        'position': position, 'color': obj['color']})

    return data
