import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time


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
        'color': 'yellow', 'GM': 3.793120749865224E+07}
]

for obj in objects:
    GM = obj['GM'] * u.km**3 / u.s**2
    # if obj['mass'] == None:
    #G in km^3 / kg * s^2
    G = 6.674e-20 * u.km**3 * u.kg**-1 * u.s**-2
    # G in km^3 / kg * s^2

    obj['mass'] = GM / G
    #(6.67408e-20 * u.au**3 / u.day**2)

print(objects)
