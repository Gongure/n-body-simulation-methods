from bruteforcealgorithm import *
from utilities import *
from treebasedalgorithm import *
from visualization import *

from astropy.time import Time
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

import time

from copy import deepcopy

# matplotlib.rcParams["toolbar"] = "toolmanager"

SKIP = False
body_number = None
theta = 0.5


def setup(start_date, total_time, time_step_size, simulation_type):
    start_julian_date = Time(start_date).jd
    # print(start_julian_date)
    end_julian_date = start_julian_date + total_time.to(u.day).value

    # print(end_julian_date)

    inital_conditions = fetch_data('on', start_julian_date, time_step_size)

    '''real_conditions = fetch_data(
        start_julian_date, end_julian_date, time_step_size)'''
    real_conditions = fetch_data('on', end_julian_date, time_step_size)

    global body_number

    # From initial conditions, only keep the first c bodies
    inital_conditions = inital_conditions[:body_number]

    if not SKIP:
        if input('Add bodies? (y/n) [This will break the accuracy prediction]: ') == 'y':
            body = input(
                'Enter: name, mass, x, y, z, vx, vy, vz, color (with xyz in AU & vxvyvz in AU/day && no spacing): ')
            body = body.split(',')
            body = {'name': body[0], 'mass': float(body[1]) * u.kg, 'position': [np.array([float(body[2]), float(body[3]), float(
                body[4])]) * u.au], 'velocity': [np.array([float(body[5]), float(body[6]), float(body[7])]) * u.au / u.day], 'color': body[8]}

            inital_conditions.append(body)
            for i in range(len(real_conditions[1]['position'])):
                real_conditions[-1]['position'][i] = body['position'][0]

    time_steps = int(total_time / time_step_size)

    results = []
    bboxes = []

    # calculate the efficiency of the simulation in terms of time
    time_result = 0

    if simulation_type == 'brute-force':
        time_result = time.time()
        results = bruteForceSimulation(
            time_steps, time_step_size, inital_conditions)
        time_result = time.time() - time_result

        bboxes = None
    elif simulation_type == 'tree-based':
        global theta
        time_result = time.time()
        results, bboxes = treeBasedAlgorithm(
            time_steps, time_step_size, inital_conditions, theta)
        time_result = time.time() - time_result

    animationtype = 0
    if not SKIP:
        animationtype = input(
            'Choose animation type:\n[1] Live3D\n[2] EndPic\n[3] Skip\n')

    if animationtype == '1':
        animate(copy.deepcopy(results), time_steps, bboxes, simulation_type,
                copy.deepcopy(real_conditions))
    elif animationtype == '2':
        visualize_planetary_motionEndPic(copy.deepcopy(results), time_steps)

    # Evaluates the results
    return evaluate_results(results, real_conditions), time_result

    # print(evaluate_results(brute_force_simulation_results, real_conditions))
    # print(evaluate_results(tree_based_simulation_results, real_conditions))


# setup('16.08.2005', 5 * u.year, 1 * u.day, 'tree-based')

# Ask user for input
if input("Use default values? (y/n): ") == 'n':
    start_date = str(input("Enter start date (YYYY-MM-DD): "))
    total_time = float(input("Enter total time (in years): "))
    time_step_size = float(
        input("Enter time step size (in hours): "))
    total_time *= u.year
    time_step_size *= u.hour

    bruteForceSimulationResults = []
    treeBasedSimulationResults = []
    bboxes = []

    theta = float(input('Input theta for Barnes-Hut algorithm: '))
else:
    start_date = '2005-08-16'
    total_time = 1 * u.year
    time_step_size = 1 * u.day


body_number = int(input(
    'Input integer for number of bodies to simulate (Planets = 0-8, Moons = 8-40: '))


start_date = str(start_date + ' 00:00:00')

a = input(
    "Choose simulation type:\n[1] Brute Force\n[2] Barnes-Hut\n[3] Both\n[4] Acccuracy + Time Analysis\n[5] Full Time Analysis\n")

brute_force_deviation = 0
tree_based_deviation = 0


if a == '1':
    brute_force_deviation = setup(
        start_date, total_time, time_step_size, 'brute-force')
if a == '2':
    theta = float(input('Input theta for Barnes-Hut algorithm: '))
    tree_based_deviation = setup(
        start_date, total_time, time_step_size, 'tree-based')

if a == '3':
    brute_force_deviation = setup(
        start_date, total_time, time_step_size, 'brute-force')
    tree_based_deviations = setup(
        start_date, total_time, time_step_size, 'tree-based')

if a != '4' and a != '5':
    if input("Show analysis? (y/n): ") == 'y':
        if brute_force_deviation != 0:
            print("Brute force deviation: ", brute_force_deviation)
        if tree_based_deviation != 0:
            print("Tree based deviation: ", tree_based_deviation)

if a == '4':
    # Full analysis mode accuracy + time
    # Plot the deviation of the brute force algorithm and the tree based algorithm
    # For a timespan of 1 year with a time step size of 1 minute, 1 hour, 1 day, 1 week, 1 month

    a = int(input('Beginning time_step_size (in hours): '))
    b = int(input('Ending time_step_size (in hours): '))
    SKIP = True

    body_number = int(input('Number of bodies to simulate: '))

    analyis = {'Time_step_size': [], 'bf_deviation': [],
               'tb_deviation': [], 'bf_time': [], 'tb_time': []}

    for i in range(a, b):
        analyis['Time_step_size'].append(i)

        bfdeviation, bf_time = setup(
            start_date, 1 * u.year, i * u.hour, 'brute-force')

        analyis['bf_deviation'].append(100 * bfdeviation)
        analyis['bf_time'].append(bf_time)

        tbdeviation, tb_time = setup(
            start_date, 1 * u.year, i * u.hour, 'tree-based')

        analyis['tb_deviation'].append(100 * tbdeviation)
        analyis['tb_time'].append(tb_time)

        print(str(i) + " / " + str(b - 1))

        # print analyis in a nice format
    print('Time_step_size\tbf_deviation\ttb_deviation\tbf_time\ttb_time')
    for i in range(len(analyis['Time_step_size'])):
        print(str(analyis['Time_step_size'][i]) + ": " + str(analyis['bf_deviation'][i]) + ", " + str(analyis['tb_deviation'][i]) + ", " +
              str(analyis['bf_time'][i]) + ", " + str(analyis['tb_time'][i]))


if a == '5':
    # Time analysis mode
    # Plot the time it takes for the brute force algorithm and the tree based algorithm
    # For a timespan of 1 year and time increment of 1s and with 2 - 500 bodies
    time_range = input('Range of bodies (n_start,n_end): ')
    #time_step_size = input('Time step size (in hours): ')
    #time_step_size = int(time_step_size) * u.hour

    time_step_size = 1 * u.day

    SKIP = True

    time_range.split(',')
    time_result = {'bodies': [], 'bf_time': [], 'tb_time': []}

    for i in range(int(time_range[0]), int(time_range[1])):

        body_number = i

        time_result['bodies'].append(i)
        d, bf_time = setup(
            start_date, 1 * u.year, time_step_size, 'brute-force')
        time_result['bf_time'].append(bf_time)

        d, tb_time = setup(
            start_date, 1 * u.year, time_step_size, 'tree-based')
        time_result['tb_time'].append(tb_time)

        print(str(i) + " / " + str(time_range[1]))

    # print time_result in a nice  format
    print('bodies\tbf_time\ttb_time')
    for i in range(len(time_result['bodies'])):
        print(str(time_result['bodies'][i]) + ': ' + '\t' +
              str(time_result['bf_time'][i]) + ', ' + '\t' + str(time_result['tb_time'][i]))
