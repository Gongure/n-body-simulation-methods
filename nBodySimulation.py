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
# matplotlib.rcParams["toolbar"] = "toolmanager"


def setup(start_date, total_time, time_step_size, simulation_type):
    start_julian_date = Time(start_date).jd
    # print(start_julian_date)
    end_julian_date = start_julian_date + total_time.to(u.day).value
    # print(end_julian_date)
    inital_conditions = fetch_data('on', start_julian_date, time_step_size)
    final_conditions = fetch_data(
        start_julian_date, end_julian_date, time_step_size)

    time_steps = int(total_time / time_step_size)

    resuls = []
    bboxes = []

    if simulation_type == 'brute-force':
        results = bruteForceSimulation(
            time_steps, time_step_size, inital_conditions)
        bboxes = None
    elif simulation_type == 'tree-based':
        results, bboxes = treeBasedAlgorithm(
            time_steps, time_step_size, inital_conditions)

    if input("Animate? (y/n): ") == 'y':
        animate(results, time_steps, bboxes, simulation_type)

    # Evaluates the results
    evaluation = evaluate_results(results, inital_conditions, final_conditions)

    # print(evaluate_results(brute_force_simulation_results, final_conditions))
    #print(evaluate_results(tree_based_simulation_results, final_conditions))

    """"
    end_date = get_end_date(start_date, total_time)
    end_conditions = fetch_data(end_date)

    start_date = julian_date(start_date)
    time_steps = int(total_time / time_step_size)

    inital_conditions = fetch_data(start_date)

    simulation_results = simulate(
        time_steps, time_step_size, inital_conditions)

    evaluate_results(simulation_results, end_conditions)

    # Animates the results
    animate(simulation_results, time_steps)
    visualize_planetary_motionEndPic(simulation_results, time_steps)

    """


#setup('16.08.2005', 5 * u.year, 1 * u.day, 'tree-based')

# Ask user for input


if input("Use default values? (y/n): ") == 'n':
    start_date = str(input("Enter start date (YYYY-MM-DD): "))
    total_time = float(input("Enter total time (in years): "))
    time_step_size = float(input("Enter time step size (in days): "))
    total_time *= u.year
    time_step_size *= u.day

    bruteForceSimulationResults = []
    treeBasedSimulationResults = []
    bboxes = []


else:
    start_date = '2005-08-16'
    total_time = 1 * u.year
    time_step_size = 1 * u.day

start_date = str(start_date + ' 00:00:00')

a = input(
    "Choose simulation type:\n[1] Brute Force\n[2] Barnes-Hut\n[3] Both\n")

if a == '1':
    setup(start_date, total_time, time_step_size, 'brute-force')
if a == '2':
    setup(start_date, total_time, time_step_size, 'tree-based')
if a == '3':
    setup(start_date, total_time, time_step_size, 'brute-force')
    setup(start_date, total_time, time_step_size, 'tree-based')
