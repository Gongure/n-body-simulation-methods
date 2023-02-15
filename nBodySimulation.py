from bruteforcealgorithm import *
from utilities import *
from treebasedalgorithm2 import *
from visualization import *

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


def setup(start_date, total_time, time_step_size):
    start_julian_date = julian_date(start_date)
    # print(start_julian_date)
    end_julian_date = str(float(start_julian_date) +
                          (total_time.to(u.day)).value)
    # print(end_julian_date)

    inital_conditions = fetch_data(start_julian_date)
    final_conditions = fetch_data(end_julian_date)

    time_steps = int(total_time / time_step_size)

    #brute_force_simulation_results = bruteForceSimulation(time_steps, time_step_size, inital_conditions)

    tree_based_simulation_results = treeBasedAlgorithm(
        time_steps, time_step_size, inital_conditions)

    visualize_planetary_motionEndPic(tree_based_simulation_results, time_steps)

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


# setup(start_date, total_time in years, time_step_size)
setup('16.08.2005', 1 * u.year, 1 * u.day)
