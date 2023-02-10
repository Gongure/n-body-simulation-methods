
from utilities import *


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
matplotlib.rcParams["toolbar"] = "toolmanager"


def setup(start_date, end_date, time_step_size):
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


# setup(start_date, end_date, time_step_size)
setup('16.08.2005', '16.08.2008', 1 * u.day)
