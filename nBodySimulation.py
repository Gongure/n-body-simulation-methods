import time

from bruteforcealgorithm import *
from treebasedalgorithm import *
from utilities import *
from visualization import *

SKIP = False
body_number = None
theta = 0.5


# Die Setup Funktion bildet die Schnittstelle zwischen der Benutzereingabe und dem eigentlichen Programm.
# Time_step_size ist die Zeit zwischen zwei Iterationen. In Bezug auf die Euler Methode entspricht das der Schrittweite.
def setup(start_date, total_time, time_step_size, simulation_type):

    # Das Datum wird in Julianische Tage umgerechnet
    start_julian_date = Time(start_date).jd
    end_julian_date = start_julian_date + total_time.to(u.day).value

    # Die Daten werden aus der Nasa Datenbank abgerufen
    # Es werden drei Parameter übergeben: das anfangs Datum, das Enddatum und die Schrittweite
    initial_conditions = fetch_data('on', start_julian_date, time_step_size)
    real_conditions = fetch_data(
        start_julian_date, end_julian_date, time_step_size)

    global body_number
    # Aus den initialen Bedingungen nur die ersten [body_number] Körper behalten
    initial_conditions = initial_conditions[:body_number]
    real_conditions = real_conditions[:body_number]

    # Falls die Benutzereingabe nicht übersprungen werden soll,
    # Gibt es die Möglichkeit, zusätzliche Körper hinzuzufügen
    if not SKIP:
        if input('Add bodies? (y/n) [This will break the accuracy prediction]: ') == 'y':
            body = input(
                'Enter: name, mass, x, y, z, vx, vy, vz, color (with xyz in AU & vxvyvz in AU/day && no spacing): ')
            body = body.split(',')
            body = {'name': body[0],
                    'mass': float(body[1]) * u.kg,
                    'position': [np.array([float(body[2]), float(body[3]), float(body[4])]) * u.au],
                    'velocity': [np.array([float(body[5]), float(body[6]), float(body[7])]) * u.au / u.day],
                    'color': body[8]}

            initial_conditions.append(body)
            for j in range(len(real_conditions[1]['position'])):
                real_conditions[-1]['position'][j] = body['position'][0]

    time_steps = int(total_time / time_step_size)

    results = []
    s_bboxes = []

    # calculate the efficiency of the simulation in terms of time
    s_time_result = 0

    if simulation_type == 'brute-force':
        s_time_result = time.time()
        results = brute_force_simulation(
            time_steps, time_step_size, initial_conditions)
        s_time_result = time.time() - s_time_result

        s_bboxes = None
    elif simulation_type == 'tree-based':
        global theta
        s_time_result = time.time()
        results, s_bboxes = tree_based_algorithm(
            time_steps, time_step_size, initial_conditions, theta)
        s_time_result = time.time() - s_time_result

    animation_type = 0
    if not SKIP:
        animation_type = input(
            'Choose animation type:\n[1] Live3D\n[2] EndPic\n[3] Skip\n')

    if animation_type == '1':
        animate(copy.deepcopy(results), time_steps, s_bboxes, simulation_type,
                copy.deepcopy(real_conditions))
    elif animation_type == '2':
        visualize_planetary_motion_end_pic(copy.deepcopy(results))

    # Evaluates the results
    return evaluate_results(results, real_conditions), s_time_result

    # print(evaluate_results(brute_force_simulation_results, real_conditions))
    # print(evaluate_results(tree_based_simulation_results, real_conditions))


# setup('16.08.2005', 5 * u.year, 1 * u.day, 'tree-based')


# Ask user for input
if input("Use default values? (y/n): ") == 'n':
    setup_start_date = str(input("Enter start date (YYYY-MM-DD): "))
    setup_total_time = float(input("Enter total time (in years): "))
    setup_time_step_size = float(
        input("Enter time step size (in hours): "))
    setup_total_time *= u.year
    setup_time_step_size *= u.hour

    bruteForceSimulationResults = []
    treeBasedSimulationResults = []
    bboxes = []

    theta = float(input('Input theta for Barnes-Hut algorithm: '))
else:
    setup_start_date = '2005-08-16'
    setup_total_time = 1 * u.year
    setup_time_step_size = 1 * u.day


body_number = int(input(
    'Input integer for number of bodies to simulate (Planets = 0-8, Moons = 8-40): '))


setup_start_date = str(setup_start_date + ' 00:00:00')

a = input(
    "Choose simulation type:\n[1] Brute Force\n[2] Barnes-Hut\n[3] Both\n[4] Accuracy + Time Analysis\n[5] Full Time "
    "Analysis\n")

brute_force_deviation = 0
tree_based_deviation = 0

if a == '1':
    brute_force_deviation = setup(
        setup_start_date, setup_total_time, setup_time_step_size, 'brute-force')
if a == '2':
    theta = float(input('Input theta for Barnes-Hut algorithm: '))
    tree_based_deviation = setup(
        setup_start_date, setup_total_time, setup_time_step_size, 'tree-based')

if a == '3':
    brute_force_deviation = setup(
        setup_start_date, setup_total_time, setup_time_step_size, 'brute-force')
    tree_based_deviations = setup(
        setup_start_date, setup_total_time, setup_time_step_size, 'tree-based')

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

    a = int(input('Beginning setup_time_step_size (in hours): '))
    b = int(input('Ending setup_time_step_size (in hours): '))
    SKIP = True

    body_number = int(input('Number of bodies to simulate: '))

    analysis = {'Time_step_size': [], 'bf_deviation': [],
                'tb_deviation': [], 'bf_time': [], 'tb_time': []}

    for i in range(a, b):
        analysis['Time_step_size'].append(i)

        bf_deviation, bf_time = setup(
            setup_start_date, 1 * u.year, i * u.hour, 'brute-force')

        analysis['bf_deviation'].append(100 * bf_deviation)
        analysis['bf_time'].append(bf_time)

        tb_deviation, tb_time = setup(
            setup_start_date, 1 * u.year, i * u.hour, 'tree-based')

        analysis['tb_deviation'].append(100 * tb_deviation)
        analysis['tb_time'].append(tb_time)

        print(str(i) + " / " + str(b - 1))

        # print analysis in a nice format
    print('Time_step_size\tbf_deviation\ttb_deviation\tbf_time\ttb_time')
    for i in range(len(analysis['Time_step_size'])):
        print(str(analysis['Time_step_size'][i]) + ": " +
              str(analysis['bf_deviation'][i]) + ", " + str(analysis['tb_deviation'][i]) + ", " +
              str(analysis['bf_time'][i]) + ", " + str(analysis['tb_time'][i]))


if a == '5':
    # Time analysis mode
    # Plot the time it takes for the brute force algorithm and the tree based algorithm
    # For a timespan of 1 year and time increment of 1s and with 2 - 500 bodies
    time_range = input('Range of bodies (n_start,n_end): ')
    setup_time_step_size = float(input('Time step size (in hours): ')) * u.hour

    SKIP = True

    time_range.split(',')
    time_result = {'bodies': [], 'bf_time': [], 'tb_time': []}

    for i in range(int(time_range[0]), int(time_range[1])):

        body_number = i

        time_result['bodies'].append(i)

        # d is not important
        d, bf_time = setup(
            setup_start_date, 1 * u.year, setup_time_step_size, 'brute-force')
        time_result['bf_time'].append(bf_time)

        d, tb_time = setup(
            setup_start_date, 1 * u.year, setup_time_step_size, 'tree-based')
        time_result['tb_time'].append(tb_time)

        print(str(i) + " / " + str(time_range[1]))

    # print time_result in a nice  format
    print('bodies\tbf_time\ttb_time')
    for i in range(len(time_result['bodies'])):
        print(str(time_result['bodies'][i]) + ': ' + '\t' +
              str(time_result['bf_time'][i]) + ', ' + '\t' + str(time_result['tb_time'][i]))
