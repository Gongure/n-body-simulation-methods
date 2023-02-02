import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def func2():

    plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True
    # if figure.autolayout was disabled what would be different
    #

    N = 50
    fps = 250
    frn = 75

    x = np.linspace(-4, 4, N + 1)
    x, y = np.meshgrid(x, x)
    zarray = np.zeros((N + 1, N + 1, frn))

    def f(x, y, sig): return 1 / np.sqrt(sig) * \
        np.exp(-(x ** 2 + y ** 2) / sig ** 2)

    for i in range(frn):
        zarray[:, :, i] = f(x, y, 1.5 + np.sin(i * 2 * np.pi / frn))

    def change_plot(frame_number, zarray, plot):
        plot[0].remove()
        plot[0] = ax.plot_surface(
            x, y, zarray[:, :, frame_number], cmap="afmhot_r")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    plot = [ax.plot_surface(x, y, zarray[:, :, 0],
                            color='0.75', rstride=1, cstride=1)]

    ax.set_zlim(0, 1.1)
    ani = animation.FuncAnimation(
        fig, change_plot, frn, fargs=(zarray, plot), interval=1000 / fps)

    ax.axis('on')

    plt.show()


def func1():
    # Fixing random state for reproducibility
    np.random.seed(19680801)

    def random_walk(num_steps, max_step=0.05):
        """Return a 3D random walk as (num_steps, 3) array."""
        start_pos = np.random.random(3)
        steps = np.random.uniform(-max_step, max_step, size=(num_steps, 3))
        walk = start_pos + np.cumsum(steps, axis=0)
        return walk

    def update_lines(num, walks, lines):
        for line, walk in zip(lines, walks):
            # NOTE: there is no .set_data() for 3 dim data...
            line.set_data(walk[:num, :2].T)

            line.set_3d_properties(walk[:num, 2])
            #print("2: " + str(walk[:num, 2]))

        return lines

    # Data: 40 random walks as (num_steps, 3) arrays
    num_steps = 40
    walks = [random_walk(num_steps) for index in range(40)]

    # Attaching 3D axis to the figure
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    # Create lines initially without data
    lines = [ax.plot([], [], [])[0] for _ in walks]

    # Setting the axes properties
    ax.set(xlim3d=(0, 1), xlabel='X')
    ax.set(ylim3d=(0, 1), ylabel='Y')
    ax.set(zlim3d=(0, 1), zlabel='Z')

    # Creating the Animation object
    ani = animation.FuncAnimation(
        fig, update_lines, num_steps, fargs=(walks, lines), interval=100)

    plt.show()


func2()
