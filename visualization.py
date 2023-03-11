import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backend_tools import ToolBase
matplotlib.rcParams["toolbar"] = "toolmanager"


def animate(data, time_steps, bboxes, name, real_data):
    # Entfernt die Einheiten, damit die Animation funktioniert
    for x in data:
        for i in range(len(x['position'])):
            x['position'][i] = x['position'][i].value

    # Tue das gleiche für die realen Daten, falls diese vorhanden sind
    if real_data is not None:
        for x in real_data:
            for i in range(len(x['position'])):
                x['position'][i] = x['position'][i].value

    # Die folgenden Variablen werden für die Animation benötigt
    PLOT_MAX = 7
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    fig.suptitle(name, fontsize=12)

    colors = [body['color'] for body in data]

    setup_buttons(fig)
    global show_axes
    show_axes = True

    # Für jede Boundary Box werden die Seiten berechnet, um diese dann zu plotten
    if bboxes is not None:
        new_bboxes = []
        for line in bboxes:
            newline = []
            for bbox in line:
                ooo = bbox[0]
                iii = bbox[1]
                ooi = [ooo[0], ooo[1], iii[2]]
                oio = [ooo[0], iii[1], ooo[2]]
                oii = [ooo[0], iii[1], iii[2]]
                ioo = [iii[0], ooo[1], ooo[2]]
                ioi = [iii[0], ooo[1], iii[2]]
                iio = [iii[0], iii[1], ooo[2]]

                # connect ooo ioo -ioi ioo- iio - iii iio - oio - oii oio - ooo, ooi ioi iii oii ooi
                array = np.array([ooo, ioo, ioi, ioo, iio, iii, iio,
                                  oio, oii, oio, ooo, ooi, ioi, iii, oii, ooi])
                newline.append(array)
            new_bboxes.append(newline)

    for i in range(time_steps):
        while is_paused:
            plt.pause(0.001)

        plt.cla()

        if not show_axes:
            plt.axis('off')
        else:
            plt.axis('on')

        # Update den Titel und die Achsenbeschriftung
        ax.set_title('Time: ' + str(i) + ' days')

        ax.set_xlabel('X in AE')
        ax.set_ylabel('Y in AE')
        ax.set_zlabel('Z in AE')

        # Setze die Limits der Achsen
        ax.set_xlim(-PLOT_MAX, PLOT_MAX)
        ax.set_ylim(-PLOT_MAX, PLOT_MAX)
        ax.set_zlim(-PLOT_MAX, PLOT_MAX)

        current_positions_x = []
        current_positions_y = []
        current_positions_z = []

        # Aus allen dem Array der Positionen für den aktuellen Zeitpunkt die x-, y- und z-Koordinaten extrahieren
        for body in data:
            current_positions_x += [body['position'][i][0]]
            current_positions_y += [body['position'][i][1]]
            current_positions_z += [body['position'][i][2]]

        # Um diese dann zu plotten
        ax.scatter(
            xs=current_positions_x, ys=current_positions_y, zs=current_positions_z, c=colors)

        # Das gleiche für die realen Daten, falls diese vorhanden sind und die Option zum Vergleich aktiviert ist
        if show_comparison:
            real_positions_x = []
            real_positions_y = []
            real_positions_z = []

            for body in real_data:
                real_positions_x += [body['position'][i][0]]
                real_positions_y += [body['position'][i][1]]
                real_positions_z += [body['position'][i][2]]

            ax.scatter(
                xs=real_positions_x, ys=real_positions_y, zs=real_positions_z, c='red')

        # Plotten der Boundary Boxes
        if bboxes is not None:
            line = new_bboxes[i]
            for bbox_points in line:
                x_points = bbox_points[:, 0]
                y_points = bbox_points[:, 1]
                z_points = bbox_points[:, 2]
                ax.plot(x_points, y_points, z_points, color='black')
                if is_slow:
                    plt.pause(0.001)

        plt.pause(0.001)

# Einmaliges anzeigen der Bahnen


def visualize_planetary_motion_end_pic(data):
    for x in data:
        for i in range(len(x['position'])):
            x['position'][i] = x['position'][i].value

    ax = plt.axes(projection="3d")

    PLOT_MAX = 5

    ax.set_xlim(-PLOT_MAX, PLOT_MAX)
    ax.set_ylim(-PLOT_MAX, PLOT_MAX)
    ax.set_zlim(-PLOT_MAX, PLOT_MAX)

    for planet in data:
        name = planet["name"]
        positions = planet["position"]

        x = [pos[0] for pos in positions]
        y = [pos[1] for pos in positions]
        z = [pos[2] for pos in positions]

        ax.plot3D(x, y, z, label=name)

    plt.legend()
    plt.show()


# Knöpfe für die Animation
class Pause(ToolBase):
    image = r"/assets/pic.png"
    description = "Pause"

    plt.show()

    def trigger(self, *args, **kwargs):
        global is_paused
        is_paused = not is_paused


class Axes(ToolBase):
    image = r"/assets/pic.png"
    description = "Show Axes"

    plt.show()

    def trigger(self, *args, **kwargs):
        global show_axes
        show_axes = not show_axes


class Slow(ToolBase):
    image = r"/assets/pic.png"
    description = "Show Tree Construction"

    plt.show()

    def trigger(self, *args, **kwargs):
        global is_slow
        is_slow = not is_slow


class Comparison(ToolBase):
    image = r"/assets/pic.png"
    description = "Show Comparison"

    plt.show()

    def trigger(self, *args, **kwargs):
        global show_comparison
        show_comparison = not show_comparison


def setup_buttons(fig):
    # add tools
    tm = fig.canvas.manager.toolmanager
    tm.add_tool("Pause", Pause)
    fig.canvas.manager.toolbar.add_tool(
        tm.get_tool("Pause"), "toolgroup")
    global is_paused
    is_paused = False

    tm = fig.canvas.manager.toolmanager
    tm.add_tool("Axes", Axes)
    fig.canvas.manager.toolbar.add_tool(
        tm.get_tool("Axes"), "toolgroup")
    global show_axes
    show_axes = False

    tm = fig.canvas.manager.toolmanager
    tm.add_tool("slow", Slow)
    fig.canvas.manager.toolbar.add_tool(
        tm.get_tool("slow"), "toolgroup")
    global is_slow
    is_slow = True

    tm = fig.canvas.manager.toolmanager
    tm.add_tool("comparison", Comparison)
    fig.canvas.manager.toolbar.add_tool(
        tm.get_tool("comparison"), "toolgroup")
    global show_comparison
    show_comparison = False
