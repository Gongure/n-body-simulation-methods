import copy
from astropy import units as u
from visualization import *

# Die Gravitationskonstante
G = 6.674e-11 * u.m**3 * u.kg**-1 * u.s**-2


# Erstelle eine Klasse für Knoten im Baum
class Node:
    children = None
    mass = None
    center_of_mass = None
    # boundary box
    bbox = None


def tree_based_algorithm(time_steps, time_step_size, initial_conditions, input_theta):
    global theta
    theta = input_theta

    boundary_boxes = []
    current_conditions = initial_conditions

    # Iterate über die Zeitschritte
    for i in range(time_steps):

        # Erstelle den Baum
        root = Node()
        # Finde die Bounding Box des gesamten Systems
        root.bbox = find_root_bbox([body['position'][-1].value
                                   for body in current_conditions])

        # clear current_boxes globally
        global current_boxes
        current_boxes = [root.bbox]

        # Körper in Baum einfügen
        for body in current_conditions:
            m = copy.deepcopy(body['mass'])
            p = copy.deepcopy(body['position'][-1])
            insert_in_tree(root, p, m)
        boundary_boxes.append(current_boxes)

        # Kraft berechnen
        for body in current_conditions:

            # Kopiere die Masse und Position des Körpers, damit diese nicht verändert werden
            m = copy.deepcopy(body['mass'])
            p = copy.deepcopy(body['position'][-1])
            resulting_force = calculate_force_for_tree(
                root, p, m)

            # Berechne die Beschleunigung des Körpers mit a = F / m
            acceleration = resulting_force / body['mass']

            # Berechne die neue Geschwindigkeit des Körpers mit v = v0 + a * dt
            body['velocity'].append(
                body['velocity'][-1] + (acceleration * time_step_size))

        # Berechne die neue Position des Körpers mit p = p0 + v * dt für jeden Körper
        for body in current_conditions:
            body['position'].append(
                body['position'][-1] + body['velocity'][-1] * time_step_size)

        print(str(i) + ' / ' +
              str(time_steps))

    return current_conditions, boundary_boxes


def insert_in_tree(node, body_position, body_mass):
    global current_boxes
    # Fall 1: Falls der Knoten x ein externer Knoten ist, der noch keinen Körper enthält, dann
    # füge den Körper b in den Knoten x ein.
    if node.mass is None:
        node.mass = copy.deepcopy(body_mass)
        node.center_of_mass = copy.deepcopy(body_position)
        return

    # Fall 2: Falls der Knoten x ein interner Knoten ist, dann aktualisiere den Mittelpunkt der Masse und die Gesamtmasse von x.
    # Füge den Körper b rekursiv in das entsprechende Octanten ein.
    elif node.children is not None:
        node.center_of_mass = copy.deepcopy((node.center_of_mass * node.mass +
                                             body_position * body_mass) / (node.mass + body_mass))
        node.mass += copy.deepcopy(body_mass)
        octant = get_octant(node.bbox, copy.deepcopy(body_position))
        if node.children[octant] is None:
            node.children[octant] = Node()
            node.children[octant].bbox = find_bbox(node.bbox, octant)
            current_boxes.append(node.children[octant].bbox)
        insert_in_tree(node.children[octant], copy.deepcopy(
            body_position), copy.deepcopy(body_mass))
        return

    # Fall 3: Falls der Knoten x ein externer Knoten ist, der bereits einen Körper enthält, dann
    # füge den Körper b in den Knoten x ein.
    # Unterteile die Region weiter in vier Kinder.
    # Füge dann beide Körper rekursiv in die entsprechenden Octanten ein.
    # Da b und c noch im selben Octanten landen können, kann es mehrere Unterteilungen während einer einzelnen Einfügung geben.
    # Aktualisiere schließlich den Mittelpunkt der Masse und die Gesamtmasse von x.
    elif node.children is None:
        node.children = [None, None, None, None, None, None, None, None]

        # Erhalte den Octanten des alten Körpers
        old_octant = get_octant(node.bbox, node.center_of_mass)
        # Erhalte den Octanten des neuen Körpers
        new_octant = get_octant(node.bbox, body_position)

        node.children[old_octant] = Node()
        node.children[old_octant].bbox = find_bbox(node.bbox, old_octant)

        current_boxes.append(node.children[old_octant].bbox)

        if new_octant != old_octant:
            node.children[new_octant] = Node()
            node.children[new_octant].bbox = find_bbox(node.bbox, new_octant)
            current_boxes.append(node.children[new_octant].bbox)

        # Füge den alten Körper in den entsprechenden Octanten ein
        insert_in_tree(node.children[old_octant], copy.deepcopy(
            node.center_of_mass), copy.deepcopy(node.mass))
        # Füge den neuen Körper in den entsprechenden Octanten ein
        insert_in_tree(node.children[new_octant], copy.deepcopy(
            body_position), copy.deepcopy(body_mass))

        # Update den Mittelpunkt der Masse und die Gesamtmasse von x
        node.center_of_mass = copy.deepcopy((node.center_of_mass * node.mass +
                                             body_position * body_mass) / (node.mass + body_mass))
        node.mass += copy.deepcopy(body_mass)
        return


def calculate_force_for_tree(node, body_position, mass):
    # Falls unser Körper b dem aktuellen Knoten x entspricht, dann
    # gibt es keine Kraft, die auf b wirkt.
    if node.mass == mass:
        return np.array([0, 0, 0]) * u.N

    # Falls der Knoten x ein externer Knoten ist,
    # berechne die Kraft, die der aktuelle Knoten auf den Körper b ausübt, und füge diese Menge zu b’s Gesamtkraft hinzu.
    elif node.children is None:
        resulting_force = calculate_gravity(
            node.center_of_mass, body_position, node.mass, mass)
        return resulting_force

    else:
        # Sonst berechne den Verhältnis s/d.
        # s = Seitenlänge des Octanten des aktuellen Knotens
        s = np.linalg.norm(node.bbox[1][0] - node.bbox[0][0])
        # d = Entfernung zwischen dem Mittelpunkt der Masse des aktuellen Knotens und dem Mittelpunkt der Masse von b
        d = np.linalg.norm((body_position - node.center_of_mass).value)
        # Falls s/d < θ, behandele diesen internen Knoten als einen einzelnen Körper und berechne die Kraft, die er auf den Körper b ausübt,
        # und füge diese Menge zu b’s Gesamtkraft hinzu.
        if s/d < theta:
            resulting_force = calculate_gravity(
                node.center_of_mass, body_position, node.mass, mass)

            return resulting_force
        else:
            # Ansonsten führe die Prozedur rekursiv auf jeden der aktuellen Knotens’ Kinder aus.
            resulting_force = np.array([0, 0, 0]) * u.N

            for child in node.children:
                if child is not None:
                    resulting_force += calculate_force_for_tree(
                        child, body_position, mass)
            return resulting_force


########################## Helper functions ##########################

# Berechne die Kraft, die der Körper a auf den Körper b ausübt
def calculate_gravity(other_body_position, body_position, other_body_mass, body_mass):
    # Berechne den Verbindungsvektor zwischen den beiden Körpern
    connection_vector = other_body_position - body_position

    # Berechne die Entfernung zwischen den beiden Körpern
    distance = np.linalg.norm(connection_vector)

    # Normalisiere den Verbindungsvektor um die Richtung zu erhalten
    direction = connection_vector / distance

    # Berechne die Kraft, die der Körper a auf den Körper b ausübt
    # F = G * (m1 * m2) / (d**2)
    force = G * (body_mass * other_body_mass) / (distance**2)

    # Berechne die gerichtete Kraft und gib sie zurück
    resulting_force = force * direction

    return resulting_force

# Finde die Bounding Box des Wurzelknotens
# Diese Bounding Box ist ein Würfel, der jeden Körper der Simulation enthält


def find_root_bbox(array_of_positions):

    # Finde die kleinste und größte x-, y- und z-Koordinate
    min_x = min(array_of_positions, key=lambda x: x[0])[0]
    max_x = max(array_of_positions, key=lambda x: x[0])[0]
    min_y = min(array_of_positions, key=lambda x: x[1])[1]
    max_y = max(array_of_positions, key=lambda x: x[1])[1]
    min_z = min(array_of_positions, key=lambda x: x[2])[2]
    max_z = max(array_of_positions, key=lambda x: x[2])[2]

    # Finde die größte Differenz zwischen den Koordinaten
    max_diff = max(max_x - min_x, max_y - min_y, max_z - min_z)

    # Finde den Mittelpunkt der Bounding Box
    center = np.array([min_x, min_y, min_z]) + \
        np.array([max_x-min_x, max_y-min_y, max_z-min_z])/2

    # Finde die Ecken der Bounding Box und gebe diese zurück
    cmin = center - np.array([max_diff, max_diff, max_diff])/2
    cmax = center + np.array([max_diff, max_diff, max_diff])/2

    return [cmin, cmax]

# Finde den Octanten, in den der Körper position passt


def get_octant(bbox, position):
    position = position.value
    cmin = bbox[0]
    cmax = bbox[1]
    a = 0
    s = cmax[0] - cmin[0]
    if position[0] > s/2 + cmin[0]:
        a = 1

    if position[1] > s/2 + cmin[1]:
        a += 2

    if position[2] > s/2 + cmin[2]:
        a += 4

    return a

# Finde die Bounding Box des Octanten octant


def find_bbox(bbox, octant):
    cmin = bbox[0]
    cmax = bbox[1]

    x = (cmax[0] - cmin[0])/2
    y = (cmax[1] - cmin[1])/2
    z = (cmax[2] - cmin[2])/2

    center = cmin + np.array([x, y, z])

    if octant == 0:
        return [cmin, center]
    elif octant == 1:
        return [np.array([center[0], cmin[1], cmin[2]]), np.array([cmax[0], center[1], center[2]])]
    elif octant == 2:
        return [np.array([cmin[0], center[1], cmin[2]]), np.array([center[0], cmax[1], center[2]])]
    elif octant == 3:
        return [np.array([center[0], center[1], cmin[2]]), np.array([cmax[0], cmax[1], center[2]])]
    elif octant == 4:
        return [np.array([cmin[0], cmin[1], center[2]]), np.array([center[0], center[1], cmax[2]])]
    elif octant == 5:
        return [np.array([center[0], cmin[1], center[2]]), np.array([cmax[0], center[1], cmax[2]])]
    elif octant == 6:
        return [np.array([cmin[0], center[1], center[2]]), np.array([center[0], cmax[1], cmax[2]])]
    elif octant == 7:
        return [center, cmax]
