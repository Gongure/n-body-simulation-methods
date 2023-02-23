import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord
from astropy import units as u

from visualization import *

import copy


# Calculate the gravitational constant
G = 6.674e-11 * u.m**3 * u.kg**-1 * u.s**-2
# in AU?


class Node:
    children = None
    mass = None
    center_of_mass = None
    # boundary box
    bbox = None


def treeBasedAlgorithm(time_steps, time_step_size, initial_conditions, input_theta):
    global theta
    theta = input_theta
    # tim_step_size defines the time between each iteration
    boundary_boxes = []
    # Setup the initial conditions
    current_conditions = initial_conditions

    # Iterate over the time steps
    for i in range(time_steps):

        # Construct the tree
        root = Node()
        # root.center_of_mass = []
        root.bbox = find_root_bbox([body['position'][-1].value
                                   for body in current_conditions])

        # clear current_boxes gloobally
        global current_boxes
        current_boxes = []
        current_boxes.append(root.bbox)

        for body in current_conditions:
            m = copy.deepcopy(body['mass'])
            p = copy.deepcopy(body['position'][-1])
            insertInTree(root, p, m)
        boundary_boxes.append(current_boxes)

        for body in current_conditions:

            m = copy.deepcopy(body['mass'])
            p = copy.deepcopy(body['position'][-1])
            resultingForce = calculateForce(
                root, p, m)

            # Calculate the acceleration of the body
            acceleration = resultingForce / body['mass']

            # Calculate the new velocity of the body
            body['velocity'].append(
                body['velocity'][-1] + (acceleration * time_step_size))

        """
        """
        for body in current_conditions:
            # Calculate the new position of the body
            body['position'].append(
                body['position'][-1] + body['velocity'][-1] * time_step_size)

        print(str(i) + ' / ' +
              str(time_steps))

    return current_conditions, boundary_boxes


def insertInTree(node, body_position, body_mass):
    global current_boxes
    # If node x does not contain a body, put the new body b here.
    if node.mass is None:
        node.mass = copy.deepcopy(body_mass)
        node.center_of_mass = copy.deepcopy(body_position)
        return

    # If node x is an internal node, update the center-of-mass and total mass of x. Recursively insert the body b in the appropriate quadrant.
    elif node.children is not None:
        node.center_of_mass = copy.deepcopy((node.center_of_mass * node.mass +
                                             body_position * body_mass) / (node.mass + body_mass))
        node.mass += copy.deepcopy(body_mass)
        octant = get_octant(node.bbox, copy.deepcopy(body_position))
        if node.children[octant] is None:
            node.children[octant] = Node()
            node.children[octant].bbox = find_bbox(node.bbox, octant)
            current_boxes.append(node.children[octant].bbox)
        insertInTree(node.children[octant], copy.deepcopy(
            body_position), copy.deepcopy(body_mass))
        return

    # If node x is an external node, say containing a body named c, then there are two bodies b and c in the same region. Subdivide the region further by creating four children. Then, recursively insert both b and c into the appropriate quadrant(s). Since b and c may still end up in the same quadrant, there may be several subdivisions during a single insertion. Finally, update the center-of-mass and total mass of x.
    elif node.children is None:
        node.children = [None, None, None, None, None, None, None, None]

        # get the octant of the current node
        old_octant = get_octant(node.bbox, node.center_of_mass)
        # get the octant of the new body
        new_octant = get_octant(node.bbox, body_position)

        node.children[old_octant] = Node()
        node.children[old_octant].bbox = find_bbox(node.bbox, old_octant)

        current_boxes.append(node.children[old_octant].bbox)

        if new_octant != old_octant:
            node.children[new_octant] = Node()
            node.children[new_octant].bbox = find_bbox(node.bbox, new_octant)
            current_boxes.append(node.children[new_octant].bbox)

        # insert the old body in the appropriate quadrant
        insertInTree(node.children[old_octant], copy.deepcopy(
            node.center_of_mass), copy.deepcopy(node.mass))
        # insert the new body in the appropriate quadrant
        insertInTree(node.children[new_octant], copy.deepcopy(
            body_position), copy.deepcopy(body_mass))

        # update the center of mass and mass of the current node
        node.center_of_mass = copy.deepcopy((node.center_of_mass * node.mass +
                                             body_position * body_mass) / (node.mass + body_mass))
        node.mass += copy.deepcopy(body_mass)
        return


def calculateForce(node, body_position, mass):
    if node.mass == mass:
        return np.array([0, 0, 0]) * u.N
    # If the current node is an external node (and it is not body b), calculate the force exerted by the current node on b, and add this amount to b’s net force.
    elif node.children is None:
        resultingForce = calculateGravity(
            node.center_of_mass, body_position, node.mass, mass)
        return resultingForce

    else:
        # Otherwise, calculate the ratio s/d.

        # size of the square boundary box
        s = np.linalg.norm(node.bbox[1][0] - node.bbox[0][0])
        # distance between the center of mass and the body
        d = np.linalg.norm(body_position - node.center_of_mass).value
        # If s/d < θ, treat this internal node as a single body, and calculate the force it exerts on body b, and add this amount to b’s net force.
        if s/d < theta:
            resultingForce = calculateGravity(
                node.center_of_mass, body_position, node.mass, mass)
            return resultingForce
        else:
            # Otherwise, run the procedure recursively on each of the curren t node’s children.
            resultingForce = np.array([0, 0, 0]) * u.N

            for child in node.children:
                if child is not None:
                    resultingForce += calculateForce(child,
                                                     body_position, mass)
            return resultingForce


########################## Helper functions##########################


def calculateGravity(other_body_position, body_position, other_body_mass, body_mass):
    # Calculate the connection vector between the bodies
    connectionVector = other_body_position - body_position

    # Calculate the length of the connection vector
    distance = np.linalg.norm(connectionVector)

    # Normalize the connection vector
    direction = connectionVector / distance

    # Calculate the gravitational force between the bodies
    force = G * (body_mass * other_body_mass) / (distance**2)

    # Calculate the resultant force
    resultingForce = force * direction

    return resultingForce


def find_root_bbox(array_of_positions):
    # Fin the smallest and largest x, y and z values

    min_x = min(array_of_positions, key=lambda x: x[0])[0]
    max_x = max(array_of_positions, key=lambda x: x[0])[0]
    min_y = min(array_of_positions, key=lambda x: x[1])[1]
    max_y = max(array_of_positions, key=lambda x: x[1])[1]
    min_z = min(array_of_positions, key=lambda x: x[2])[2]
    max_z = max(array_of_positions, key=lambda x: x[2])[2]

    # Find the largest difference between the smallest and largest values
    max_diff = max(max_x - min_x, max_y - min_y, max_z - min_z)

    # Find the center of the root bounding box
    center = np.array([min_x, min_y, min_z]) + \
        np.array([max_x-min_x, max_y-min_y, max_z-min_z])/2

    # Fin the bottom left corner of the root bounding box
    cmin = center - np.array([max_diff, max_diff, max_diff])/2
    cmax = center + np.array([max_diff, max_diff, max_diff])/2

    return [cmin, cmax]

    '''
    ooo = np.array([min_x, min_y, min_z])
    ooi = np.array([min_x, min_y, max_z])
    oio = np.array([min_x, max_y, min_z])
    oii = np.array([min_x, max_y, max_z])
    ioo = np.array([max_x, min_y, min_z])
    ioi = np.array([max_x, min_y, max_z])
    iio = np.array([max_x, max_y, min_z])
    iii = np.array([max_x, max_y, max_z])


    return [ooo, ooi, oio, oii, ioo, ioi, iio, iii]

    '''

    '''
    cmin = np.array([min_x, min_y, min_z])
    cmax = np.array([max_x, max_y, max_z])
    '''

    # smallest cube boundary that contains all the bodies


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

        # boundary of the octant cube
