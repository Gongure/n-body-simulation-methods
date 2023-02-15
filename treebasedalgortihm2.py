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

theta = 0.5

# Calculate the gravitational constant
G = 6.674e-11 * u.m**3 * u.kg**-1 * u.s**-2


class Node:
    children = None
    mass = None
    center_of_mass = None
    # boundary box
    bbox = None


def insertInTree(node, body_position, body_mass):
    # If node x does not contain a body, put the new body b here.
    if node.mass is None:
        node.mass = body_mass
        node.center_of_mass = body_position
        return
    # If node x is an internal node, update the center-of-mass and total mass of x. Recursively insert the body b in the appropriate quadrant.
    elif node.children is not None:
        node.mass += body_mass
        node.center_of_mass = (node.center_of_mass * node.mass +
                               body_position * body_mass) / (node.mass + body_mass)
        octand = get_octand(node.bbox, body_position)
        insertInTree(node.children[octand], body_position, body_mass)
        return
    # If node x is an external node, say containing a body named c, then there are two bodies b and c in the same region. Subdivide the region further by creating four children. Then, recursively insert both b and c into the appropriate quadrant(s). Since b and c may still end up in the same quadrant, there may be several subdivisions during a single insertion. Finally, update the center-of-mass and total mass of x.
    elif node.children is None:
        node.children = [Node(), Node(), Node(), Node(),
                         Node(), Node(), Node(), Node()]

        # get the octand of the current node
        old_octand = get_octand(node.bbox, node.center_of_mass)
        # get the octand of the new body
        new_octand = get_octand(node.bbox, body_position)

        # update the center of mass and mass of the current node
        node.mass += body_mass
        node.center_of_mass = (node.center_of_mass * node.mass +
                               body_position * body_mass) / (node.mass + body_mass)

        # insert the old body in the appropriate quadrant
        insertInTree(node.children[old_octand], node.center_of_mass, node.mass)
        # insert the new body in the appropriate quadrant
        insertInTree(node.children[new_octand], body_position, body_mass)
        return


def calculateForce(node, body_position, mass):
    # If the current node is an external node (and it is not body b), calculate the force exerted by the current node on b, and add this amount to b’s net force.
    if node.children is None and node.mass != mass:

        resultingForce = calculateForce(
            node.center_of_mass, body_position, node.mass, mass)
        return resultingForce

    else:
        # Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal node as a single body, and calculate the force it exerts on body b, and add this amount to b’s net force.
        # size of the square boundary box
        s = np.linalg.norm(node.bbox[1] - node.bbox[0])
        # distance between the center of mass and the body
        d = np.linalg.norm(body_position - node.center_of_mass)
        if s/d < theta:
            resultingForce = calculateForce(
                node.center_of_mass, body_position, node.mass, mass)
            return resultingForce
        else:
            # Otherwise, run the procedure recursively on each of the current node’s children.
            resultingForce = np.array([0, 0, 0]) * u.N
            for child in node.children:
                resultingForce += calculateForce(child, body_position, mass)
            return resultingForce


def treeBasedAlgorithm(time_steps, time_step_size, initial_conditions):
    # tim_step_size defines the time between each iteration

    # Setup the initial conditions
    current_conditions = initial_conditions

    # Iterate over the time steps
    for i in range(time_steps):

        # Construct the tree
        root = Node()
        # root.center_of_mass = []
        root.bbox = find_root_bbox([object['position'][-1]
                                   for body in current_conditions])

        for body in current_conditions:
            insertInTree(root, body['position'][-1], body['mass'])

        for body in current_conditions:

            resultingForce = np.array([0, 0, 0]) * u.N

            reultingForce = calculateForce(
                root, body['position'][-1], body['mass'])

            # Calculate the acceleration of the body
            acceleration = resultingForce / body['mass']

            # Calculate the new velocity of the body
            body['velocity'].append(
                body['velocity'][-1] + (acceleration * time_step_size))

        for body in current_conditions:
            # Calculate the new position of the body
            body['position'].append(
                body['position'][-1] + body['velocity'][-1] * time_step_size)

        print(str(i) + ' / ' +
              str(time_steps))

    return current_conditions


########################## Helper functions##########################

def calculateForce(other_body_position, body_position, other_body_mass, body_mass):
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
