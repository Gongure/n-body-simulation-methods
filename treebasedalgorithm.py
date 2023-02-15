from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import random
import astropy.units as u

theta = 0.5
AU = (149.6e6 * 1000)     # 149.6 million km, in meters.
G = 6.674e-11 * u.m**3 * u.kg**-1 * u.s**-2
fig1 = plt.figure()
sim = fig1.add_subplot(111, aspect='equal')
fig2 = plt.figure()
quadt = fig2.add_subplot(111, aspect='equal')


class Node:
    children = None
    mass = None
    center_of_mass = np.array([0, 0, 0]) * u.au
    bbox = None
    velocity = np.array([0, 0, 0]) * u.au / u.day


def oct_insert(root, position, mass):  # octree
    if root.mass is None:  # when the root is empty, add the first particle
        root.mass = mass
        root.center_of_mass = position
        return
    elif root.children is None:
        root.children = [None, None, None, None, None, None, None, None]
        old_quadrant = octant_of_particle(
            root.bbox, root.center_of_mass)
        if root.children[old_quadrant] is None:
            # can root.children[old_quadrant] not be None?
            # if root.children is None, then root.children[old_quadrant] is Non
            # the if statement is redundant
            root.children[old_quadrant] = Node()
            root.children[old_quadrant].bbox = octant_bbox(
                root.bbox, old_quadrant)
        oct_insert(root.children[old_quadrant],
                   root.center_of_mass, root.mass)
        new_quadrant = octant_of_particle(root.bbox, position)
        if root.children[new_quadrant] is None:
            root.children[new_quadrant] = Node()
            root.children[new_quadrant].bbox = octant_bbox(
                root.bbox, new_quadrant)
        oct_insert(root.children[new_quadrant], position, mass)
        root.center_of_mass = (root.center_of_mass *
                               root.mass + position * mass) / (root.mass + mass)
        root.mass = root.mass + mass
    else:
        new_quadrant = octant_of_particle(root.bbox, position)
        if root.children[new_quadrant] is None:
            root.children[new_quadrant] = Node()
            root.children[new_quadrant].bbox = octant_bbox(
                root.bbox, new_quadrant)
        oct_insert(root.children[new_quadrant], position, mass)
        root.center_of_mass = (root.center_of_mass *
                               root.mass + position * mass) / (root.mass + mass)
        root.mass = root.mass + mass


def display(root):
    if root.mass is None:
        return
    if root.children is not None:
        x = (root.bbox[0] + root.bbox[1]) / 2
        y = (root.bbox[2] + root.bbox[3]) / 2
        width = x-root.bbox[0]
        plt_node(root.bbox[0], root.bbox[2], width)
        plt_node(root.bbox[0], y, width)
        plt_node(x, root.bbox[2], width)
        plt_node(x, y, width)
        for i in range(4):
            if root.children[i] is not None:
                display(root.children[i])
    else:
        quadt.scatter(root.center_of_mass[0], root.center_of_mass[1])


def integrate(time_steps, time_step_size, bodies):
    results = bodies

    # do a version where results isnt overwrittens
    for i in range(time_steps):
        particles_force = {}
        root = Node()
        root.center_of_mass = []
        root.bbox = find_root_bbox(bodies)
        for body in bodies:
            oct_insert(root, body['position'][-1], body['mass'])
        for body in bodies:
            total_force = compute_force(
                root, body['position'][-1], body['mass'])
            particles_force[body['name']] = total_force
        for body in bodies:
            force = particles_force[body['name']]
            acceleration = force / body['mass']

            body['velocity'].append(
                body['velocity'][-1] + (acceleration * time_step_size))

            body['position'].append(
                body['position'][-1] + body['velocity'][-1] * time_step_size)
        print(str(i) + " / " + str(time_steps))
    return results


def compute_force(root, position, m):
    if root.mass is None:
        return np.array([0, 0, 0]) * u.N
    if root.center_of_mass.value.all() == position.value.all() and root.mass == m:
        return np.array([0, 0, 0]) * u.N
    d = root.bbox[1] - root.bbox[0]

    r = position - root.center_of_mass
    r = np.linalg.norm(r)

    if d/r < theta or root.children is None:
        return force(m, position, root.mass, root.center_of_mass)
    else:
        f = np.array([0, 0, 0]) * u.N

        for i in range(8):
            if root.children[i] is not None:
                f += compute_force(root.children[i], position, m)
        return f

################################################# SUPPORTING FUNCTION ##############################################################


def force(m, position, mcm, pcm):
    d = np.linalg.norm(position - pcm)
    direction = (position - pcm) / d
    f = (G * m * mcm) / (d**2)
    f = f * direction.value

    return f


def plt_node(x, y, width):
    quadt.add_patch(patches.Rectangle((x, y), width, width, fill=False))


def find_root_bbox(array):
    """ Create a suitable square boundary box for the input particles
    """
    if len(array) == 0 or len(array) == 1:
        return None

    # create a search algorithm to find the min and max of x, y, and z stored in the array as a np.array([x,y,z]
    # (hint: use np.min and np.max)
    xmin = xmax = ymin = ymax = zmin = zmax = 0

    for body in array:
        if body['position'][-1][0] > xmax:
            xmax = body['position'][-1][0]
        if body['position'][-1][0] < xmin:
            xmin = body['position'][-1][0]
        if body['position'][-1][1] > ymax:
            ymax = body['position'][-1][1]
        if body['position'][-1][1] < ymin:
            ymin = body['position'][-1][1]
        if body['position'][-1][2] > zmax:
            zmax = body['position'][-1][2]
        if body['position'][-1][2] < zmin:
            zmin = body['position'][-1][2]
    # make sure the boundary box is a square
    if xmax - xmin == ymax - ymin == zmax - zmin:
        return xmin, xmax, ymin, ymax, zmin, zmax
    elif xmax - xmin > ymax - ymin == zmax - zmin:
        return xmin, xmax, ymin, ymax+(xmax-xmin-ymax+ymin), zmin, zmax+(xmax-xmin-zmax+zmin)
    elif xmax - xmin == ymax - ymin > zmax - zmin:
        return xmin, xmax+(ymax-ymin-xmax+xmin), ymin, ymax, zmin, zmax+(ymax-ymin-zmax+zmin)
    elif xmax - xmin == zmax - zmin > ymax - ymin:
        return xmin, xmax+(zmax-zmin-xmax+xmin), ymin, ymax, zmin, zmax
    elif ymax - ymin == zmax - zmin > xmax - xmin:
        return xmin, xmax, ymin, ymax+(zmax-zmin-ymax+ymin), zmin, zmax
    elif xmax - xmin > ymax - ymin and xmax - xmin > zmax - zmin:
        return xmin, xmax, ymin, ymax+(xmax-xmin-ymax+ymin), zmin, zmax+(xmax-xmin-zmax+zmin)
    elif ymax - ymin > xmax - xmin and ymax - ymin > zmax - zmin:
        return xmin, xmax+(ymax-ymin-xmax+xmin), ymin, ymax, zmin, zmax+(ymax-ymin-zmax+zmin)
    else:
        return xmin, xmax+(zmax-zmin-xmax+xmin), ymin, ymax+(zmax-zmin-ymax+ymin), zmin, zmax


def octant_of_particle(bbox, position):
    # retun the octant of the particle position = np.array([x,y,z])+
    x = position[0]
    y = position[1]
    z = position[2]

    if z >= (bbox[5] + bbox[4])/2:
        if y >= (bbox[3] + bbox[2])/2:
            if x <= (bbox[1] + bbox[0])/2:
                return 0
            else:
                return 1
        else:
            if x >= (bbox[1] + bbox[0])/2:
                return 2
            else:
                return 3
    else:
        if y >= (bbox[3] + bbox[2])/2:
            if x <= (bbox[1] + bbox[0])/2:
                return 4
            else:
                return 5
        else:
            if x >= (bbox[1] + bbox[0])/2:
                return 6
            else:
                return 7


def octant_bbox(bbox, octant):
    """Return the coordinate of the octant
    """
    x = (bbox[0] + bbox[1]) / 2
    y = (bbox[2] + bbox[3]) / 2
    z = (bbox[4] + bbox[5]) / 2

    # Octant 0: (xmin, x, y, ymax, zmin, z)
    if octant == 0:
        return bbox[0], x, y, bbox[3], bbox[4], z
    # Octant 1: (x, xmax, y, ymax, zmin, z)
    elif octant == 1:
        return x, bbox[1], y, bbox[3], bbox[4], z
    # Octant 2: (x, xmax, ymin, y, zmin, z)
    elif octant == 2:
        return x, bbox[1], bbox[2], y, bbox[4], z
    # Octant 3: (xmin, x, ymin, y, zmin, z)
    elif octant == 3:
        return bbox[0], x, bbox[2], y, bbox[4], z
    # Octant 4: (xmin, x, y, ymax, z, zmax)
    elif octant == 4:
        return bbox[0], x, y, bbox[3], z, bbox[5]
    # Octant 5: (x, xmax, y, ymax, z, zmax)
    elif octant == 5:
        return x, bbox[1], y, bbox[3], z, bbox[5]
    # Octant 6: (x, xmax, ymin, y, z, zmax)
    elif octant == 6:
        return x, bbox[1], bbox[2], y, z, bbox[5]
    # Octant 7: (xmin, x, ymin, y, z, zmax)
    elif octant == 7:
        return bbox[0], x, bbox[2], y, z, bbox[5]


def treeBasedAlgorithm(time_steps, time_step_size, data):
    results = integrate(time_steps, time_step_size, data)
    return results
