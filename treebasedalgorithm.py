import numpy as np


class Node:
    def __init__(self, position, mass, size, children=None):
        self.position = position
        self.mass = mass
        self.size = size
        self.children = children or []

    def is_leaf(self):
        return not self.children


def simulate_n_bodies(bodies, time_steps, time_step_size, theta=0.5):
    def build_tree(bodies, depth=0):
        n = len(bodies)
        if n <= 1:
            return bodies[0] if n else None

        if depth % 3 == 0:
            axis = 0
        elif depth % 3 == 1:
            axis = 1
        else:
            axis = 2

        bodies.sort(key=lambda b: b['position'][0][axis])
        mid = n // 2
        left = build_tree(bodies[:mid], depth + 1)
        right = build_tree(bodies[mid:], depth + 1)
        mid_body = bodies[mid]
        size = np.abs(left.position - right.position).max()
        position = (left.mass * left.position + right.mass *
                    right.position) / (left.mass + right.mass)

        return Node(position, left.mass + right.mass, size, [left, right])

    def update_force(node, body, theta, force):
        if node.is_leaf() or node.size / np.linalg.norm(node.position - body['position'][0]) < theta:
            displacement = node.position - body['position'][0]
            force += node.mass * displacement / np.linalg.norm(displacement)**3
        else:
            for child in node.children:
                update_force(child, body, theta, force)

    def update_forces(root, bodies, theta):
        for body in bodies:
            force = np.zeros(3)
            update_force(root, body, theta, force)
            body['velocity'].append(
                body['velocity'][-1] + force * time_step_size)
            body['position'].append(
                body['position'][-1] + body['velocity'][-1] * time_step_size)

    for t in range(time_steps):
        root = build_tree(bodies)
        update_forces(root, bodies, theta)
    return bodies
