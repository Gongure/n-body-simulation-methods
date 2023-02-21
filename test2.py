
class Node:
    mass = None


def myfunc(c):
    b = Node()
    b.mass = c
    b.mass += 1


a = [1, 4, 2]
print(a)
myfunc(a[0])
print(a)
