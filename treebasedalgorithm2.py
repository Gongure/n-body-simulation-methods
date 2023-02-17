from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
x, y, z = [1, 1.5, 2, 3], [1, 2.4, 4, 2], [3.4, 1.4, 2, 3]
ax.scatter(x, y, z, c='red', s=100)
ax.plot(x, y, z, color='black')
plt.show()
