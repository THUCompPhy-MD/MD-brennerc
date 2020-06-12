import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# matplotlib.use("Qt5Agg")

index = 0

with open("xmol.d") as coord:
    while coord.readable():
        points = []
        head = coord.readline()
        if head == "":
            break
        basic_data = coord.readline()
        NP, _, _, _ = basic_data.split()
        coord.readline()
        cube = coord.readline()
        cubex, cubey, cubez = cube.split()
        for i in range(0, int(NP)):
            line = coord.readline()
            _, _, x, y, z = line.split()
            points.append([float(x), float(y), float(z)])

        points = np.array(points)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(points[:, 0], points[:, 1], points[:, 2])
        ax.set_xlim([-float(cubex) / 2, float(cubex) / 2])
        ax.set_ylim([-float(cubey) / 2, float(cubey) / 2])
        ax.set_zlim([-float(cubez) / 2, float(cubez) / 2])
        plt.savefig("fig/sim2_{}".format(index))
        plt.close(fig)
        index = index + 1
