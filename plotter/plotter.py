from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations

def plot_room(config:dict):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_aspect("auto")

    # PLOT THE ROOM

    points = []

    x_size = config["room"]["x"]
    y_size = config["room"]["y"]
    z_size = config["room"]["z"]

    ax.set_xlim([0, config["room"]["x"]])
    ax.set_ylim([0, config["room"]["y"]])
    ax.set_zlim([0, config["room"]["z"]])

    for x in [0, x_size]:
        for y in [0, y_size]:
            for z in [0, z_size]:
                points.append([x,y,z])

    # Convert to array
    points = np.array(points)

    # Plot
    ax.scatter3D(points[:, 0], points[:, 1], points[:, 2],c="black")
    for s, e in combinations(points, 2):
        # 2 coordinates should be the same for straight lines
        diff = list(s - e)
        if diff.count(0) == 2:
            ax.plot3D(*zip(s, e), color="k")

    # PLOT THE STRIPES
    for stripe in config["radio_stripes"]:
        print(f"{len(config['radio_stripes'])} Stripes found")
        units = []
        for unit in stripe:
            unit = unit["radio_unit"]
            if "x" in unit:
                units.append([unit["x"], unit["y"], unit["z"]])
        units = np.array(units)
        print(units)
        ax.scatter3D(units[:, 0], units[:, 1], units[:, 2], c="red")
        ax.plot3D(units[:, 0], units[:, 1], units[:, 2], color="red")

    # PLOT THE UEs
    ues = []
    for ue_pos in config["ue_positions"]:
        print(f"{len(config['ue_positions'])} UEs found")
        if "x" in ue_pos:
            ues.append([ue_pos["x"], ue_pos["y"], ue_pos["z"]])
    ues = np.array(ues)
    print(ues)
    ax.scatter3D(ues[:, 0], ues[:, 1], ues[:, 2], c="blue", alpha=0.2)

    plt.show()
