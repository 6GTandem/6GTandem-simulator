import os

import matplotlib.pyplot as plt
import numpy as np
import yaml
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations

dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")


class Simulator:
    def __init__(self):
        self.env = None
        self.waveform = None
        self.config = None

    def load_env(self, env_file):
        self.env = Environment(env_file)

    def load_config(self, config_file):
        config_path = os.path.join(dir_path, "configurations")

        # Read the YAML file
        with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
            self.config = yaml.safe_load(file)

    def plot_env(self):
        self.env.plot()

class Environment:
    def __init__(self, config_file):
        config_path = os.path.join(dir_path, "environments")

        # Read the YAML file
        with open(os.path.join(config_path, config_file), "r", encoding="utf8") as file:
            self.config = yaml.safe_load(file)

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_aspect("auto")

        # PLOT THE ROOM

        points = []

        x_size = self.config["room"]["x"]
        y_size = self.config["room"]["y"]
        z_size = self.config["room"]["z"]

        for x in [0, x_size]:
            for y in [0, y_size]:
                for z in [0, z_size]:
                    points.append([x, y, z])

        # Convert to array
        points = np.array(points)

        # Plot
        ax.scatter3D(points[:, 0], points[:, 1], points[:, 2], c="black")
        for s, e in combinations(points, 2):
            # 2 coordinates should be the same for straight lines
            diff = list(s - e)
            if diff.count(0) == 2:
                ax.plot3D(*zip(s, e), color="k")

        # PLOT THE STRIPES
        for stripe in self.config["radio_stripes"]:
            units = []
            for unit in stripe:
                if "x" in unit:
                    units.append([unit["x"], unit["y"], unit["z"]])
            units = np.array(units)
            ax.scatter3D(units[:, 0], units[:, 1], units[:, 2], c="red")
            ax.plot3D(units[:, 0], units[:, 1], units[:, 2], color="red", label="Stripe")

        # PLOT THE USERS
        for user in self.config["users"]:
            ax.scatter3D(user["x"], user["y"], user["z"], c="blue")

        # TODO add labels
        plt.show()
