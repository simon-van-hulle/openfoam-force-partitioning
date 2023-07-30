#!/usr/bin/env python3

"""
This plots the analytic solution for the harmonic function (Zhang, 2015). This
Auxiliary potential should be calculated from
.. math::
    \nabla^2 \phi^{(i)} = 0 with
"""

import matplotlib.pyplot as plt
import numpy as np

def harm_func(x, y, R=1):
    r = np.sqrt(x*x + y*y)
    theta = np.arctan2(y, x) + np.pi / 2

    return (r > R) * -1 * R*R*R * np.cos(theta) / (2 * r*r)


def main():
    fig, ax = plt.subplots()
    # ax.add_patch(plt.Circle((0, 0), 1))

    N = 1000
    x = np.linspace(-4, 4, N)
    y = np.linspace(-4, 4, N)
    xx, yy = np.meshgrid(x, y)
    PhiI = harm_func(xx, yy)
    #  print(PhiI)
    print("Plotting PhiI analytic solution")

    plt.pcolor(xx, yy, PhiI, cmap="BrBG")


    ax.axis("equal")
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    main()
