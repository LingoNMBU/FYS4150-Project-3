# !python
# -*- coding: utf-8 -*

__author__ = 'Erling Ween Eriksen'
__email__ = 'erlinge@nmbu.no'

import pandas as pd
import matplotlib.pyplot as plt

pos_ana = pd.read_csv('positions_ana', header=None)
pos_eu = pd.read_csv('positions_euler', header=None)
pos_rk = pd.read_csv('positions_rk4', header=None)

print(pos_eu)
print(pos_ana)
print(pos_rk)

pos_ana = pos_ana.transpose()
pos_eu = pos_eu.transpose()
pos_rk = pos_rk.transpose()
pos_ana.columns = ['x', 'y', 'z']
pos_eu.columns = ['x', 'y', 'z']
pos_rk.columns = ['x', 'y', 'z']

print(pos_eu)

plot = True

if plot:

    plt.plot(pos_ana.x, label='x analytical', marker='^')
    plt.plot(pos_rk.x, label='x runge kutta 4', marker='o')
    plt.plot(pos_eu.x, label='x euler', marker='o')
    plt.legend()
    plt.show()

    plt.plot(pos_rk.y, label='y runge kutta 4', marker='o')
    plt.plot(pos_eu.y, label='y euler', marker='o')
    plt.plot(pos_ana.y, label='y analytical', marker='^')
    plt.legend()
    plt.show()

    plt.plot(pos_eu.z, label='z euler', marker='o')
    plt.plot(pos_ana.z, label='z analytical', marker='^')
    plt.plot(pos_rk.z, label='z runge kutta 4', marker='o')
    plt.legend()
    plt.show()

    orb = plt.figure()
    axes = orb.add_axes([0.1, 0.1, 0.8, 0.8])
    axes.scatter(pos_eu.x, pos_eu.y, label='euler', marker='o')
    axes.scatter(pos_ana.x, pos_ana.y, label='analytical', marker='^')
    axes.scatter(pos_rk.x, pos_rk.y, label='runge kutta 4', marker='o')
    # axes.scatter(0, 0, label='0,0')
    plt.axis('equal')
    plt.legend()
    plt.xlabel('x-coordinate')
    plt.ylabel('y-coordinate')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(pos_eu.x, pos_eu.y, pos_eu.z, marker='o')
    ax.scatter(pos_ana.x, pos_ana.y, pos_ana.z, marker='o')

    # Set the axis labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # Rotate the axes and update
    #ax.view_init(90, 90, 90)
    plt.show()

sum_err = (pos_ana.x - pos_eu.x) + (pos_ana.y - pos_eu.y) + (pos_ana.z - pos_eu.z)

print(sum_err)
