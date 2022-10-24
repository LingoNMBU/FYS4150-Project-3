# !python
# -*- coding: utf-8 -*

__author__ = 'Erling Ween Eriksen'
__email__ = 'erlinge@nmbu.no'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def process_data(path):
    df = pd.read_csv(path, header=None)
    df = df.transpose()
    df.columns = ['x', 'y', 'z']
    df.index = 50 * df.index / len(df.index)
    return df


def plot(x, y):
    return 0


plot1 = False

if plot1:
    pos_rk = process_data('positions_single_rk4')

    f1, ax1 = plt.subplots(figsize=(5, 4))
    ax1.set_ylim(-25, +22)
    ax1.plot(pos_rk.z, label='z runge kutta 4', marker='o')
    plt.xlabel('time [microseconds]')
    plt.ylabel(f'z [micrometers]')
    #plt.savefig('single_particle_z.pdf')
    ax1.legend()
    plt.show()

# Two particles no particle interactions
plot2 = False

if plot2:
    pos_rk_d1 = process_data('r_double_rk4_1')
    pos_rk_d2 = process_data('r_double_rk4_2')

    fig, ax = plt.subplots(figsize=(5, 4))
    # ax.set_ylim(-25, +22)
    ax.scatter(pos_rk_d1.x, pos_rk_d1.y, label='particle 1', marker='o')
    ax.scatter(pos_rk_d2.x, pos_rk_d2.y, label='particle 2', marker='o')
    plt.xlabel('x [micrometers]')
    plt.ylabel(f'y [micrometers]')
    #plt.savefig('two_particles_no_int_xy.pdf')
    ax.legend()
    plt.show()

    pos_rk_d1 = process_data('r_double_rk4_1_int')
    pos_rk_d2 = process_data('r_double_rk4_2_int')

    fig2, ax = plt.subplots(figsize=(5, 4))
    # ax.set_ylim(-25, +22)
    ax.scatter(pos_rk_d1.x, pos_rk_d1.y, label='particle 1', marker='o')
    ax.scatter(pos_rk_d2.x, pos_rk_d2.y, label='particle 2', marker='o')
    plt.xlabel('x [micrometers]')
    plt.ylabel(f'y [micrometers]')
    #plt.savefig('two_particles_with_int_xy.pdf')
    ax.legend()
    plt.show()

# plot4 = False
#
# if plot4:
#     pos_rk_d1 = process_data('r_double_rk4_1_int')
#     pos_rk_d2 = process_data('r_double_rk4_2_int')
#     vel_rk_d1 = process_data('v_double_rk4_1_int')
#     vel_rk_d2 = process_data('v_double_rk4_2_int')
#
#     fig, ax = plt.subplots(figsize=(5, 4))
#     # ax.set_ylim(-25, +22)
#     ax.scatter(pos_rk_d2.x, vel_rk_d2.x, label='particle 2 x ', marker='o')
#     ax.scatter(pos_rk_d2.z, vel_rk_d2.z, label='particle 2 z ', marker='o')
#     plt.xlabel('r [micrometers]')
#     plt.ylabel(f'v [m/s]')
#     plt.savefig('particle2_with_int_rv.pdf')
#     ax.legend()
#     plt.show()
#
#     fig, ax = plt.subplots(figsize=(5, 4))
#     # ax.set_ylim(-25, +22)
#     ax.scatter(pos_rk_d1.x, vel_rk_d1.x, label='particle 1 x ', marker='^')
#     ax.scatter(pos_rk_d1.z, vel_rk_d1.z, label='particle 1 z ', marker='^')
#     plt.xlabel('r [micrometers]')
#     plt.ylabel(f'v [m/s]')
#     plt.savefig('particle1_with_int_rv.pdf')
#     ax.legend()
#     plt.show()

plot5 = False

if plot5:
    pos_rk_d1 = process_data('r_double_rk4_1')
    pos_rk_d2 = process_data('r_double_rk4_2')
    vel_rk_d1 = process_data('v_double_rk4_1')
    vel_rk_d2 = process_data('v_double_rk4_2')
    pos_rk_d1_int = process_data('r_double_rk4_1_int')
    pos_rk_d2_int = process_data('r_double_rk4_2_int')
    vel_rk_d1_int = process_data('v_double_rk4_1_int')
    vel_rk_d2_int = process_data('v_double_rk4_2_int')

    fig, ax = plt.subplots(figsize=(5, 4))
    # ax.set_ylim(-25, +22)
    ax.scatter(pos_rk_d1_int.x, vel_rk_d1_int.x, label='particle 1 with interactions ', marker='o')
    ax.scatter(pos_rk_d2_int.x, vel_rk_d2_int.x, label='particle 2 with interactions ', marker='o')
    plt.xlabel('x [micrometers]')
    plt.ylabel(f'v_x [m/s]')
    #plt.savefig('particle2_with_int_rv.pdf')
    ax.legend()
    plt.show()

    fig, ax = plt.subplots(figsize=(5, 4))
    # ax.set_ylim(-25, +22)
    ax.scatter(pos_rk_d1_int.z, vel_rk_d1_int.z, label='particle 1 with interactions ', marker='^')
    ax.scatter(pos_rk_d2_int.z, vel_rk_d2_int.z, label='particle 2 with interactions ', marker='^')
    plt.xlabel('z [micrometers]')
    plt.ylabel(f'v_z [m/s]')
    #plt.savefig('particle1_no_int_rv.pdf')
    ax.legend()
    plt.show()


    fig, ax = plt.subplots(figsize=(5, 4))
    # ax.set_ylim(-25, +22)
    ax.scatter(pos_rk_d1.x, vel_rk_d1.x, label='particle 1 ', marker='o')
    ax.scatter(pos_rk_d2.x, vel_rk_d2.x, label='particle 2 ', marker='o')
    plt.xlabel('x [micrometers]')
    plt.ylabel(f'v_x [m/s]')
    #plt.savefig('particle2_no_int_rv.pdf')
    ax.legend()
    plt.show()

    fig, ax = plt.subplots(figsize=(5, 4))
    # ax.set_ylim(-25, +22)
    ax.scatter(pos_rk_d1.z, vel_rk_d1.z, label='particle 1  ', marker='^')
    ax.scatter(pos_rk_d2.z, vel_rk_d2.z, label='particle 2 ', marker='^')
    plt.xlabel('z [micrometers]')
    plt.ylabel(f'v_z [m/s]')
    #plt.savefig('particle1_no_int_rv.pdf')
    ax.legend()
    plt.show()

plot6 = False

if plot6:
    pos_rk_d1 = process_data('r_double_rk4_1')
    pos_rk_d2 = process_data('r_double_rk4_2')

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(pos_rk_d1.x, pos_rk_d1.y, pos_rk_d1.z, marker='o', label='Particle 1')
    ax.scatter(pos_rk_d2.x, pos_rk_d2.y, pos_rk_d2.z, marker='o', label='Particle 2')

    # Set the axis labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # Rotate the axes and update
    # ax.view_init(90, 90, 90)
    plt.show()

plot7 = False

if plot7:
    pos_rk_d1 = process_data('r_double_rk4_1_int')
    pos_rk_d2 = process_data('r_double_rk4_2_int')

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(pos_rk_d1.x, pos_rk_d1.y, pos_rk_d1.z, marker='o', label='Particle 1')
    ax.scatter(pos_rk_d2.x, pos_rk_d2.y, pos_rk_d2.z, marker='o', label='Particle 2')

    # Set the axis labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.legend()

    # Rotate the axes and update
    # ax.view_init(90, 90, 90)
    plt.show()

plot7 = False

if plot7:
    pos_eu = process_data('rs_eu2_n_32000')
    pos_ana = process_data('rs_ana2_n_32000')
    pos_rk = process_data('rs_rk2_n_32000')

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.plot(pos_rk.index, pos_ana.x, label='analytical', marker='o')
    ax.plot(pos_rk.index, pos_eu.x, label='forward euler', marker='o')
    ax.plot(pos_rk.index, pos_rk.x, label='runge kutta 4', marker='o')

    # Set the axis labels
    ax.set_xlabel('t')
    ax.set_ylabel('x')
    plt.legend()

    # Rotate the axes and update
    # ax.view_init(90, 90, 90)
    plt.show()

plot8 = True
if plot8:
    particles_df = process_data('trapped1')


    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot()
    n_steps = round((2.5 - 0.2) / 0.02)
    trap_index = np.linspace(0.2, 2.5, n_steps)
    ax.plot(trap_index, particles_df.x,label='f = 0.1', marker='o')
    ax.plot(trap_index, particles_df.y, label='f = 0.4', marker='o')
    ax.plot(trap_index, particles_df.z, label='f = 0.7', marker='o')




    # Set the axis labels
    ax.set_xlabel('w_v')
    ax.set_ylabel('Trapped Particles')
    plt.legend()
    plt.savefig('trapped_particles_wide2.pdf')

    # Rotate the axes and update
    # ax.view_init(90, 90, 90)
    plt.show()

plot8 = True
if plot8:
    particles_df = process_data('trapped1')


    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot()
    n_steps = round((2.5 - 0.2) / 0.02)
    trap_index = np.linspace(0.2, 2.5, n_steps)
    ax.plot(trap_index, particles_df.x,label='f = 0.1', marker='o')
    ax.plot(trap_index, particles_df.y, label='f = 0.4', marker='o')
    ax.plot(trap_index, particles_df.z, label='f = 0.7', marker='o')




    # Set the axis labels
    ax.set_xlabel('w_v')
    ax.set_ylabel('Trapped Particles')
    plt.legend()
    plt.savefig('trapped_particles_wide2.pdf')

    # Rotate the axes and update
    # ax.view_init(90, 90, 90)
    plt.show()

