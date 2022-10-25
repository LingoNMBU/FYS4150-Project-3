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

def rel_error(u, v):
    """
    Calculates relative error
    :param u: exact value
    :param v: model value
    :return:
    """
    errors = np.zeros(len(u))
    k = 0
    for i, j in zip(u, v):
        if i - j == 0:
            errors[k] = 0
            k += 1
        else:
            relres = (i - j) / i
            errors[k] = np.log(abs((relres)))[:]
            k += 1
    return errors



plot1 = False

if plot1:
    pos_rk = process_data('positions_single_rk4')

    f1, ax1 = plt.subplots(figsize=(4.5,3))
    ax1.set_ylim(-25, +22)
    ax1.plot(pos_rk.z, label='z runge kutta 4', color='g', linestyle='dashed')
    plt.xlabel('time [μs]')
    plt.ylabel(f'z\n[μm]', rotation=0)
    plt.savefig('single_particle_z.pdf', bbox_inches="tight")
    #ax1.YAxis.set_units('μm')
    #ax1.XAxis.set_units('μs')
    #ax1.legend()
    plt.tight_layout()
    plt.show()

# Two particles no particle interactions
plot2 = False

if plot2:
    pos_rk_d1 = process_data('r_double_rk4_1')
    pos_rk_d2 = process_data('r_double_rk4_2')
    pos_rk_d1_int = process_data('r_double_rk4_1_int')
    pos_rk_d2_int = process_data('r_double_rk4_2_int')

    #fig, ax = plt.subplots(nrows=1, ncols=2, sharex='none', sharey='none', squeeze=True, subplot_kw=None, gridspec_kw=None)
    fig = plt.figure(figsize=(5.5,3))
    #plt.xlim(60, 60)
    #plt.ylim(60, 60)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, sharey=ax1, sharex=ax1)
    plt.setp([ax1, ax2], xlabel=f'x [μm]')
    plt.setp([ax2], ylabel=f'y[μm]')
    # ax.set_ylim(-25, +22)
    ax1.scatter(pos_rk_d1.x, pos_rk_d1.y, marker='o', s=0)
    ax1.scatter(pos_rk_d2.x, pos_rk_d2.y, marker='o', s=0)
    ax1.plot(pos_rk_d1.x, pos_rk_d1.y, label='particle 1', linestyle='solid')
    ax1.plot(pos_rk_d2.x, pos_rk_d2.y, label='particle 2', linestyle='solid')
    ax1.yaxis.set_label(f'y\n[μm]')
    ax1.xaxis.set_label(f'x [μm]')
    #plt.xlabel(f'x[μm]', rotation=0)

    ax2.scatter(pos_rk_d1_int.x, pos_rk_d1_int.y, marker='o', s=0, vmin = 55, vmax = 55)
    ax2.scatter(pos_rk_d2_int.x, pos_rk_d2_int.y, marker='o', s=0, vmin = 55, vmax = 55)
    ax2.plot(pos_rk_d1_int.x, pos_rk_d1_int.y, label='particle 1', linestyle='solid')
    ax2.plot(pos_rk_d2_int.x, pos_rk_d2_int.y, label='particle 2', linestyle='solid')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    #fig.text(0.5, 0.04, 'x [μm]', ha='center')
    plt.ylabel(f'y\n[μm]', rotation=0)
    #plt.xlabel(f'x[μm]', rotation=0)
    plt.savefig('two_particles_xy.pdf', bbox_inches="tight")
    #ax2.legend()

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

    fig = plt.figure(figsize=(6, 6))
    # plt.xlim(60, 60)
    # plt.ylim(60, 60)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222, sharey=ax1, sharex=ax1)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224, sharey=ax3, sharex=ax3)
    axs = [ax1, ax2, ax3, ax4]
    plt.setp([ax2, ax1], xlabel=f'x [μm]', ylabel = f'v_x')
    plt.setp([ax3, ax4], xlabel=f'z [μm]', ylabel = f'v_z')
    ax1.xaxis.set_label(f'x [μm]')
    ax2.yaxis.tick_right()
    ax4.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax4.yaxis.set_label_position("right")

    ax2.scatter(pos_rk_d1_int.x, vel_rk_d1_int.x, label='particle 1 with interactions ', marker='o', s=0.1)
    ax2.scatter(pos_rk_d2_int.x, vel_rk_d2_int.x, label='particle 2 with interactions ', marker='o', s = 0.1)

    ax4.scatter(pos_rk_d1_int.z, vel_rk_d1_int.z, label='particle 1 with interactions ', marker='^', s= 0.1)
    ax4.scatter(pos_rk_d2_int.z, vel_rk_d2_int.z, label='particle 2 with interactions ', marker='^', s = 0.1)

    ax1.scatter(pos_rk_d1.x, vel_rk_d1.x, label='particle 1 ', marker='o', s=0.1)
    ax1.scatter(pos_rk_d2.x, vel_rk_d2.x, label='particle 2 ', marker='o', s=0.1)

    ax3.scatter(pos_rk_d1.z, vel_rk_d1.z, label='particle 1  ', marker='^', s=0.1)
    ax3.scatter(pos_rk_d2.z, vel_rk_d2.z, label='particle 2 ', marker='^', s=0.1)
    plt.savefig('two_particles_vx_vz.pdf', bbox_inches="tight")

    plt.show()

plot3D = True

if plot3D:
    pos_rk_d1 = process_data('r_double_rk4_1')
    pos_rk_d2 = process_data('r_double_rk4_2')
    pos_rk_d1_int = process_data('r_double_rk4_1_int')
    pos_rk_d2_int = process_data('r_double_rk4_2_int')

    # fig, ax = plt.subplots(nrows=1, ncols=2, sharex='none', sharey='none', squeeze=True, subplot_kw=None, gridspec_kw=None)
    fig = plt.figure(figsize=(3, 6))
    # plt.xlim(60, 60)
    # plt.ylim(60, 60)
    ax1 = fig.add_subplot(211, projection= '3d')
    ax2 = fig.add_subplot(212, sharey=ax1, sharex=ax1, projection= '3d')

    # ax.set_ylim(-25, +22)
    ax1.scatter(pos_rk_d1.x, pos_rk_d1.y, pos_rk_d1.z, marker='o', label='Particle 1', s = 0.1)
    ax1.scatter(pos_rk_d2.x, pos_rk_d2.y, pos_rk_d2.z, marker='o', label='Particle 2', s = 0.1)

    ax2.scatter(pos_rk_d1_int.x, pos_rk_d1_int.y, pos_rk_d1_int.z, marker='o', label='Particle 1', s = 0.1)
    ax2.scatter(pos_rk_d2_int.x, pos_rk_d2_int.y, pos_rk_d2_int.z, marker='o', label='Particle 2', s = 0.1)

    # Set the axis labels
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    # plt.xlabel(f'x[μm]', rotation=0)

    plt.savefig('two_particles_xyz.pdf', bbox_inches="tight")

    #ax2.yaxis.tick_right()
    #ax2.yaxis.set_label_position("right")
    # fig.text(0.5, 0.04, 'x [μm]', ha='center')
    #plt.ylabel(f'y\n[μm]', rotation=0)
    # plt.xlabel(f'x[μm]', rotation=0)

    # ax2.legend()

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
    ns = [4000, 8000, 16000, 32000]
    # fig, ax = plt.subplots(nrows=1, ncols=2, sharex='none', sharey='none', squeeze=True, subplot_kw=None, gridspec_kw=None)
    fig = plt.figure(figsize=(6, 6))
    # plt.xlim(60, 60)
    # plt.ylim(60, 60)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222, sharey=ax1, sharex=ax1)
    ax3 = fig.add_subplot(223, sharey=ax1, sharex=ax1)
    ax4 = fig.add_subplot(224, sharey=ax1, sharex=ax1)
    axs = [ax1, ax2, ax3, ax4]
    plt.yscale('log', base=10)
    plt.setp(axs, xlabel=f't [μs]')
    ax1.xaxis.set_label(f'x [μm]')
    ax2.yaxis.tick_right()
    ax4.yaxis.tick_right()

    # plt.xlabel(f'x[μm]', rotation=0)



    #ax2.yaxis.set_label_position("right")
    # fig.text(0.5, 0.04, 'x [μm]', ha='center')
    #plt.ylabel(f'y\n[μm]', rotation=0)
    # plt.xlabel(f'x[μm]', rotation=0)
    #plt.savefig('two_particles_xy.pdf', bbox_inches="tight")
    # ax2.legend()
    for i in range(4):



        err_rk = pd.read_csv(f'rk_errors_n_{ns[i]}', header=None, names=['error_rk'])
        err_eu = pd.read_csv(f'eu_errors_n_{ns[i]}', header=None, names=['error_eu'])
        err_rk.index = 50 * err_rk.index / len(err_rk.index)
        err_eu.index = err_rk.index

        axs[i].plot(abs(err_eu.error_eu), label='E', linestyle='solid')
        axs[i].plot(abs(err_rk.error_rk), label=f'RK', linestyle='solid')
        axs[i].text(25, 2000, f'{ns[i]} points', horizontalalignment='center')

        axs[i].legend()
    plt.tight_layout()

    plt.savefig('errors_n.pdf', bbox_inches = 'tight')
    plt.show()

plot8 = False
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

plot8 = False
if plot8:
    particles_df = process_data('trapped1')


    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot()
    n_steps = round((2.5 - 0.2) / 0.02)
    trap_index = np.linspace(0.2, 2.5, n_steps)
    ax.plot(trap_index, particles_df.x,label='f = 0.1', marker='o')
    ax.plot(trap_index, particles_df.y, label='f = 0.4', marker='o')
    ax.plot(trap_index, particles_df.z, label='f = 0.7', marker='o')


plot9 = False
if plot9:
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot()

    trapped_fine_df = pd.read_csv('trapped_particles_fine2', header=None)
    trapped_fine_df = trapped_fine_df.transpose()
    trapped_fine_df.columns = ['x', 'y']
    n_steps = 99
    trap_fine_index = np.linspace(0.8, 0.9, n_steps)

    ax.plot(trap_fine_index, trapped_fine_df.x,label='f = No particle interactions', marker='o')
    ax.plot(trap_fine_index, trapped_fine_df.y, label='f = With particle interactions', marker='o')
    # Set the axis labels
    ax.set_xlabel('w_v')
    ax.set_ylabel('Trapped particles')
    plt.legend()
    plt.savefig('trapped_particles_fine.pdf')

    # Rotate the axes and update
    # ax.view_init(90, 90, 90)
    plt.show()

