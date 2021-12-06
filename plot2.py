import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from plane import plane2xy

def interp(pointx, pointy, values):
    x = np.linspace(-1, 1, 500)
    y =  np.linspace(-1, 1, 500)
    X, Y = np.meshgrid(x, y)
    data =  griddata((pointx, pointy), values, (X, Y), method='linear')
    return X, Y, data

def draw_stereonet(x, y, z, ax, fig, title, stress_state, colormap=None, levels=None, directions=None):
    sigma1 = stress_state.orientation.sigma1
    sigma2 = stress_state.orientation.sigma2
    sigma3 = stress_state.orientation.sigma3

    p = ax.contourf(*interp(x, y, z), cmap=colormap or 'Oranges', levels=levels)
    ax.set_aspect('equal', 'box')
    ax.set_title(title)
    ax.axis('off')
    stereonet_border = plt.Circle((0, 0), 1, color='black', fill=False)
    ax.scatter(*plane2xy(sigma1), marker="s", color='black')
    ax.scatter(*plane2xy(sigma2), marker="o", color='black')
    ax.scatter(*plane2xy(sigma3), marker="^", color='black')
    ax.text(*plane2xy(sigma1), r'$\sigma1$', fontsize=10)
    ax.text(*plane2xy(sigma3), r'$\sigma3$', fontsize=10)
    if directions:
        for x1, y1, x2, y2 in zip(*directions):
            ax.plot((x1, x2), (y1, y2), color='black')
            ax.scatter(x1, y1, marker="o", color='black')


    ax.add_patch(stereonet_border)
    fig.colorbar(p, ax=ax)

def prepare_lists_of_data(stresses_on_plane, stress_state, mu_fracture):
    xx = []
    yy = []
    taus = []
    snns = []
    taus_reduced = []
    snns_reduced = []
    slip_tendency = []
    dilation_tendency = []
    fracture_susceptibility = []

    s1 = stress_state.values.sigma1
    s2 = stress_state.values.sigma2
    s3 = stress_state.values.sigma3

    for stress_on_plane in stresses_on_plane:
        x, y = plane2xy(stress_on_plane.plane)
        s_nn = stress_on_plane.s_nn
        s_nm = stress_on_plane.s_nm
        s_nt = stress_on_plane.s_nt
        tau_n = stress_on_plane.tau_n
        xx.append(x)
        yy.append(y)
        taus.append(stress_on_plane.tau_n)
        taus_reduced.append(stress_on_plane.tau_n_reduced)
        snns_reduced.append(stress_on_plane.s_nn_reduced)
        snns.append(stress_on_plane.s_nn)
        slip_tendency.append(tau_n/s_nn)
        dilation_tendency.append((s1-s_nn)/(s1-s3))
        fracture_susceptibility.append(s_nn - tau_n/mu_fracture)
    return xx, yy, taus, snns, taus_reduced, snns_reduced, slip_tendency, dilation_tendency, fracture_susceptibility


def plot(stresses_on_plane, stress_state):
    x, y, taus, snns, taus_reduced, snns_reduced, slip_tendency, dilation_tendency, fracture_susceptibility = prepare_lists_of_data(stresses_on_plane, stress_state, mu_fracture=0.5)

    fig, axs = plt.subplots(2, 3)

    # axs[0, 0].scatter(snns_reduced, taus_reduced)
    axs[0, 0].scatter(snns, taus)
    axs[0, 0].set_aspect('equal', 'box')
    axs[0, 0].set_title("Morh diagram")
    axs[0, 0].set_xlabel(r"$\sigma_{nn}$")
    axs[0, 0].set_ylabel(r"$\tau_{n}$")

    # axs[2, 1].scatter(snns_reduced, taus_reduced)
    # axs[2, 1].set_aspect('equal', 'box')

     
    draw_stereonet(x, y, snns_reduced, ax=axs[0, 1], fig=fig, title=r'$\sigma_{nn}$', stress_state=stress_state)
    draw_stereonet(x, y, taus_reduced, ax=axs[0, 2], fig=fig, title=r'$\tau_{n}$', stress_state=stress_state)
    draw_stereonet(x, y, slip_tendency, ax=axs[1, 0], fig=fig, title='Slip Tendency', stress_state=stress_state)
    draw_stereonet(x, y, dilation_tendency, ax=axs[1, 1], fig=fig, title='Dilation Tendency', stress_state=stress_state)
    draw_stereonet(x, y, fracture_susceptibility, ax=axs[1, 2], fig=fig, title='Fracture Susceptibility', stress_state=stress_state)
    # draw_stereonet(x, y, phis, axs[2, 0], r'$\phi$', colormap='hsv', levels=180)


    fig.suptitle(r'$\mu_s$ = %.2f, $\phi$ = %.2f' % (stress_state.values.mu_s, stress_state.values.phi), fontsize=16)
    fig.tight_layout()

    plt.show()

