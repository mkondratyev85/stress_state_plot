import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from plane import plane2xy
from calculate import fracture_criteria_reduced

from matplotlib.widgets import Slider, Button

def interp(pointx, pointy, values):
    x = np.linspace(-1, 1, 500)
    y =  np.linspace(-1, 1, 500)
    X, Y = np.meshgrid(x, y)
    data =  griddata((pointx, pointy), values, (X, Y), method='linear')
    # data = np.nan_to_num(data)
    return X, Y, data

class Plot:

    output = None
    output_morh = None
    gui_update_func = None

    def __call__(self, stresses_on_plane, stress_state, output=None, stresses_on_fractures=None, gui=False, gui_update_func=None):
        if output:
            self.output = output
            self.output_morh = output[:-4] + "_morh.jpg"
        self.replot(stresses_on_plane, stresses_on_fractures, stress_state)
        if gui:
            self.gui_update_func = gui_update_func
            self.prepare_controls()
            plt.plot()


    def replot(self, stresses_on_plane, stresses_on_fractures, stress_state, init_figs=True):
        self.prepare_everything(stresses_on_plane, stresses_on_fractures, stress_state)
        self.plot_stereo(init_figs)
        self.plot_morh(init_figs)

    def prepare_everything(self, stresses_on_plane, stresses_on_fractures, stress_state):
        self.stresses_on_plane = stresses_on_plane
        self.stresses_on_fractures = stresses_on_fractures
        self.stress_state = stress_state
        (self.xx, self.yy, self.taus, self.snns, self.taus_reduced, self.snns_reduced, 
            self.slip_tendency, self.dilation_tendency, self.fracture_susceptibility, self.fracture_criteria) = \
            self.prepare_lists_of_data(mu_fracture=0.5, planes=self.stresses_on_plane)

        self.fractures_xx = None
        self.fractures_yy = None

        if self.stresses_on_fractures:
            (self.fractures_xx, self.fractures_yy, self.fracture_taus, self.fracture_snns,
                self.fracture_taus_reduced, self.fracture_snns_reduced, _, _, _, _) = \
                self.prepare_lists_of_data(mu_fracture=0.5, planes=self.stresses_on_fractures)


    def prepare_lists_of_data(self, mu_fracture, planes):
        xx = []
        yy = []
        taus = []
        snns = []
        taus_reduced = []
        snns_reduced = []
        slip_tendency = []
        dilation_tendency = []
        fracture_susceptibility = []
        fracture_criteria = []

        s1 = self.stress_state.values.sigma1
        s2 = self.stress_state.values.sigma2
        s3 = self.stress_state.values.sigma3

        self.sigma1 = self.stress_state.orientation.sigma1
        self.sigma2 = self.stress_state.orientation.sigma2
        self.sigma3 = self.stress_state.orientation.sigma3

        for stress_on_plane in planes:
            x, y = plane2xy(stress_on_plane.plane)
            s_nn = stress_on_plane.s_nn
            s_nm = stress_on_plane.s_nm
            s_nt = stress_on_plane.s_nt
            tau_n = stress_on_plane.tau_n
            xx.append(x)
            yy.append(y)
            taus.append(stress_on_plane.tau_n)
            snns.append(stress_on_plane.s_nn)
            taus_reduced.append(stress_on_plane.tau_n_reduced)
            snns_reduced.append(stress_on_plane.s_nn_reduced)
            fracture_criteria.append(fracture_criteria_reduced(stress_on_plane, tau_f=0, k_f=0.5))
            slip_tendency.append(tau_n/s_nn)
            dilation_tendency.append((s1-s_nn)/(s1-s3))
            fracture_susceptibility.append(s_nn - tau_n/mu_fracture)
        return xx, yy, taus, snns, taus_reduced, snns_reduced, slip_tendency, dilation_tendency, fracture_susceptibility, fracture_criteria

    def prepare_controls(self):
        self.fig_controls, axs_controls = plt.subplots(2, 3, dpi=300, figsize=(12,6))

        self.slider_mu_sigma = Slider(
            ax=axs_controls[0,0],
            label="Amplitude",
            valmin=-1,
            valmax=1,
            valinit=0,
            orientation="horizontal"
        )
        self.slider_mu_sigma.on_changed(self.update_plot)

        plt.show()

    def update_plot(self, val):
        mu_s = self.slider_mu_sigma.val
        stress_state, stresses_on_plane = self.gui_update_func(mu_s)
        self.replot(stresses_on_plane, stresses_on_fractures=None, stress_state=stress_state, init_figs=False)

        self.fig.canvas.draw_idle()

    def plot_morh(self, init_figs=True):
        if init_figs:
            self.fig_morh, self.morh_axs = plt.subplots(1, 2, dpi=300, figsize=(12,6))

        self.morh_axs[0].scatter(self.snns, self.taus)
        self.morh_axs[0].set_aspect('equal', 'box')
        self.morh_axs[0].set_title("Morh diagram")
        self.morh_axs[0].set_xlabel(r"$\sigma_{nn}$")
        self.morh_axs[0].set_ylabel(r"$\tau_{n}$")
        if self.fractures_xx and self.fractures_yy:
            self.morh_axs[0].scatter(self.fracture_snns, self.fracture_taus, marker="+", color='black')

        tau_f = 0
        k_f = 0.8
        y1, y2 = 0, 1
        x1 = tau_f
        x2 = (y2 - tau_f) / k_f


        self.morh_axs[1].scatter(self.snns_reduced, self.taus_reduced)
        self.morh_axs[1].set_aspect('equal', 'box')
        self.morh_axs[1].set_title("Morh diagram in reduced stress")
        self.morh_axs[1].set_xlabel(r"$\sigma_{nn}$")
        # self.morh_axs[1].set_ylabel(r"$\tau_{n}$")
        if self.fractures_xx and self.fractures_yy:
            self.morh_axs[1].scatter(self.fracture_snns_reduced, self.fracture_taus_reduced, marker="+", color='black')
        self.morh_axs[1].plot((x1, x2),(y1, y2), color='red')

        self.fig_morh.suptitle(r'$\mu_s$ = %.2f, $\phi$ = %.2f' % (self.stress_state.values.mu_s, self.stress_state.values.phi), fontsize=16)
        self.fig_morh.tight_layout()

        if self.output_morh:
            plt.savefig(self.output_morh)


    def plot_stereo(self, init_figs=True):
        if init_figs:
            self.fig, self.fig_axs = plt.subplots(2, 3, dpi=300, figsize=(12,6))

        self.draw_stereonet(self.snns_reduced, ax=self.fig_axs[0, 1], title=r'$\sigma_{nn}$')
        self.draw_stereonet(self.taus_reduced, ax=self.fig_axs[0, 2], title=r'$\tau_{n}$')
        self.draw_stereonet(self.slip_tendency, ax=self.fig_axs[1, 0], title='Slip Tendency')
        self.draw_stereonet(self.dilation_tendency, ax=self.fig_axs[1, 1], title='Dilation Tendency')
        self.draw_stereonet(self.fracture_susceptibility, ax=self.fig_axs[1, 2], title='Fracture Susceptibility')
        self.draw_stereonet(self.fracture_criteria, ax=self.fig_axs[0, 0], title='Fracture Susceptibility')
        # draw_stereonet(x, y, phis, axs[2, 0], r'$\phi$', colormap='hsv', levels=180)


        self.fig.suptitle(r'$\mu_s$ = %.2f, $\phi$ = %.2f' % (self.stress_state.values.mu_s, self.stress_state.values.phi), fontsize=16)
        self.fig.tight_layout()

        if self.output:
            plt.savefig(self.output)

    def draw_stereonet(self, z, ax, title, colormap=None, levels=None, directions=None):

        p = ax.contourf(*interp(self.xx, self.yy, z), cmap=colormap or 'Oranges', levels=levels)
        ax.set_aspect('equal', 'box')
        ax.set_title(title)
        ax.axis('off')
        stereonet_border = plt.Circle((0, 0), 1, color='black', fill=False)
        ax.scatter(*plane2xy(self.sigma1), marker="s", color='red')
        ax.scatter(*plane2xy(self.sigma2), marker="o", color='red')
        ax.scatter(*plane2xy(self.sigma3), marker="^", color='red')
        ax.text(*plane2xy(self.sigma1), r'$\sigma1$', fontsize=10)
        ax.text(*plane2xy(self.sigma3), r'$\sigma3$', fontsize=10)
        if self.fractures_xx and self.fractures_yy:
            ax.scatter(self.fractures_xx, self.fractures_yy, marker="+", color='black')

        if directions:
            for x1, y1, x2, y2 in zip(*directions):
                ax.plot((x1, x2), (y1, y2), color='black')
                ax.scatter(x1, y1, marker="o", color='black')

        ax.add_patch(stereonet_border)
        # self.fig.colorbar(p, ax=ax, spacing='proportional')

plot = Plot()
