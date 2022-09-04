import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.widgets import Slider, Button
import numpy as np
from scipy.interpolate import griddata

from .plane import plane2xy
from .calculate import fracture_criteria_reduced
from .entities import StressState
import matplotlib.colors as colors


def interp(pointx, pointy, values):
    x = np.linspace(-1, 1, 500)
    y = np.linspace(-1, 1, 500)
    X, Y = np.meshgrid(x, y)
    data = griddata((pointx, pointy), values, (X, Y), method="linear")
    return X, Y, data


class Plot:

    output = None
    output_morh = None
    gui_update_func = None
    gui_array_resolution = None
    fractures = None

    def __call__(
        self,
        stresses_on_plane,
        stress_state,
        k_f,
        tau_f,
        output=None,
        stresses_on_fractures=None,
        gui=False,
        gui_update_func=None,
    ):
        if output:
            self.output = output
            self.output_morh = output[:-4] + "_morh.png"
        if gui:
            self.gui_update_func = gui_update_func

        self.k_f = k_f
        self.tau_f = tau_f

        shape = int(len(stresses_on_plane) ** 0.5)
        self.gui_array_resolution = (shape, shape)
        self.make_plot(stresses_on_plane, stresses_on_fractures, stress_state)
        self.planes = [p.plane if p else None for p in stresses_on_plane]
        if stresses_on_fractures:
            self.fractures = [p.plane if p else None for p in stresses_on_fractures]
        if gui:
            self.prepare_controls(stress_state)
            plt.plot()

    def make_plot(self, stresses_on_plane, stresses_on_fractures, stress_state):
        self.prepare_everything(stresses_on_plane, stresses_on_fractures, stress_state)
        if self.gui_update_func:
            self.plot_gui()
        else:
            self.plot_stereo()
            self.plot_morh()

    def prepare_everything(
        self, stresses_on_plane, stresses_on_fractures, stress_state
    ):
        self.stresses_on_plane = stresses_on_plane
        self.stresses_on_fractures = stresses_on_fractures
        self.stress_state = stress_state
        (
            self.xx,
            self.yy,
            self.taus,
            self.snns,
            self.taus_reduced,
            self.snns_reduced,
            self.slip_tendency,
            self.dilation_tendency,
            self.fracture_susceptibility,
            self.fracture_criteria,
        ) = self.prepare_lists_of_data(mu_fracture=self.k_f, planes=self.stresses_on_plane)

        self.fractures_xx = None
        self.fractures_yy = None

        if self.stresses_on_fractures:
            (
                self.fractures_xx,
                self.fractures_yy,
                self.fracture_taus,
                self.fracture_snns,
                self.fracture_taus_reduced,
                self.fracture_snns_reduced,
                _,
                _,
                _,
                _,
            ) = self.prepare_lists_of_data(
                mu_fracture=self.k_f,
                planes=self.stresses_on_fractures,
                fractures=True,
            )

    def prepare_lists_of_data(self, mu_fracture, planes, fractures=False):
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
            x, y = None, None
            if stress_on_plane:
                if fractures:
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
                fracture_criteria.append(
                    fracture_criteria_reduced(stress_on_plane, tau_f=self.tau_f, k_f=self.k_f)
                )
                slip_tendency.append(tau_n / s_nn)
                dilation_tendency.append((s1 - s_nn) / (s1 - s3))
                fracture_susceptibility.append(s_nn - tau_n / mu_fracture)
            else:
                xx.append(np.nan)
                yy.append(np.nan)
                taus.append(np.nan)
                snns.append(np.nan)
                taus_reduced.append(np.nan)
                snns_reduced.append(np.nan)
                fracture_criteria.append(np.nan)
                slip_tendency.append(np.nan)
                dilation_tendency.append(np.nan)
                fracture_susceptibility.append(np.nan)

        return (
            xx,
            yy,
            taus,
            snns,
            taus_reduced,
            snns_reduced,
            slip_tendency,
            dilation_tendency,
            fracture_susceptibility,
            fracture_criteria,
        )

    def prepare_controls(self, stress_state):
        slider_width = 0.25
        slider_height = 0.03

        self.slider_sigma1_dir = Slider(
            ax=plt.axes([0.05, 0.6, slider_width, slider_height]),
            label=r"$\sigma1 \alpha$",
            valmin=0,
            valmax=360,
            valinit=stress_state.orientation.sigma1.dir,
            orientation="horizontal",
        )
        self.slider_sigma1_dir.on_changed(self.update_plot)

        self.slider_sigma1_dip = Slider(
            ax=plt.axes([0.05, 0.52, slider_width, slider_height]),
            label=r"$\sigma1 \beta$",
            valmin=0,
            valmax=360,
            valinit=stress_state.orientation.sigma1.dip,
            orientation="horizontal",
        )
        self.slider_sigma1_dip.on_changed(self.update_plot)

        self.slider_sigma3_dir = Slider(
            ax=plt.axes([0.05, 0.44, slider_width, slider_height]),
            label=r"$\sigma3 \alpha$",
            valmin=0,
            valmax=360,
            valinit=stress_state.orientation.sigma3.dir,
            orientation="horizontal",
        )
        self.slider_sigma3_dir.on_changed(self.update_plot)

        self.slider_mu_sigma = Slider(
            ax=plt.axes([0.05, 0.3, slider_width, slider_height]),
            label=r"$\mu_{\sigma}$",
            valmin=-1,
            valmax=1,
            valinit=stress_state.values.mu_s,
            orientation="horizontal",
        )
        self.slider_mu_sigma.on_changed(self.update_plot)

        self.slider_tau = Slider(
            ax=plt.axes([0.05, 0.22, slider_width, slider_height]),
            label=r"$\tau$",
            valmin=0,
            valmax=100,
            valinit=stress_state.values.tau,
            orientation="horizontal",
        )
        self.slider_tau.on_changed(self.update_plot)

        self.slider_p = Slider(
            ax=plt.axes([0.05, 0.14, slider_width, slider_height]),
            label=r"$p$",
            valmin=-100,
            valmax=0,
            valinit=stress_state.values.p,
            orientation="horizontal",
        )
        self.slider_p.on_changed(self.update_plot)

        self.slider_tau_f = Slider(
            ax=plt.axes([0.05, 0.78, slider_width, slider_height]),
            label=r"$\tau_f$",
            valmin=0,
            valmax=1,
            valinit=self.tau_f,
            orientation="horizontal",
        )
        self.slider_tau_f.on_changed(self.update_plot)

        self.slider_k_f = Slider(
            ax=plt.axes([0.05, 0.86, slider_width, slider_height]),
            label=r"$k_f$",
            valmin=0.01,
            valmax=1,
            valinit=self.k_f,
            orientation="horizontal",
        )
        self.slider_k_f.on_changed(self.update_plot)

        plt.show()

    def update_plot(self, val):
        self.k_f = self.slider_k_f.val
        self.tau_f = self.slider_tau_f.val
        stress_state = StressState(
            orientations={
                "sigma1_sigma3": (
                    self.slider_sigma1_dir.val,
                    self.slider_sigma1_dip.val,
                    self.slider_sigma3_dir.val,
                )
            },
            values={
                "tau": self.slider_tau.val,
                "mu_s": self.slider_mu_sigma.val,
                "p": self.slider_p.val,
            },
        )
        stresses_on_plane = self.gui_update_func(
            stress_state=stress_state,
            planes=self.planes,
        )
        if self.fractures:
            stresses_on_fractures = self.gui_update_func(
                stress_state=stress_state,
                planes=self.fractures,
            )
        else:
            stresses_on_fractures = None
        self.prepare_everything(stresses_on_plane, stresses_on_fractures, stress_state)
        self.update_morh()
        self.update_stereo(stress_state)

        self.fig.canvas.draw_idle()

    def reshape_for_imshow(self, z):
        z = np.array(z).reshape(self.gui_array_resolution)
        return np.fliplr(np.flipud(z.T))

    def update_countourf(self, z, im, stress_state):
        im, s1, s2, s3, t1, t3 = im
        im.set_data(self.reshape_for_imshow(z))

        x, y = plane2xy(stress_state.orientation.sigma1)
        s1.set_offsets(np.c_[x, y])
        t1.set_position((x, y))

        x, y = plane2xy(stress_state.orientation.sigma2)
        s2.set_offsets(np.c_[x, y])

        x, y = plane2xy(stress_state.orientation.sigma3)
        s3.set_offsets(np.c_[x, y])
        t3.set_position((x, y))

    def update_stereo(self, stress_state):
        self.update_countourf(self.snns_reduced, self.snn_cont, stress_state)
        self.update_countourf(self.taus_reduced, self.taus_cont, stress_state)
        self.update_countourf(self.fracture_criteria, self.fc, stress_state)

    def update_morh(self):
        self.morh_reduced_scatter.set_offsets(
            np.c_[self.snns_reduced, self.taus_reduced]
        )

        self.fracture_reduced_scatter.set_offsets(
            np.c_[
                self.fracture_snns_reduced,
                self.fracture_taus_reduced,
            ]
        )

        x1, x2, y1, y2 = self.get_morh_line_coordinates()
        self.morh_line.set_ydata((y1, y2))
        self.morh_line.set_xdata((x1, x2))

    def get_morh_line_coordinates(self):
        tau_f = self.tau_f
        k_f = self.k_f
        y1, y2 = 0, 1
        x1 = (y1 - tau_f)/k_f
        x2 = (y2 - tau_f)/k_f
        return x1, x2, y1, y2
 
    def plot_single_morh(
        self, ax, x, y, title="Morh diagram in reduced stress", plot_line=False
    ):
        scatter = ax.scatter(self.snns_reduced, self.taus_reduced)
        ax.set_aspect("equal", "box")
        ax.set_title(title)
        ax.set_xlabel(r"$\sigma_{nn}$")
        ax.set_ylabel(r"$\tau_{n}$")
        ax.set_xlim((-1.5, 1.5))

        fracture_scatter = None
        if self.fractures_xx and self.fractures_yy:
            fracture_scatter = ax.scatter(
                self.fracture_snns_reduced,
                self.fracture_taus_reduced,
                marker="+",
                color="black",
            )
        if plot_line:
            x1, x2, y1, y2 = self.get_morh_line_coordinates()
            self.morh_line, = ax.plot((x1, x2), (y1, y2), color="red")
        return scatter, fracture_scatter

    def plot_morh(self, init_figs=True):
        if init_figs:
            self.fig_morh, self.morh_axs = plt.subplots(1, 2, dpi=300, figsize=(12, 6))

        self.morh_scatter = self.plot_single_morh(
            self.morh_axs[0], self.snns, self.taus, title="Morh diagram"
        )
        self.morh_reduced_scatter = self.plot_single_morh(
            self.morh_axs[1],
            self.snns_reduced,
            self.taus_reduced,
            title="Morh diagram in reduced stress",
            plot_line=True,
        )

        self.fig_morh.suptitle(
            r"$\mu_s$ = %.2f, $\phi$ = %.2f"
            % (self.stress_state.values.mu_s, self.stress_state.values.phi),
            fontsize=16,
        )
        self.fig_morh.tight_layout()

        if self.output_morh:
            plt.savefig(self.output_morh)

    def plot_gui(self):
        self.fig, self.fig_axs = plt.subplots(2, 3, dpi=100, figsize=(18, 9))

        self.snn_cont = self.draw_stereonet(
            self.snns_reduced,
            ax=self.fig_axs[0, 1],
            title=r"$\sigma_{nn}$",
            colormap="PuOr",
            norm_zero=True,
        )
        self.taus_cont = self.draw_stereonet(
            self.taus_reduced,
            ax=self.fig_axs[0, 2],
            title=r"$\tau_{n}$",
            colormap="Oranges",
            norm_zero=False,
        )

        self.fc = self.draw_stereonet(
            self.fracture_criteria,
            ax=self.fig_axs[1, 2],
            title="Fracture Criteria",
            colormap="Binary",
            norm_zero=False,
        )

        (
            self.morh_reduced_scatter,
            self.fracture_reduced_scatter,
        ) = self.plot_single_morh(
            self.fig_axs[1, 1],
            self.snns_reduced,
            self.taus_reduced,
            title="Morh diagram in reduced stress",
            plot_line=True,
        )

        self.fig.delaxes(self.fig_axs[0, 0])
        self.fig.delaxes(self.fig_axs[1, 0])

        self.fig.tight_layout()

    def plot_stereo(self):
        self.fig, self.fig_axs = plt.subplots(2, 3, dpi=300, figsize=(12, 6))

        self.snn_cont = self.draw_stereonet(
            self.snns_reduced, ax=self.fig_axs[0, 1], title=r"$\sigma_{nn}$",
            use_contourf=True,
        )
        self.taus_cont = self.draw_stereonet(
            self.taus_reduced,
            ax=self.fig_axs[0, 2],
            title=r"$\tau_{n}$",
            colormap="Oranges",
            norm_zero=False,
            use_contourf=True,
        )
        self.slip_cont = self.draw_stereonet(
            self.slip_tendency,
            ax=self.fig_axs[1, 0],
            title="Slip Tendency",
            colormap="Oranges",
            norm_zero=False,
            use_contourf=True,
        )
        self.dilation_cont = self.draw_stereonet(
            self.dilation_tendency,
            ax=self.fig_axs[1, 1],
            title="Dilation Tendency",
            colormap="Oranges",
            norm_zero=False,
            use_contourf=True,
        )
        self.fs_cont = self.draw_stereonet(
            self.fracture_susceptibility,
            ax=self.fig_axs[1, 2],
            title="Fracture Susceptibility",
            colormap="Oranges",
            norm_zero=False,
            use_contourf=True,
        )
        self.fc = self.draw_stereonet(
            self.fracture_criteria,
            ax=self.fig_axs[0, 0],
            title="Fracture Criteria",
            colormap="Binary",
            norm_zero=False,
            use_contourf=True,
        )
        # draw_stereonet(x, y, phis, axs[2, 0], r'$\phi$', colormap='hsv', levels=180)

        self.fig.suptitle(
            r"$\mu_s$ = %.2f, $\phi$ = %.2f"
            % (self.stress_state.values.mu_s, self.stress_state.values.phi),
            fontsize=16,
        )

        self.fig.tight_layout()

        if self.output:
            plt.savefig(self.output)

    def draw_stereonet(
        self,
        z,
        ax,
        title,
        colormap=None,
        levels=None,
        directions=None,
        norm_zero=True,
        use_contourf=False,
        show_cbar=True,
    ):

        cmap = colormap or "PuOr"
        array = self.reshape_for_imshow(z)

        if cmap == "Binary":
            cmap = matplotlib.colors.ListedColormap(['white', 'orange'])
            show_cbar = False

        if use_contourf:
            splot = ax.contourf(
                array,
                cmap=cmap,
                levels=levels,
                norm=None if not norm_zero else colors.CenteredNorm(),
                extent=[-1, 1, -1, 1],
                origin="upper",
            )
        else:
            splot = ax.imshow(
                array,
                cmap,
                norm=None if not norm_zero else colors.CenteredNorm(),
                extent=[-1, 1, -1, 1],
            )
        ax.set_aspect("equal", "box")
        ax.set_title(title)
        ax.axis("off")
        stereonet_border = plt.Circle((0, 0), 1, color="black", fill=False)
        s1 = ax.scatter(*plane2xy(self.sigma1), marker="s", color="red")
        s2 = ax.scatter(*plane2xy(self.sigma2), marker="o", color="red")
        s3 = ax.scatter(*plane2xy(self.sigma3), marker="^", color="red")
        t1 = ax.text(*plane2xy(self.sigma1), r"$\sigma1$", fontsize=10)
        t3 = ax.text(*plane2xy(self.sigma3), r"$\sigma3$", fontsize=10)
        if self.fractures_xx and self.fractures_yy:
            ax.scatter(self.fractures_xx, self.fractures_yy, marker="+", color="black")

        if directions:
            for x1, y1, x2, y2 in zip(*directions):
                ax.plot((x1, x2), (y1, y2), color="black")
                ax.scatter(x1, y1, marker="o", color="black")

        ax.add_patch(stereonet_border)
        if show_cbar:
            self.fig.colorbar(splot, ax=ax, spacing="proportional")
        return splot, s1, s2, s3, t1, t3


plot = Plot()
