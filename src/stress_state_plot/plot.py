"""Module for visualisation of stress state via matplotlib.
"""

import numpy as np
from numpy.typing import ArrayLike

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import colors
from matplotlib.text import Text
from matplotlib.image import AxesImage
from matplotlib.collections import PathCollection


from .plane import plane2xy
from .calculate import fracture_criteria as fracture_criteria_func
from .entities import StressState, FrictionState


class Morhplot:
    def __init__(
        self,
        ax,
        snn,
        taus,
        friction_state: FrictionState,
        stress_state: StressState,
        fractures_snns=None,
        fractures_taus=None,
        title: str = "Morh diagram in reduced stress",
        plot_line: bool = False,
    ):

        self.friction_state = friction_state
        self.stress_state = stress_state
        self.ax = ax

        self.init_plot(title)
        self.scatter = self.ax.scatter(snn, taus)

        if fractures_snns and fractures_taus:
            self.fracture_scatter = ax.scatter(
                fractures_snns,
                fractures_taus,
                marker="+",
                color="black",
            )
        if plot_line:
            self.plot_lines()
            self.update_scale()

    def init_plot(self, title: str):
        self.ax.set_aspect("equal", "box")
        self.ax.set_title(title)
        self.ax.set_xlabel(r"$\sigma_{nn}$")
        self.ax.set_ylabel(r"$\tau_{n}$")

    def plot_lines(self):
        x1, x2, y1, y2 = self.brittle_friction_line_coords
        (self.morh_line,) = self.ax.plot((x1, x2), (y1, y2), color="red")
        (self.morh_line_brittle,) = self.ax.plot((x1, x2), (y1, y2), color="red", linestyle="dashed")
        x1, x2, y1, y2 = self.limit_friction_line_coords
        self.morh_line.set_ydata((y1, y2))
        self.morh_line.set_xdata((x1, x2))
        self.ax.plot((0, 0), (0, 50), color="red")

    def update_scale(self):
        x_min, x_max = self.stress_state.snn_limits
        tau = self.stress_state.values.tau
        self.ax.set_xlim((x_min - 0.2 * tau, x_max + 0.2 * tau))
        self.ax.set_ylim(0, 1.2 * tau)

    def update(
        self, snns, taus, fractures_snns, fractures_taus, friction_state: FrictionState, stress_state: StressState,
    ):
        self.scatter.set_offsets(np.c_[snns, taus])
        self.friction_state = friction_state
        self.stress_state = stress_state

        if fractures_snns and fractures_taus:
            self.fracture_scatter.set_offsets(np.c_[fractures_snns, fractures_taus])

        x1, x2, y1, y2 = self.limit_friction_line_coords
        self.morh_line.set_ydata((y1, y2))
        self.morh_line.set_xdata((x1, x2))

        x1, x2, y1, y2 = self.brittle_friction_line_coords
        self.morh_line_brittle.set_ydata((y1, y2))
        self.morh_line_brittle.set_xdata((x1, x2))

        self.update_scale()

    @property
    def limit_friction_line_coords(self):
        tau_f = self.friction_state.tau_f
        k_f = self.friction_state.k_f
        y1, y2 = 0, 2 * self.stress_state.values.tau
        x1 = (y1 - tau_f) / k_f
        x2 = (y2 - tau_f) / k_f
        return x1, x2, y1, y2

    @property
    def brittle_friction_line_coords(self):
        k_f = self.friction_state.k_f
        x1, y1 = 0, 0
        y2 = 2 * self.stress_state.values.tau
        x1 = 0
        x2 = y2 / k_f
        return x1, x2, y1, y2


class Stereonet:
    """Result of drawing stereonet."""

    plot: AxesImage
    s1: PathCollection
    s2: PathCollection
    s3: PathCollection
    t1: Text
    t3: Text

    def __init__(
        self,
        array: ArrayLike,
        stress_state: StressState,
        ax,
        title: str,
        fractures_xx=None,
        fractures_yy=None,
        colormap=None,
        levels=None,
        directions=None,
        norm_zero: bool = True,
        use_contourf: bool = False,
        show_cbar: bool = True,
    ):
        """Draw single stereonet plot."""

        sigma1 = stress_state.orientation.sigma1
        sigma2 = stress_state.orientation.sigma2
        sigma3 = stress_state.orientation.sigma3

        if use_contourf:
            self.plot = ax.contourf(
                array,
                cmap=colormap,
                levels=levels,
                norm=None if not norm_zero else colors.CenteredNorm(),
                extent=[-1, 1, -1, 1],
                origin="upper",
            )
        else:
            self.plot = ax.imshow(
                array,
                cmap=colormap,
                norm=None if not norm_zero else colors.CenteredNorm(),
                extent=[-1, 1, -1, 1],
            )

        ax.set_aspect("equal", "box")
        ax.set_title(title)
        ax.axis("off")
        stereonet_border = plt.Circle((0, 0), 1, color="black", fill=False)

        self.s1 = ax.scatter(*plane2xy(sigma1), marker="s", color="red")
        self.s2 = ax.scatter(*plane2xy(sigma2), marker="o", color="red")
        self.s3 = ax.scatter(*plane2xy(sigma3), marker="^", color="red")
        self.t1 = ax.text(*plane2xy(sigma1), r"$\sigma1$", fontsize=10)
        self.t3 = ax.text(*plane2xy(sigma3), r"$\sigma3$", fontsize=10)

        if fractures_xx and fractures_yy:
            ax.scatter(fractures_xx, fractures_yy, marker="+", color="black")

        if directions:
            for x1, y1, x2, y2 in zip(*directions):
                ax.plot((x1, x2), (y1, y2), color="black")
                ax.scatter(x1, y1, marker="o", color="black")

        ax.add_patch(stereonet_border)
        if show_cbar:
            plt.colorbar(self.plot, ax=ax, spacing="proportional")

    def update(self, array: ArrayLike, stress_state: StressState):
        """Update plot for given stress_state and background array"""
        self.plot.set_data(array)

        x, y = plane2xy(stress_state.orientation.sigma1)
        self.s1.set_offsets(np.c_[x, y])
        self.t1.set_position((x, y))

        x, y = plane2xy(stress_state.orientation.sigma2)
        self.s2.set_offsets(np.c_[x, y])

        x, y = plane2xy(stress_state.orientation.sigma3)
        self.s3.set_offsets(np.c_[x, y])
        self.t3.set_position((x, y))


class Plot:
    """
    Visualisation of stress state in interactive and non-interactive mode.
    """

    output = None
    output_morh = None
    gui_update_func = None
    gui_array_resolution = None
    fractures = None
    k_f = None
    tau_f = None

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
        ) = self.prepare_lists_of_data(
            mu_fracture=self.k_f, planes=self.stresses_on_plane
        )

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
        # s2 = self.stress_state.values.sigma2
        s3 = self.stress_state.values.sigma3

        for stress_on_plane in planes:
            x, y = None, None
            if stress_on_plane:
                if fractures:
                    x, y = plane2xy(stress_on_plane.plane)
                s_nn = stress_on_plane.s_nn
                # s_nm = stress_on_plane.s_nm
                # s_nt = stress_on_plane.s_nt
                tau_n = stress_on_plane.tau_n
                xx.append(x)
                yy.append(y)
                taus.append(stress_on_plane.tau_n)
                snns.append(stress_on_plane.s_nn)
                taus_reduced.append(stress_on_plane.tau_n_reduced)
                snns_reduced.append(stress_on_plane.s_nn_reduced)
                fracture_criteria.append(
                    fracture_criteria_func(
                        stress_on_plane, tau_f=self.tau_f, k_f=self.k_f
                    )
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
            valmax=90,
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
            valmin=0.01,
            valmax=10,
            valinit=stress_state.values.tau,
            orientation="horizontal",
        )
        self.slider_tau.on_changed(self.update_plot)

        self.slider_p = Slider(
            ax=plt.axes([0.05, 0.14, slider_width, slider_height]),
            label=r"$p$",
            valmin=-10,
            valmax=10,
            valinit=stress_state.values.p,
            orientation="horizontal",
        )
        self.slider_p.on_changed(self.update_plot)

        self.slider_tau_f = Slider(
            ax=plt.axes([0.05, 0.78, slider_width, slider_height]),
            label=r"$\tau_f$",
            valmin=0.01,
            valmax=5,
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

        self.snn_cont.update(self.reshape_for_imshow(self.snns_reduced), stress_state)
        self.taus_cont.update(self.reshape_for_imshow(self.taus_reduced), stress_state)
        self.fc.update(self.reshape_for_imshow(self.fracture_criteria), stress_state)

        self.morh_reduced.update(
            snns=self.snns,
            taus=self.taus,
            fractures_snns=self.fracture_snns,
            fractures_taus=self.fracture_taus,
            friction_state=FrictionState(k_f=self.k_f, tau_f=self.tau_f),
            stress_state=self.stress_state,
        )

        self.fig.canvas.draw_idle()

    def reshape_for_imshow(self, z):
        z = np.array(z).reshape(self.gui_array_resolution)
        return np.fliplr(np.flipud(z.T))

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

        self.morh_reduced = Morhplot(
            self.fig_axs[1, 1],
            self.snns,
            self.taus,
            friction_state=FrictionState(k_f=self.k_f, tau_f=self.tau_f),
            stress_state=self.stress_state,
            fractures_snns=self.fracture_snns,
            fractures_taus=self.fracture_taus,
            title="Morh diagram",
            plot_line=True,
        )

        self.fig.delaxes(self.fig_axs[0, 0])
        self.fig.delaxes(self.fig_axs[1, 0])

        self.fig.tight_layout()

    def plot_morh(self, init_figs=True):
        if init_figs:
            self.fig_morh, self.morh_axs = plt.subplots(1, 2, dpi=300, figsize=(12, 6))

        self.morh = Morhplot(
            self.morh_axs[0],
            self.snns,
            self.taus,
            friction_state=FrictionState(k_f=self.k_f, tau_f=self.tau_f),
            stress_state=self.stress_state,
            fractures_snns=self.fracture_snns,
            fractures_taus=self.fracture_taus,
            title="Morh diagram",
            plot_line=False,
        )

        self.morh_reduced = Morhplot(
            self.morh_axs[1],
            self.snns_reduced,
            self.taus_reduced,
            friction_state=FrictionState(k_f=self.k_f, tau_f=self.tau_f),
            stress_state=self.stress_state,
            fractures_snns=self.fracture_snns_reduced,
            fractures_taus=self.fracture_taus_reduced,
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

    def plot_stereo(self):
        self.fig, self.fig_axs = plt.subplots(2, 3, dpi=300, figsize=(12, 6))

        self.snn_cont = self.draw_stereonet(
            self.snns_reduced,
            ax=self.fig_axs[0, 1],
            title=r"$\sigma_{nn}$",
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
            rf"$\mu_s$ = {self.stress_state.values.mu_s:.2f}, $\phi$ = {self.stress_state.values.phi:.2f}",
            fontsize=16,
        )

        self.fig.tight_layout()

        if self.output:
            plt.savefig(self.output)

    def draw_stereonet(
        self,
        z: ArrayLike,
        ax,
        title,
        colormap=None,
        levels=None,
        directions=None,
        norm_zero: bool = True,
        use_contourf: bool = False,
        show_cbar: bool = True,
    ) -> Stereonet:
        """Draw single stereonet plot."""

        cmap = colormap or "PuOr"
        if cmap == "Binary":
            cmap = colors.ListedColormap(["grey", "white", "orange"])
            show_cbar = False

        array = self.reshape_for_imshow(z)

        return Stereonet(
            array,
            stress_state=self.stress_state,
            ax=ax,
            title=title,
            fractures_xx=self.fractures_xx,
            fractures_yy=self.fractures_yy,
            colormap=cmap,
            levels=levels,
            directions=directions,
            norm_zero=norm_zero,
            use_contourf=use_contourf,
            show_cbar=show_cbar,
        )


plot = Plot()
