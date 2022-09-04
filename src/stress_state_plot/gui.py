"""Gui support for stress state
"""

from .entities import StressState
from .calculate import calculate_stress_on_planes, calculate_stresses_on_fractures
from .plot2 import plot


def recalculate_stress(stress_state, planes=None):
    stresses_on_plane = calculate_stress_on_planes(
        stress_state, resolution=50, planes=planes
    )
    return stresses_on_plane


def gui(stress_state: StressState, k_f, tau_f, stresses_on_fractures=None):
    stresses_on_plane = recalculate_stress(stress_state)

    plot(
        stresses_on_plane,
        stress_state,
        stresses_on_fractures=stresses_on_fractures,
        gui=True,
        k_f=k_f,
        tau_f=tau_f,
        gui_update_func=recalculate_stress,
    )
