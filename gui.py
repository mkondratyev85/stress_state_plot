"""Gui support for stress state
"""

from entities import StressState
from calculate import calculate_stress_on_planes, calculate_stresses_on_fractures
from plot2 import plot


def recalculate_stress(stress_state, planes=None):
    stresses_on_plane = calculate_stress_on_planes(
        stress_state, resolution=30, planes=planes
    )
    return stresses_on_plane


def gui(stress_state: StressState, fractures=None):
    stresses_on_plane = recalculate_stress(stress_state)
    if fractures:
        stresses_on_fractures = recalculate_stress(stress_state, planes=fractures)
    else:
        stresses_on_fractures = None
    plot(
        stresses_on_plane,
        stress_state,
        stresses_on_fractures=stresses_on_fractures,
        gui=True,
        gui_update_func=recalculate_stress,
    )
