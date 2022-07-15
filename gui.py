"""Gui support for stress state
"""

from entities import StressState
from calculate import calculate_stress_on_planes, calculate_stresses_on_fractures
from plot2 import plot


def recalculate_stress(mu_sigma: float = 0):
    stress_state = StressState(
        orientations = {
            'sigma1_sigma3': (125, 11, 30),
            },
        values = {
            'tau': 10,
            'mu_s': mu_sigma,
            'p': 40,
            },
    )
    stresses_on_plane = calculate_stress_on_planes(stress_state, resolution=30)
    return stress_state, stresses_on_plane

def gui(stress_state: StressState):
    stress_state, stresses_on_plane = recalculate_stress()
    plot(stresses_on_plane, stress_state, gui=True, gui_update_func=recalculate_stress)
