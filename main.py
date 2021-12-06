from plane import Plane
from calculate import calculate_stress_on_planes
from plot2 import plot
from entities import StressState




class Morh:

    def __call__(self, stress_state):
        self.stress_state = stress_state
        self.stresses_on_plane = calculate_stress_on_planes(stress_state)
        plot(self.stresses_on_plane, self.stress_state, output="/tmp/fig1.png")
        print(stress_state)

f = Morh()

if __name__ == '__main__':
    stress_state = StressState(
        orientations = {
            'sigma1': (144,20), 
            'sigma3': Plane(300, 49),
            },
        values = {
            'sigma3': 52, 
            'sigma2': 60,
            'sigma1': 64, 
            },
    )
    f(stress_state=stress_state)

