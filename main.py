import click

from plane import Plane
from calculate import calculate_stress_on_planes, calculate_stresses_on_fractures
from plot2 import plot
from entities import StressState
from load_from_file import load_fractures



def parse_sigma_orientation(string_value):
    if string_value is None:
        return None
    values = string_value.split(" ")
    if len(values) != 2:
        raise ValueError("Sigma orientation must be set as two floats separated by a space\n"
                f"{string_value} was giving instead"
                )
    dir = float(values[0])
    dip = float(values[1])
    return dir, dip



class Morh:

    def __call__(self, stress_state):

        fractures = load_fractures("fractures.txt")

        self.stress_state = stress_state
        self.stresses_on_plane = calculate_stress_on_planes(stress_state)
        self.stresses_on_fractures = calculate_stresses_on_fractures(stress_state, fractures)
        plot(self.stresses_on_plane, self.stress_state, stresses_on_fractures=self.stresses_on_fractures, output="/tmp/fig1.png")
        print(stress_state)

f = Morh()

@click.command()
@click.option('--pressure', help='Value of pressure.', type=float)
@click.option('--tau', help='Value of tau.', type=float)
@click.option('--mu_sigma', help='Value of mu_sigma.', type=float)
@click.option('--phi', help='Phi.', type=float)
@click.option('--sigma1_value', help='Value of sigma1.')
@click.option('--sigma2_value', help='Value of sigma2.')
@click.option('--sigma3_value', help='Value of sigma3.')
@click.option('--sigma1_orientation', help='Orientation of sigma1.')
@click.option('--sigma2_orientation', help='Orientation of sigma2.')
@click.option('--sigma3_orientation', help='Orientation of sigma3.')
def run(pressure, tau, mu_sigma, phi, sigma1_value, sigma2_value, sigma3_value, sigma1_orientation, sigma2_orientation, sigma3_orientation):
    try:
        stress_state = StressState(
            orientations = {
                'sigma1': parse_sigma_orientation(sigma1_orientation),
                'sigma2': parse_sigma_orientation(sigma2_orientation),
                'sigma3': parse_sigma_orientation(sigma3_orientation),
                },
            values = {
                'sigma1': sigma1_value,
                'sigma2': sigma2_value,
                'sigma3': sigma3_value, 
                'tau': tau,
                'mu_s': mu_sigma,
                'phi': phi, 'p': pressure,
                },
        )
    except ValueError as e:
        click.echo("An error occured while setting stress state:", err=True)
        click.echo(e, err=True)
        return

    f(stress_state=stress_state)


if __name__ == '__main__':
    run()
    
