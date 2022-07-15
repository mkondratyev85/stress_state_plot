import click

from calculate import calculate_stress_on_planes, calculate_stresses_on_fractures
from plot2 import plot
from entities import StressState
from load_from_file import load_fractures
from save_to_xlsx import save_to_xlsx
from gui import gui



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

def parse_sigma_and_direction(sigma_and_direction):
    if sigma_and_direction is None:
        return None
    try:
        sigma_dir, sigma_dip, direction = sigma_and_direction.split(" ")
        direction = float(direction)
        sigma_dir = float(sigma_dir)
        sigma_dip = float(sigma_dip)
        return sigma_dir, sigma_dip, direction
    except ValueError:
        raise ValueError(f"Error while parsing sigma and direction argument {sigma_and_direction}")



class Morh:

    def __call__(self, stress_state, file_with_fractures, xlsx_path, png_path, gui_flag=False):

        if file_with_fractures:
            fractures = load_fractures(file_with_fractures)
            stresses_on_fractures = calculate_stresses_on_fractures(stress_state, fractures)
        else:
            fractures = None
            stresses_on_fractures = None
        print(stress_state)

        if gui_flag:
            gui(stress_state)
            return

        stresses_on_plane = calculate_stress_on_planes(stress_state, resolution = 15)


        if png_path:
            plot(stresses_on_plane, stress_state, stresses_on_fractures=stresses_on_fractures, output=png_path)
            print(f"Plots saved to {png_path}")

        if xlsx_path:
            save_to_xlsx(xlsx_path=xlsx_path, stresses_on_plane=stresses_on_fractures)
            print(f"Report saved to {xlsx_path}")


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
@click.option('--sigma1_orientation_and_sigma3_direction', help='Set orientations only via single sigma1_orientationa and direction of sigma3')
@click.option('--sigma3_orientation_and_sigma1_direction', help='Set orientations only via single sigma3_orientationa and direction of sigma1')
@click.option('--fractures', help='Path to the input file with fractures.')
@click.option('--xlsx_path', help='Path to the output xlsx file with report.')
@click.option('--png_path', help='Path to the output xlsx file with report.')
@click.option('--gui', help='Plots stress in interactive GUI.', is_flag=True)
def run(pressure, tau, mu_sigma, phi, sigma1_value, sigma2_value, sigma3_value, 
        sigma1_orientation, sigma2_orientation, sigma3_orientation, 
        sigma1_orientation_and_sigma3_direction, sigma3_orientation_and_sigma1_direction, 
        fractures, xlsx_path, png_path, gui):
    try:
        stress_state = StressState(
            orientations = {
                'sigma1': parse_sigma_orientation(sigma1_orientation),
                'sigma2': parse_sigma_orientation(sigma2_orientation),
                'sigma3': parse_sigma_orientation(sigma3_orientation),
                'sigma1_sigma3': parse_sigma_and_direction(sigma1_orientation_and_sigma3_direction),
                'sigma3_sigma1': parse_sigma_and_direction(sigma3_orientation_and_sigma1_direction),
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

    f(stress_state=stress_state, file_with_fractures=fractures, xlsx_path=xlsx_path, png_path=png_path, gui_flag=gui)


if __name__ == '__main__':
    run()
    
