import argparse
def get_argument_parser():
    parser = argparse.ArgumentParser(
        prog='ssf_figure_export',
        description= 'Export figures from .xdmf file using ParaView'
    )
    # POSITIONAL ARGUMENTS (mandatory)
    parser.add_argument('simulation_name', help='Simulation name',
                        type=str)
    parser.add_argument('--input_dir',
                        type=str, default='./results/')
    parser.add_argument('--output_dir',
                        type=str, default='./output/')
    parser.add_argument('-o', '--output', help='Output format',
                        type=str, default='svg')

    # OPTIONAL ARGUMENTS
    parser.add_argument('-T', '--final_time',
                        type=float, default=1E-2)
    parser.add_argument('-N', '--time_steps',
                        type=int, default=1)

    parser.add_argument('-maxP', '--maximum_particle_concentration',
                        type=float, default=40.0)
    parser.add_argument('-maxL', '--maximum_fluid_concentration',
                        type=float, default=998.0)
    parser.add_argument('-f', '--axes_font',
                        type=int, default=24)
    parser.add_argument('-n', '--axes_ticks',
                        type=int, default=3)
    # SWITCHES
    parser.add_argument('--contour', action='store_true')
    parser.add_argument('--colorbar', action='store_true')
    parser.add_argument('--axes', action='store_true')
    parser.add_argument('--quiver', action='store_true')
    parser.add_argument('--all', action='store_true')
    parser.add_argument('--compress', action='store_true')
    return parser
