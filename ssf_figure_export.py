from getpass import getuser
from numpy import linspace
import paraview_biofilm
from pathlib import Path

# GET ARGUMENTS FROM COMMAND LINE
from cli_utilities import get_argument_parser
parser = get_argument_parser()
args = parser.parse_args()

DENSITY_PARTICLES = 1117
DENSITY_LIQUIDS = 998
PHI_STAR = 1E-2

if args.all:
    print_quiver = print_contour = print_axes = print_colorbar = [1, 0]
else:
    print_quiver = [int(args.quiver)]
    print_contour = [int(args.contour)]
    print_axes = [int(args.axes)]
    print_colorbar = [int(args.colorbar)]

if args.time_steps == -1:
    export_times = [0.]
elif args.time_steps == 0:
    export_times = [args.final_time]
else:
    export_times = linspace(0., args.final_time, args.time_steps + 1)

parameters = {
    "simulation_name" : args.sim_name,
    "input_dir" : args.input_dir,
    "rho_P" : DENSITY_PARTICLES,
    "rho_L" : DENSITY_LIQUIDS,
    "u_star" : DENSITY_PARTICLES * PHI_STAR,
    "max_u" : args.maximum_particle_concentration,
    "max_c" : args.maximum_particle_concentration,
    "max_s" : args.maximum_fluid_concentration,
    "axes_font_size" : args.font_size,
    "n_ticks" : args.number_of_ticks,
    "print_contour" : print_contour,
    "print_axes" : print_axes,
    "print_bars" : print_colorbar,
    "print_quiver" : print_quiver,
    "export_path" : args.output_dir + args.simulation_name + '/',
    "export_times" : export_times,
    "export_extension" : args.output,
    "export_params" : {
        'Plottitle' : 'ParaView GL2PS Export',
        'Compressoutputfile' : args.compress,
        'Linewidthscalingfactor' : 1.0,
        'Rasterize3Dgeometry' : 0,
        'Pointsizescalingfactor' : 1.0
    }
}
paraview_biofilm.export_results(**parameters)
