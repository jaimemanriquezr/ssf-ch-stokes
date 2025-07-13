import argparse
import numpy as np

def get_argument_parser():
    parser = argparse.ArgumentParser(
        description = "Run SSF experiment"
    )
    parser.add_argument("--name", type=str)
    parser.add_argument("--output_dir", type=str,
                        default="./results/")
    parser.add_argument("--preset_name", type=str,
                        default="blobs")
    parser.add_argument("--preset_parameter", type=float,
                        default=.05)

    parser.add_argument("--filter_length", type=float,
                        default=1.0)
    parser.add_argument("--filter_area", type=float,
                        default=1.0)
    parser.add_argument("--inlet_length", type=float,
                        default=0.1)
    parser.add_argument("--inlet_area", type=float,
                        default=0.1)
    parser.add_argument("--inlet_center", type=float,
                        default=0.2)
    parser.add_argument("--inflow_velocity", type=float,
                        default=0.0)
    parser.add_argument("--inflow_profile", type=str,
                        default="parabolic")
    parser.add_argument("--inflow_concentration", type=float,
                        default=0.0)

    parser.add_argument("--density_particles", type=float,
                        default=1117)
    parser.add_argument("--density_liquids", type=float,
                        default=998)

    parser.add_argument("--phi_star", type=float,
                        default=1E-2)
    parser.add_argument("--potential_name", type=str,
                        default="klapper")

    parser.add_argument("--kappa", type=float,
                        default=1E-8/2)
    parser.add_argument("--mobility_factor", type=float,
                        default=200.0)
    parser.add_argument("--mobility_coefficient", type=int,
                        default=0)
    parser.add_argument("--bulk_potential_factor", type=float,
                        default=1.0)

    parser.add_argument("--c_1", type=float,
                        default=1.0)
    parser.add_argument("--c_2", type=float,
                        default=1.0)
    parser.add_argument("--k_11", type=float,
                        default=2.00E-02)
    parser.add_argument("--k_12", type=float,
                        default=1.00E-02)
    parser.add_argument("--k_21", type=float,
                        default=2.00E-02)
    parser.add_argument("--k_22", type=float,
                        default=4.00E-02)

    parser.add_argument("--viscosity_liquid", type=float,
                        default=1E-3)
    parser.add_argument("--viscosity_biofilm", type=float,
                        default=1.0)

    parser.add_argument("--friction_factor", type=float,
                        default=1E-12)

    parser.add_argument("--n_steps", type=int,
                        default=20)
    parser.add_argument("--time_step_snap", type=float,
                        default=1E-7)
    parser.add_argument("--time_step_size", type=float,
                        default=1E-7)
    parser.add_argument("--n_mesh", type=int,
                        default=30)

    parser.add_argument("--gravity", type=float,
                        default=9.81)

    parser.add_argument("--is_mass_lumped", 
                        action='store_true')

    parser.add_argument("--abs_tol", type=float,
                        default=1E-9)
    parser.add_argument("--rel_tol", type=float,
                        default=1E-8)
    parser.add_argument("--max_iterations", type=int,
                        default=150)
    return parser
