from dolfin import *
from mshr import Polygon, generate_mesh
import numpy as np
import argparse
from scipy.io import savemat
from pathlib import Path
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True
IS_MASS_LUMPED = True
#----------------------------------------------------------------------#
from ssf.fem import *
from ssf.initial_conditions import get_preset
from ssf.geometry import get_filter_geometry
from ssf.functions import get_bulk_potential
from ssf.cli import get_argument_parser
##### ============================================================ #####
##### ==================== INPUT ARGUMENTS ======================= #####
##### ============================================================ #####
parser = get_argument_parser()
args = parser.parse_args()

OUTPUT_DIR = args.output_dir
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

SIMULATION_NAME = args.name
PRESET_NAME = args.preset_name
PRESET_PARAMETER = args.preset_parameter

FILTER_LENGTH = args.filter_length
FILTER_AREA = args.filter_area
INLET_LENGTH = args.inlet_length
INLET_AREA = args.inlet_area
INLET_CENTER = args.inlet_center
INFLOW_VELOCITY = -args.inflow_velocity
INFLOW_PROFILE = args.inflow_profile
INFLOW_CONCENTRATION = args.inflow_concentration
DENSITY_PARTICLES = args.density_particles
DENSITY_LIQUIDS = args.density_liquids
PHI_STAR = args.phi_star
POTENTIAL_NAME = args.potential_name
KAPPA = args.kappa
LAMBDA = args.mobility_factor
GAMMA = args.mobility_coefficient
PSI_0 = args.bulk_potential_factor
C_1 = args.c_1
C_2 = args.c_2
K_11 = args.k_11
K_12 = args.k_12
K_21 = args.k_21
K_22 = args.k_22
VISCOSITY_LIQUID = args.viscosity_liquid
VISCOSITY_BIOFILM = args.viscosity_biofilm
ALPHA = args.friction_factor
GRAVITY_CONSTANT = args.gravity
TIME_STEP_NUMBER = args.n_steps
TIME_STEP_SNAP = args.time_step_snap
TIME_STEP_SIZE = args.time_step_size
N_MESH = args.n_mesh
NL_ABS_TOL = args.abs_tol
NL_REL_TOL = args.rel_tol
NL_MAX_IT = args.max_iterations
#----------------------------------------------------------------------#
dt = TIME_STEP_SIZE
t_final = TIME_STEP_NUMBER * dt
inc_mod = round(TIME_STEP_SNAP / TIME_STEP_SIZE)
#----------------------------------------------------------------------#
u_init, c_1_init, s_1_init = get_preset(PRESET_NAME, DENSITY_PARTICLES,
                                        DENSITY_LIQUIDS, KAPPA,
                                        PRESET_PARAMETER)
#----------------------------------------------------------------------#
psi, dpsi, dpsi_i, dpsi_e = get_bulk_potential(POTENTIAL_NAME, PSI_0,
                                               DENSITY_PARTICLES, 
                                               PHI_STAR)
##### ============================================================ #####

##### ============================================================ #####
##### ================= GEOMETRY AND B.C. ======================== #####
##### ============================================================ #####
mesh, boundaries, ds = get_filter_geometry(FILTER_LENGTH, FILTER_AREA,
                                           INLET_LENGTH, INLET_AREA,
                                           INLET_CENTER, N_MESH,
                                           quads=args.quadrilaterals)
XDMFFile(OUTPUT_DIR + SIMULATION_NAME + "_mesh.xdmf").write(mesh)
#----------------------------------------------------------------------#
q_null = Constant((0., 0.))
q_in = INFLOW_VELOCITY
q_out = q_in * INLET_AREA/ FILTER_AREA
match INFLOW_PROFILE:
    case "constant":
        q_inflow = Constant((0., q_in))
        q_outflow = Constant((0., q_out))
    case "parabolic":
        _expr = ("0.0", "-4*q / pow(b - a, 2) * (x[0] - a)*(x[0] - b)")
        inflow_params = {'q': q_in, 
                         'a': INLET_CENTER - INLET_AREA/2,
                         'b': INLET_CENTER + INLET_AREA/2,
                         'degree': 2}
        q_inflow = Expression(_expr, **inflow_params)
        outflow_params = {'q': q_out,
                          'a': 0.0,
                          'b': FILTER_AREA,
                          'degree': 2}
        q_outflow = Expression(_expr, **outflow_params)
##### ============================================================ #####

##### ============================================================ #####
##### ========== FINITE ELEMENT SPACES AND FUNCTIONS ============= #####
##### ============================================================ #####
biofilm_space, stokes_space = get_fem_spaces(mesh)
#----------------------------------------------------------------------#
## The following functions will hold the solutions of the scheme
## at each time step. In other words, they correspond to u^{n+1}
biofilm_solution  = Function(biofilm_space)
biofilm_solution.rename("U_h", "")
u, mu, w, c_1, s_1, u_hat, c_1_hat, s_1_hat = split(biofilm_solution)
stokes_solution  = Function(stokes_space)
stokes_solution.rename("Q_h", "")
q, p, xi = split(stokes_solution)
#----------------------------------------------------------------------#
## The following functions hold the data of the initial data
## of each time step iteration, i.e. they correspond to u^{n}.
## They are used to couple the schemes and are updated after each step.
## Therefore, each is treated as an independent scalar/vector function,
## whereas the ones corresponding to u^{n+1} belong to a mixed space.
u_0 = Function(get_sub_space(u))
u_0.rename("u", "concentration of biofilm")
w_0 = Function(get_sub_space(w))
w_0.rename("u_tilde", "concentration of biofilm")
mu_0 = Function(get_sub_space(mu))
mu_0.rename("mu", "chemical potential")
c_1_0 = Function(get_sub_space(c_1))
c_1_0.rename("c_1", "concentration of biofilm component 1")
s_1_0 = Function(get_sub_space(s_1))
s_1_0.rename("s_1", "concentration of liquid component 1")
c_2_0 = Function(get_sub_space(c_1))
c_2_0.rename("c_2", "concentration of biofilm component 2")
s_2_0 = Function(get_sub_space(s_1))
s_2_0.rename("s_2", "concentration of liquid component 2")
q_0 = Function(get_sub_space(q[0]))
q_0.rename("q", "Stokes velocity of fluid")
#----------------------------------------------------------------------#
## INITIAL CONDITIONS
u_0.assign(u_init)
c_1_0.assign(c_1_init)
s_1_0.assign(s_1_init)

_c_2 = u_0 - c_1_0
proj_assign(c_2_0, _c_2)
_s_2 = DENSITY_LIQUIDS * (1. - (u_0 / DENSITY_PARTICLES)) - s_1_0
proj_assign(s_2_0, _s_2)

u_00 = interpolate(u_init, FunctionSpace(mesh, "CG", 2))
k_rho = KAPPA / DENSITY_PARTICLES
proj_assign(mu_0, -k_rho*(div(grad(u_00))) + dpsi(u_00))

__dw = TrialFunction(w_0.function_space())
__w = TestFunction(w_0.function_space())
__mass_form = __dw * __w * dx
__mass = assemble(__mass_form)
if IS_MASS_LUMPED:
    __mass.zero()
    __mass.set_diagonal(assemble(action(__mass_form, Constant(1.))))
solve(__mass, w_0.vector(), assemble(u_0 * __w * dx), 'mumps')

## At first, both the u^{n+1} and u^{n} objects are initialized to u^0
assign(biofilm_solution, [u_0, mu_0, w_0,
                          c_1_0, s_1_0,
                          u_0, c_1_0, s_1_0])
##### ============================================================ #####

##### ============================================================ #####
##### ================= VARIATIONAL PROBLEMS ===================== #####
##### ============================================================ #####
### STOKES
phi_0 = u_0 / DENSITY_PARTICLES
# Stokes functions:
def strain(v): return 0.5 * (grad(v) + grad(v).T)
nu = VISCOSITY_LIQUID * (1 - phi_0) + VISCOSITY_BIOFILM * (phi_0)

g = Constant((0., GRAVITY_CONSTANT))
g *= (DENSITY_PARTICLES - DENSITY_LIQUIDS) * phi_0
#----------------------------------------------------------------------#
dq, dp, dxi = TrialFunctions(stokes_space)
_q, _p, _xi = TestFunctions(stokes_space) 
#----------------------------------------------------------------------#
stokes_LHS = (nu*inner(strain(dq),strain(_q))*dx - dp*(div(_q) + _xi)*dx
              + _p * (div(dq) + dxi) * dx)
stokes_RHS = -dot(g, _q)*dx - ALPHA*dot(mu_0*grad(w_0),_q)*dx
#----------------------------------------------------------------------#
stokes_problem = (stokes_LHS == stokes_RHS)
_space = stokes_space.sub(0)
if INLET_AREA > 0:
    stokes_bcs = [
        DirichletBC(_space, q_inflow, boundaries["inlet"]),
        DirichletBC(_space, q_null, boundaries["wall"]),
        DirichletBC(_space, q_outflow, boundaries["outlet"])
    ]
else:
    stokes_bcs = [
        DirichletBC(_space, q_null, boundaries["whole"])
    ]
stokes_parameters = {"linear_solver" : "mumps"}
#----------------------------------------------------------------------#
### BIOFILM
_lt = lambda v : lt(v, DENSITY_PARTICLES / (2 + GAMMA))
u_m = DENSITY_PARTICLES / (2 + GAMMA)

def u_M(v):
    _phi = v / DENSITY_PARTICLES
    _M = LAMBDA * v * pow(1. - _phi, 1 + GAMMA)
    return (_M + abs(_M)) / 2
def u_M_up(v):
    return conditional(_lt(v), u_M(v), u_M(u_m))
def u_M_dn(v):
    return conditional(_lt(v), Constant(0.), u_M(v) - u_M(u_m))

def c_M(u):
    _phi = u / DENSITY_PARTICLES
    _M = LAMBDA * pow(1. - _phi, 1 + GAMMA)
    return (_M + abs(_M)) / 2
def c_M_up(u):
    return conditional(_lt(u), c_M(u), u_M(u_m) / u)
def c_M_dn(u):
    return conditional(_lt(u), Constant(0.), c_M(u) - u_M(u_m) / u)

def s_M(v):
    _phi = v / DENSITY_PARTICLES
    _M = LAMBDA * _phi * pow(1. - _phi, GAMMA)
    return (_M + abs(_M)) / 2

n_K = FacetNormal(mesh)('+')
ndot_pos = lambda b : (abs(dot(b, n_K)) + dot(b, n_K)) / 2
ndot_neg = lambda b : (abs(dot(b, n_K)) - dot(b, n_K)) / 2
def a_convection(b, p, w):
    f_in = ndot_pos(avg(b)) * p('+')
    f_out = ndot_neg(avg(b)) * p('-')
    return (f_in - f_out) * jump(w) * dS

def a_diffusion_u(b, u, w):
    f_in = ndot_pos(avg(b)) * (u_M_up(u('+')) 
                               + u_M_dn(u('-')))
    f_out = ndot_neg(avg(b)) * (u_M_up(u('-')) 
                                + u_M_dn(u('+')))
    return (f_in - f_out) * jump(w) * dS

def a_diffusion_c(b, u, c, w):
    f_in = ndot_pos(avg(b)) * (c_M_up(u('+'))*c('+') 
                               + c_M_dn(u('-'))*c('-'))
    f_out = ndot_neg(avg(b)) * (c_M_up(u('-'))*c('-') 
                                + c_M_dn(u('+'))*c('+'))
    return (f_in - f_out) * jump(w) * dS

def a_diffusion_s(b, u, s, w):
    f_in = -ndot_pos(avg(b)) * s_M(u('+')) * s('+')
    f_out = -ndot_neg(avg(b)) * s_M(u('-')) * s('-')
    return (f_in - f_out) * jump(w) * dS
#----------------------------------------------------------------------#
dU_h = TrialFunction(biofilm_space)
du, dmu, dw, dc_1, ds_1, du_hat, dc_1_hat, ds_1_hat = split(dU_h) 
_U_h = TestFunction(biofilm_space)
_u, _mu, _w, _c_1, _s_1, _u_hat, _c_1_hat, _s_1_hat = split(_U_h)
#----------------------------------------------------------------------#
u_flux = (a_convection(q_0, u_hat, _u_hat)
          + LAMBDA * a_diffusion_u(-grad(mu), u_hat, _u_hat))
eq_u_hat = (u_hat - u_0)*_u_hat*dx + dt*u_flux

c_flux = (a_convection(q_0, c_1_hat, _c_1_hat)
          + LAMBDA * a_diffusion_c(-grad(mu), u_hat, c_1_hat, _c_1_hat))
eq_c_hat = (c_1_hat - c_1_0)*_c_1_hat*dx  + dt*c_flux

s_flux = (a_convection(q_0, s_1_hat, _s_1_hat)
          + LAMBDA * a_diffusion_s(-grad(mu), u_hat, s_1_hat, _s_1_hat))
eq_s_hat = (s_1_hat - s_1_0)*_s_1_hat*dx  + dt*s_flux

if np.abs(INFLOW_VELOCITY) > 0:
    q_n = dot(q_inflow, FacetNormal(mesh))
    s_in = INFLOW_CONCENTRATION
    eq_s_hat += dt * (q_n*s_in*_s_1_hat*ds(1) 
                      + q_n*s_1_hat*_s_1_hat*ds(2))
#----------------------------------------------------------------------#
r_1 = C_1 * s_1_0 / (s_1_0 + K_11) * c_1_0 * c_2_0 / (c_2_0 + K_12)
r_2 = C_2 * s_1_0 / (s_1_0 + K_21) * c_1_0 / (c_1_0 + K_22) * c_2_0

r_u = ( 0.5*r_1 + 0.5*r_2 )
r_c = ( 1.0*r_1 - 0.5*r_2 )
r_s = (-0.5*r_1 - 0.5*r_2 )

eq_u = ((u - u_hat) - dt*r_u) * _u  * dx
eq_c = ((c_1 - c_1_hat) - dt*r_c) * _c_1 * dx
eq_s = ((s_1 - s_1_hat) - dt*r_s) * _s_1 * dx
#----------------------------------------------------------------------#
dpsi_ey = PSI_0 * (dpsi_i(u) + dpsi_e(u_0))
eq_mu = (mu*_mu*dx - (k_rho*dot(grad(w),grad(_mu))*dx + dpsi_ey*_mu*dx))

## !!! u_bio = u_hat OR u_bio = u
u_bio = u_hat
## w * wbar * dx is implemented separately to account for mass lumping
# eq_w = (w - u_bio) * wbar * dx
eq_w  = -u_bio * _w * dx
if IS_MASS_LUMPED:
    mass_form = ( du*_u*dx + dmu*_mu*dx + dw*_w*dx 
                  + dc_1*_c_1*dx + ds_1*_s_1*dx 
                  + du_hat*_u_hat*dx 
                  + dc_1_hat*_c_1_hat*dx + ds_1_hat*_s_1_hat*dx )
    biofilm_mass = assemble(mass_form)
    biofilm_mass.zero()
    _1 = Constant((0., 0., 1., 0., 0., 0., 0., 0.))
    biofilm_mass.set_diagonal(assemble(action(mass_form, _1)))
else:
    mass_form = dw*_w*dx
    biofilm_mass = assemble(mass_form)
#----------------------------------------------------------------------#
biofilm_F = ( eq_u + eq_mu + eq_w 
              + eq_c + eq_s 
              + eq_u_hat + eq_c_hat + eq_s_hat )
biofilm_jacobian  = derivative(biofilm_F, biofilm_solution, dU_h)
biofilm_problem = CahnHilliardEquation(f = biofilm_F, 
                                       df = biofilm_jacobian, 
                                       mass = biofilm_mass)
##### ============================================================ #####

##### ============================================================ #####
##### ==================== OUTPUT RESULTS ======================== #####
##### ============================================================ #####
var_file = open(OUTPUT_DIR + SIMULATION_NAME + "_parameters.txt", "w")
var_file.write("rho_b, rho_f = %.2e, %.2e\n" 
                % (DENSITY_PARTICLES, DENSITY_LIQUIDS))
var_file.write("phi_star, = %.2e, psi_0 = %.2e\n" 
                % (PHI_STAR, PSI_0))
var_file.write("kappa, lambda, ALPHA, GAMMA = %.2e, %.2e, %.2e, %i\n" 
                % (KAPPA, LAMBDA, ALPHA, GAMMA))
var_file.write("nu_0, nu_1 = %.2e, %.2e\n" 
                % (VISCOSITY_LIQUID, VISCOSITY_BIOFILM))
var_file.write("q_in, s_in = %.2e, %.2e\n" 
                % (INFLOW_VELOCITY, INFLOW_CONCENTRATION))
var_file.write("C_1, C_2 = %.2e, %.2e\n" 
                % (C_1, C_2))
var_file.write("K_11, K_12, K_21, K_22 = %.2e, %.2e, %.2e, %.2e\n" 
                % (K_11, K_12, K_21, K_22))
var_file.close()
#----------------------------------------------------------------------#
tracked_quantities = []
def add_tracker(name, form):
    _list = []
    _dict = {"Form" : form, "List" : _list, "Name" : name}
    tracked_quantities.append(_dict)
def update_trackers():
    for quantity in tracked_quantities:
        quantity["List"].append(assemble(quantity["Form"]))
l2sq = lambda v : dot(v, v)
add_tracker("CahnHilliardEnergy", (k_rho/2*l2sq(grad(u_0)) + psi(u_0))*dx)
add_tracker("BiofilmMass", u_0 * dx)
add_tracker("StokesEnergy", .5 * l2sq(q_0) * dx)
add_tracker("MuGradientEnergy", l2sq(grad(mu_0)) * dx)
add_tracker("ImbalanceU", (u_hat - u_0) * dx)
add_tracker("ImbalanceC", (c_1_hat - c_1_0) * dx)
add_tracker("ImbalanceS", (s_1_hat - s_1_0) * dx)
if (C_1**2 + C_2**2) > 0:
    add_tracker("ReactionsEnergy", r_u * mu_0 * dx)
export_vars = {qty["Name"] : qty["List"] for qty in tracked_quantities}
#----------------------------------------------------------------------#
output_file = XDMFFile(OUTPUT_DIR + SIMULATION_NAME + ".xdmf")
output_file.parameters['rewrite_function_mesh'] = False
output_file.parameters["functions_share_mesh"]  = True
output_file.parameters["flush_output"]          = True
#----------------------------------------------------------------------#
time_list = []
iteration_list = []
def write_results(t, n_it):
    print("t = %.8e" % t)
    for fx in [q_0, u_0, mu_0, w_0, c_1_0, c_2_0, s_1_0, s_2_0]:
        output_file.write(fx, t)
    update_trackers()
    time_list.append(t)
    iteration_list.append(n_it)
##### ============================================================ #####

##### ============================================================ #####
##### =================== TIME INTEGRATION ======================= #####
##### ============================================================ #####
NLsolver = NewtonSolver()
NLsolver.parameters['absolute_tolerance'] = NL_ABS_TOL
NLsolver.parameters['relative_tolerance'] = NL_REL_TOL
NLsolver.parameters['maximum_iterations'] = NL_MAX_IT
NLsolver.parameters['linear_solver']      = "mumps"
t, inc, n_it = 0., 0, 0
while (t < t_final):
    solve(stokes_problem, stokes_solution, stokes_bcs, 
          solver_parameters=stokes_parameters)
    assign(q_0, stokes_solution.sub(0))
    #------------------------------------------------------------------#
    if (inc % inc_mod) == 0:
        write_results(t, n_it)
    #------------------------------------------------------------------#
    n_it, _ = NLsolver.solve(biofilm_problem, biofilm_solution.vector())
    #------------------------------------------------------------------#
    t   += dt
    inc += 1
    # print("Iteration: %d/%d" % (inc, TIME_STEP_NUMBER))
    #------------------------------------------------------------------#
    assign(u_0,  biofilm_solution.sub(0))
    assign(mu_0, biofilm_solution.sub(1))
    assign(w_0,  biofilm_solution.sub(2))
    assign(c_1_0, biofilm_solution.sub(3))
    assign(s_1_0, biofilm_solution.sub(4))
    _assign(c_2_0, _c_2)
    _assign(s_2_0, _s_2)
write_results(t, n_it)
#----------------------------------------------------------------------#
print("Saving tracked quantities...")
export_vars |= {"Time" : time_list, "NewtonIterations" : iteration_list}
savemat(OUTPUT_DIR + SIMULATION_NAME + ".mat", export_vars)
print("Done!")
with XDMFFile(OUTPUT_DIR + SIMULATION_NAME + "_checkpoints.xdmf") as file:
    file.write_checkpoint(u_0, "u", t)
    file.write_checkpoint(q_0, "q", t, append=True)
    file.write_checkpoint(mu_0, "mu", t, append=True)
    file.write_checkpoint(c_1_0, "c_1", t, append=True)
    file.write_checkpoint(s_1_0, "s_1", t, append=True)
##### ============================================================ #####
