from dolfin import *
from mshr import Polygon, generate_mesh
import numpy as np
import argparse
from scipy.io import savemat
from pathlib import Path
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True
class CahnHilliardEquation(NonlinearProblem):
    def __init__(self, **kwargs):
        NonlinearProblem.__init__(self)
        self.__dict__.update(kwargs)
    def F(self, b, x):
        assemble(self.f, tensor = b)
        # Mass lumping: F(x) = f(x) +  M*x
        b += self.mass * x
    def J(self, A, x):
        assemble(self.df, tensor = A)
        # Mass lumping: dF(x) = df + M
        A += self.mass
PARTICLES_N = 2
LIQUIDS_N = 2
IS_MASS_LUMPED = True
#---------------------------------------------------------------------#
from geometry import *
from initial_conditions import *
from cli_utilities import get_argument_parser
def _assign(u, v):
    _v = project(v, u.function_space())
    u.assign(_v)
##### ============================================================ #####
##### ==================== INPUT ARGUMENTS ======================= #####
##### ============================================================ #####
parser = get_argument_parser()
args = parser.parse_args()

SIMULATION_NAME = args.name
OUTPUT_DIR = args.output_dir
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
PRESET_NAME = args.preset_name
PRESET_PARAMETER = args.preset_parameter

FILTER_LENGTH = args.filter_length
FILTER_AREA = args.filter_area
INLET_LENGTH = args.inlet_length
INLET_AREA = args.inlet_area
INLET_CENTER = args.inlet_center
INFLOW_VELOCITY = args.inflow_velocity
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

TIME_STEP_NUMBER = args.n_steps
TIME_STEP_SNAP = args.time_step_snap
TIME_STEP_SIZE = args.time_step_size
N_MESH = args.n_mesh

GRAVITY_CONSTANT = args.gravity

NL_ABS_TOL = args.abs_tol
NL_REL_TOL = args.rel_tol
NL_MAX_IT = args.max_iterations
#----------------------------------------------------------------------#
dt = TIME_STEP_SIZE
t_final = TIME_STEP_NUMBER * dt
inc_mod = round(TIME_STEP_SNAP / TIME_STEP_SIZE)
#----------------------------------------------------------------------#
u_init, c_1_init, s_1_init = get_preset_conditions(PRESET_NAME,
                                                   DENSITY_PARTICLES,
                                                   DENSITY_LIQUIDS,
                                                   KAPPA,
                                                   PRESET_PARAMETER)
#----------------------------------------------------------------------#
match POTENTIAL_NAME.lower():
    case "klapper":
        def F(u):
            phi = u / DENSITY_PARTICLES
            f0 = DENSITY_PARTICLES * PSI_0
            return f0 * pow(phi, 3) * (phi - 4/3 * PHI_STAR)
        def f(u):
            phi = u / DENSITY_PARTICLES
            f0 = 4 * PSI_0
            return f0 * pow(phi, 2) * (phi - PHI_STAR)
        def f_i(u):
            phi = u / DENSITY_PARTICLES
            f0 = 4 * PSI_0
            return f0 * (3 - 2*PHI_STAR) * phi
        def f_e(u):
            phi = u / DENSITY_PARTICLES
            f0 = 4 * PSI_0
            f_convex = pow(phi, 2) * (phi - PHI_STAR)
            f_linear = (3 - 2*PHI_STAR) * phi
            return f0 * (f_convex - f_linear)
    case "double-well":
        def F(u):
            phi = u / DENSITY_PARTICLES
            f0 = DENSITY_PARTICLES * PSI_0
            return f0 * pow(phi, 2) * pow(1. - phi, 2)
        def f(u):
            phi = u / DENSITY_PARTICLES
            f0 = PSI_0 / 2
            return f0 * phi * (2*pow(phi, 2) - 3*phi + 1.)
        def f_i(u):
            phi = u / DENSITY_PARTICLES
            f0 = 3/4 * PSI_0
            return f0 * phi
        def f_e(u):
            phi = u / DENSITY_PARTICLES
            f0 = 1/4 * PSI_0
            return f0 * (4*pow(phi, 3) - 6*pow(phi, 2) - pow(phi, 2))
##### ============================================================ #####

##### ============================================================ #####
##### ================= GEOMETRY AND B.C. ======================== #####
##### ============================================================ #####
filter_args = {"filter_length": FILTER_LENGTH,
               "filter_area" : FILTER_AREA,
               "inlet_length" : INLET_LENGTH,
               "inlet_area" : INLET_AREA,
               "inlet_center" : INLET_CENTER,
               "inflow_velocity" : INFLOW_VELOCITY,
               "inflow_profile" : INFLOW_PROFILE}
mesh, bdy_fncs, ds = get_filter_geometry(mesh_N=N_MESH, **filter_args)

q_null = Constant((0., 0.))
q_in = INFLOW_VELOCITY
q_out = q_in * INLET_AREA/ FILTER_AREA
match INFLOW_PROFILE:
    case "constant":
        q_inflow = Constant((0., q_in))
        q_outflow = Constant((0., q_out))
    case "parabolic":
        _expr = ("0.0",
                 "-4*q / pow(b - a, 2) * (x[0] - a)*(x[0] - b)")
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
DG_0 = FiniteElement("DG", mesh.ufl_cell(), 0)
CG_1 = FiniteElement("CG", mesh.ufl_cell(), 1)
V_CG_2 = VectorElement("CG", mesh.ufl_cell(), 2)
R = FiniteElement("Real", mesh.ufl_cell(), 0)

u_list   = [DG_0, CG_1, CG_1]
c_s_list   = ((PARTICLES_N - 1) + (LIQUIDS_N - 1)) * [DG_0]
hat_list = [DG_0] + c_s_list

biofilm_elements  = u_list + c_s_list + hat_list
biofilm_space    = FunctionSpace(mesh, MixedElement(biofilm_elements))

stokes_elements = [V_CG_2, DG_0, R]
stokes_space    = FunctionSpace(mesh, MixedElement(stokes_elements))
#----------------------------------------------------------------------#
## The following functions will hold the solutions of the scheme
## at each time step. In other words, they correspond to u^{n+1}
biofilm_solution  = Function(biofilm_space)
u, mu, w, c_1, s_1, u_hat, c_1_hat, s_1_hat = split(biofilm_solution)

stokes_solution  = Function(stokes_space)
q, p, xi = split(stokes_solution)

## The following functions hold the data of the initial data
## of each time step iteration, i.e. they correspond to u^{n}.
## They are used to couple the schemes and are updated after each step.
## Therefore, each is treated as an independent scalar/vector function,
## whereas the ones corresponding to u^{n+1} belong to a mixed space.
def get_sub_space(v):
    mixed_function = v.ufl_operands[0]
    mixed_space = mixed_function.function_space()
    index = v.ufl_operands[1][0]._value
    return mixed_space.extract_sub_space([index]).collapse()
u_0 = Function(get_sub_space(u))
u_0.rename("u", "concentration of biofilm")
w_0 = Function(get_sub_space(w))
w_0.rename("u_tilde", "concentration of biofilm")
mu_0 = Function(get_sub_space(mu))
mu_0.rename("mu", "chemical potential")
c_1_0 = Function(get_sub_space(c_1))
c_1_0.rename("c_1", "concentration of biofilm component 1")
c_2_0 = Function(get_sub_space(c_1))
c_2_0.rename("c_2", "concentration of biofilm component 2")
s_1_0 = Function(get_sub_space(s_1))
s_1_0.rename("s_1", "concentration of liquid component 1")
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
_assign(c_2_0, _c_2)
_s_2 = DENSITY_LIQUIDS * (1. - (u_0 / DENSITY_PARTICLES)) - s_1_0
_assign(s_2_0, _s_2)

u_00 = interpolate(u_init, FunctionSpace(mesh, "CG", 2))
k_rho = KAPPA / DENSITY_PARTICLES
_assign(mu_0, -k_rho*(div(grad(u_00))) + f(u_00))

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
phi_b_0 = u_0 / DENSITY_PARTICLES
# Stokes functions:
def strain(v): return 0.5 * (grad(v) + grad(v).T)
nu = VISCOSITY_LIQUID * (1 - phi_b_0) + VISCOSITY_BIOFILM * (phi_b_0)

g = Constant((0., GRAVITY_CONSTANT))
g *= (DENSITY_PARTICLES - DENSITY_LIQUIDS) * phi_b_0
#----------------------------------------------------------------------#
dq, dp, dxi = TrialFunctions(stokes_space)
_q, _p, _xi = TestFunctions(stokes_space) 

stokes_LHS = (nu*inner(strain(dq),strain(_q))*dx - dp*(div(_q) + _xi)*dx
              + _p * (div(dq) + dxi) * dx)
stokes_RHS = -dot(g, _q)*dx - ALPHA*dot(mu_0*grad(w_0),_q)*dx
#----------------------------------------------------------------------#
stokes_problem = (stokes_LHS == stokes_RHS)
_space = stokes_space.sub(0)
if INLET_AREA > 0:
    stokes_bcs = [
        DirichletBC(_space, q_inflow, bdy_fncs["inlet"]),
        DirichletBC(_space, q_null, bdy_fncs["wall"]),
        DirichletBC(_space, q_outflow, bdy_fncs["outlet"])
    ]
else:
    def on_filter_boundary(x, on_boundary): return on_boundary
    stokes_bcs = [
        DirichletBC(_space, q_null, on_filter_boundary)
    ]
stokes_parameters = {"linear_solver" : "mumps"}
#----------------------------------------------------------------------#
### BIOFILM
_lt = lambda v : lt(v, DENSITY_PARTICLES / (2 + GAMMA))
u_m = DENSITY_PARTICLES / (2 + GAMMA)

def u_M(v):
    _M = LAMBDA * v * pow(1. - v/DENSITY_PARTICLES, 1 + GAMMA)
    return (_M + abs(_M)) / 2
def u_M_up(v):
    return conditional(_lt(v), u_M(v), u_M(u_m))
def u_M_dn(v):
    return conditional(_lt(v), Constant(0.), u_M(v) - u_M(u_m))

def c_M(u):
    _M = LAMBDA * pow(1. - u/DENSITY_PARTICLES, 1 + GAMMA)
    return (_M + abs(_M)) / 2
def c_M_up(u):
    return conditional(_lt(u), c_M(u), u_M(u_m) / u)
def c_M_dn(u):
    return conditional(_lt(u), Constant(0.), c_M(u) - u_M(u_m) / u)

def s_M(v):
    _M = LAMBDA * (v / DENSITY_PARTICLES) * pow(1. - v/DENSITY_PARTICLES, GAMMA)
    return (_M + abs(_M)) / 2

n_K = FacetNormal(mesh)('+')
ndot_pos = lambda b : (abs(dot(b, n_K)) + dot(b, n_K)) / 2
ndot_neg = lambda b : (abs(dot(b, n_K)) - dot(b, n_K)) / 2
def a_convection(b, p, w):
    f_in = ndot_pos(avg(b)) * p('+')
    f_out = ndot_neg(avg(b)) * p('-')
    return (f_in - f_out) * jump(w) * dS

def a_diffusion_u(b, u, w):
    f_in = ndot_pos(avg(b))*(u_M_up(u('+')) + u_M_dn(u('-')))
    f_out = ndot_neg(avg(b))*(u_M_up(u('-')) + u_M_dn(u('+')))
    return (f_in - f_out) * jump(w) * dS

def a_diffusion_c(b, u, c, w):
    f_in = ndot_pos(avg(b))*(c_M_up(u('+'))*c('+') + c_M_dn(u('-'))*c('-'))
    f_out = ndot_neg(avg(b))*(c_M_up(u('-'))*c('-') + c_M_dn(u('+'))*c('+'))
    return (f_in - f_out) * jump(w) * dS

def a_diffusion_s(b, u, s, w):
    f_in = -ndot_pos(avg(b))*s_M(u('+'))*s('+')
    f_out = -ndot_neg(avg(b))*s_M(u('-'))*s('-')
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
f = PSI_0 * (f_i(u) + f_e(u_0))
eq_mu = (mu*_mu*dx - (k_rho*dot(grad(w), grad(_mu))*dx + f*_mu*dx))

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
output_file = XDMFFile(OUTPUT_DIR + SIMULATION_NAME + ".xdmf")
output_file.parameters['rewrite_function_mesh'] = False
output_file.parameters["functions_share_mesh"]  = True
output_file.parameters["flush_output"]          = True
def save_results(t):
    for fx in [u_0, mu_0, w_0, c_1_0, c_2_0, s_1_0, s_2_0]:
        output_file.write(fx, t)
#---------------------------------------------------------------------#
var_file = open(OUTPUT_DIR + SIMULATION_NAME + "_parameters.txt", "w")
var_file.write("k_b, k_f = %i, %i\n" 
                % (PARTICLES_N, LIQUIDS_N))
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
#---------------------------------------------------------------------#
tracked_quantities = []
ch_energy = []
energy_CH = {"Form" : (k_rho/2*(dot(grad(u_0),grad(u_0))) + F(u_0)) * dx,
             "List" : ch_energy,
             "Name" : "CahnHilliardEnergy"}
tracked_quantities.append(energy_CH)
bf_mass = []
mass_CH = {"Form" : u_0 * dx,
           "List" : bf_mass,
           "Name" : "BiofilmMass"}
tracked_quantities.append(mass_CH)
st_energy = []
energy_ST = {"Form" : .5 * dot(q_0, q_0) * dx,
             "List" : st_energy,
             "Name" : "StokesEnergy"}
tracked_quantities.append(energy_ST)
if (C_1**2 + C_2**2) > 0:
    rx_energy = []
    energy_RX = {"Form" : r_u * mu_0 * dx,
                 "List" : rx_energy,
                 "Name" : "ReactionsEnergy"}
    tracked_quantities.append(energy_RX)
#---------------------------------------------------------------------#
variables_to_save = {}
for qty in tracked_quantities:
    variables_to_save |= {qty["Name"] : qty["List"]}
##### ============================================================ #####

##### ============================================================ #####
##### =================== TIME INTEGRATION ======================= #####
##### ============================================================ #####
NLsolver = NewtonSolver()
NLsolver.parameters['absolute_tolerance'] = NL_ABS_TOL
NLsolver.parameters['relative_tolerance'] = NL_REL_TOL
NLsolver.parameters['maximum_iterations'] = NL_MAX_IT
NLsolver.parameters['linear_solver']      = "mumps"
t, inc = 0., 0
save_results(0.)
while (t < t_final):
    solve(stokes_problem, stokes_solution, stokes_bcs, 
          solver_parameters=stokes_parameters)
    assign(q_0, stokes_solution.sub(0))
    #-----------------------------------------------------------------#
    NLsolver.solve(biofilm_problem, biofilm_solution.vector())
    assign(u_0,  biofilm_solution.sub(0))
    assign(mu_0, biofilm_solution.sub(1))
    assign(w_0,  biofilm_solution.sub(2))
    assign(c_1_0, biofilm_solution.sub(3))
    _assign(c_2_0, _c_2)
    assign(s_1_0, biofilm_solution.sub(4))
    _assign(s_2_0, _s_2)
    #-----------------------------------------------------------------#
    t   += dt
    inc += 1
    # print("Iteration: %d/%d" % (inc, TIME_STEP_NUMBER))
    #-----------------------------------------------------------------#
    for qty in tracked_quantities:
        qty["List"].append(assemble(qty["Form"]))
    if (inc % inc_mod) == 0:
        save_results(t)
print("Saving tracked quantities...")
savemat(OUTPUT_DIR + SIMULATION_NAME + "_ENERGY.mat", variables_to_save)
print("Done!")
##### ============================================================ #####
