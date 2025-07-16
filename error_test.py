from dolfin import *
import numpy as np
from ssf.fem import get_fem_spaces, get_sub_space

N_0 = 50
INPUT_DIR = "./results/"
SIM_NAME = "norm_rx_TOL9_blob_far"

def load_data(n):
    _mesh = Mesh()
    filename = INPUT_DIR + SIM_NAME + ("_N_%i_dt_1E-07" % n)
    XDMFFile(filename + "_mesh.xdmf").read(_mesh)
    h = _mesh.hmax()

    _bf_space, _st_space = get_fem_spaces(_mesh)
    _U = Function(_bf_space)
    _u, _mu, _w, _c_1, _s_1, *_ = split(_U)
    _Q = Function(_st_space)
    _q, *_ = split(_Q)

    u = Function(get_sub_space(_u))
    c_1 = Function(get_sub_space(_c_1))
    s_1 = Function(get_sub_space(_s_1))
    q = Function(get_sub_space(_q[0]))

    data_file = XDMFFile(filename + "_checkpoints.xdmf")
    data_file.read_checkpoint(u, "u")
    data_file.read_checkpoint(c_1, "c_1")
    data_file.read_checkpoint(s_1, "s_1")
    data_file.read_checkpoint(q, "q")

    return _mesh, u, c_1, s_1, q

mesh_0, u_0, c_1_0, s_1_0, q_0 = load_data(N_0)

mesh_list = []
err_u = []
err_c_1 = []
err_s_1 = []
err_q = []

def compute_error(v, v_ref):
    proj_v = project(v_ref, v.function_space())
    return np.sqrt(assemble(((proj_v - v) ** 2) * dx))
def error_rate(err_list, mesh_list):
    err = np.array(err_list)
    h = np.array([m.hmax() for m in mesh_list])
    return np.log(err[1:] / err[:-1]) / np.log(h[1:] / h[:-1])

for n in [25, 30, 35]:
    mesh, u, c_1, s_1, q = load_data(n)

    mesh_list.append(mesh)
    err_u.append(compute_error(u, u_0))
    err_c_1.append(compute_error(c_1, c_1_0))
    err_s_1.append(compute_error(s_1, s_1_0))
    err_q.append(compute_error(q, q_0))

r_u = error_rate(err_u, mesh_list)
r_c_1 = error_rate(err_c_1, mesh_list)
r_s_1 = error_rate(err_s_1, mesh_list)
r_q = error_rate(err_q, mesh_list)
