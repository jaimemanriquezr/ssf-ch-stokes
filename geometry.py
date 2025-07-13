from dolfin import *
from mshr import Polygon, generate_mesh
import numpy as np
def get_filter_geometry(filter_length, filter_area,
                        inlet_length, inlet_area, inlet_center,
                        mesh_N, **_):
    outline_filter = [Point(0., filter_length),
                      Point(0., 0.),
                      Point(filter_area, 0.),
                      Point(filter_area, filter_length)]

    if inlet_area > 0:
        inlet_left = inlet_center - inlet_area/2
        inlet_right = inlet_center + inlet_area/2
        total_length = filter_length + inlet_length
        outline_inlet = [Point(inlet_right, filter_length), 
                         Point(inlet_right, total_length),
                         Point(inlet_left, total_length), 
                         Point(inlet_left, filter_length)]
    else:
        outline_inlet = []

    def on_filter_inlet(x, on_boundary):
        _radius = inlet_area/2 + DOLFIN_EPS
        on_inlet_0 = np.abs(x[0] - inlet_center) < _radius

        _length = filter_length + inlet_length - DOLFIN_EPS
        on_inlet_1 = x[1] > _length

        on_inlet = on_inlet_0 and on_inlet_1
        return on_boundary and on_inlet

    def on_filter_outlet(x, on_boundary):
        on_bottom = (x[1] < DOLFIN_EPS)
        return on_boundary and on_bottom

    def on_filter_wall(x, on_boundary):
        on_left = (x[0] < DOLFIN_EPS)
        on_right = (x[0] > filter_area - DOLFIN_EPS)
        on_filter_sides = on_left or on_right

        _length_up = filter_length + inlet_length - DOLFIN_EPS
        _length_down = filter_length - DOLFIN_EPS
        on_top = (_length_down < x[1]) and (x[1] < _length_up)
        _radius = inlet_area/2 - DOLFIN_EPS
        not_inlet_0 = np.abs(x[0] - inlet_center) > _radius
        # includes inlet sides + filter top (two L shapes)
        on_filter_top = on_top and not_inlet_0

        on_wall = on_filter_top or on_filter_sides
        return on_boundary and on_wall

    filter_domain = Polygon(outline_filter + outline_inlet)
    mesh = generate_mesh(filter_domain, mesh_N)

    boundary_functions = {"inlet" : on_filter_inlet,
                          "outlet" : on_filter_outlet,
                          "wall" : on_filter_wall}
    ds = get_boundary_measure(mesh, boundary_functions)
    return mesh, boundary_functions, ds

def get_boundary_measure(mesh, boundary_functions):
    #https://sites.math.rutgers.edu/~falk/math575/Boundary-conditions.html
    facet_dim = mesh.ufl_cell().topological_dimension() - 1
    mf = MeshFunction("size_t", mesh, facet_dim)
    mf.set_all(0)

    on_inlet = boundary_functions["inlet"]
    class InletBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_inlet(x, on_boundary)
    filter_inlet = InletBoundary()
    filter_inlet.mark(mf, 1)

    on_outlet = boundary_functions["outlet"]
    class OutletBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_outlet(x, on_boundary)
    filter_outlet = OutletBoundary()
    filter_outlet.mark(mf, 2)

    on_wall = boundary_functions["wall"]
    class WallBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_wall(x, on_boundary)
    filter_inlet = WallBoundary()
    filter_inlet.mark(mf, 3)

    _ds = ds(subdomain_data = mf)
    return _ds
