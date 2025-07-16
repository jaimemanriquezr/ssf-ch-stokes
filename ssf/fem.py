from dolfin import (
    NonlinearProblem, assemble,
    FiniteElement, VectorElement, MixedElement,
    FunctionSpace, Expression, project
)

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

def get_fem_spaces(mesh, k_particles=2, k_liquids=2):
    DG_0 = FiniteElement("DG", mesh.ufl_cell(), 0)
    CG_1 = FiniteElement("CG", mesh.ufl_cell(), 1)
    R = FiniteElement("Real", mesh.ufl_cell(), 0)

    u_list   = [DG_0, CG_1, CG_1]
    c_s_list   = ((k_particles - 1) + (k_liquids - 1)) * [DG_0]
    hat_list = [DG_0] + c_s_list

    biofilm_elements = u_list + c_s_list + hat_list
    biofilm_space = FunctionSpace(mesh, MixedElement(biofilm_elements))

    q_element = VectorElement("CG", mesh.ufl_cell(), 4)
    p_element = FiniteElement("DG", mesh.ufl_cell(), 3)
    stokes_elements = [q_element, p_element, R]
    stokes_space = FunctionSpace(mesh, MixedElement(stokes_elements))
    return biofilm_space, stokes_space

def proj_assign(u, v):
    _v = project(v, u.function_space())
    u.assign(_v)

def get_sub_space(v):
    mixed_function = v.ufl_operands[0]
    mixed_space = mixed_function.function_space()
    index = v.ufl_operands[1][0]._value
    return mixed_space.extract_sub_space([index]).collapse()
