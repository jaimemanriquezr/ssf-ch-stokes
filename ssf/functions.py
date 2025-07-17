def get_bulk_potential(name, psi_0, rho_b, phi_star):
    match name.lower():
        case "klapper":
            def F(u):
                phi = u / rho_b
                f0 = rho_b * psi_0
                return f0 * pow(phi, 3) * (phi - 4/3 * phi_star)
            def f(u):
                phi = u / rho_b
                f0 = 4 * psi_0
                return f0 * pow(phi, 2) * (phi - phi_star)
            def f_i(u):
                phi = u / rho_b
                f0 = 4 * psi_0
                return f0 * (3 - 2*phi_star) * phi
            def f_e(u):
                phi = u / rho_b
                f0 = 4 * psi_0
                f_convex = pow(phi, 2) * (phi - phi_star)
                f_linear = (3 - 2*phi_star) * phi
                return f0 * (f_convex - f_linear)
        case "double-well":
            def F(u):
                phi = u / rho_b
                f0 = rho_b * psi_0
                return f0 * pow(phi, 2) * pow(1. - phi, 2)
            def f(u):
                phi = u / rho_b
                f0 = psi_0 / 2
                return f0 * phi * (2*pow(phi, 2) - 3*phi + 1.)
            def f_i(u):
                phi = u / rho_b
                f0 = 3/4 * psi_0
                return f0 * phi
            def f_e(u):
                phi = u / rho_b
                f0 = 1/4 * psi_0
                return f0*(4*pow(phi, 3) - 6*pow(phi, 2) - pow(phi, 2))
    return F, f, f_i, f_e

