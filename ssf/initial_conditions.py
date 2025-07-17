from dolfin import Expression
def CircleIC(x0, y0, r0, z0, kappa, deg = 2):
    _norm = "sqrt(pow(x[0] - X, 2) + pow(x[1] - Y, 2))"
    text = "Z * 1/2 * ( tanh( (R - %s) / sqrt(2*K) ) + 1)" % _norm
    params = {"X" : x0, "Y" : y0, "R" : r0, "Z" : z0, "K" : kappa}
    return Expression(text, **params, degree=deg)

def RectangleIC(y0=0.0, z0=0.0, location="upper", deg = 2):
    params = {"Y" : y0, "Z" : z0}
    match location:
        case "upper":
            text = "x[1] > Y ? Z : 0."
        case "lower":
            text = "x[1] < Y ? Z : 0."
    return Expression(text, **params, degree=deg)

def get_preset(name, rhob, rhof, kappa, dist=.05, deg=2):
    match name:
        case "blobs":
            _norm1 = "sqrt(pow(x[0] - X1, 2) + pow(x[1] - Y1, 2))"
            text1 = "Z1/2 * (tanh( (R1 - %s)/sqrt(2*K) ) + 1)" % _norm1

            _norm2 = "sqrt(pow(x[0] - X2, 2) + pow(x[1] - Y2, 2))"
            text2 = "Z2/2 * (tanh( (R2 - %s)/sqrt(2*K) ) + 1)" % _norm2

            text = text1 + " + " + text2

            params = {"X1" : (.5 - dist), "Y1" : .75, "R1": .2,
                      "X2" : (.5 + dist), "Y2" : .35, "R2": .2,
                      "K": kappa, "degree" : deg}

            u_init = Expression(text, **params, 
                                Z1=2E-2*rhob, Z2=.8E-2*rhob)
            c_1_init = Expression(text, **params, 
                                  Z1=1E-2*rhob, Z2=.2E-2*rhob)
            s_1_init = CircleIC(.5, .5, .25, 20E-2*rhof, kappa)
        case "filter":
            u_init = RectangleIC(y0=dist, location="lower",
                                 z0=.05*rhob)
            c_1_init = RectangleIC(y0=dist, location="lower",
                                  z0=.025*rhob)
            s_1_init = RectangleIC(y0=.0, location="upper", 
                                   z0=.0)
    return u_init, c_1_init, s_1_init
