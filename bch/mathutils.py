import numpy as np
from sympy.abc import x, alpha
from sympy.polys.galoistools import gf_div
from sympy import ZZ, GF, Poly, Pow, Add
import logging

log = logging.getLogger("bchcodegenerator")


def order(x, p):
    tx = int(x)
    for i in range(p):
        if tx == 1:
            return i + 1
        tx = (tx * x) % p
    return None


def minimal_poly(i, n, q, irr_poly):
    ti = int(i)
    checked = np.zeros(n, dtype=bool)
    checked[ti] = True
    poly = Poly(x - alpha ** ti, x)
    for k in range(n):
        ti = (ti * q) % n
        if checked[ti]:
            polys = [(Poly(c, alpha) % irr_poly).trunc(q) for c in poly.all_coeffs()]
            for p in polys:
                if p.degree() > 0:
                    raise Exception("Couldn't find minimal polynomial")
            coeffs = [p.nth(0) for p in polys]
            return Poly(coeffs, x)
        checked[ti] = True
        poly = poly * Poly(x - alpha ** ti, x)
    return None


def flatten_frac(muls, m, p):
    log.debug("Dividing: {}".format(muls))
    if len(muls.args) != 2:
        raise Exception("Wrong case")
    inv = muls.args[0]
    add = muls.args[1]
    if type(add) == Pow:
        inv, add = add, inv
    if type(add) != Add:
        raise Exception("Wrong case")
    if (inv.args[1] > 0):
        print(inv.args)
        raise Exception("Wrong case")
    inv_poly = (Poly(inv.args[0] ** -inv.args[1])).set_domain(GF(p))
    add_poly = (Poly(add)).set_domain(GF(p))
    result = Poly([ee.numerator for ee in
                   gf_div([e.p for e in add_poly.all_coeffs()], [e.p for e in inv_poly.all_coeffs()], p, ZZ)[
                       0]],
                  alpha).trunc(p)
    log.debug("Dividing result: {}".format(result))
    return result
