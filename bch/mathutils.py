import numpy as np
from sympy.abc import x, alpha
from sympy import GF, Poly, Pow, Add, Symbol
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


def power_dict(n, irr, p):
    result = dict()
    for i in range(1, n + 1):
        test_poly = (Poly(alpha ** i, alpha) % irr).set_domain(GF(p))
        result[tuple(test_poly.all_coeffs())] = i
    return result


def flatten_frac(muls, m, p, pow_dict):
    if len(muls.args) == 0:
        return (Poly(muls, alpha) % m).set_domain(GF(p))
    log.debug("Dividing: {}".format(muls))
    if len(muls.args) != 2:
        raise Exception("Wrong case")
    inv = muls.args[0]
    add = muls.args[1]
    if type(add) == Pow:
        inv, add = add, inv
    log.debug("num: {}; denum: {}".format(add, inv))
    if type(add) != Add and type(add) != Symbol and type(add) != Pow:
        log.debug(type(add))
        add = int(add)
        if add < 0:
            inv_poly = (Poly(muls.args[0] ** -add) % m).set_domain(GF(p))
            if inv_poly.is_zero:
                raise Exception("Dividing by 0")
            result_pow = pow_dict[tuple(inv_poly.all_coeffs())]
            if result_pow < 0:
                result_pow += len(pow_dict)
            result = Poly(alpha ** result_pow, alpha).set_domain(GF(p))
            log.debug("Dividing result: {}".format(result))
            return (result % m).set_domain(GF(p))
    if (inv.args[1] > 0):
        print(inv.args)
        raise Exception("Wrong case")
    add_poly = (Poly(add) % m).set_domain(GF(p))
    if add_poly.is_zero:
        return add_poly
    i = pow_dict[tuple(add_poly.all_coeffs())]
    inv_poly = (Poly(inv.args[0] ** -inv.args[1]) % m).set_domain(GF(p))
    if inv_poly.is_zero:
        raise Exception("Dividing by 0")
    j = pow_dict[tuple(inv_poly.all_coeffs())]
    result_pow = i - j
    if result_pow < 0:
        result_pow += len(pow_dict)
    result = Poly(alpha ** result_pow, alpha).set_domain(GF(p))
    log.debug("Dividing result: {}".format(result))
    return (result % m).set_domain(GF(p))
