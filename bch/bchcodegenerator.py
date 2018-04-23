from bch.mathutils import *
from sympy.polys.galoistools import gf_irreducible, gf_irreducible_p
from sympy import lcm, ZZ
import logging

log = logging.getLogger("bchcodegenerator")


class BchCodeGenerator:

    def __init__(self, n, b, d):
        self.n = n
        self.b = b
        self.d = d
        self.q = 2
        self.m = order(self.q, self.n)
        log.info("BchCodeGenerator(n={},q={},m={},b={},d={}) initiated"
                 .format(self.n, self.q, self.m, self.b, self.d))

    def gen(self):
        irr_poly = Poly(alpha ** self.m + alpha + 1, alpha).set_domain(GF(self.q))
        if gf_irreducible_p([int(c) for c in irr_poly.all_coeffs()], self.q, ZZ):
            quotient_size = len(power_dict(self.n, irr_poly, self.q))
        else:
            quotient_size = 0
        log.info("irr(q_size: {}): {}".format(quotient_size, irr_poly))
        while quotient_size < self.n:
            irr_poly = Poly([int(c.numerator) for c in gf_irreducible(self.m, self.q, ZZ)], alpha)
            quotient_size = len(power_dict(self.n, irr_poly, self.q))
            log.info("irr(q_size: {}): {}".format(quotient_size, irr_poly))
        g_poly = None
        for i in range(self.b, self.b + self.d - 1):
            if g_poly is None:
                g_poly = minimal_poly(i, self.n, self.q, irr_poly)
            else:
                g_poly = lcm(g_poly, minimal_poly(i, self.n, self.q, irr_poly))
        g_poly = g_poly.trunc(self.q)
        log.info("g(x)={}".format(g_poly))
        return irr_poly, g_poly
