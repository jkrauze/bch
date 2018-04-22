from bch.mathutils import *
from sympy.polys.galoistools import gf_irreducible
from sympy import lcm
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
        irr_poly = Poly([int(c.numerator) for c in gf_irreducible(self.m, self.q, ZZ)], alpha)
        log.info("irr: {}".format(irr_poly))
        g_poly = None
        for i in range(self.b, self.b + self.d - 1):
            if g_poly is None:
                g_poly = minimal_poly(i, self.n, self.q, irr_poly)
            else:
                g_poly = lcm(g_poly, minimal_poly(i, self.n, self.q, irr_poly))
        g_poly = g_poly.trunc(self.q)
        log.info("g(x)={}".format(g_poly))
        return irr_poly, g_poly
