from bch.mathutils import *
from sympy.abc import x, alpha
from sympy import Matrix
import logging

log = logging.getLogger("bchcoder")


class BchCoder:

    def __init__(self, n, b, d, r_poly, g_poly):
        self.n = n
        self.b = b
        self.d = d
        self.r_poly = r_poly.set_domain(GF(2))
        self.g_poly = g_poly
        self.q = 2
        self.m = order(self.q, self.n)
        self.k = n - g_poly.degree()
        self.t = (self.d - 1) // 2
        log.info("BchCoder(n={},k={},q={},m={},b={},d={}, g={}) initiated"
                 .format(self.n, self.k, self.q, self.m, self.b, self.d, self.g_poly))

    def encode(self, msg_poly):
        shift_m_poly = msg_poly * Poly(x ** (self.n - self.k), x)
        log.info("shift_m: {}".format(shift_m_poly))
        r_poly = (shift_m_poly % self.g_poly).trunc(self.q)
        log.info("r: {}".format(r_poly))
        return (shift_m_poly - r_poly).trunc(self.q)

    def decode(self, msg_poly):
        log.debug("msg: {}".format(msg_poly))
        s = []
        for i in range(self.b, self.b + self.d - 1):
            s.append((Poly(msg_poly.eval(alpha ** i), alpha)).set_domain(GF(2)))

        log.debug("s: {}".format(s))

        error = Poly(0, alpha)
        for p in s:
            error += p
        if error.is_zero:
            return msg_poly.all_coeffs()[:self.k]

        S = Matrix(self.t, self.t, lambda i, j: s[i + j])
        S_det = (S.det()%self.r_poly).set_domain(GF(2))
        log.debug("S_det: {}".format(S_det))

        while S_det.is_zero:
            S = S[:-1, :-1]
            S_det = (S.det()%self.r_poly).set_domain(GF(2))
            log.debug("S_det: {} (t={})".format(S_det, S.shape[0]))

        log.debug("S: {}".format(S))

        t = S.shape[0]
        log.debug("t: {}".format(t))
        C = Matrix(t, 1, lambda i, j: s[i + t])
        log.debug("C: {}".format(C))

        S = S.col_insert(t, C)
        log.debug("S|C: {}".format(S.rref()))
        L = S.rref()[0].col(t)
        log.debug("L: {}".format(L))

        l_poly = Poly(1, x).set_domain(GF(2))
        log.debug("l(0): {}".format(l_poly))
        for i, p in enumerate(L[::-1], start=1):
            l_poly += Poly(flatten_frac(p, self.q).as_expr() * x ** i, x, alpha).set_domain(GF(2))
            log.debug("l({}): {}".format(i, l_poly))

        result = msg_poly.all_coeffs()
        for i in range(1, self.n + 1):
            test_poly = (Poly(Poly(l_poly, x).eval(alpha ** i), alpha) % self.r_poly).set_domain(GF(2))
            log.debug("testing: {}".format(test_poly))
            if test_poly.is_zero:
                log.info("REPAIRED ERROR ON {}th POSITION".format(i))
                result[i - 1] = 1 if result[i - 1] == 0 else 0

        log.debug("Message polynom after decoding: {}".format(result))
        return result[:self.k]
