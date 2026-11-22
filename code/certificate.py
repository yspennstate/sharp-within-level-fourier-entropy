"""Exact certificate for the finite ranges of Theorem A and Proposition 4.1.

For 2 <= n <= NMAX and every 0 <= s <= (n-1)/2 it computes, in rational arithmetic,

    M(s) = 2 pi_<=(s) pi_>(s) / pi'(s),
    D(n, s) = (sum_{r=R+1}^{K} e_r w_r + 1/64) - M(s) (sum_{r=R+1}^{K} w_r + 1/(2K)),

with R = floor(n/2), K = 4095, e_r = 2r binom(2r, r) / 4^r and w_r = 1/(r(2r-1)). By Proposition 2.7 and the tail
bounds of Section 3, 2 - T(thr_s) >= D(n, s). The script checks D(n, s) > 0 and n D(n, s)^2 >= 0.805^2.
D is decreasing in M(s), so it suffices to test the largest M(s) for each n.
It also prints the largest ratio M(s)/Llow(n), where Llow(n) = A/B below is a lower bound for
L(n) = sum_{r>R} e_r w_r / sum_{r>R} w_r (Remark 3.1).
Control: the same test with 0.805 replaced by 1 must fail.

usage: python certificate.py [NMAX]     (default 1500)
"""
import math
import platform
import sys
import time
from fractions import Fraction
from math import comb

K = 4095
G = Fraction(805, 1000)


def max_M(n):
    two, cum, best = 2 ** n, 0, Fraction(0)
    for s in range((n - 1) // 2 + 1):
        cum += comb(n, s)                               # 2^n pi_<=(s)
        # pi'(s) = binom(n-1, s) / 2^(n-1) = 2 (n-s) binom(n, s) / (n 2^n)
        M = Fraction(n * cum * (two - cum), comb(n, s) * two * (n - s))
        best = max(best, M)
    return best


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    t0 = time.time()
    print("python %s, %s" % (platform.python_version(), platform.platform()))
    assert (K + 1) == 64 ** 2                            # sum_{r>K} e_r w_r >= (K+1)^(-1/2) = 1/64
    Pe, Pw = [Fraction(0)], [Fraction(0)]
    for r in range(1, K + 1):
        w = Fraction(1, r * (2 * r - 1))
        Pe.append(Pe[-1] + Fraction(2 * r * comb(2 * r, r), 4 ** r) * w)
        Pw.append(Pw[-1] + w)
    fails, ctrl_fails, worst = [], 0, (float("inf"), None)
    ratio = (0.0, None)
    for n in range(2, nmax + 1):
        R = n // 2
        A = Pe[K] - Pe[R] + Fraction(1, 64)
        B = Pw[K] - Pw[R] + Fraction(1, 2 * K)
        Mn = max_M(n)
        D = A - Mn * B
        rt = float(Mn * B / A)
        if rt > ratio[0]:
            ratio = (rt, n)
        if not (D > 0 and n * D * D >= G * G):
            fails.append(n)
        if not (D > 0 and n * D * D >= 1):
            ctrl_fails += 1
        q = math.sqrt(n) * float(D)
        if q < worst[0]:
            worst = (q, n)
        if n in (2, 3, 4, 5, 10, 50, 100, 150, 500, 1000) or n == nmax:
            print("n %5d  sqrt(n) D = %.6f" % (n, q), flush=True)
    print("n = 2..%d: %d failures %s" % (nmax, len(fails), fails[:10]))
    print("smallest sqrt(n) D = %.6f at n = %d" % worst)
    print("largest max_s M(s)/Llow(n) = %.6f at n = %d" % ratio)
    print("control (1 in place of 0.805): fails at %d of %d values of n" % (ctrl_fails, nmax - 1))
    print("time %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
