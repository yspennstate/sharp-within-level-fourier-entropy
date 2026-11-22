"""Checks that do not use the identities of Section 2.

(1) For every non-constant symmetric f on n <= NB variables: Fourier coefficients and alpha(S) straight from the
    definitions (restrictions to the first k coordinates), then H, A, Gamma, B, I and T, and the inequalities
    H = 2A + 2Gamma, |f^(S)| <= alpha(S) <= I_1 (S nonempty), B <= (2/ln 2) I, T <= 2 - 0.805/sqrt(n), Gamma <= I/ln 2.
(2) For majority on n = 3, 5, ..., 13 variables: alpha(S) = I_1 for every nonempty S = [k].
(3) For every threshold thr_s on n <= NS variables: T computed from the exact level weights (Krawtchouk sums in
    integers) lies in [Tlow, Tlow + M(s)/(2KA)], where Tlow = sum_{r <= KA} w_r p_2r/p_1 with p_j from the walk
    on Hamming levels in integers (Lemma 2.1 truncated, tail bounded by Lemma 2.5).

usage: python checks.py [NB] [NS]      (defaults 12, 40)
"""
import math
import sys
import time
from fractions import Fraction
from math import comb, log, log2, sqrt

import numpy as np

LN2 = log(2)


def h(p):
    return 0.0 if p <= 0 or p >= 1 else -p * log(p) - (1 - p) * log(1 - p)


def coeffs(fvec, n, k):
    """c_[k](z) for all z (coordinates k..n-1), from f on {0,1}^n (bit i set = coordinate i is -1)."""
    y = np.arange(2 ** k)
    chi = (1 - 2 * (np.array([bin(v).count("1") for v in y]) % 2)).astype(np.int64)
    return (fvec.reshape(2 ** (n - k), 2 ** k) @ chi) / 2.0 ** k     # integer sums, exact


def part1(NB):
    worst_T, worst_G, worst_B, err_id, viol = -9.0, 0.0, 0.0, 0.0, []
    count = 0
    for n in range(1, NB + 1):
        wt = np.array([bin(v).count("1") for v in range(2 ** n)])
        for code in range(1, 2 ** (n + 1) - 1):                 # F(w) = -1 iff bit w of code is set
            F = np.array([-1.0 if (code >> w) & 1 else 1.0 for w in range(n + 1)])
            fvec = F[wt].astype(np.int64)
            fh, al = np.zeros(n + 1), np.ones(n + 1)
            for k in range(n + 1):
                c = coeffs(fvec, n, k)
                fh[k], al[k] = c.mean(), np.abs(c).mean()
            W = np.array([comb(n, k) for k in range(n + 1)]) * fh ** 2
            I = float(sum(k * W[k] for k in range(n + 1)))
            Inf1 = I / n
            nz = np.abs(fh) > 1e-12
            H = float(sum(W[k] * log2(1 / fh[k] ** 2) for k in range(n + 1) if nz[k]))
            A = float(sum(W[k] * log2(1 / al[k]) for k in range(n + 1) if nz[k]))
            G = float(sum(W[k] * log2(al[k] / abs(fh[k])) for k in range(n + 1) if nz[k]))
            B = float(sum(W[k] * log2(comb(n, k)) for k in range(n + 1)))
            T = n * float(sum(W[k] * h(k / n) for k in range(n + 1))) / I
            err_id = max(err_id, abs(H - 2 * A - 2 * G))
            count += 1
            ok = (abs(W.sum() - 1) < 1e-9 and all(al[k] >= abs(fh[k]) - 1e-12 for k in range(n + 1))
                  and all(al[k] <= Inf1 + 1e-12 for k in range(1, n + 1)) and B <= 2 / LN2 * I + 1e-9
                  and G <= I / LN2 + 1e-9 and (n == 1 or T <= 2 - 0.805 / sqrt(n) + 1e-12))
            if not ok:
                viol.append((n, code))
            worst_T = max(worst_T, sqrt(n) * (T - 2) if n > 1 else -9)
            worst_G, worst_B = max(worst_G, G * LN2 / I), max(worst_B, B * LN2 / (2 * I))
        print("  n %2d done" % n, flush=True)
    print("(1) %d functions, n <= %d: violations %s; max |H - 2A - 2Gamma| = %.2e" % (count, NB, viol[:5] or "none", err_id))
    print("    max sqrt(n)(T - 2) = %.6f (must be <= -0.805); max Gamma ln2/I = %.6f; max B ln2/(2I) = %.6f" % (worst_T, worst_G, worst_B))


def part2():
    for n in range(3, 14, 2):
        m = (n - 1) // 2
        wt = np.array([bin(v).count("1") for v in range(2 ** n)])
        fvec = np.where(wt > m, -1, 1).astype(np.int64)
        inf1 = comb(n - 1, m) / 2 ** (n - 1)
        dev = max(abs(np.abs(coeffs(fvec, n, k)).mean() - inf1) for k in range(1, n + 1))
        print("(2) Maj_%d: max_k |alpha([k]) - I_1| = %.2e  (I_1 = %.6f)" % (n, dev, inf1))


def T_levels(n, s):
    """T(thr_s) from the exact level weights."""
    two = 2 ** n
    tot = 0.0
    I = Fraction(0)
    for k in range(n + 1):
        # sum_{|x| = w} chi_[k](x) = coefficient of t^w in (1-t)^k (1+t)^(n-k)
        poly = [1]
        for _ in range(k):
            poly = [a - b for a, b in zip(poly + [0], [0] + poly)]
        for _ in range(n - k):
            poly = [a + b for a, b in zip(poly + [0], [0] + poly)]
        num = sum((poly[w] if w <= s else -poly[w]) for w in range(n + 1))
        Wk = Fraction(comb(n, k) * num * num, two * two)
        I += k * Wk
        tot += float(Wk) * h(k / n)
    return n * tot / float(I)


def part3(NS, KA=200, K=4095):
    """Tests Lemma 2.1 (series), Lemma 2.4 (p_2r/p_1 <= e_r for 2r <= n), Lemma 2.5 (p_2r/p_1 <= M(s)) and the
    certificate's inequality 2 - T >= D(n, s) against T computed from the exact level weights."""
    worst, nfail, ncase = -9.0, {"interval": 0, "walk": 0, "cap": 0, "D": 0}, 0
    wr = [None] + [Fraction(1, r * (2 * r - 1)) for r in range(1, K + 1)]
    er = [None] + [Fraction(2 * r * comb(2 * r, r), 4 ** r) for r in range(1, K + 1)]
    Pe, Pw = [Fraction(0)], [Fraction(0)]
    for r in range(1, K + 1):
        Pe.append(Pe[-1] + er[r] * wr[r])
        Pw.append(Pw[-1] + wr[r])
    for n in range(2, NS + 1):
        two = 2 ** n
        R = n // 2
        for s in range((n - 1) // 2 + 1):
            ncase += 1
            cum = sum(comb(n, k) for k in range(s + 1))
            M = Fraction(n * cum * (two - cum), comb(n, s) * two * (n - s))
            Wv = [1 if k > s else 0 for k in range(n + 1)]
            wts = [comb(n, k) for k in range(s + 1)]
            Tlow, num1 = Fraction(0), None
            for j in range(1, 2 * KA + 1):
                Wv = [(k * Wv[k - 1] if k > 0 else 0) + ((n - k) * Wv[k + 1] if k < n else 0) for k in range(n + 1)]
                if j == 1:
                    num1 = sum(wts[k] * Wv[k] for k in range(s + 1))
                elif j % 2 == 0:
                    r = j // 2
                    q = Fraction(sum(wts[k] * Wv[k] for k in range(s + 1)), num1 * n ** (j - 1))   # p_2r / p_1
                    if q > M:
                        nfail["cap"] += 1
                    if 2 * r <= n and q > er[r]:
                        nfail["walk"] += 1
                    Tlow += q * wr[r]
            Td = T_levels(n, s)
            lo, hi = float(Tlow), float(Tlow + M / (2 * KA))
            if not (lo - 1e-9 <= Td <= hi + 1e-9):
                nfail["interval"] += 1
                print("  !! n %d s %d: T %.10f not in [%.10f, %.10f]" % (n, s, Td, lo, hi))
            D = (Pe[K] - Pe[R] + Fraction(1, 64)) - M * (Pw[K] - Pw[R] + Fraction(1, 2 * K))
            if float(D) > 2 - Td + 1e-9:
                nfail["D"] += 1
                print("  !! n %d s %d: D %.10f > 2 - T %.10f" % (n, s, float(D), 2 - Td))
            worst = max(worst, sqrt(n) * (Td - 2))
    print("(3) thresholds n <= %d, %d cases: failures %s; max sqrt(n)(T - 2) = %.6f" % (NS, ncase, nfail, worst))


if __name__ == "__main__":
    NB = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    NS = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    import platform
    print("python %s, numpy %s, %s" % (platform.python_version(), np.__version__, platform.platform()))
    t0 = time.time()
    part2()
    part3(NS)
    part1(NB)
    print("time %.0f s" % (time.time() - t0))
