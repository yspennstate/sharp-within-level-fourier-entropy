"""The numerical constants in Sections 3 and 4, recomputed with 40 digits (mpmath)."""
import mpmath as mp
from mpmath import mpf, sqrt, log, exp, pi, erfc, findroot

mp.mp.dps = 40
import platform
print("python %s, mpmath %s, %s" % (platform.python_version(), mp.__version__, platform.platform()))
phi = lambda z: exp(-z * z / 2) / sqrt(2 * pi)
Phibar = lambda z: erfc(z / sqrt(2)) / 2

# Section 3: the two ratios at n = 100 and their limits
for n in (mpf(100), mpf(1000), mpf(10) ** 6):
    r1 = mpf("1.061") * sqrt(n * (n + 1)) / (sqrt(2) * (n - 1))
    r2 = sqrt(3) / 2 * (n + 2) * sqrt(n + 1) / ((n - 1) * sqrt(n - 2))
    print("n=%s  ratio(b) %s  ratio(a) %s" % (mp.nstr(n, 8), mp.nstr(r1, 8), mp.nstr(r2, 8)))
print("1.5 sqrt(1/2) =", mp.nstr(mpf(3) / 2 / sqrt(2), 10), "< 1.061")

# Lemma 4.2
N = mpf(999)
y_over = mpf("2.048") / (sqrt(N) * (1 - mpf("1.6") / sqrt(N)))
a = N / 2
kap = sqrt(a / (a + 1)) * exp(-y_over)
print("exponent bound %s, sqrt(a/(a+1)) %s, kappa >= %s" % (mp.nstr(y_over, 8), mp.nstr(sqrt(a / (a + 1)), 8), mp.nstr(kap, 8)))
c0 = sqrt(2 * pi) / 4
print("sqrt(2 pi)/4 = %s; times (2 - 0.933) = %s" % (mp.nstr(c0, 10), mp.nstr(c0 * (2 - mpf("0.933")), 10)))
print("case z>=1.6: %s" % mp.nstr(3 / sqrt(1000) + 1 / mpf("1.6"), 8))
print("case middle: %s" % mp.nstr(3 / sqrt(1000) + mpf("0.66865"), 8))
sig = sqrt(N) / 2
print("1.5 phi(0)/sigma = %s" % mp.nstr(mpf("1.5") * phi(0) / sig, 8))
z0 = findroot(lambda z: Phibar(z) - (mpf(1) / 2 - mpf("0.0379")), 0.1)
print("Phibar(z) = 0.4621 at z = %s" % mp.nstr(z0, 8))
print("case third: %s" % mp.nstr(1 / (4 * mpf("0.933") * phi(mpf("0.0952"))), 8))
# monotone checks on a grid for the three-case functions used with sqrt(N) <= sqrt(n)
zs = [mpf(k) / 1000 for k in range(0, 1601)]
print("max Phi Phibar / phi on [0,1.6] = %s" % mp.nstr(max((1 - Phibar(z)) * Phibar(z) / phi(z) for z in zs), 10))
print("max Phibar^2 / phi on [0,1.6]   = %s" % mp.nstr(max(Phibar(z) ** 2 / phi(z) for z in zs), 10))

# Proposition 4.1 for n >= 1000
c = 2 * sqrt(2) / sqrt(pi)
print("even n=1000: %s" % mp.nstr(c * sqrt(mpf(1000) / 1003) - mpf("0.76352"), 8))
print("odd  n=1001: %s" % mp.nstr(c * sqrt(mpf(1001) / 1003) - mpf("0.76352") * 1001 / 1000, 8))

# Theorem B: g* and psi(t1)
gs = findroot(lambda t: mp.diff(lambda u: log(1 + u * u) / u, t), 2)
print("g* = max ln(1+t^2)/t = %s at t = %s" % (mp.nstr(log(1 + gs ** 2) / gs, 12), mp.nstr(gs, 10)))
g = mpf("0.805")
t1 = (1 + sqrt(1 - g * g)) / g
t0 = (1 - sqrt(1 - g * g)) / g
print("t0 = %s, t1 = %s, psi(t1) = %s" % (mp.nstr(t0, 8), mp.nstr(t1, 8), mp.nstr(g * t1 - log(1 + t1 ** 2), 8)))
print("g*/(2 ln 2) = %s" % mp.nstr(log(1 + gs ** 2) / gs / (2 * log(2)), 8))
print("2/ln2 = %s, 2 + 1/ln2 = %s, 2/ln2 + log2(4+e) = %s" % (mp.nstr(2 / log(2), 8), mp.nstr(2 + 1 / log(2), 8), mp.nstr(2 / log(2) + log(4 + mp.e, 2), 8)))
print("OWZ C1 = %s" % mp.nstr((1 + 4 * sqrt(2 * mp.e / pi)) / log(2), 8))
