# A sharp within-level Fourier entropy bound for symmetric Boolean functions

Yitzchak Shmalo

[Paper (PDF)](paper/main.pdf) · [LaTeX source](paper/main.tex) · [Declarations](DECLARATIONS.md)

For a symmetric Boolean function the Fourier weight on level $k$ is spread evenly over the $\binom nk$ sets of size $k$, so the Fourier entropy splits as $H[f]=B[f]+L[f]$, with the within-level part $B[f]=\sum_k W_k\log_2\binom nk$ and the level part $L[f]=\sum_k W_k\log_2(1/W_k)$. O'Donnell, Wright and Zhou proved $B[f]\le 9.04\,I[f]$. The paper finds the exact constant.

## Main results

- **Theorem A.** For every non-constant symmetric $f:\{-1,1\}^n\to\{-1,1\}$,

  $$n\sum_k W_k\,h(k/n)\le 2\,I[f],\qquad\text{hence}\qquad B[f]\le\frac{2}{\ln 2}\,I[f]=2.885\ldots\,I[f],$$

  where $h$ is the binary entropy in nats. Majority shows that $2/\ln 2$ cannot be lowered.

- **Theorem B.** Write $H[f]=2A[f]+2\Gamma[f]$, where $A[f]=\sum_S\hat f(S)^2\log_2(1/\alpha(S))$ and $\alpha(S)$ is the average absolute value of the coefficient of $S$ in the restrictions of $f$ that fix the coordinates outside $S$. Then $\Gamma[f]\le I[f]/\ln 2$ for every non-constant symmetric $f$, and majority shows that $1/\ln 2$ cannot be lowered.

- **Rate.** $n\sum_k W_k h(k/n)\le(2-0.805\,n^{-1/2})\,I[f]$ for every $n\ge2$.

The proof writes $n\sum_k W_kh(k/n)/I[f]$ as a series in the return probabilities of the random walk on the cube, reduces to threshold functions, and compares the walk on Hamming levels with simple random walk on the integers.

## Computations

Python 3 with `numpy` and `mpmath`. Run from the `code` directory.

```
python3 certificate.py 1500
python3 checks.py 12 40
python3 constants.py
```

`certificate.py` works in exact rational arithmetic. For $2\le n\le1500$ and every threshold it bounds $2-T$ from below by the quantity $D(n,s)$ of Section 6 and checks $D>0$ (Theorem A) and $\sqrt n\,D\ge0.805$ (the rate). The smallest value of $\sqrt n\,D$ is $0.891458$, at $n=1500$; the same test with $1$ in place of $0.805$ fails for 1495 of the 1499 values of $n$. Analytic proofs cover $n\ge100$ (Theorem A) and $n\ge1000$ (the rate).

`checks.py` computes, from the definitions, the Fourier coefficients, $\alpha(S)$, $H$, $A$, $\Gamma$ and $T$ for all 16,356 non-constant symmetric functions on at most 12 variables and checks $H=2A+2\Gamma$, the rate and Theorem B. It checks $\alpha(S)=I_1$ for majority on 3 to 13 variables. For the 419 threshold functions $\mathrm{thr}_s$ with $2\le n\le40$ and $s\le(n-1)/2$, which by symmetry cover all thresholds, it compares $T$ from the exact level weights with the series identity, and checks the random-walk comparison, the cap and the certificate's bound.

`constants.py` recomputes the numerical constants of Sections 3 and 4 in 40-digit arithmetic.

The outputs of the runs used in the paper are `code/out_*.txt` (Python 3.13.2, numpy 2.2.3, mpmath 1.3.0, macOS 15.5). `SHA256SUMS` lists the hashes of the files.
