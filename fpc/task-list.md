
# FPC Porting Task List

## Conventions
- `fpc/_shared/` → Reusable Pascal units (algorithms)
- `fpc/<topic>/src/` → CLI demo project (.lpr + .lpi)
- `fpc/<topic>/gui/src/` → GUI (LCL/Lazarus) demo (optional)
- Unit naming: `jpm<topic>.pas`
- Base type: `Float = Double` in `jpmtypes.pas`
- Callbacks: `TFunction1`, `TFunction2`
- Binaries: always go to `fpc/bin/`

## Done
- `jpmtypes.pas`, `jpminterpolation.pas` (Akima, Lagrange), `jpmintegration.pas` (Romberg)
- `jpmstats.pas` (NormalDist, InvNormalDist, ProbST, StudNT)
- `jpmappointment.pas` (Hungarian job assignment)
- [x] **CAT 3 ALL** — `jpmroots.pas`: Bisection, Newton, Secant, RegulaFalsi, Brent, Mueller, Steffensen, Aitken, Bairstow ✅ COMMITTED
- CLI demos: akima, lagrange, romberg, normal, student, prime, primes, appointment, roots
- GUI demos: normal, appointment

---

## Remaining Tasks

### CAT 1 — Statistics (`stat/`)
- [ ] 1.1 Chi-Square: `ChiSquareDist`, `InvChiSquare` → `jpmstats.pas`; demo `fpc/chi2/src/`
- [ ] 1.2 F-Distribution: `FDist`, `InvFDist` → `jpmstats.pas`; demo `fpc/fdist/src/`
- [ ] 1.3 Gamma/Beta: `GammaFunc`, `LnGamma`, `BetaFunc`, `IncompleteBeta` → `jpmspecial.pas`; demo `fpc/gamma/src/`
- [ ] 1.4 Distributions (Binomial, Poisson, etc.) → `jpmstats.pas`; demo `fpc/distributions/src/`
- [ ] 1.5 Moments: Mean, Variance, StdDev, Skewness, Kurtosis, Median → `jpmstats.pas`; demo `fpc/moments/src/`

### CAT 2 — Integration & Derivatives (`functions1/`)
- [ ] 2.1 Derivatives: `Derivative`, `SecondDerivative` → `jpmderivative.pas`; demo `fpc/derivative/src/`
- [ ] 2.2 Simpson, Gauss-Legendre → `jpmintegration.pas`; demo `fpc/integration/src/`
- [ ] 2.3 Chebyshev: Fit, Eval, Deriv, Integ → `jpmchebyshev.pas`; demo `fpc/chebyshev/src/`
- [ ] 2.4 Continued Fractions → `jpmcontinued.pas`; demo `fpc/confract/src/`

### CAT 3 — Root Finding (`roots/`) ✅ DONE

### CAT 4 — Optimization (`functions1/`, `functions2/`)
- [ ] 4.1 Golden Section → `jpmoptimize.pas`; demo `fpc/optimize/src/`
- [ ] 4.2 Brent Minimization → `jpmoptimize.pas`
- [ ] 4.3 Nelder-Mead (Amoeba) → `jpmoptimize.pas`
- [ ] 4.4 Powell's Method → `jpmoptimize.pas`
- [ ] 4.5 Steepest Descent → `jpmoptimize.pas`

### CAT 5 — Linear Algebra / Matrices (`matrices/`) ← **NEXT**
- [ ] 5.1 LU Decompose/Solve/Inverse → `jpmmatrices.pas`; demo `fpc/matrices/src/`
- [ ] 5.2 Gauss-Seidel iteration → `jpmmatrices.pas`
- [ ] 5.3 Cholesky decomposition → `jpmmatrices.pas`
- [ ] 5.4 Eigenvalues (Jacobi, HQR) → `jpmmatrices.pas`
- [ ] 5.5 Tridiagonal systems → `jpmmatrices.pas`
- [ ] 5.6 Determinant → `jpmmatrices.pas`

### CAT 6 — Least Squares (`lstsqr/`)
- [ ] 6.1 Polynomial Fit → `jpmlstsqr.pas`; demo `fpc/lstsqr/src/`
- [ ] 6.2 Linear/Multiple Regression → `jpmlstsqr.pas`
- [ ] 6.3 Levenberg-Marquardt → `jpmlstsqr.pas`

### CAT 7 — Differential Equations (`diffeqa/`)
- [ ] 7.1 Runge-Kutta 4 → `jpmdiffeq.pas`; demo `fpc/diffeq/src/`
- [ ] 7.2 RKF45 → `jpmdiffeq.pas`
- [ ] 7.3 Adams-Bashforth/Moulton → `jpmdiffeq.pas`
- [ ] 7.4 Gear (stiff ODEs) → `jpmdiffeq.pas`
- [ ] 7.5 Bulirsch-Stoer → `jpmdiffeq.pas`

### CAT 8 — Polynomials (`polynomials/`)
- [ ] 8.1 Horner Eval → `jpmpolynomials.pas`; demo `fpc/polynomials/src/`
- [ ] 8.2 Poly Arithmetic (Mul, Div, Deriv) → `jpmpolynomials.pas`
- [ ] 8.3 Poly GCD → `jpmpolynomials.pas`
- [ ] 8.4 Partial Fractions → `jpmpolynomials.pas`

### CAT 9 — Special Functions (`functions2/`, `series/`, `bessel/`)
- [ ] 9.1 Legendre P, Q → `jpmspecial.pas`; demo `fpc/special/src/`
- [ ] 9.2 Hermite → `jpmspecial.pas`
- [ ] 9.3 Laguerre → `jpmspecial.pas`
- [ ] 9.4 Airy Ai, Bi → `jpmspecial.pas`
- [ ] 9.5 Elliptic K, E → `jpmspecial.pas`
- [ ] 9.6 Bessel J, Y, I, K → `jpmbessel.pas`; demo `fpc/bessel/src/`
- [ ] 9.7 Hypergeometric → `jpmspecial.pas`

### CAT 10 — Sorting & Searching (`sorting/`)
- [ ] 10.1 BubbleSort, MergeSort, QuickSort → `jpmsorting.pas`; demo `fpc/sorting/src/`
- [ ] 10.2 BinarySearch, LinearSearch → `jpmsorting.pas`

### CAT 11 — Signal Processing (`signal/`)
- [ ] 11.1 FFT, InverseFFT → `jpmsignal.pas`; demo `fpc/fft/src/`
- [ ] 11.2 Digital Filters (Low/High/Band) → `jpmsignal.pas`
- [ ] 11.3 Convolve, Deconvolve → `jpmsignal.pas`
- [ ] 11.4 Savitzky-Golay smoothing → `jpmsignal.pas`

### CAT 12 — Arithmetic / Number Theory (`arith/`)
- [ ] 12.1 Primes — verify existing `fpc/arith/prime/` and `fpc/arith/primes/` demos
- [ ] 12.2 GCD, LCM → `jpmarith.pas`; demo `fpc/arith/src/`
- [ ] 12.3 Combinatorics: Factorial, Binomial, Permutation → `jpmarith.pas`
- [ ] 12.4 Fractions → `jpmarith.pas`
- [ ] 12.5 Base Conversion → `jpmarith.pas`
- [ ] 12.6 Diophantine equations → `jpmarith.pas`

### CAT 13 — Geometry (`geometry/`)
- [ ] 13.1 Triangle solver → `jpmgeometry.pas`; demo `fpc/geometry/src/`
- [ ] 13.2 Conic sections → `jpmgeometry.pas`
- [ ] 13.3 Circle fitting → `jpmgeometry.pas`

### CAT 14 — Linear Programming (`linearprog/`)
- [ ] 14.1 Simplex → `jpmsimplex.pas`; demo `fpc/simplex/src/`
- [ ] 14.2 Transportation problem → `jpmsimplex.pas`
- [ ] 14.3 Simulated Annealing → `jpmanneal.pas`; demo `fpc/anneal/src/`

### CAT 15 — Series (`series/`)
- [ ] 15.1 Chebyshev series → `jpmchebyshev.pas`
- [ ] 15.2 Verify `InvNormalDist` completeness in `jpmstats.pas`
- [ ] 15.3 Asymptotic error function `AsymErfc` → `jpmspecial.pas`

### CAT 16 — Miscellaneous
- [ ] 16.1 Complex numbers: `TComplex` + ops → `jpmcomplex.pas`
- [ ] 16.2 Physics demos (rel_mass, refract) → `fpc/physics/src/`
- [ ] 16.3 Morse Code demo → `fpc/morse/src/`

---

## Cross-Cutting Tasks
- [ ] X.1 Unit test framework for all `_shared` units
- [ ] X.2 GUI demos (priority: roots, matrices, diffeq, signal)
- [ ] X.3 `README.md` for each `fpc/<topic>/` folder
- [ ] X.4 API review: consistent naming, error handling, `EJPMError` exception

---

## Suggested Priority Order
1. ✅ `jpmroots.pas`
2. `jpmmatrices.pas` ← **NOW**
3. `jpmspecial.pas` (Gamma/Beta — needed by stats)
4. `jpmstats.pas` extensions (Chi2, F-dist, distributions, moments)
5. `jpmdiffeq.pas`
6. `jpmpolynomials.pas`
7. `jpmoptimize.pas`
8. `jpmsignal.pas`
9. `jpmlstsqr.pas`
10. `jpmbessel.pas`
11. `jpmcomplex.pas`
12. Remaining categories
