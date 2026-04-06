
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
- [x] **CAT 1 ALL** — `jpmstats.pas` extended: ChiSquareDist/InvChiSquare, FDist/InvFDist, BinomialProb/CDF, PoissonProb/CDF, Mean/Variance/StdDev/Skewness/Kurtosis/Median ✅ COMMITTED (c158916)
- [x] **CAT 1.3** — `jpmspecial.pas`: GammaFunc, LnGamma, BetaFunc, LnBeta, IncompleteBeta, IncompleteGammaP/Q, ErfFunc, ErfcFunc ✅ COMMITTED (13ec193)
- [x] **CAT 2.1** — `jpmderivative.pas`: Deriv1/2/3/4, DerivRichardson, DerivRomberg, Gradient, Jacobian ✅ COMMITTED (290ca5d)
- [x] **CAT 2.2** — `jpmintegration.pas` extended: TrapezoidIntegral, SimpsonIntegral, AdaptiveSimpson, GaussLegendre ✅ COMMITTED (a49ca43)
- [x] **CAT 2.3** — `jpmchebyshev.pas`: ChebFit, ChebEval, ChebInteg, ChebDeriv ✅ COMMITTED (9173acc)
- [x] **CAT 2.4** — `jpmcontinued.pas`: FloatToCF, CFToFloat, CFConvergents, GenCFEval, SqrtCF, BestRational ✅ COMMITTED (71f08a6)
- [x] **CAT 4 ALL** — `jpmoptimize.pas`: GoldenSection, BrentMin, BracketMin, NelderMead, SteepestDescent, PowellMin ✅ COMMITTED (a10f1d8, fa82776)
- [x] **CAT 5 ALL** — `jpmmatrices.pas`: LU decompose/solve/inverse/det, Gauss-Seidel, Cholesky, Thomas tridiagonal, Jacobi eigenvalues ✅ COMMITTED (451132e)
- [x] **CAT 6 ALL** — `jpmlstsqr.pas`: LeastSquares, PolyRegression, LinearRegression, MultipleRegression, ChiSquareFit, LevenbergMarquardt ✅ COMMITTED (6074b62, 948cb2d)
- [x] **CAT 7 ALL** — `jpmdiffeq.pas`: RK4, RKF45, AdamsBashforth4, AdamsMoulton4, GearStiff (Rosenbrock), BulirschStoer ✅ COMMITTED (273f1f5, 1e698e0, b8c5274)
- [x] **CAT 8 ALL** — `jpmpolynomials.pas`: PolyEval/Add/Sub/Mul/Div/Deriv/Integ/GCD, TPolyFrac, PolyFracPartial ✅ COMMITTED (654b419, 102eab5)
- [x] **CAT 9.1–9.5,9.7** — `jpmspecial.pas` extended: LegendrePn/Qn, HermiteHn, LaguerreLn, AiryAi/Bi, EllipticK/E, Hypergeometric2F1 ✅ COMMITTED (e86c7c0)
- [x] **CAT 9.6** — `jpmbessel.pas`: BesselJ0/J1/Y0/Y1/I0/I1/K0/K1, BesselJn, BesselJ0Zero/J1Zero ✅ COMMITTED (e3dd57d)
- [x] **CAT 10 ALL** — `jpmsort.pas`: BubbleSort, InsertionSort, SelectionSort, QuickSort, MergeSort, HeapSort, ShellSort, BinarySearch, RankSort ✅ COMMITTED (7ef74cb)
- [x] **CAT 11 ALL** — `jpmsignal.pas`: DFT, FFT (Cooley-Tukey), PowerSpectrum, MovingAverage, SavitzkyGolay5, FourierCoeff ✅ COMMITTED (099ea83)
- [x] **CAT 12 ALL** — `jpmarith.pas`: GCD/LCM, ExtGCD, Factorial/BinomCoeff/Permutation, TFraction, IntToBase/BaseToInt, IsPrime, PrimeSieve, Factorize, SolveDiophantine ✅ COMMITTED (e4ad8eb)
- [x] **CAT 13 ALL** — `jpmgeometry.pas`: SolveTriangle, SolveArcCircle, PointDistance2D/3D, TriangleArea, PolygonArea, CircleArea/Perimeter/Chord/ArcLen, EllipseArea/Perimeter, HyperbolaAsymptoteAngle ✅ COMMITTED (8fef706)
- [x] **CAT 14.1** — `jpmsimplex.pas`: SimplexMaximize, SimplexMinimize ✅ COMMITTED (fcb64d3)
- [x] **CAT 14.3** — `jpmanneal.pas`: SimulatedAnnealing (N-dim continuous, geometric cooling) ✅ COMMITTED (679e39d)
- [x] **CAT 16.1** — `jpmcomplex.pas`: TComplex record + full arithmetic, transcendentals, trig ✅ COMMITTED (381fcbd)

---

## Remaining Tasks

### CAT 14 — Linear Programming (`linearprog/`)
- [ ] 14.2 Transportation problem → `jpmsimplex.pas`

### CAT 15 — Series (`series/`)
- [ ] 15.2 Verify `InvNormalDist` completeness in `jpmstats.pas`
- [ ] 15.3 Asymptotic error function `AsymErfc` → `jpmspecial.pas`

### CAT 16 — Miscellaneous
- [ ] 16.2 Physics demos (rel_mass, refract) → `fpc/physics/src/`
- [ ] 16.3 Morse Code demo → `fpc/morse/src/`

---

## Cross-Cutting Tasks
- [ ] X.1 Unit test framework for all `_shared` units
- [ ] X.2 GUI demos (priority: roots, matrices, diffeq, signal)
- [ ] X.3 `README.md` for each `fpc/<topic>/` folder
- [ ] X.4 API review: consistent naming, error handling, `EJPMError` exception
