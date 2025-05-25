# jpmMath Library: Complex Numbers

This folder contains a collection of Free Pascal source code modules dedicated to performing various mathematical operations and algorithms in the complex number domain. It provides functionalities ranging from elementary complex arithmetic to advanced topics such as complex function evaluation, polynomial root finding, and linear algebra operations on complex matrices.

## 1. Core Complex Number Utilities

The foundation of this library lies in robust complex number representation and basic arithmetic operations. The `COMPLEX1` and `COMPLEX2` units provide the fundamental tools for handling complex numbers.

### 1.1. Complex Number Representation

Complex numbers are typically represented using a `Record` type, encompassing both Cartesian (real `R`/`x`, imaginary `I`/`y`) and polar (modulus `r`, argument `t`) forms:

```pascal
Type
  Complex = Record
    R, I: Real { or Double in some units }
    { Optional: r, t: Real for polar form, maintained by conversion procedures }
  End;
```

### 1.2. Basic Operations (`COMPLEX1.pas`, `COMPLEX2.pas`, `ticgt.txt`, `cdetmat.pas`)

These units offer a comprehensive set of procedures and functions for fundamental complex number arithmetic and conversions:

*   **`RectPol(VAR n: COMPLEX)` / `AssignXY(VAR Z:Complex; x,y:DOUBLE)`:** Converts a complex number from Cartesian to polar form.
*   **`PolRect(VAR n: COMPLEX)` / `AssignRT(VAR Z:Complex; r,t:DOUBLE)`:** Converts a complex number from polar to Cartesian form.
*   **`ZSum(z1,z2: COMPLEX; Var z: COMPLEX)` / `CADD(Z1,Z2:Complex; Var Z:Complex)`:** Adds two complex numbers.
*   **`ZMinus(z1,z2: COMPLEX; Var z: COMPLEX)` / `CDIF(Z1,Z2:Complex; Var Z:Complex)`:** Subtracts two complex numbers.
*   **`ZMult(z1,z2: COMPLEX; Var z: COMPLEX)` / `CMUL(Z1,Z2:Complex; Var Z:Complex)`:** Multiplies two complex numbers.
*   **`ZDiv(z1,z2: COMPLEX; Var z: COMPLEX)` / `CDIV(Z1,Z2:Complex; Var Z:Complex)`:** Divides two complex numbers. Includes zero-division check.
*   **`CABS(Z:Complex): Real`:** Calculates the absolute value (modulus) of a complex number.
*   **`ZSqr(z1: COMPLEX; Var z: COMPLEX)`:** Computes the square of a complex number.
*   **`ZSqrt(z1: COMPLEX; Var z: COMPLEX)` / `CSQRT(z:COMPLEX; Var z1:COMPLEX)`:** Computes the square root of a complex number.
*   **`ZInv(z1: COMPLEX; Var z: COMPLEX)`:** Computes the inverse (1/Z) of a complex number.
*   **`ZPower(z1,z2: COMPLEX; Var z: COMPLEX)` / `PowerComplex(VAR Z2:Complex; Z,Z1:Complex)`:** Computes Z1 raised to the power of Z2 (Z1^Z2).
*   **`ZExp(z1: COMPLEX; Var z: COMPLEX)` / `ExpComplex(VAR Z1:Complex; Z:Complex)`:** Computes the complex exponential `e^Z`.
*   **`ZLn(z1: COMPLEX; Var z: COMPLEX)` / `LnComplex(VAR Z1:Complex; Z:Complex)`:** Computes the complex natural logarithm `Ln(Z)`.

### 1.3. Trigonometric and Hyperbolic Functions (`COMPLEX1.pas`, `COMPLEX2.pas`)

The library also includes functions for trigonometric and hyperbolic operations in the complex domain:

*   **`ZSin(z1: COMPLEX; Var z: COMPLEX)` / `SinComplex(VAR Z1:Complex; Z:Complex)`:** Sine of a complex number.
*   **`ZCos(z1: COMPLEX; Var z: COMPLEX)` / `CosComplex(VAR Z1:Complex; Z:Complex)`:** Cosine of a complex number.
*   **`ZTan(z1: COMPLEX; Var z: COMPLEX)` / `TanComplex(VAR Z1:Complex; Z:Complex)`:** Tangent of a complex number.
*   **`ZSh(z1: COMPLEX; Var z: COMPLEX)` / `ShComplex(VAR Z1:Complex; Z:Complex)`:** Hyperbolic sine of a complex number.
*   **`ZCh(z1: COMPLEX; Var z: COMPLEX)` / `ChComplex(VAR Z1:Complex; Z:Complex)`:** Hyperbolic cosine of a complex number.
*   **`ZTh(z1: COMPLEX; Var z: COMPLEX)` / `ThComplex(VAR Z1:Complex; Z:Complex)`:** Hyperbolic tangent of a complex number.
*   **`ZArcsin(z1: COMPLEX; Var z: COMPLEX)` / `ArcSinComplex(VAR Z1:Complex; Z:Complex)`:** Inverse sine of a complex number.
*   **`ZArccos(z1: COMPLEX; Var z: COMPLEX)` / `ArcCosComplex(VAR Z1:Complex; Z:Complex)`:** Inverse cosine of a complex number.
*   **`ZArctan(z1: COMPLEX; Var z: COMPLEX)` / `ArcTanComplex(VAR Z1:Complex; Z:Complex)`:** Inverse tangent of a complex number.
*   **`ZArgsh(z1: COMPLEX; Var z: COMPLEX)` / `ArgShComplex(VAR Z1:Complex; Z:Complex)`:** Inverse hyperbolic sine of a complex number.
*   **`ZArgch(z1: COMPLEX; Var z: COMPLEX)` / `ArgChComplex(VAR Z1:Complex; Z:Complex)`:** Inverse hyperbolic cosine of a complex number.
*   **`ZArgth(z1: COMPLEX; Var z: COMPLEX)` / `ArgThComplex(VAR Z1:Complex; Z:Complex)`:** Inverse hyperbolic tangent of a complex number.

### 1.4. `ucomplex.pas` - Complex Numbers Calculator

The `ucomplex.pas` program provides an interactive, stack-based calculator for complex numbers. It leverages the functions from `COMPLEX1.pas` to perform operations. This serves as a practical demonstration of the `COMPLEX1` unit's capabilities.

**Real-World Applications:** Core complex number utilities are fundamental in electrical engineering (AC circuit analysis), signal processing (Fourier transforms), quantum mechanics, control systems, and fluid dynamics, where quantities are inherently complex or benefit from complex representation.

## 2. Specialized Complex Functions

Beyond elementary operations, the library implements several specialized mathematical functions for complex arguments.

### 2.1. Gamma Function (`mcgama.pas`)

Computes the Gamma function `Γ(z)` or `Ln[Γ(z)]` for a complex argument `z = x + iy`.

*   **Procedure:** `CGAMA(X,Y:Double; KF:Integer; Var GR,GI:Double)`
*   **Inputs:** `X` (real part of z), `Y` (imaginary part of z), `KF` (function code: `0` for `Ln[Γ(z)]`, `1` for `Γ(z)`).
*   **Outputs:** `GR` (real part of result), `GI` (imaginary part of result).
*   **Application Idea:** Useful in probability, statistics, and various areas of physics and engineering where generalized factorials and integrals are involved.

### 2.2. Psi (Digamma) Function (`mcpsi.pas`)

Calculates the Psi (or Digamma) function `ψ(z)` for a complex argument `z = x + iy`.

*   **Procedure:** `CPSI(X,Y:Double; Var PSR,PSI:Double)`
*   **Inputs:** `X` (real part of z), `Y` (imaginary part of z).
*   **Outputs:** `PSR` (real part of `ψ(z)`), `PSI` (imaginary part of `ψ(z)`).
*   **Note:** `ψ(z) = Γ'(z) / Γ(z)`.
*   **Application Idea:** Applied in statistical mechanics, quantum field theory, and special function theory.

### 2.3. Exponential Integral E1(z) (`me1z.pas`)

Computes the complex exponential integral `E1(z)`.

*   **Procedure:** `E1Z(Z:COMPLEX; Var CE1:COMPLEX)`
*   **Input:** `Z` (complex argument).
*   **Output:** `CE1` (complex result `E1(z)`).
*   **Application Idea:** Found in radiative transfer, neutron diffusion, and heat conduction problems.

### 2.4. Legendre Polynomials (`mclpn.pas`, `mclpmn.pas`, `mclqn.pas`, `mclqmn.pas`)

This section provides implementations for various Legendre polynomials and associated functions for complex arguments.

#### 2.4.1. Legendre Polynomials of First Kind (`mclpn.pas`)

Computes the Legendre polynomials `Pn(z)` and their derivatives `Pn'(z)` for a complex argument `z = x + iy`.

*   **Procedure:** `CLPN(N:Integer; X,Y: Double; Var CPN,CPD:pCVEC)`
*   **Inputs:** `N` (degree), `X` (real part of z), `Y` (imaginary part of z).
*   **Outputs:** `CPN` (array of `Pn(z)` values), `CPD` (array of `Pn'(z)` values).
*   **Application Idea:** Used in electrostatics, quantum mechanics (angular momentum), and numerical analysis for approximations.

#### 2.4.2. Associated Legendre Functions of First Kind (`mclpmn.pas`)

Computes the associated Legendre functions `Pmn(z)` and their first derivatives `Pmn'(z)` for a complex argument `z = x + iy`.

*   **Procedure:** `CLPMN(M,N:Integer; X,Y:Double; Var CPM,CPD:pCTab)`
*   **Inputs:** `M` (order, `0 <= m <= n`), `N` (degree), `X` (real part of z), `Y` (imaginary part of z).
*   **Outputs:** `CPM` (table of `Pmn(z)` values), `CPD` (table of `Pmn'(z)` values).
*   **Application Idea:** Crucial in solving partial differential equations in spherical coordinates, common in physics and engineering (e.g., multipole expansions, quantum mechanics).

#### 2.4.3. Legendre Polynomials of Second Kind (`mclqn.pas`)

Computes the Legendre polynomials of the second kind `Qn(z)` and their derivatives `Qn'(z)` for a complex argument `z = x + iy`.

*   **Procedure:** `CLQN(N:Integer; X,Y: Double; CQN,CQD:pCVEC)`
*   **Inputs:** `N` (degree), `X` (real part of z), `Y` (imaginary part of z).
*   **Outputs:** `CQN` (array of `Qn(z)` values), `CQD` (array of `Qn'(z)` values).
*   **Application Idea:** Often appear in solutions to boundary value problems in potential theory and wave propagation, complementing the first kind Legendre polynomials.

#### 2.4.4. Associated Legendre Functions of Second Kind (`mclqmn.pas`)

Computes the associated Legendre functions of the second kind `Qmn(z)` and their first derivatives `Qmn'(z)` for a complex argument `z = x + iy`.

*   **Procedure:** `CLQMN(M,N:Integer; X,Y: Double; Var CQM,CQD:pCTab)`
*   **Inputs:** `M` (order), `N` (degree), `X` (real part of z), `Y` (imaginary part of z).
*   **Outputs:** `CQM` (table of `Qmn(z)` values), `CQD` (table of `Qmn'(z)` values).
*   **Application Idea:** Similar to `Pmn(z)`, these are used in more complex boundary value problems and advanced physics applications involving spherical harmonics.

## 3. Polynomial Operations

This section covers various algorithms for evaluating complex polynomials and finding their roots.

### 3.1. Complex Polynomial Evaluation (Horner's Rule - `chorner.pas`)

Evaluates a complex polynomial `P(Z)` at a given complex argument `Z` using Horner's rule, which is an efficient algorithm for polynomial evaluation.

*   **Procedure:** `CHorner(A:pVec; N: integer; X0: COMPLEX; Var Y0:COMPLEX)`
*   **Inputs:** `A` (pointer to complex coefficients array, `A^[1]` to `A^[N+1]` where `A^[N+1]` is the highest degree coefficient), `N` (order of polynomial), `X0` (complex argument Z).
*   **Output:** `Y0` (complex result `P(Z)`).
*   **Application Idea:** Fundamental in many numerical methods, especially iterative root-finding algorithms, where polynomial values need to be computed repeatedly.

### 3.2. Complex Root Counting (`rootnum.pas`)

Calculates the number of complex roots within a specified circle by counting (u,v) transitions around its circumference. It implements a graphical method for root localization.

*   **Procedure:** `RootNum` (uses global variables for input/output).
*   **Inputs (global):** `w` (radius of search circle), `x0,y0` (center of the circle), `m` (evaluation points per quadrant).
*   **Outputs (global):** `nn` (number of complete cycles/roots found), `a` (residual, non-zero indicates a failure).
*   **`Eval(x,y:DOUBLE;VAR u,v:DOUBLE)`:** An example complex function `f(z) = z^2 + 1` used for demonstration.
*   **Application Idea:** Useful for initial root localization before applying more precise root-finding algorithms, especially for functions with many roots.

### 3.3. Root Finding Algorithms

#### 3.3.1. Newton's Method (`znewton.txt`, `znewton.pas`, `tnewton.pas`)

Implements Newton's iterative method in the complex domain for finding roots of analytic functions `f(z) = µ(x,y) + i v(x,y)`. It requires the function and its partial derivatives.

*   **Procedure:** `ZNewton` (demonstrated in `znewton.pas`).
*   **Inputs (global):** `x0, y0` (initial guess), `e` (convergence criterion), `n` (maximum iterations).
*   **Outputs (global):** `x, y` (calculated root), `k` (actual iterations).
*   **`Eval(x,y:double;VAR u,v,u1,v1,u2,v2:double)`:** An example function `f(z) = z^2 + 1` providing `u,v` and their partial derivatives `du/dx, dv/dx, du/dy, dv/dy` (referred as `u1, v1, u2, v2` in code comments).
*   **`CNEWTON (N:Integer; A,B: pVEC; Var X,Y:pVEC; Var IER:Integer; G,H,T,F:pVEC)` (`tnewton.pas`):** Finds all roots of a complex polynomial using Newton's iterative formulation.
*   **Inputs:** `N` (order of polynomial), `A, B` (real and imaginary parts of coefficients in decreasing powers).
*   **Outputs:** `X, Y` (real and imaginary parts of roots), `IER` (error code).
*   **Application Idea:** Fast convergence for well-behaved functions when a good initial guess is available. Widely used in optimization, numerical analysis, and solving non-linear equations.

#### 3.3.2. Mueller's Method (`zmueller.txt`, `zmueller.pas`)

Applies Mueller's one-dimensional parabolic method directly to the complex plane to search for roots. It's robust for functions with distinct roots.

*   **Procedure:** `ZMueller` (demonstrated in `zmueller.pas`).
*   **Inputs (global):** `x0, y0` (initial guess), `b1, b2` (bounds on error), `e` (convergence criterion), `n` (maximum iterations).
*   **Outputs (global):** `x, y` (calculated zero), `k` (iterations performed).
*   **`Z(x,y:DOUBLE; VAR u,v:DOUBLE)`:** An example function `f(z) = z^2 + 1` to find roots.
*   **Application Idea:** Effective for finding roots of analytic and non-analytic functions, particularly when Newton's method struggles (e.g., near local extrema or multiple roots).

#### 3.3.3. Laguerre's Method (`tclague.pas`)

Finds all roots of a complex polynomial using Laguerre's iterative formulation, known for its good convergence properties.

*   **Procedure:** `CLAGUE(N: Integer; A: CTab; ITMAX: Integer; EPS,EPS2: Double; Var IMP: Integer; Var X: CTab; B, C, D: CTab)`
*   **Inputs:** `N` (order of polynomial), `A` (complex coefficients in decreasing power), `ITMAX` (max iterations per root), `EPS` (min relative error), `EPS2` (max relative error).
*   **Outputs:** `IMP` (flag: `0` for convergence, `1` for no convergence), `X` (array of found complex roots).
*   **Application Idea:** A reliable method for finding all roots of polynomials, suitable for problems in control theory, digital filter design, and signal processing.

#### 3.3.4. Bauhuber's Method (`fbauhube.pas`, `tbauhube.pas`)

Uses Bauhuber's Method to find all real or complex roots of a polynomial of degree `n`, including complex coefficients. It includes polynomial scaling options for stability and performance.

*   **Function:** `bauhub(real0, scale, n: integer; ar, ai: VEC; Var rootr, rooti, absf: VEC):Integer`
*   **Inputs:** `real0` (0: complex coeffs, <>0: real coeffs), `scale` (0: no scaling, <>0: scaling), `n` (degree), `ar, ai` (real/imaginary parts of coefficients).
*   **Outputs:** `rootr, rooti` (real/imaginary parts of roots), `absf` (absolute value of function at roots).
*   **Return Value:** `0` (ok), `1` (invalid input), `2` (leading coefficient zero), `3` (iteration max exceeded).
*   **Dependencies:** `scpoly`, `bauroot`, `chorner`, `polydiv` (internal procedures within the unit).
*   **Application Idea:** Highly versatile for finding roots of any polynomial, crucial for stability analysis in engineering systems or solving characteristic equations.

### 3.4. Roots of Specific Degree Equations

#### 3.4.1. Second Degree Equation (`tequa2.pas`)

Solves quadratic equations `az^2 + bz + c = 0` with complex coefficients.

*   **Procedure:** `equa2c(a,b,c: COMPLEX; VAR s1, s2: COMPLEX)`
*   **Inputs:** `a, b, c` (complex coefficients).
*   **Outputs:** `s1, s2` (complex roots).
*   **Application Idea:** Essential for basic circuit analysis, mechanical vibration problems, and finding critical points in complex functions.

#### 3.4.2. Third Degree Equation (`croot3.pas`)

Solves cubic equations `az^3 + bz^2 + cz + d = 0` with complex coefficients using Cardano's method adapted for complex numbers.

*   **Procedure:** `Croot3(a,b,c,d:Complex; Var z1,z2,z3:Complex)`
*   **Inputs:** `a, b, c, d` (complex coefficients).
*   **Outputs:** `z1, z2, z3` (complex roots).
*   **Application Idea:** Used in advanced mechanics, quantum chemistry, and solving equations of state in thermodynamics.

#### 3.4.3. Fourth Degree Equation (`tcroot4.pas`)

Solves quartic equations `az^4 + bz^3 + cz^2 + dz + e = 0` with complex coefficients. It reduces the quartic to a cubic equation and then solves two quadratic equations.

*   **Procedure:** `Croot4(a,b,c,d,e:Complex; Var z1,z2,z3,z4:Complex)`
*   **Inputs:** `a, b, c, d, e` (complex coefficients).
*   **Outputs:** `z1, z2, z3, z4` (complex roots).
*   **Application Idea:** Relevant in kinematics, signal processing, and control system design where higher-order characteristic equations may arise.

## 4. Linear Algebra with Complex Matrices

This section provides routines for fundamental linear algebra operations involving complex matrices.

### 4.1. Matrix Inversion (Gauss Method - `ticgt.txt`)

Calculates the inverse of a complex square matrix using the Gauss Method with full pivoting. The sample run demonstrates its use and verifies the result by multiplying the inverse with the original matrix.

*   **Procedure:** `ICGT(eps:Real; N:integer; A:Matc; Var it:integer; Var AM1:Matc)`
*   **Inputs:** `eps` (required precision), `N` (size of matrix), `A` (input complex matrix).
*   **Outputs:** `it` (`0` if singular, `1` if regular), `AM1` (inverse complex matrix).
*   **Dependencies:** `TSCGT` (triangularization) and `BSCGT` (back substitution) procedures within the same file.
*   **Application Idea:** Crucial for solving systems of linear equations, calculating determinants, and in control systems where matrix inversion is frequently required.

### 4.2. Determinant of a Complex Square Matrix

#### 4.2.1. By Gauss Method (`cdetmat.pas`)

Computes the complex determinant of a square matrix using Gauss elimination with full pivoting.

*   **Procedure:** `DCGT(eps:Real; N:integer; A:Matc; Var det:Complex)`
*   **Inputs:** `eps` (required precision), `N` (size of matrix), `A` (complex matrix).
*   **Output:** `det` (complex determinant).
*   **Dependencies:** `TSCGT` (triangularization) procedure within the same file.
*   **Application Idea:** Essential for checking matrix invertibility, analyzing system stability, and in tensor analysis.

#### 4.2.2. Using Function `CFindDet` (`cfinddet.pas`)

Calculates the determinant of a complex square matrix by converting it to upper triangular form through row operations and then multiplying the diagonal elements.

*   **Procedure:** `CFindDet(matrix: Matc; n: integer; Var cdet: COMPLEX)`
*   **Inputs:** `matrix` (complex matrix), `n` (size).
*   **Output:** `cdet` (complex determinant).
*   **Application Idea:** Provides an alternative method for determinant calculation, useful for verifying results or where a specific algorithm preference exists.

### 4.3. Solving Complex Linear Systems

#### 4.3.1. By Gauss-Jordan Method (`csysmat.pas`)

Solves a complex linear matrix system `AX=B` using the Gauss-Jordan method with full pivoting. It also calculates the inverse of A and its determinant.

*   **Procedure:** `CMATINV(N,M:integer; VAR AA:MATC; VAR BB:MATC1;VAR DET:Complex)`
*   **Inputs:** `N` (size of A), `M` (number of right-hand side vectors/columns in B), `AA` (complex matrix A, modified in-place), `BB` (complex matrix B, modified in-place to X).
*   **Outputs:** `AA` (inverse of A), `DET` (complex determinant of A), `BB` (solution matrix X).
*   **Application Idea:** Widely used in engineering and scientific computing for structural analysis, circuit simulation, and solving systems of coupled equations.

#### 4.3.2. Homogeneous Linear Systems (`rshcgt.pas`)

Solves a complex homogeneous linear system `AX = 0` using Gauss Method with full pivoting.

*   **Procedure:** `RSHCGT(eps:real; N:integer; A:Matc; Var R0,M0:integer; Var VX:Matc)`
*   **Inputs:** `eps` (required precision), `N` (size of A), `A` (complex matrix, modified in-place).
*   **Outputs:** `R0` (rank of A), `M0` (dimension of solution space), `VX` (solution vectors, stored in columns).
*   **Application Idea:** Essential for finding null spaces of matrices, which has applications in determining linear dependencies, solving underdetermined systems, and analyzing network flows.

#### 4.3.3. By LU Decomposition (`test_clu.pas`, `clu.pas`)

Solves a complex linear system `AX = B` using LU decomposition. The `CLU` unit provides the `LUDCMP` (LU decomposition) and `LUBKSB` (back substitution) routines.

*   **Unit:** `CLU`
*   **Procedure:** `LUDCMP(VAR A:pCVEC; N:INTEGER; VAR INDX:pIVEC; VAR D:INTEGER; VAR CODE:INTEGER)`:
    *   **Inputs:** `A` (complex matrix stored as vector), `N` (size).
    *   **Outputs:** `A` (LU decomposition), `INDX` (row permutation), `D` (sign of determinant), `CODE` (error code).
*   **Procedure:** `LUBKSB(VAR A:pCVEC; N:INTEGER; INDX: pIVEC; VAR B:pCVEC)`:
    *   **Inputs:** `A` (LU decomposed matrix), `N` (size), `INDX` (row permutation), `B` (right-hand side vector).
    *   **Output:** `B` (solution vector X).
*   **Application Idea:** Efficient for solving multiple linear systems with the same matrix A but different right-hand sides B, commonly found in time-dependent simulations and iterative numerical methods.

#### 4.3.4. Matrix Inversion by LU Decomposition (`inv_lu.pas`)

Performs the inversion of a complex square matrix using LU decomposition, building upon the `CLU` unit.

*   **Program:** `Inversion_LU`
*   **Inputs:** Reads matrix A from a data file.
*   **Outputs:** Prints the inverted matrix Y and a verification `A*Y = I`.
*   **Dependencies:** Uses `CLu` unit for `LUDCMP` and `LUBKSB`.
*   **Application Idea:** Alternative to Gauss-Jordan for matrix inversion, often preferred for its numerical stability and efficiency for large matrices.

#### 4.3.5. Solving Linear Systems by Conversion to Real System (`zmatsys.pas`)

Solves a complex linear system `AZ=B` by converting it into an equivalent real linear system of twice the dimension, then solving the real system using Gauss-Jordan elimination.

*   **Program:** `Zmatsys`
*   **Input:** Complex matrix A and vector B (hardcoded in example).
*   **Output:** Complex solution vector Z and determinant of the real system.
*   **Dependencies:** `MATINV` (Gauss-Jordan for real matrices, internal to the program).
*   **Application Idea:** Provides a conceptual alternative for complex linear system solution, potentially useful for environments optimized for real matrix operations.

### 4.4. Eigenvalues and Eigenvectors of Complex Matrices

#### 4.4.1. Using QR Algorithm (`tceigen.pas`)

Calculates the eigenvalues and eigenvectors of a general complex square matrix using the QR algorithm. The method involves reduction to Hessenberg form by unitary transformations, followed by the QR algorithm.

*   **Procedure:** `CEIGEN(Var AR,AI:pMAT;NM,N:Integer;AVEC:Boolean; Var WR,WI:pVEC; Var ZR,ZI:pMAT; Var WORK:pVEC; Var IERR: Integer)`
*   **Inputs:** `AR, AI` (real/imaginary parts of matrix A), `NM` (1st dimension), `N` (2nd dimension), `AVEC` (true for eigenvectors), `WORK` (working space).
*   **Outputs:** `WR, WI` (real/imaginary parts of eigenvalues), `ZR, ZI` (real/imaginary parts of eigenvectors), `IERR` (error code).
*   **Dependencies:** `CBALAN`, `CUNITH`, `COMQR2`, `CBALBK`, `NORMAL` (internal procedures).
*   **Application Idea:** Fundamental in stability analysis of dynamical systems, quantum mechanics, and principal component analysis for complex data.

#### 4.4.2. Using Jacobi Method (`tcomeig.pas`)

Calculates the eigenvalues and eigenvectors of a complex matrix using the Jacobi method, which is also recommended for real matrices with complex eigenvalues.

*   **Procedure:** `COMEIG (NDIM,N,NN:Integer; A,Z,T,U:pMAT; Var IER:Integer; EN:pVEC)`
*   **Inputs:** `NDIM, N` (dimensions), `NN` (max iterations), `A, Z` (real/imaginary parts of matrix, modified in-place).
*   **Outputs:** `A(J,J), Z(J,J)` (eigenvalues), `T, U` (eigenvectors in columns), `IER` (error code), `EN` (working zone).
*   **Application Idea:** Provides an alternative to QR for smaller, dense matrices, often preferred for its simplicity and robustness in certain scenarios like solving problems in solid-state physics and quantum computing.

## 5. File Manifest

This section lists the primary source files and their corresponding functionalities within the library.

-   **`_info.txt`**: Program descriptions and overview.
-   **`ticgt.txt`**: Pascal program demonstrating the calculation of a complex square matrix inverse using Gauss Method with full pivoting. (Note: Although a .txt file, it contains Pascal code and sample run output).
-   **`zmueller.txt`**: Explanation file for the `ZMueller` program, detailing Mueller's Method in the complex plane. (Contains documentation, not executable code).
-   **`znewton.txt`**: Explanation file for the `ZNewton` program, detailing Newton's Method in the complex domain. (Contains documentation, not executable code).
-   **`cdetmat.pas`**: Program to calculate the determinant of a complex square matrix by Gauss Method with full pivoting.
-   **`cfinddet.pas`**: Program to calculate the determinant of a complex square matrix using the `CFindDet` function (upper triangularization method).
-   **`chorner.pas`**: Program to evaluate a complex polynomial by Horner's rule.
-   **`clu.pas`**: Pascal unit containing LU decomposition routines (`LUDCMP`, `LUBKSB`) for complex matrices, used by `test_clu.pas` and `inv_lu.pas`.
-   **`complex1.pas`**: Pascal unit providing basic complex number operations and functions (e.g., arithmetic, exp, ln, power, trig, hyperbolic functions). Used by `ucomplex.pas`.
-   **`complex2.pas`**: Pascal unit providing elementary operations with complex numbers, similar to `complex1.pas` but with slightly different procedure names. Used by `tcomplex.pas`.
-   **`croot3.pas`**: Program to solve a 3rd degree equation with complex coefficients.
-   **`csysmat.pas`**: Program to solve a complex linear system `AX=B` by Gauss-Jordan Method.
-   **`fbauhube.pas`**: Pascal unit implementing Bauhuber's Method to find all real or complex roots of a polynomial of degree n. Used by `tbauhube.pas`.
-   **`inv_lu.pas`**: Program for the inversion of a complex square matrix by LU decomposition.
-   **`mcgama.pas`**: Program to calculate the Gamma Function with a complex argument.
-   **`mclpmn.pas`**: Program to calculate the Associated Legendre Functions of the First Kind and their First Derivatives for a Complex Argument.
-   **`mclpn.pas`**: Program to calculate the Legendre Polynomials of the First Kind for a Complex Argument.
-   **`mclqmn.pas`**: Program to calculate the Associated Legendre Functions of the Second Kind and their First Derivatives for a Complex Argument.
-   **`mclqn.pas`**: Program to calculate the Legendre Polynomials of the Second Kind for a Complex Argument.
-   **`mcpsi.pas`**: Program to calculate the Psi function for a complex argument.
-   **`me1z.pas`**: Program to calculate the complex exponential integral E1(z) with a complex argument.
-   **`rootnum.pas`**: Program to demonstrate the complex root counting subroutine.
-   **`rshcgt.pas`**: Program to solve a complex homogeneous linear system by Gauss Method with full pivoting.
-   **`tbauhube.pas`**: Test program for Bauhuber's method for finding polynomial roots.
-   **`tceigen.pas`**: Test program for eigenvalues and eigenvectors of a general complex square matrix using the QR algorithm.
-   **`tclague.pas`**: Test program for finding all roots of a complex polynomial using Laguerre formulation.
-   **`tcomeig.pas`**: Test program for eigenvalues/eigenvectors of a complex matrix by the Jacobi method.
-   **`tcomplex.pas`**: Test program for elementary operations on complex numbers using `complex2.pas`.
-   **`tcroot4.pas`**: Program to solve a 4th degree equation with complex coefficients.
-   **`tequa2.pas`**: Program to find roots of a second degree equation with complex coefficients.
-   **`test_clu.pas`**: Test program for solving a complex linear system by LU decomposition.
-   **`tnewton.pas`**: Test program for finding all roots of a complex polynomial using Newton's iterative formulation.
-   **`trslcgtc.pas`**: Program to solve a complex linear system by Gauss Method with full pivoting and correction process.
-   **`ucomplex.pas`**: Complex numbers calculator (interactive program).
-   **`zcircle.pas`**: Program to demonstrate the zero searching algorithm (using a shrinking circle method).
-   **`zmatsys.pas`**: Program for solving a complex linear system by converting it to a real system.
-   **`zmueller.pas`**: Program to demonstrate the complex domain Mueller's subroutine.
-   **`znewton.pas`**: Program to demonstrate the Newton root subroutine in the complex domain.

## 6. References

The algorithms and implementations in this library are based on various numerical methods and texts, as referenced in the source code comments:

-   `[BIBLI 01]`: BASIC Scientific Subroutines, Vol. II by F.R. Ruckdeschel, BYTE/McGRAW-HILL, 1981.
-   `[BIBLI 03]`: Mathématiques en Turbo-Pascal By M. Ducamp and A. Reverchon (vol 2), Eyrolles, Paris, 1988.
-   `[BIBLI 05]`: Mathématiques en Turbo-Pascal By M. Ducamp and A. Reverchon (vol 2), Eyrolles, Paris, 1988. (Note: Same as [BIBLI 03], referenced by different files).
-   `[BIBLI 11]`: Numerical algorithms with C, By Gisela Engeln-Muellges and Frank Uhlig, Springer-Verlag, 1996.
-   `[BIBLI 18]`: Numath Library By Tuan Dang Trong in Fortran 77.
-   `[Others]`: Algèbre, Algorithmes et programmes en Pascal By Jean-Louis Jardrin, DUNOD Paris, 1988; Fortran Routines for Computation of Special Functions (jin.ece.uiuc.edu/routines/routines.html).

