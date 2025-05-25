
# jpmMath Library - A Comprehensive Math Toolkit by Jean-Pierre Moreau

## Introduction
The `jpmMath` library is a robust and extensive collection of numerical and scientific computing routines implemented in Pascal. Designed for engineers, scientists, researchers, and students, this library provides a foundational toolkit for solving a wide array of mathematical, physics, and engineering problems. It emphasizes direct mathematical implementation, numerical accuracy, and modularity, making it an excellent resource for both practical applications and educational purposes in understanding classical algorithms.

## License ##
Jean-Pierre Moreau does not explicitely specify any conditions how his code can be used. After inquiry about this (https://www.lazarusforum.de/viewtopic.php?p=150255#p150255), his son responded in a mail:

```
Hello,

My father, now very old, who created this website to capitalize and share his passion for
developing numerical computing algorithms, would be very happy if his work could benefit other people.

You can recover and reuse his work with great pleasure.

Kind regards,
François-Xavier Moreau
```

## Core Principles and Philosophy

The `jpmMath` library adheres to the following core principles:

*   **Direct Mathematical Implementation**: Algorithms are typically implemented directly from their mathematical formulations, fostering a deep understanding of their inner workings.
*   **Numerical Accuracy and Stability**: Focuses on robust numerical methods, often incorporating techniques for error control, precision management, and stability enhancement (e.g., pivoting, iterative refinement, appropriate choice of algorithms for specific problem types).
*   **Modularity and Reusability**: The library is organized into self-contained units and programs, promoting reusability of core functionalities across different applications and simplifying maintenance.
*   **Educational Value**: Serves as a practical learning resource for numerical analysis, linear algebra, differential equations, statistics, and various scientific computing disciplines. The clear structure and detailed comments in the source code, coupled with comprehensive documentation, aim to facilitate understanding and further development.
*   **Problem-Oriented Design**: Each module and program is designed to address specific classes of real-world problems, from basic arithmetic and complex numbers to advanced topics like finite element analysis and shock spectrum calculations.

## Library Modules Overview

The `jpmMath` library is logically organized into several modules, each focusing on a distinct area of numerical and scientific computing. Click on each module name to access its detailed documentation.

### [arith - Basic Arithmetic and Number Theory](arith/README.md)

This module provides fundamental tools for tasks such as numerical base conversions, high-precision arithmetic, equation solving (Diophantine), combinatorial analysis, and prime number operations. It emphasizes direct mathematical implementation within the Pascal environment.

### [bessel - Bessel Functions](bessel/README.md)

This module focuses on the computation of various Bessel functions (real, spherical, and complex), their derivatives, and their roots. These functions are critical in solving problems across wave propagation, heat conduction, vibrations, and quantum mechanics, particularly those involving cylindrical or spherical symmetry.

### [complex - Complex Numbers and Linear Algebra](complex/README.md)

This module offers a comprehensive set of functionalities for complex number arithmetic, evaluation of specialized complex functions (Gamma, Psi, Exponential Integral), polynomial operations (evaluation, root finding using various methods), and linear algebra operations on complex matrices (inversion, determinants, solving systems, eigenvalues).

### [diffeqa - Differential Equations](diffeqa/README.md)

This module provides a collection of numerical methods for solving ordinary differential equations (ODEs) as initial value problems (IVPs) and boundary value problems (BVPs), as well as partial differential equations (PDEs). It includes methods like Runge-Kutta, Adams-Bashforth, Bulirsch-Stoer, and solvers for stiff systems, emphasizing accuracy and adaptive step-size control.

### [functions1 - Mathematical Functions (Part 1)](functions1/README.md)

This module covers a broad range of numerical routines including various interpolation methods (Akima Spline, Chebyshev, Lagrange, Newton, Cubic Spline, Rational, Trigonometric), numerical differentiation (formal, Lagrange, Nth derivative, Romberg), numerical integration (Gauss, Simpson, Romberg, Clenshaw-Curtis, Primitive), multidimensional integration (cubature), and optimization algorithms (Bracketing, Golden Section, Brent, Nelder-Mead Simplex, Steepest Descent). It also includes tools for Taylor series expansion and function plotting.

### [functions2 - Mathematical Functions (Part 2)](functions2/README.md)

This module offers implementations for a variety of specialized mathematical functions and further optimization techniques. It includes calculations for Bernoulli and Euler numbers, Airy functions and their zeros, Beta function, Cosine and Sine Integrals, Mathieu functions, Confluent Hypergeometric functions, Struve functions (H and L series), hyperbolic and inverse hyperbolic functions, and complete elliptic integrals. It also provides methods for Hermite, Lagrange, and Legendre polynomial coefficients.

### [geometry - Geometric and Astronomical Calculations](geometry/README.md)

This module provides programs for solving specific geometrical and astronomical problems. It includes tools for calculating arc parameters, reducing conic equations, finding conics from five points, demonstrating Euler's circle, determining planetary positions, computing polygon surface areas, and resolving triangles.

### [linearprog - Linear Programming and Operational Research](linearprog/README.md)

This module focuses on operational research and linear programming methods. It includes implementations for the Appointment Method (assignment problems), Dantzig's Model (optimal path finding), Simplex Method (basic and general forms), Time P.E.R.T. Model (project scheduling), Transport Model (resource dispatching), and the Traveling Salesman Problem using Simulated Annealing.

### [lstsqr - Least Squares Approximation](lstsqr/README.md)

This module offers a comprehensive collection of programs for least squares approximation and regression analysis. It supports linear, parabolic, and Nth-order polynomial regression, multidimensional least squares, orthogonal polynomial fitting, parametric (non-linear) least squares, and iterated regression for error reduction. It also includes Chi-square statistic calculation and multi-dimensional curve fitting via the Simplex method.

### [matrices - Numerical Linear Algebra](matrices/README.md)

This module provides a versatile collection of routines for various numerical linear algebra operations. It includes solvers for general dense, iterative, symmetric, and specialized linear systems (banded, tridiagonal, Vandermonde, Toeplitz). It also covers matrix inversion (LU, Householder, Cholesky), determinant calculation, and a variety of algorithms for computing eigenvalues and eigenvectors (Power methods, Jacobi, QL, QR algorithms, Characteristic Polynomials).

### [mechanics - Mechanical Engineering Simulations](mechanics/README.md)

This module contains programs designed to solve various problems in mechanical engineering, including structural analysis, vibrations, and composite materials. It covers 1-DOF and N-DOF system responses, beam analysis, composite layer stress/stiffness, bouncing ball simulations, and a Finite Element Method (FEM) demonstration for 3D structures.

### [miscellaneous - Diverse Computational and Algorithmic Concepts](miscellaneous/README.md)

This module is a collection of standalone Pascal programs covering a diverse range of computational and algorithmic concepts. It includes fundamental physics calculations (relativistic mass, Snell's Law), astronomical phenomena (stellar magnitude), classical cryptography techniques (transposition, substitution ciphers, Morse code), and demonstrations of cellular automata (Sierpinski Triangle).

### [polynomials - Algebraic Operations on Polynomials](polynomials/README.md)

This module provides a collection of routines for performing various algebraic operations on polynomials and polynomial fractions. It includes basic arithmetic (addition, multiplication), differentiation, division (Euclidian, increasing powers), greatest common divisor (GCD), smallest common multiple (SCM), substitution, and symbolic expansion. It handles real, integer, and fractional coefficients.

### [roots - Root Finding Methods](roots/README.md)

This module offers a comprehensive collection of numerical methods for finding roots of functions. It includes single-variable root-finding algorithms (Bisection, Newton, Secant, Brent, Mueller, Pegasus, Lagrange, Aitken Acceleration, Zeroin), polynomial root-finding methods (Bairstow, Bernouilli, Lin, direct solvers for degree 2,3,4), and solvers for systems of non-linear equations (Brown, Levenberg-Marquardt, NNES, 2D Mueller).

### [series - Series Expansions and Special Functions](series/README.md)

This module focuses on numerical methods related to series expansions and special functions. It includes calculations for asymptotic error functions, Chebyshev economization and series coefficients, Chi-square distribution functions (PDF and CDF), logarithm of factorial, complex series evaluation, Horner's rule for polynomial shifting, inverse normal distribution, polynomial inversion, and series reversion.

### [signal - Signal Processing Applications](signal/README.md)

This module provides a robust collection of tools for various signal processing applications. It covers Fourier analysis (Fourier series, FFT), numerical filtering (Butterworth), signal smoothing (FFT-based, Savitzky-Golay), single-degree-of-freedom oscillator response, shock spectrum analysis, and deconvolution.

### [sorting - Sorting and Searching Algorithms](sorting/README.md)

This module contains implementations of various classical sorting algorithms (Bubble Sort, Straight Insertion Sort, Shell Sort, Heapsort, Quicksort, Merge Sort) and searching algorithms (Linear Search, Binary Search). It provides fundamental tools for data organization and retrieval, emphasizing time and space complexity characteristics.

### [stat - Statistical Functions](stat/README.md)

This module offers a suite of Pascal programs for performing statistical computations. It includes implementations for special mathematical functions (Gamma, Incomplete Beta), various probability distributions (Binomial, Poisson, Normal, Chi-square, Student's T, F-distribution), descriptive statistics and moments (mean, deviation, skewness, kurtosis), median calculation, and simple linear regression minimizing Sum of Absolute Deviations (SAD). It also includes a parser for evaluating probability expressions.

### [utilit - General Utilities for Graphics and Timing](utilit/README.md)

This module houses a collection of foundational Pascal units designed to provide essential functionalities for graphical applications, time measurement, and screen interaction within a Windows environment. It includes units for CRT window management, 2D graphing, screen save/load utilities, timing, and printer interfacing.

## General Usage Notes

*   **Programming Language**: All components of the `jpmMath` library are implemented in **Turbo Pascal** or **Borland Pascal for Windows**. (In an on-going activity, routines are being be ported to modern **Free Pascal** - see folder `fpc`).
*   **Data Types**: The library primarily utilizes `DOUBLE` (aliased as `REAL_AR` in many units) for floating-point numbers to ensure high precision in numerical computations. Integers (`INTEGER`, `LONGINT`) are used for degrees, counts, and indices.
*   **Input/Output**: Most programs interact with the user via the console for input and display results directly to the console. Some advanced modules, particularly those handling large datasets or providing graphical output, may read from structured text files (e.g., `.dat` files) and write comprehensive reports to output files (e.g., `.lst` or `.txt` files).
*   **Global Variables**: Many procedures and functions within the library's units frequently use global variables for inputs, outputs, and internal state. When integrating these routines into larger applications, users should be aware of this scope and manage variable assignments carefully to avoid unintended side effects.
*   **Numerical Accuracy and Error Handling**: The documentation for individual modules (and inline code comments) often discusses the numerical accuracy, stability characteristics, and potential error conditions for each algorithm. Return codes (e.g., `0` for success, non-zero for specific error types) are commonly used for robust error handling. Users are encouraged to consult these notes to understand appropriate use cases and limitations.
*   **Modularity**: The library is designed with modularity in mind. Many complex algorithms are broken down into smaller, reusable procedures and functions encapsulated within Pascal `UNIT`s. This allows developers to import and utilize specific functionalities without needing to understand the entire codebase.

## Contribution Guidelines

Contributions aimed at improving the clarity, robustness, or extending the functionality of this library are welcome. When contributing:

*   **Documentation Focus**: All new or modified code should be thoroughly documented, following the existing style of explanatory files (`.txt` or markdown) and inline comments in `.pas` files.
*   **Technical Language**: Maintain a technical and precise language, avoiding commercial jargon or promotional statements.
*   **Focus on Existing Features**: Ensure documentation solely reflects currently implemented features. Avoid documenting future or "to-be-developed" functionalities.
*   **Referenced Files/Functions**: Verify that all references to files, procedures, and functions within documentation accurately reflect their existence and naming in the codebase to prevent broken links or misleading information.
*   **Educational Value**: Consider how contributions can help future developers understand the mathematical and computational aspects of the methods, fostering a deeper understanding.
*   **Code Structure**: While some duplication of utility functions exists for demonstration purposes, consider proposing consolidation into shared units (e.g., Pascal `unit` files) for larger enhancements, if appropriate for Pascal.

## References

The theoretical foundations and implementations within the `jpmMath` library are drawn from established numerical analysis, scientific computing, and statistical texts.

*   **[BIBLI 01]**: Ruckdeschel, F.R. *BASIC Scientific Subroutines, Vol. II*. BYTE/McGRAW-HILL, 1981.
*   **[BIBLI 03]**: Ducamp, M., & Reverchon, A. (1991). *Analyse en Turbo Pascal versions 5.5 et 6.0*. Eyrolles, Paris.
*   **[BIBLI 04]**: Nowakowski, C. (1984). *Méthode de calcul numérique, Tome 2 - Programmes en Basic et en Pascal*. Edition du P.S.I.
*   **[BIBLI 05]**: Ducamp, M., & Reverchon, A. (1988). *Mathématiques en Turbo-Pascal (vol 2)*. Eyrolles, Paris.
*   **[BIBLI 08]**: Press, W.H., Flannery, B.P., Teukolsky, S.A., & Vetterling, W.T. (1986). *Numerical Recipes*. Cambridge University Press.
*   **[BIBLI 11]**: Engeln-Muellges, G., & Uhlig, F. (1996). *Numerical Algorithms with C*. Springer-Verlag.
*   **[BIBLI 12]**: Dony, R. (1990). *Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0*. MASSON.
*   **[BIBLI 13]**: Haut, H. (1981). *Mathematiques et statistiques*. PSI Editions, France.
*   **[BIBLI 16]**: Lalanne, M., Berthier, P., & Der Hagopian, J. (1980). *Mécanique des vibrations linéaires*. Masson, Paris.
*   **[BIBLI 18]**: Dang Trong, T. (n.d.). *Numath Library*. Fortran 77.

*Additional references cited within individual module documentations:*

*   Abramowitz, M., & Stegun, I.A. (1964). *Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables*. National Bureau of Standards, Applied Mathematics Series 55.
*   Adams, A.G. (1969). Areas under the Normal Curve, Algorithm 39. *Computer J.*, *12*, 197-198.
*   Amos, D.E. (1983, 1985, 1986). Sandia National Laboratories reports for Bessel functions.
*   Ball, L.A. (n.d.). *Algorithms for RPN calculators*. Wiley and sons.
*   Becket, R., & Hurt, J. (n.d.). *Numerical Calculaiions and Algorithms*.
*   Brent, R.P. (1973). *Algorithms for Minimization without Derivatives*. Prentice-Hall.
*   Brown, K.M. (1969). A quadratically convergent Newton-like method based upon Gaussian elimination. *SIAM J. Numer. Anal.*, *6*, 560.
*   Clenshaw, C.W. (n.d.). *Chebyshev Series for Mathematical Functions*. (Mentioned in `tbessi.pas`, `tbessj.pas`, `tbessk.pas`).
*   CRC Standard Mathematical Tables, 24th Edition.
*   Fike, C.T. (1968). *Computer Evaluation of Mathematical Functions*. Englewood Cliffs, NJ: Prentice-Hall.
*   Hart, J.F. et al. (1968). *Computer Approximations*. Wiley.
*   Hewlett-Packard statistics programs, 1974.
*   Hill, I.D., & Pike, M.C. (1967). Algorithm 299. *Collected Algorithms for the CACM*, p. 243.
*   Horner, W.G. (1819). A new method of solving numerical equations of all orders, by continuous approximation. *Philosophical Transactions of the Royal Society of London*, *109*, 308-335.
*   Jardrin, J.-L. (1988). *Algèbre, Algorithmes et programmes en Pascal*. DUNOD Paris.
*   Jin.ece.uiuc.edu/routines/routines.html (Fortran Routines for Computation of Special Functions).
*   Journal of Applied Statistics (1968) vol.17, p.189.
*   Journal of Applied Statistics (1972) vol.21, p.226.
*   Journal of Applied Statistics (1973) vol.22 no.3.
*   Journal of Applied Statistics (1978) vol.27 no.3.
*   Lafara, R.L. (n.d.). *Computer Methods for Science and Engineering*.
*   Lehning, H., & Jakubowicz, D. (1982). *Mathematiques par l'informatique individuelle - Programmes en basic*. MASSON Paris.
*   McMahon, J. (n.d.). *Tables of Bessel Functions*. (Referenced for zeros in `trootj.pas`).
*   More, J.J., Garbow, B.S., & Hillstrom, K.E. (1980). *User Guide for MINPACK-1*. Argonne National Laboratory.
*   Peirce, B.O. (1957). *A Short Table of Integrals*. Ginn and Company.
*   Shoup, T.E. (n.d.). *A Practical Guide to Computer Methods*.
*   Texas Instruments SR-51 Owners Manual, 1974.
*   Walker, J. (n.d.). CHI2 and Inverse CHI2 Functions. (Adapted from Gary Perlman). Online resource: `http://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html`.

