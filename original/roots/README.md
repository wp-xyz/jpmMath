
# jpmMath - roots - Comprehensive Collection of Numerical Methods

## Table of Contents

1.  [Single Variable Root Finding Methods](#single-variable-root-finding-methods)
    *   [Aitken Acceleration](#aitken-acceleration)
    *   [Aitken-Steffenson Iteration](#aitken-steffenson-iteration)
    *   [Bisection Method](#bisection-method)
    *   [Brent's Method (Van Wijngaarden-Dekker-Brent Method)](#brents-method-van-wijngaarden-dekker-brent-method)
    *   [Lagrange Method](#lagrange-method)
    *   [Modified False Position (Regula Falsi) Method](#modified-false-position-regula-falsi-method)
    *   [Mueller's Method (One-Dimensional)](#muellers-method-one-dimensional)
    *   [Newton's Method](#newtons-method)
    *   [NextRoot Subroutine](#nextroot-subroutine)
    *   [Pegasus Method](#pegasus-method)
    *   [Quartile Method](#quartile-method)
    *   [Secant Method](#secant-method)
    *   [Zeroin Method](#zeroin-method)
2.  [Polynomial Root Finding Methods and Operations](#polynomial-root-finding-methods-and-operations)
    *   [Allroot Subroutine (Complex Domain Mueller's)](#allroot-subroutine-complex-domain-muellers)
    *   [Bairstow's Method](#bairstows-method)
    *   [Bernouilli's Method](#bernouillis-method)
    *   [Lin's Method](#lins-method)
    *   [Roots of Algebraic Equations (Degree 2, 3, 4)](#roots-of-algebraic-equations-degree-2-3-4)
    *   [Polynomial Root Testing](#polynomial-root-testing)
    *   [Synthetic Division of Polynomials](#synthetic-division-of-polynomials)
3.  [Systems of Non-Linear Equations Solvers](#systems-of-non-linear-equations-solvers)
    *   [Brown's Method](#browns-method)
    *   [Levenberg-Marquardt (LM) Algorithm](#levenberg-marquardt-lm-algorithm)
    *   [Mueller's Method (Two-Dimensional)](#muellers-method-two-dimensional)
    *   [NNES (Nonlinear Equation Solver)](#nnes-nonlinear-equation-solver)
    *   [Nonlinear System of Two Variables](#nonlinear-system-of-two-variables)
4.  [Real-World Applications and Problem Solving](#real-world-applications-and-problem-solving)
5.  [Project Structure](#project-structure)
6.  [Utility Modules](#utility-modules)
    *   [Utils (`utils.pas`)](#utils-utilspas)
    *   [Function Definitions (`fonction.pas`, `fcn.pas`)](#function-definitions-fonctionpas-fcnpas)
7.  [Further Information](#further-information)

---

## Single Variable Root Finding Methods

These methods are designed to find a single real root `x` of a univariate function `f(x) = 0`.

### Aitken Acceleration

*   **Description**: Aitken acceleration is a technique to improve the convergence rate of an iterative sequence `x_n+1 = g(x_n)`. It employs three successive terms of the sequence (`x_n`, `x_n-1`, `x_n-2`) to parabolically predict an improved value `x'_n`, which then serves as the starting point for a new sequence of three iterations. This method effectively accounts for the curvature of the function to speed up convergence.
*   **Mathematical Principle**: The acceleration formula is given by:
    `x'_n = x_n - (x_n - x_{n-1})^2 / (x_n - 2x_{n-1} + x_{n-2})`.
    This technique is, in essence, based on accounting for the curvature of the function as the root is approached.
*   **Key Parameters**:
    *   `X0` (double): Initial guess.
    *   `C` (double): Convergence factor. A constant value (e.g., -1) is often used, but its optimal choice can affect convergence. If divergence occurs, smaller and/or positive values should be tried.
    *   `E` (double): Error criterion. The iteration returns when the iterative change in `x` is less than `E`.
    *   `M` (integer): Maximum number of iterations to prevent infinite loops.
*   **Outputs**:
    *   `X` (double): The calculated estimate of the position of the zero.
    *   `Y(X)` (double): The associated function value at `X`.
    *   `N` (integer): The number of iterations used to obtain that estimate.
*   **Characteristics**:
    *   **Strengths**: Can significantly accelerate the convergence of linearly convergent sequences (fixed-point iterations). It is known for its stable convergence properties, often requiring less accurate initial guesses compared to methods like Newton's or Secant. The number of iterations tends to be reasonably constant.
    *   **Limitations**:
        *   Requires three previous iteration steps before the acceleration can be performed.
        *   The choice of `C` (convergence factor) is empirical and can impact the convergence rate or lead to convergence on different roots if too large.
        *   The method is susceptible to round-off error if the denominator (`x_n - 2x_{n-1} + x_{n-2}`) approaches zero, potentially leading to a reduction in efficiency or a divide-by-zero overflow (though the implementation guards against overflow).
        *   The returned accuracy is not necessarily `E`; `E` serves as a rough error estimate.
        *   The formulation of `y(x)` can affect accuracy, potentially introducing "false roots" due to round-off errors.
*   **Real-world Applications**:
    *   **Iterative Solutions of Differential Equations**: Accelerating Picard iterations for solving initial value problems.
    *   **Numerical Simulations**: Improving the efficiency of simulations where quantities are iteratively calculated until convergence (e.g., in fluid dynamics or heat transfer).
    *   **Optimization**: Speeding up iterative algorithms for finding optimal parameters in engineering models.
*   **References**: F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `aitken.pas`, `aitken.txt`

### Aitken-Steffenson Iteration

*   **Description**: Aitken-Steffenson iteration is a direct variation of the Aitken acceleration scheme, specifically tailored for fixed-point iterations `x = g(x)`. It directly incorporates the acceleration into each iteration step, aiming for quadratic convergence without explicitly requiring the derivative of `g(x)`.
*   **Mathematical Principle**: The iteration formula is:
    `x'_n = x_n - [g(x_n) - x_n]^2 / (g[g(x_n)] - 2g(x_n) + x_n)`.
    The calculated `x'_n` is then used to directly compute the next iterate `x_{n+1} = g(x'_n)`. This differs from the original Aitken acceleration which requires three prior values for the acceleration calculation.
*   **Key Parameters**:
    *   `X0` (double): Initial guess.
    *   `C` (double): Convergence factor.
    *   `E` (double): Convergence criterion (accuracy).
    *   `M` (integer): Maximum number of iterations.
*   **Outputs**:
    *   `X` (double): The calculated zero.
    *   `Y(X)` (double): The associated function value at `X`.
    *   `N` (integer): The number of iterations used.
*   **Characteristics**:
    *   **Strengths**: Offers very good (quadratic) convergence properties if the initial guess is in the vicinity of the root. It is often more stable than some other algorithms, avoiding divergence even when very close to the root, and achieves quadratic convergence without requiring an analytical derivative.
    *   **Limitations**:
        *   Can be somewhat slower than the classic Aitken acceleration in terms of the total number of iterations required for a given precision.
        *   A significant portion of examples may fail to terminate, likely due to round-off error problems in the calculations, especially when the denominator approaches zero.
        *   Convergence can be slow if the chosen `C` is too small. Conversely, if `C` is too large or has the wrong sign, it can lead to "jumps" in the iteration, potentially converging to a different root.
        *   Poor convergence can occur if the derivative `|dg/dx|` is close to 1 near the root, particularly if `dy/dx` is close to zero.
*   **Real-world Applications**:
    *   **Fixed-Point Problems**: Solving implicit equations (e.g., equations of state in thermodynamics, or in fluid dynamics where solutions are found by repeated substitution).
    *   **Numerical Analysis Research**: Exploring advanced convergence acceleration techniques for various iterative schemes.
*   **References**: F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `steffen.pas`, `steffen.txt`

### Bisection Method

*   **Description**: The Bisection method is a robust and guaranteed-convergent root-finding algorithm. It works by repeatedly halving an interval `[a, b]` where a root is known to exist (i.e., `f(a)` and `f(b)` have opposite signs). The process continues until the interval size falls below a specified tolerance.
*   **Mathematical Principle**: Starting with an interval `[x0, x1]` where `Y(x0) * Y(x1) < 0`, the midpoint `x_mid = (x0 + x1) / 2` is calculated. The interval is then updated to `[x0, x_mid]` or `[x_mid, x1]` based on which sub-interval still contains a sign change of the function, thereby guaranteeing a root within the new, smaller interval.
*   **Key Parameters**:
    *   `X0`, `X1` (double): The initial range `(X0, X1)`. Must bracket the root.
    *   `E` (double): The convergence criterion (desired interval width).
*   **Outputs**:
    *   `X` (double): The calculated zero (midpoint of the final interval).
    *   `Y(X)` (double): The associated function value at `X`.
    *   `M` (integer): The number of steps (iterations) performed.
*   **Characteristics**:
    *   **Strengths**: Guaranteed convergence if the function is continuous and a root is bracketed within the initial interval. It is simple to implement and very robust, not sensitive to the quality of the initial guess (beyond bracketing).
    *   **Limitations**:
        *   Relatively slow convergence rate (linear), as the error is halved at each step.
        *   Requires an initial interval where the function changes sign.
        *   Cannot find multiple roots within the same interval without additional external logic.
        *   Does not work for roots where the function touches the x-axis without crossing (e.g., `x^2` at `x=0`).
*   **Real-world Applications**:
    *   **Initial Root Localization**: For finding the approximate location of a root when a reliable, albeit slow, first approximation is sufficient.
    *   **Financial Analysis**: Determining the break-even point for a project or the internal rate of return (IRR) where a cost function crosses a revenue function, and a guaranteed solution is needed.
    *   **Control Systems**: Finding a stable operating point for a system where monotonicity and bracketing can be ensured.
*   **Associated Files**: `bisect.pas`, `tquart.pas`

### Brent's Method (Van Wijngaarden-Dekker-Brent Method)

*   **Description**: Brent's method is a highly robust and efficient root-finding algorithm for one-dimensional real functions. It intelligently combines the guaranteed convergence of the Bisection method with the faster convergence rates of the Secant method and inverse quadratic interpolation. This hybrid approach ensures reliability while often achieving superlinear convergence.
*   **Mathematical Principle**: The algorithm maintains a bracketing interval `[a, b]` where a root is known to exist. It attempts to use inverse quadratic interpolation, which is faster, but if the interpolated point falls outside the bracket or if convergence is too slow, it falls back to a bisection step. This adaptive strategy ensures both guaranteed convergence and competitive speed.
*   **Key Parameters**:
    *   `x1`, `x2` (double): The interval `(X1, X2)` that must bracket the root (`f(x1) * f(x2) <= 0`).
    *   `Tolerance` (double): The desired accuracy for the root.
    *   `maxIterations` (integer): The maximum number of iterations to perform.
*   **Outputs**:
    *   `result` (double): The estimated root value.
    *   `valueAtRoot` (double): The function value at the estimated root.
    *   `niter` (integer): The number of iterations performed.
    *   `error` (integer): An error code (0: all OK, 1: no root found in interval, 2: max iterations reached).
*   **Characteristics**:
    *   **Strengths**: Widely regarded as the method of choice for general one-dimensional root finding when only function values are available (derivative-free). It is guaranteed to converge if a root is bracketed, and it combines the sureness of bisection with the speed of higher-order methods when appropriate.
    *   **Limitations**: Can still exhibit slower convergence for pathological functions or if the second derivative changes sharply near the root. Requires an initial interval where the function changes sign.
*   **Real-world Applications**:
    *   **Financial Modeling**: Calculating implied volatility from the Black-Scholes model, where an implicit non-linear equation needs to be solved reliably.
    *   **Chemical Reaction Engineering**: Determining equilibrium constants or conversion rates for complex reactions.
    *   **Optimization**: Finding roots of functions without explicit derivatives, as often arises in unconstrained optimization problems.
    *   **Scientific Computing**: A default choice in many numerical libraries due to its balance of reliability and efficiency for a wide range of root-finding problems.
*   **References**: R.P. Brent, "Algorithms for Minimization without Derivatives", Prentice-Hall, 1973. (Also references "BORLAND MATHEMATICAL LIBRARY").
*   **Associated Files**: `zbrent.pas`, `zbrent.txt`

### Lagrange Method

*   **Description**: The Lagrange method, as implemented here, is a variant of the False Position (Regula Falsi) method for finding a root of a non-linear function. It iteratively approximates the function as a straight line between two points that bracket the root.
*   **Mathematical Principle**: Given two points `(X1, Y1)` and `(X2, Y2)` such that `Y1` and `Y2` have opposite signs, the method calculates a new estimate `XB` where the secant line connecting `(X1, Y1)` and `(X2, Y2)` intersects the x-axis: `XB = X2 - (X2 - X1) * Y2 / (Y2 - Y1)`. The interval is then updated by replacing either `X1` or `X2` with `XB` to ensure the root remains bracketed.
*   **Key Parameters**:
    *   `X1`, `X2` (double): Starting and ending X values for the interval. Must bracket the root.
    *   `E` (double): Maximum error (convergence criterion).
    *   `MAXITER` (integer): Maximum number of iterations.
*   **Outputs**:
    *   `XB` (double): The estimated root.
    *   `K` (integer): The number of iterations performed.
    *   `FUNC(XB)` (double): The function value at the estimated root.
*   **Characteristics**:
    *   **Strengths**: Guaranteed to converge if a root is bracketed and the function is continuous. Does not require the function's derivative.
    *   **Limitations**: Can be slower than Newton's method or even Secant method, especially if one endpoint gets "stuck" (i.e., doesn't update frequently), leading to many iterations with a fixed point. Requires an initial interval where the function changes sign.
*   **Real-world Applications**:
    *   **Implicit Equation Solving**: For equations in engineering or physics where analytical inversion is difficult, and numerical derivatives are not preferred.
    *   **Economic Models**: Finding market equilibrium points or other parameters from implicit functions.
    *   **Empirical Modeling**: When approximating roots from experimental data points where the analytical form of the function might not be known or is complex.
*   **References**: Claude Nowakowski, "Méthode de calcul numérique - Programmes en Basic", PSI Editions, France, 1981.
*   **Associated Files**: `lagrange1.pas`

### Modified False Position (Regula Falsi) Method

*   **Description**: This subroutine calculates the zeroes of a function `Y(x)` using Hamming's modification to the classical False Position method. The False Position method itself is a bracketing method that uses a secant line to determine the next approximation. Hamming's modification addresses the issue where one endpoint remains fixed for many iterations, significantly slowing convergence.
*   **Mathematical Principle**: The method starts with two points `(x0, Y(x0))` and `(x1, Y(x1))` that bracket the root. A new estimate `x` is found using the secant line. If `Y(x)` has the same sign as `Y(x0)`, then `x0` is replaced by `x`, and `Y(x1)` is halved. If `Y(x)` has the same sign as `Y(x1)`, then `x1` is replaced by `x`, and `Y(x0)` is halved. This "halving" trick prevents one endpoint from becoming permanently "stuck", forcing the interval to shrink more symmetrically.
*   **Key Parameters**:
    *   `X0`, `X1` (double): Two initial guesses that must bracket the root.
    *   `E` (double): Convergence factor (tolerance). The iteration stops if `|X1 - X| < E`.
    *   `M` (integer): Maximum number of iterations.
*   **Outputs**:
    *   `X` (double): The calculated zero.
    *   `Y(X)` (double): The associated function value at `X`.
    *   `N` (integer): The actual number of iterations used.
*   **Characteristics**:
    *   **Strengths**: Guaranteed convergence because it always maintains a bracketing interval around the root. Hamming's modification significantly improves its convergence speed compared to the classical False Position method. Does not require derivative information.
    *   **Limitations**: Still requires an initial interval where the root is bracketed. While faster than Bisection, it generally does not achieve the quadratic convergence of Newton's method.
*   **Real-world Applications**:
    *   **Chemical Engineering**: Solving implicit equations for properties like activity coefficients or fugacity.
    *   **Environmental Modeling**: Finding the concentration of a pollutant given a non-linear decay model.
    *   **Engineering Design**: Finding critical dimensions or operating points where a design criterion is met, and a reliable bracketing method is needed.
*   **References**: F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `regula.pas`

### Mueller's Method (One-Dimensional)

*   **Description**: Mueller's method is a powerful root-finding algorithm that approximates the function with a parabola passing through three given points. It then uses the roots of this quadratic equation as the next approximation. This method can find both real and complex roots and is known for its robustness and fast convergence.
*   **Mathematical Principle**: An extension of the secant method, Mueller's method fits a parabola through three points: `(x1, Y(x1))`, `(x2, Y(x2))`, and `(x3, Y(x3))`. It determines the roots of this parabolic fit. The root closest to `x3` is chosen as the new `x3` for the next iteration, and `x1` is dropped, `x2` becomes `x1`, `x3` becomes `x2`. This iterative refinement continues until a convergence criterion is met.
*   **Key Parameters**:
    *   `X0` (double): Initial guess.
    *   `D` (double): A bound on the error of this initial guess, used to set the initial three evaluation points (`X0`, `X0-D`, `X0+D`).
    *   `E` (double): Convergence criterion. The iteration stops when `|X_new - X_old| < E`.
    *   `N` (integer): Maximum number of iterations.
*   **Outputs**:
    *   `X` (double): The calculated zero.
    *   `Y(X)` (double): The associated function value at `X`.
    *   `K` (integer): The number of iterations performed.
*   **Characteristics**:
    *   **Strengths**: Offers a very powerful convergence rate (superlinear, generally better than Secant and approaching quadratic) and good stability. It can naturally find complex roots of real-valued functions. It is more robust than Newton's method in some challenging cases, especially when the initial guess is far from a root or if the function exhibits high curvature.
    *   **Limitations**:
        *   The algebraic complexity is higher compared to linear approximation methods.
        *   Requires three initial evaluation points.
        *   Convergence can be slow if the initial guess is very far from any root, as the parabolic approximation may not accurately represent the function globally.
        *   The algorithm contains checks to handle cases where the parabolic fit does not cross the X-axis (leading to imaginary roots), substituting artificial intersections to keep the iteration going.
        *   The final accuracy is not strictly guaranteed to be less than the error criterion `E`; `E` serves as a rough estimate.
*   **Real-world Applications**:
    *   **Electrical Engineering**: Finding resonance frequencies or poles in circuit analysis, which might involve complex numbers.
    *   **Signal Processing**: Designing filters or analyzing system stability by finding roots of characteristic equations.
    *   **Chemical Kinetics**: Solving for reaction rate constants that involve non-linear functions with potentially complex roots.
    *   **Advanced Control Systems**: Analyzing system stability from characteristic equations whose roots may be complex.
*   **References**: R. Becket and J. Hurt, "Numerical Calculaiions and Algorithms". F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `mueller.pas`, `mueller.txt`

### Newton's Method

*   **Description**: Newton's method (also known as the Newton-Raphson method) is a powerful and typically fast root-finding algorithm that requires knowledge of both the function `f(x)` and its derivative `f'(x)`. It uses the tangent line at the current approximation to find the next, typically converging very rapidly if the initial guess is sufficiently close to the root.
*   **Mathematical Principle**: The method is based on the idea that if a point `x_n` is close to the root, the tangent line at `(x_n, Y(x_n))` will intersect the x-axis at a point `x_{n+1}` that is even closer to the root. The iterative formula is:
    `x_{n+1} = x_n - Y(x_n) / Y'(x_n)`.
    The `Y` function also returns its derivative `Y1`.
*   **Key Parameters**:
    *   `X0` (double): Initial guess.
    *   `E` (double): Convergence factor (tolerance). The iteration stops when `|Y(x0) / Y1|` is less than `E`.
    *   `M` (integer): Maximum number of iterations.
*   **Outputs**:
    *   `X0` (double): The calculated zero (updated in place).
    *   `YY` (double): The associated function value at `X0`.
    *   `N` (integer): The number of steps (iterations) performed.
*   **Characteristics**:
    *   **Strengths**: Exhibits quadratic convergence when sufficiently close to a simple root, meaning the number of accurate digits approximately doubles with each iteration. It is extremely fast when it converges.
    *   **Limitations**:
        *   Requires the analytical derivative `Y'(x)` to be explicitly known and computable. This can be a significant drawback for complex or implicitly defined functions.
        *   Highly sensitive to the quality of the initial guess. If the initial guess is poor, the method may diverge, oscillate, or converge to a different root.
        *   Problems arise if the derivative `Y'(x)` is zero or very small near the root or the current iteration point, potentially leading to division by zero or large jumps.
        *   It can be difficult to guarantee convergence globally.
*   **Real-world Applications**:
    *   **Optimization Problems**: Finding local minima or maxima by setting the gradient of an objective function to zero.
    *   **Numerical Solutions of Differential Equations**: Used in implicit numerical schemes (e.g., backward Euler method) for solving non-linear algebraic equations that arise at each time step.
    *   **Robotics**: Solving inverse kinematics problems for robot manipulators, where precise joint angles are required.
    *   **Computer Graphics**: Calculating intersections of rays with implicitly defined surfaces (e.g., spheres, tori).
    *   **Financial Engineering**: Valuing complex derivatives that involve implicit equations (e.g., interest rate swaps).
*   **References**: F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `newton1.pas`, `newton1.txt`

### NextRoot Subroutine

*   **Description**: The `NextRoot` subroutine is designed to find additional roots of a function `Y(x)` given a set of already established roots. It utilizes a modified Newton-Raphson iteration to effectively "deflate" the function by the known roots, allowing the search for other zeros.
*   **Mathematical Principle**: When some roots `A_i` of a function `Y(x)` are known, the function can be implicitly deflated by considering the problem of finding roots of `Y(x) / Product(x - A_i)`. The Newton-Raphson iteration is applied to this modified function. The derivative of `ln(Y(x) / Product(x - A_i))` is `Y'(x)/Y(x) - Sum(1/(x - A_i))`. This allows the algorithm to avoid converging to previously found roots.
*   **Key Parameters**:
    *   `L` (integer): The number of roots that have already been determined.
    *   `A` (array of double): An array containing the values of the `L` known roots.
    *   `X0` (double): The initial guess for the *next* root to be found.
    *   `E` (double): Accuracy criteria. The iteration stops when `|X_new - X_old| < E`.
    *   `M` (integer): Maximum number of iterations.
*   **Outputs**:
    *   `X` (double): The calculated new zero.
    *   `Y(X)` (double): The associated function value at `X`.
    *   `N` (integer): The number of steps (iterations) performed.
*   **Characteristics**:
    *   **Strengths**: Enables finding multiple roots of a single function systematically. It is particularly useful for functions with many roots, preventing Newton's method from repeatedly finding the same root. Leverages the fast convergence properties of Newton-Raphson.
    *   **Limitations**: The accuracy of subsequently found roots can be negatively affected by accumulated errors from imprecise previously found roots.
*   **Real-world Applications**:
    *   **Modal Analysis**: In mechanical or structural engineering, finding all natural frequencies (modes) of a vibrating system by solving its characteristic equation, where each frequency corresponds to a root.
    *   **Control Systems Design**: Identifying all poles and zeros of a complex transfer function, which can be critical for stability and performance analysis.
    *   **Chemical Reaction Engineering**: Finding multiple equilibrium points in complex reacting systems.
*   **Associated Files**: `nextroot.pas`

### Pegasus Method

*   **Description**: The Pegasus method is an efficient root-finding algorithm for continuous functions, known for its guaranteed convergence. It improves upon the False Position method by introducing a modification that ensures the bracketing interval always shrinks effectively, preventing cases where one endpoint remains fixed for many iterations.
*   **Mathematical Principle**: The method starts with an interval `[x1, x2]` where `fct(x1)` and `fct(x2)` have opposite signs. A new approximation `x3` is found using the secant line. If `fct(x2)` and `fct(x3)` have opposite signs, the interval `[x1, x2]` is updated to `[x2, x3]` (standard False Position). However, if `fct(x1)` and `fct(x3)` have opposite signs, the method updates `x2` to `x3` and `fct(x1)` is scaled down (multiplied by `fct(x2) / (fct(x2) + fct(x3))`) to ensure that the next secant line is "pulled" towards the root, thus guaranteeing faster convergence.
*   **Key Parameters**:
    *   `x1`, `x2` (double): Starting values for the interval. Must satisfy `fct(x1) * fct(x2) <= 0`.
*   **Outputs**:
    *   `x2` (double): The computed approximation for the root.
    *   `f2` (double): The function value at the approximate root (should be nearly zero).
    *   `iter` (integer): The number of iterations performed.
    *   **Return Code (integer)**:
        *   `-1`: No inclusion (initial `fct(x1) * fct(x2) > 0`).
        *   `0`: Root has been found (`ABS(f2) < FCTERR`).
        *   `1`: Break-off due to small step size (`ABS(x_new - x_old) < ABSERR + x_new * RELERR`).
        *   `2`: Iteration limit reached.
*   **Characteristics**:
    *   **Strengths**: Guaranteed convergence (unlike Secant or Newton's methods) because it continuously maintains a bracketing interval around the root. Often converges faster than the classical False Position method and Bisection.
    *   **Limitations**: Requires an initial interval where the function changes sign.
*   **Real-world Applications**:
    *   **Process Control**: In systems requiring highly reliable root finding for real-time adjustments, such as determining optimal setpoints in chemical reactors.
    *   **Numerical Optimization**: For finding optima of unimodal functions where a bracketing method is suitable, and a relatively fast convergence is desired.
    *   **Any application where guaranteed convergence with improved speed over Bisection is crucial**.
*   **References**: Gisela Engeln-Muellges and Frank Uhlig, "Numerical Algorithms with C", Springer-Verlag, 1996.
*   **Associated Files**: `fpegasus.pas`, `fonction.pas` (defines test functions), `tpegasus.pas` (test program)

### Quartile Method

*   **Description**: The Quartile method is an interval-reducing root-finding algorithm, conceptually similar to the Bisection method. Instead of simply halving the interval at each step, it divides the interval into quarters (or other proportional segments) and then selects the sub-interval that still brackets the root. This can sometimes lead to slightly faster convergence or improved accuracy compared to pure bisection.
*   **Mathematical Principle**: Given an interval `[a, b]` where `Y(a) * Y(b) < 0`, the method calculates two interior points, often at `a + 0.25*(b-a)` and `a + 0.75*(b-a)`. It then evaluates the function at these points and identifies the new sub-interval that still contains a root (i.e., where the function changes sign). The specific implementation here appears to choose between `a + co*(b-a)` and `a + (1-co)*(b-a)` based on `ABS(Y(a))` vs `ABS(Y(b))`, picking the point closer to the smaller magnitude function value.
*   **Key Parameters**:
    *   `a`, `b` (double): The initial range (interval) where the root is sought. Must bracket the root.
    *   `tol` (double): The tolerance for the interval width (`abs(a-b) < tol`).
*   **Outputs**:
    *   `root` (double): The found root (midpoint of the final interval).
    *   `step` (integer): The number of steps (iterations) used.
*   **Characteristics**:
    *   **Strengths**: Like the Bisection method, it guarantees convergence if a root is bracketed. It is generally slightly faster and can be marginally more accurate than the classical Bisection method in practice.
    *   **Limitations**: It still exhibits linear convergence, meaning it's slower than superlinear methods like Secant or Newton's. Requires an initial interval where the function changes sign.
*   **Real-world Applications**:
    *   **Robust Search**: As a refined alternative to Bisection for applications where bracketing and guaranteed convergence are paramount, and slight speed improvements over traditional bisection are beneficial without the risks of divergence from open methods.
    *   **Educational Tool**: Useful for demonstrating variations of interval-reducing methods and their comparative performance.
    *   **Simple Optimization**: Finding extrema of unimodal functions within a given range by applying the root-finding method to the derivative.
*   **Associated Files**: `tquart.pas` (demonstrates both Bisection and Quartile for comparison)

### Secant Method

*   **Description**: The Secant method is a widely used root-finding algorithm that approximates the derivative of the function using a secant line between two points. It is a faster alternative to the Bisection method but does not guarantee convergence.
*   **Mathematical Principle**: Starting with two initial guesses `x0` and `x1`, the method generates successive approximations `x` by finding the intersection of the x-axis with the secant line passing through `(x0, Y(x0))` and `(x1, Y(x1))`. The formula is:
    `x = (x0 * Y(x1) - x1 * Y(x0)) / (Y(x1) - Y(x0))`.
    For the next iteration, `x0` is replaced by `x1`, and `x1` is replaced by `x`.
*   **Key Parameters**:
    *   `X0`, `X1` (double): Two initial guesses. These do not necessarily need to bracket the root.
    *   `E` (double): Convergence factor (tolerance). The iteration stops if `|X1 - X0| < E`.
    *   `M` (integer): Maximum number of iterations.
*   **Outputs**:
    *   `X0` (double): The estimated root (updated in place).
    *   `Y(X0)` (double): The associated function value at `X0`.
    *   `N` (integer): The number of iterations performed.
*   **Characteristics**:
    *   **Strengths**: Does not require the computation of the analytical derivative, making it easier to apply than Newton's method for some functions. It converges faster than the Bisection or False Position methods, exhibiting superlinear convergence (order approximately 1.618).
    *   **Limitations**:
        *   Does not guarantee convergence. If the initial guesses are far from the root or if the function has problematic curvature, the method may diverge.
        *   Does not guarantee bracketing of the root, meaning the iteration can "jump" away from the actual root.
        *   Can suffer from division by zero or large numerical errors if `Y(x1)` is very close to `Y(x0)`.
*   **Real-world Applications**:
    *   **Chemical Engineering**: Solving implicit equations for fluid properties (e.g., vapor pressure or compressibility factor from complex equations of state).
    *   **Control Systems**: Finding roots of characteristic equations in control systems where analytical derivatives are difficult to obtain.
    *   **Optimization**: Approximating roots of the gradient of a function for unconstrained optimization problems.
    *   **Physics Simulations**: Solving transcendental equations where explicit analytical solutions are not feasible.
*   **References**: F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `secant.pas`

### Zeroin Method

*   **Description**: The Zeroin method is a highly robust and efficient root-finding algorithm that combines the guaranteed convergence of the Bisection method with the faster convergence rates of the Secant method and inverse quadratic interpolation. It adaptively chooses the best method at each step to ensure both reliability and speed.
*   **Mathematical Principle**: The algorithm operates on an interval `[a, b]` where `f(a) * f(b) < 0`. It intelligently decides whether to use a bisection step (guaranteed to reduce the interval by half), a secant step (using the two latest points), or an inverse quadratic interpolation step (using three latest points). The choice is based on criteria that ensure rapid convergence when possible while falling back to bisection if other methods perform poorly or would lead outside the bracketing interval, thus guaranteeing convergence.
*   **Key Parameters**:
    *   `abserr` (REAL): Absolute error bound.
    *   `relerr` (REAL): Relative error bound. These bounds ensure termination criteria.
    *   `fmax` (integer): Maximum number of calls for the function `Fkt(x)`.
    *   `protnam` (string): Name of the log file for intermediate results (can be an empty string for no logging).
    *   `a`, `b` (REAL): Endpoints of the interval. Must include a root (i.e., `Fkt(a) * Fkt(b) < 0`).
*   **Outputs**:
    *   `b` (REAL): The approximate root.
    *   `fb` (REAL): The function value at the root `b`.
    *   `fanz` (integer): The number of actual function calls made.
    *   `rc` (integer): Error code.
        *   `-2`: `abserr` or `relerr` negative, or both zero, or `fmax < 1`.
        *   `-1`: Initial `Fkt(a) * Fkt(b)` is not negative (no inclusion).
        *   `0`: Desired accuracy has been reached.
        *   `1`: `b` on output is a root where `Fkt(b) = 0`.
        *   `2`: Either `a` or `b` was a root at input.
        *   `3`: After `fmax` calls, no root found.
        *   `4`: Error opening intermediate result file.
*   **Characteristics**:
    *   **Strengths**: Considered one of the most reliable and efficient general-purpose root-finding algorithms. It combines the guaranteed convergence of bisection with the faster rates of secant and inverse quadratic interpolation, switching strategies dynamically for optimal performance.
    *   **Limitations**: Requires an initial interval where the function changes sign (i.e., a bracketed root).
*   **Real-world Applications**:
    *   **Numerical Libraries**: Often implemented as the standard root-finding routine in mathematical software due to its robustness.
    *   **Financial Engineering**: Calculating implied volatility or yield to maturity, which are implicit equations.
    *   **Physics and Engineering Simulations**: For any application requiring highly reliable root finding where functions can be expensive to evaluate, and derivatives are not available.
*   **References**: Gisela Engeln-Muellges and Frank Uhlig, "Numerical Algorithms with C", Springer-Verlag, 1996.
*   **Associated Files**: `fzeroin.pas`, `tzeroin.pas` (test program)

## Polynomial Root Finding Methods and Operations

These methods are specialized for finding roots of polynomial equations `P(x) = 0` or performing operations on polynomials.

### Allroot Subroutine (Complex Domain Mueller's)

*   **Description**: This routine attempts to calculate multiple roots, including complex roots, of a given function or polynomial. It achieves this by repeatedly using the `ZMueller` subroutine (an implementation of Mueller's method adapted for the complex domain) and then "removing" (deflating) the roots already found by division.
*   **Mathematical Principle**: The core idea is to find one root using `ZMueller`, then divide the original function/polynomial by the factor `(z - root)` (for a single root) or `(z^2 + pz + q)` (for a quadratic factor representing a complex conjugate pair). This process, known as polynomial deflation, reduces the degree of the polynomial, allowing the method to search for the next root in the remaining polynomial.
*   **Key Parameters**:
    *   `X0, Y0` (double): Initial guesses for the root in the complex plane (`X0` real part, `Y0` imaginary part).
    *   `b1, b2` (double): Bounds on the error in these initial guesses for X and Y, respectively.
    *   `e` (double): Convergence criterion.
    *   `n` (integer): Maximum number of iterations for each root.
    *   `n2` (integer): The total number of roots to be sought.
    *   `n3` (integer): A flag: `1` for a generic function `F(Z)` (defined in `Eval`), `2` for a polynomial series (defined in `Series`).
*   **Outputs**:
    *   `X(i), Y(i)` (arrays of double): The `n2` estimated roots, where `X(i)` is the real part and `Y(i)` is the imaginary part of the i-th root.
    *   `k` (integer): The last number of iterations performed for finding one root.
*   **Characteristics**:
    *   **Strengths**: Capable of finding multiple roots, including complex roots, of functions and polynomials effectively. It leverages the robustness of Mueller's method and systematically deflates the polynomial to find subsequent roots.
    *   **Limitations**: The overall success and accuracy depend on the robustness of the underlying `ZMueller` method and the numerical stability of the deflation process, especially when roots are very close or when errors accumulate during repeated deflation.
*   **Real-world Applications**:
    *   **Electrical Engineering**: Analyzing circuits with complex impedances, where poles and zeros of transfer functions can be complex.
    *   **Quantum Mechanics**: Solving for energy eigenvalues in systems where the characteristic equations can have complex roots.
    *   **Control Systems**: Comprehensive analysis of system stability and response by identifying all poles and zeros of complex transfer functions.
    *   **Applied Physics**: Finding complex roots of dispersion relations in wave propagation phenomena.
*   **Associated Files**: `allroot.pas`

### Bairstow's Method

*   **Description**: Bairstow's method is a powerful iterative algorithm for finding both real and complex roots of a polynomial `f(x)`. It works by iteratively extracting quadratic factors from the polynomial, incorporating a Newton-Raphson step to accelerate convergence.
*   **Mathematical Principle**: The method assumes a polynomial `f(x)` can be divided by a quadratic factor `x^2 + px + q`, leaving a remainder `rx + s`. Bairstow's method uses synthetic division to perform this division and obtain `r` and `s` as functions of `p` and `q`. It then applies a 2x2 Newton's method to find `p` and `q` such that `r(p,q) = 0` and `s(p,q) = 0`. Once such `p` and `q` are found, the roots of `x^2 + px + q = 0` are calculated (which can be real or a complex conjugate pair). The process is then repeated on the deflated polynomial (the quotient) until all roots are found or the degree is reduced to 1 or 2.
*   **Key Parameters**:
    *   `A` (array of double): The coefficients of the polynomial.
    *   `M` (integer): The order (degree) of the polynomial.
    *   `e` (double): Convergence factor (tolerance for the changes in `p` and `q`).
    *   `n` (integer): Maximum number of iterations for finding each quadratic factor.
    *   `aa`, `bb` (double): Initial guesses for `p` and `q`.
*   **Outputs**:
    *   `x1, y1`, `x2, y2` (double): The real and imaginary parts of the two roots found for the extracted quadratic factor. (If `y1=y2=0`, they are real roots).
    *   `k` (integer): The number of iterations used to find the last factor.
    *   When finding *all* roots (`bairsto1.pas`): A list of all roots and an approximation of their error.
*   **Characteristics**:
    *   **Strengths**:
        *   Can find both real and complex roots directly without the need for complex arithmetic (by working with real coefficients `p` and `q`).
        *   Exhibits quadratic convergence when sufficiently close to the roots, making it significantly faster than Lin's method.
        *   Robust enough to handle polynomials with multiple roots.
    *   **Limitations**:
        *   Can be sensitive to the initial guess for `p` and `q`, particularly for higher-degree polynomials or if roots are very close.
        *   Precision of subsequently found roots can be affected by accumulated round-off errors from repeatedly deflating the polynomial.
*   **Real-world Applications**:
    *   **Control Systems**: Finding poles and zeros of transfer functions in continuous and discrete-time systems, which are critical for stability analysis and controller design.
    *   **Electrical Engineering**: Analyzing RLC circuits, filter design, and stability of feedback systems, where characteristic equations are polynomials.
    *   **Mechanical Engineering**: Vibration analysis of multi-degree-of-freedom systems, where natural frequencies and damping ratios are derived from polynomial roots.
    *   **Chemical Kinetics**: Solving complex systems of non-linear algebraic equations derived from reaction mechanisms.
*   **References**:
    *   H. Lehning and D. Jakubowicz, "Mathematiques par l'informatique individuelle - Programmes en basic", MASSON Paris, 1982.
    *   R.L. Lafara, "Computer Methods for Science and Engineering".
    *   R. Becket and J. Hurt, "Numerical Calculaiions and Algorithms".
*   **Associated Files**: `bairstow.pas` (two complex roots), `bairsto1.pas` (all roots), `bairstow.txt`

### Bernouilli's Method

*   **Description**: Bernouilli's method is an iterative algorithm primarily used to find the root with the largest absolute value (the dominant root) of a polynomial equation. Once found, the polynomial can be deflated (divided by a linear factor corresponding to the root), and the process can be repeated to find other roots in descending order of magnitude.
*   **Mathematical Principle**: For a polynomial `P_n(x) = a_0 x^n + a_1 x^{n-1} + ... + a_n`, it considers a linear recurrence relation `a_0 y_k + a_1 y_{k-1} + ... + a_n y_{k-n} = 0`. The ratio of successive terms `y_{k+1} / y_k` in a sequence generated by this recurrence relation (with specific starting values for `y`) converges to the root of `P_n(x)` with the largest absolute value, provided such a unique dominant root exists.
*   **Key Parameters**:
    *   `nd` (integer): The order (degree) of the polynomial.
    *   `A` (array of double): The coefficients of the polynomial.
*   **Outputs**:
    *   `X` (double): The calculated root. The program iteratively finds and prints each root.
*   **Characteristics**:
    *   **Strengths**: Conceptually straightforward for identifying the dominant root of a polynomial. Can be extended to find all roots through repeated deflation.
    *   **Limitations**:
        *   Only guaranteed to find the largest root in magnitude.
        *   Convergence can be slow if the dominant root is not well-separated in magnitude from other roots.
        *   Difficulties arise with multiple roots of the same largest magnitude or with complex conjugate roots of the largest magnitude (which might not be unique).
        *   Requires `A[0]` (coefficient of `x^n`) to be non-zero.
        *   Accuracy can degrade for subsequent roots due to accumulated errors from polynomial deflation.
*   **Real-world Applications**:
    *   **Structural Analysis**: Finding the dominant natural frequency of a vibrating system, which corresponds to the largest eigenvalue (and thus the largest root of the characteristic polynomial).
    *   **Population Dynamics**: Determining the long-term growth rate of a population modeled by a linear recurrence relation.
    *   **Time Series Analysis**: Identifying dominant trends or components in time series data.
    *   **Control Systems**: Preliminary analysis to find dominant poles that affect system stability and response.
*   **References**: Claude Nowakowski, "Methodes de calcul numerique", PSI Editions, France, 1981.
*   **Associated Files**: `bernou.pas`, `bernou.txt`

### Lin's Method

*   **Description**: Lin's method is an iterative procedure for extracting quadratic factors from a polynomial. It capitalizes on the fact that complex roots of real-valued polynomials always appear in conjugate pairs, meaning a quadratic factor with real coefficients can be extracted. This allows for finding complex roots using only real-number operations.
*   **Mathematical Principle**: The method assumes a quadratic factor of the form `z^2 + Az + B`. It iteratively performs synthetic division of the polynomial `P(z)` by this factor to obtain a remainder. The goal is to determine `A` and `B` such that the remainder coefficients are zero. The current implementation uses an update rule for `A` and `B` based on the remainder terms.
*   **Key Parameters**:
    *   `A` (array of double): Polynomial coefficients.
    *   `M` (integer): Order (degree) of the polynomial. Should be `>=4` for meaningful quadratic factor extraction.
    *   `aa`, `bb` (double): Initial guesses for `A` and `B` coefficients.
    *   `e` (double): Convergence factor (tolerance). The iteration stops if `ABS(aa - a1) + ABS(bb - b1) < e^2`.
    *   `n` (integer): Maximum number of iterations.
*   **Outputs**:
    *   `x1, y1`, `x2, y2` (double): The real and imaginary parts of the two complex conjugate roots (or real roots if `y1=y2=0`) found for the extracted quadratic factor.
    *   `k` (integer): The number of iterations performed.
*   **Characteristics**:
    *   **Strengths**: Uses only real-number operations for finding complex conjugate root pairs, simplifying implementation. Applicable to polynomials of degree four or higher.
    *   **Limitations**:
        *   Often suffers from slow convergence.
        *   Convergence is not always guaranteed, especially when multiple roots are present or if the initial guesses are poor.
        *   Small relative errors in the final estimates of `A` and `B` can result in much larger relative errors in the finally derived roots.
        *   Sensitive to the choice of initial values for `A` and `B`.
        *   May encounter divide-by-zero issues if the quotient coefficients become zero.
*   **Real-world Applications**:
    *   **Preliminary Factorization**: Used for initial factorization of high-degree polynomials in control systems design or network analysis, providing preliminary factors that can be refined by faster methods like Bairstow's.
    *   **Historical/Educational Context**: Understanding the development of polynomial root-finding algorithms, as Bairstow's method is a significant improvement building upon Lin's.
*   **References**: T.E. Shoup, "A Practical Guide to Computer Methods". F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `lin.pas`, `lin.txt`

### Roots of Algebraic Equations (Degree 2, 3, 4)

*   **Description**: This program provides direct and specialized methods for calculating real or complex roots of algebraic equations (polynomials) of degree two, three, and four. It leverages analytical formulas where available and uses a bisection method for the cubic case.
*   **Mathematical Principle**:
    *   **Degree 2 (Quadratic)**: `x^2 + bx + c = 0`. Solved using the quadratic formula `x = [-b +/- sqrt(b^2 - 4c)] / 2`. The nature of the roots (real or complex conjugate) depends on the sign of the discriminant `b^2 - 4c`.
    *   **Degree 3 (Cubic)**: `x^3 + a2 x^2 + a1 x + a0 = 0`. First, one real root (`xc`) is found. If `a0 = 0`, `xc = 0` is a root. Otherwise, a bisection method is used to find a real root within a defined interval (e.g., `[-(1+lambda), 0]` if `a0 > 0`, `[0, (1+lambda)]` if `a0 < 0`, where `lambda` is the greatest absolute value of `a_i` coefficients). Once `xc` is found, the cubic equation is divided by `(x - xc)` to yield a quadratic equation, which is then solved using the degree 2 method.
    *   **Degree 4 (Quartic)**: `x^4 + a x^3 + b x^2 + c x + d = 0`. A substitution `x = y - a/4` transforms the equation into a "depressed" quartic. This can then be factored into a product of two quadratic equations: `(y^2 + 2ky + l)(y^2 - 2ky + m) = 0`. The parameter `k^2` is a real root of an auxiliary cubic equation. Solving this cubic equation (using the degree 3 method) provides `k`, which then allows for `l` and `m` to be determined, leading to the solution of the two quadratic factors.
*   **Key Parameters**:
    *   `n` (integer): The degree of the equation (must be 2, 3, or 4).
    *   `a` (matrix): A matrix storing the coefficients of the polynomial (e.g., `a[n,i]` for `x^i`).
*   **Outputs**:
    *   `r` (array): An array storing the real parts of the roots.
    *   `im`, `ii` (double): Variables storing the imaginary parts of the roots where applicable (for complex conjugate pairs).
*   **Characteristics**:
    *   **Strengths**: Provides direct (analytical) solutions for quadratic, cubic, and quartic equations, avoiding the iterative convergence issues of general methods for these lower-degree cases once the auxiliary equations are solved. Highly efficient and accurate for their respective degrees.
    *   **Limitations**: Strictly applicable only for polynomials of degree 2, 3, or 4. The bisection method used for the auxiliary cubic equation still requires bracketing.
*   **Real-world Applications**:
    *   **Engineering Design**: Rapidly solving common polynomial equations arising in mechanical engineering (e.g., beam deflection, gear ratios), electrical engineering (e.g., resonant circuits), or civil engineering (e.g., structural analysis).
    *   **Physics**: Solving problems in kinematics, dynamics, and optics that often reduce to quadratic or cubic equations.
    *   **Geometry**: Finding intersection points of curves and surfaces, or solving problems involving volumes and areas that lead to low-degree polynomials.
*   **References**: "Mathématiques et statistiques - Programmes en BASIC", Editions du P.S.I., 1981.
*   **Associated Files**: `root4.pas`, `root4.txt`

### Polynomial Root Testing

*   **Description**: This program assists in the pre-analysis of polynomial root-finding by providing insights into the characteristics and potential locations of roots without necessarily finding their exact values. It can estimate the maximum magnitude of the largest root and provide information about the number of positive and negative real roots using Descartes' Rule of Signs.
*   **Mathematical Principle**:
    *   **Bounds on Roots**: It estimates an upper bound for the magnitude of the largest root based on the coefficients of the polynomial.
    *   **Descartes' Rule of Signs**: This rule states that the number of positive real roots of a polynomial `P(x)` is either equal to the number of sign changes between consecutive non-zero coefficients, or is less than it by an even number. For negative real roots, the rule is applied to `P(-x)`.
*   **Key Parameters**:
    *   `N` (integer): The degree of the polynomial.
    *   `A` (array of double): The coefficients of the polynomial (`A(0)` to `A(N)`).
*   **Outputs**:
    *   The total number of roots (`N`).
    *   An upper bound on the magnitude of the largest root.
    *   Statements indicating the presence of at least one negative real root (if applicable).
    *   Statements regarding the maximum number of positive and negative real roots (based on Descartes' rule).
    *   Inferences about the existence of complex roots.
*   **Characteristics**:
    *   **Strengths**: Serves as a valuable preliminary analysis tool. It provides quick, qualitative heuristics that can guide the choice of specific iterative root-finding algorithms and help in selecting appropriate initial guesses or bracketing intervals, thereby improving the efficiency and reliability of subsequent computations.
    *   **Limitations**: Provides only qualitative or bounding information, not exact root locations. Descartes' Rule of Signs gives an upper bound, not an exact count, and only applies to real roots.
*   **Real-world Applications**:
    *   **Algorithm Selection**: Helping to decide which root-finding method (e.g., real-only vs. complex-capable) is most appropriate for a given polynomial.
    *   **Initial Guess Generation**: Guiding the process of selecting starting points for iterative solvers by indicating regions where roots are likely to be found.
    *   **Robustness Enhancement**: Used as a preliminary step in a multi-stage root-finding process to improve overall reliability by pre-conditioning the search.
*   **References**: F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `roottest.pas`

### Synthetic Division of Polynomials

*   **Description**: This subroutine performs synthetic division of polynomials, specifically calculating `A(x) = C(x) / B(x)`, where `C(x)` is the dividend polynomial and `B(x)` is the divisor polynomial. Synthetic division is a simplified and efficient method for dividing a polynomial, particularly useful when the divisor is a linear factor `(x - r)` or a quadratic factor.
*   **Mathematical Principle**: The procedure calculates the coefficients of the quotient polynomial `A(x)` (of maximum order `n1 - n2`) by iteratively applying the synthetic division algorithm. It essentially "divides" the coefficients of `C(x)` by the coefficients of `B(x)` in a streamlined manner.
*   **Key Parameters**:
    *   `n1` (integer): The degree of the first polynomial `C(x)`.
    *   `n2` (integer): The degree of the second polynomial `B(x)`.
    *   `C` (array of double): Coefficients of the first polynomial (dividend).
    *   `B` (array of double): Coefficients of the second polynomial (divisor).
*   **Outputs**:
    *   `A` (array of double): The coefficients of the resulting quotient polynomial.
*   **Characteristics**:
    *   **Strengths**: Highly efficient for polynomial division, especially when `n2` is small (linear or quadratic factors). It is a fundamental operation for polynomial deflation in many root-finding algorithms (e.g., Bairstow's method, Bernouilli's method, NextRoot subroutine).
    *   **Limitations**: The primary function is division; it does not directly find roots. Assumes real coefficients. The implementation typically calculates the quotient and implicitly discards the remainder (or assumes exact division).
*   **Real-world Applications**:
    *   **Polynomial Deflation**: After a root `r` is found for a polynomial `P(x)`, synthetic division by `(x - r)` yields a new polynomial `Q(x)` of lower degree, from which remaining roots can be found. This is a core step in many polynomial root-finding algorithms.
    *   **Root Testing**: In some contexts, synthetic division can be used to quickly evaluate a polynomial at a given point (the remainder will be `P(r)`).
    *   **Control Systems**: Used in algorithms for polynomial factorization, which is relevant for filter design and system identification.
*   **References**: F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `rsyndiv.pas`

## Systems of Non-Linear Equations Solvers

These methods are designed to find the roots of a system of `n` non-linear equations in `n` variables, `F(x) = 0`.

### Brown's Method

*   **Description**: Brown's method is a quadratically convergent Newton-like method for solving systems of non-linear equations `F(x) = 0`. It is based on Gaussian elimination and offers robust convergence properties by iteratively linearizing the system.
*   **Mathematical Principle**: For a system of `n` non-linear equations, Brown's method proceeds by sequentially linearizing each equation. It involves a modified form of Gaussian elimination to solve the system for each iteration, updating the solution vector `x`. This approach avoids explicitly forming and inverting the full Jacobian matrix at each step, although it still requires partial derivative information (computed here via finite differences).
*   **Key Parameters**:
    *   `n` (integer): The number of equations and variables in the system.
    *   `x0` (vector): The starting value for the iteration (initial guess vector).
    *   `eps` (REAL): The error bound (desired accuracy for convergence).
    *   `prot` (integer): Protocol switch (0: no output, 1: detailed output of intermediate results).
    *   `maxit` (integer): Maximal number of steps (iterations).
*   **Outputs**:
    *   `x1` (vector): The solution vector.
    *   `itanz` (integer): The actual number of steps (iterations) performed.
    *   `rc` (integer): Error code.
        *   `0`: All OK (successful iteration).
        *   `1`: Desired accuracy not achieved after `maxit` iterations.
        *   `2`: System matrix (linearized Jacobian) singular.
        *   `3`: Lack of memory.
        *   `4`: Wrong input parameters (`n < 1` or `maxit < 1`).
        *   `5`: Error in evaluating the function.
*   **Characteristics**:
    *   **Strengths**: Exhibits quadratic convergence when sufficiently close to the solution, offering fast convergence rates. It is a Newton-like method that can be more numerically stable than a direct Newton's method for some problems, and it avoids explicit matrix inversion.
    *   **Limitations**: Can be sensitive to the quality of the initial guess `x0`. May encounter issues if the linearized system (Jacobian) becomes singular or ill-conditioned. Computational cost can still be significant for very large systems.
*   **Real-world Applications**:
    *   **Chemical Process Design**: Solving complex chemical reaction networks, material and energy balances in multi-stage reactors.
    *   **Economic Modeling**: Finding equilibrium points in general equilibrium models with many interacting non-linear relationships.
    *   **Electrical Engineering**: Solving power flow equations in electrical grids or analyzing complex non-linear circuits.
    *   **Robotics**: Inverse kinematics for complex robotic manipulators with multiple degrees of freedom.
*   **References**:
    *   K.M. Brown, "A quadratically convergent Newton-like method based upon Gaussian elimination", SIAM J. Numer. Anal. Vol. 6, (1969), p. 560.
    *   G. Engeln-Mueller and F. Uhlig, "Numerical Algorithms with C", Springer-Verlag, 1996.
*   **Associated Files**: `brown.pas` (unit), `brownts1.pas` (test program)

### Levenberg-Marquardt (LM) Algorithm

*   **Description**: The Levenberg-Marquardt algorithm is a widely used optimization algorithm for solving non-linear least squares problems. It is particularly effective for data fitting or when finding the "best fit" solution to an overdetermined system of non-linear equations (more equations than variables). It smoothly interpolates between the Gauss-Newton algorithm and the method of gradient descent.
*   **Mathematical Principle**: The algorithm minimizes the sum of squares of `m` non-linear functions `F(x)` by iteratively updating the parameter vector `x`. It solves a modified normal equation at each step: `(J^T J + λI)Δx = J^T F(x)`, where `J` is the Jacobian matrix, `λ` is the damping parameter, `I` is the identity matrix, and `Δx` is the update step. When `λ` is small, it approximates Gauss-Newton (fast, but can diverge if far from minimum). When `λ` is large, it approximates gradient descent (slow, but robust). The algorithm adaptively adjusts `λ` to ensure convergence, typically decreasing `λ` on successful steps and increasing it on unsuccessful ones. The Jacobian can be provided analytically or approximated by finite differences.
*   **Key Parameters**:
    *   `m` (integer): Number of functions (equations).
    *   `n` (integer): Number of variables ($n \le m$).
    *   `x` (pVec): On input, contains an initial estimate of the solution vector. On output, contains the final estimate.
    *   `fvec` (pVec): Output array containing the functions evaluated at the output `x`.
    *   `tol` (double): Nonnegative input variable. Termination occurs when the algorithm estimates either that the relative error in the sum of squares is at most `tol` or that the relative error between `x` and the solution is at most `tol`.
*   **Outputs**:
    *   `info` (integer): Integer output variable indicating termination reason.
        *   `0`: Improper input parameters.
        *   `1`: Both actual and predicted relative reductions in the sum of squares are at most `tol`.
        *   `2`: Relative error between two consecutive iterates is at most `tol`.
        *   `3`: Conditions for info = 1 and info = 2 both hold.
        *   `4`: `fvec` is orthogonal to the columns of the Jacobian to machine precision.
        *   `5`: Number of calls to function evaluations has reached or exceeded `maxfev`.
        *   `6`: `tol` is too small; no further reduction in sum of squares is possible.
        *   `7`: `tol` is too small; no further improvement in `x` is possible.
    *   `nfev` (integer): Number of function evaluations.
    *   `njev` (integer): Number of Jacobian evaluations (if applicable).
*   **Characteristics**:
    *   **Strengths**: Very robust and widely used for non-linear least squares problems. Its hybrid nature ensures convergence from a wider range of initial guesses compared to pure Gauss-Newton, while retaining good convergence speed near the minimum. Effectively handles cases where the Jacobian is ill-conditioned or rank-deficient.
    *   **Limitations**: Primarily designed for problems that can be formulated as minimizing a sum of squares. Can be computationally intensive for very large systems due to the need to compute or approximate the Jacobian matrix.
*   **Real-world Applications**:
    *   **Data Fitting/Nonlinear Regression**: Fitting experimental data to complex non-linear models in all scientific and engineering fields (e.g., kinetic models in chemistry, dose-response curves in biology, material constitutive laws in engineering).
    *   **Parameter Estimation**: Calibrating physical models by adjusting parameters to best match observed data.
    *   **Computer Vision**: Image registration, camera calibration, and 3D reconstruction from 2D images.
    *   **Robotics**: Robot calibration and inverse kinematics solutions where the objective is to minimize the error between desired and actual end-effector positions.
*   **References**: Argonne National Laboratory, MINPACK Project, March 1980 (Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More).
*   **Associated Files**: `lm.pas` (unit), `tlm.pas` (test program), `fcn.pas` (defines test functions)

### Mueller's Method (Two-Dimensional)

*   **Description**: This is an extension of the one-dimensional Mueller's method adapted to find roots of a two-dimensional function `W(x, y) = 0`. It combines the parabolic fitting technique with successive substitution to iteratively refine the `x` and `y` coordinates of the root.
*   **Mathematical Principle**: The algorithm proceeds iteratively. In each pass, it first applies the 1D Mueller's method to find a better estimate for the `X` coordinate (treating `Y` as constant). Then, it uses this newly found `X` to apply the 1D Mueller's method to find a better estimate for the `Y` coordinate (treating `X` as constant). This sequential refinement process continues until the sum of the absolute changes in `X` and `Y` falls below a chosen error criterion.
*   **Key Parameters**:
    *   `X0`, `Y0` (double): Initial guesses for the root location `(X0, Y0)`.
    *   `b1`, `b2` (double): Bounds on the error in the initial guess for the X and Y directions, respectively. These are used to set the initial three evaluation points for the 1D Mueller's parabolic fit in each dimension.
    *   `e` (double): Convergence criterion. The iteration stops when `|delta_X| + |delta_Y| < E`.
    *   `n` (integer): Maximum number of iterations.
*   **Outputs**:
    *   `X`, `Y` (double): The estimated root location `(X, Y)`.
    *   `W(X,Y)` (double): The associated function value at the estimated root.
    *   `K` (integer): The number of iterations performed.
*   **Characteristics**:
    *   **Strengths**: Extends a powerful 1D root-finding technique to higher dimensions, making it suitable for systems of two non-linear equations. It can effectively handle functions exhibiting high curvature due to its parabolic approximation.
    *   **Limitations**:
        *   Limited to solving systems of two variables.
        *   Convergence can be slow, especially with highly curved functions or if the initial guess is far from the true root.
        *   The method is vulnerable to divergence if the function has no unique root (e.g., `x^2 - y^2 + 1 = 0`) or if the iteration starts near saddle points or local minima. The provided implementation does not explicitly guard against non-unique roots.
*   **Real-world Applications**:
    *   **Engineering Design**: Solving coupled non-linear equations in mechanical design (e.g., finding the intersection of two complex curves or surfaces), or in fluid dynamics (e.g., two coupled flow equations).
    *   **Chemical Engineering**: Solving for two unknown concentrations or temperatures in coupled material and energy balance equations.
    *   **Physics**: Finding equilibrium points in simple multi-particle systems or in systems described by two coupled non-linear equations.
    *   **Geometry**: Calculating the intersection points of two implicitly defined curves in a plane.
*   **References**: R. Becket and J. Hurt, "Numerical Calculaiions and Algorithms". F.R. Ruckdeschel, "BASIC Scientific Subroutines, Vol. II", BYTE/McGRAW-HILL, 1981.
*   **Associated Files**: `mueller2.pas`, `mueller2.txt`

### NNES (Nonlinear Equation Solver)

*   **Description**: NNES (Nonmonotonic Nonlinear Equation Solver) is a comprehensive and highly configurable module for solving sets of non-linear equations. It implements advanced techniques such as nonmonotonic line searches, trust region methods (single and double dogleg), and various Jacobian approximation strategies including Broyden updates and Lee and Lee updates. It includes extensive options for adaptive scaling, convergence criteria, and detailed output diagnostics.
*   **Mathematical Principle**: NNES is built upon robust optimization theory. It typically seeks to minimize an objective function, usually the sum of squares of the residuals (`F(x)^T F(x)`). It offers two main globalization strategies:
    *   **Line Search Methods**: Determine a search direction (e.g., Newton step) and then find an optimal step length along that direction to reduce the objective function.
    *   **Trust Region Methods**: Define a region around the current iterate where a quadratic model of the function is deemed trustworthy, and then find a step within this region that minimizes the model. This includes "dogleg" (single and double) strategies.
    It supports:
    *   **Jacobian Evaluation**: Analytical (user-supplied) or finite differences (forward, backward, central).
    *   **Quasi-Newton Updates**: Broyden's method and Lee and Lee method to iteratively update the Jacobian approximation, reducing the need for expensive re-evaluations.
    *   **Adaptive Scaling**: Dynamically adjusts scaling factors for variables and functions to improve the conditioning of the problem and numerical stability.
    *   **Nonmonotonic Step Acceptance**: Allows for steps that temporarily increase the objective function value, which can help in escaping local minima and improving global convergence.
    *   **Cholesky Decomposition**: Used for solving linear systems and perturbing the Jacobian for ill-conditioned problems.
*   **Key Parameters**: NNES is highly configurable. Key parameters are often set via the `setup` procedure in the test program and are detailed in `nneshelp.txt`.
    *   **Logical Flags**: `absnew`, `cauchy`, `deuflh`, `geoms`, `linesr`, `newton`, `overch`, `matsup`, `wrnsup` (control specific algorithm choices like absolute Newton, Cauchy point for trust region, Deuflhard initialization, geometric search, line search vs. trust region, pure Newton, overflow checking, matrix printing suppression, warning suppression).
    *   **Integer Parameters**: `acptcr`, `itsclf`, `itsclx`, `jactyp`, `jupdm`, `maxexp`, `maxit`, `maxns`, `maxqns`, `mgll`, `minqns`, `narmij`, `niejev`, `njacch`, `output`, `qnupdm`, `stopcr`, `supprs`, `trupdm` (control acceptance criteria, adaptive scaling start iterations, Jacobian type, Jacobian update method, maximum exponent, maximum iterations, maximum line search/trust region steps, nonmonotonic history length, minimum quasi-Newton steps, Armijo steps, initial non-quasi-Newton steps, Jacobian check iterations, output detail, QR update method, stopping criteria, output suppression, trust region update method).
    *   **Real Parameters**: `alpha`, `confac`, `delta`, `delfac`, `epsmch`, `etafac`, `fdtolj`, `ftol`, `lam0`, `mstpf`, `nsttol`, `omega`, `ratiof`, `sigma`, `stptol` (control Armijo constant, constraint factor, initial trust region size, trust region update factors, machine precision, dogleg shape, Jacobian finite difference tolerance, function tolerance, initial relaxation factor, maximum step size factor, Newton step tolerance, Lee & Lee omega factor, ratio test for quasi-Newton updates, geometric reduction factor, step tolerance).
    *   **Vector Parameters**: `boundl`, `boundu`, `scalef`, `scalex` (arrays for lower/upper bounds, function scaling factors, variable scaling factors).
*   **Outputs**:
    *   `xplus` (pVec): The final solution vector.
    *   `fcnnew` (double): The final objective function value (sum of squares of residuals).
    *   `itnum` (integer): Total number of iterations.
    *   `nfunc` (integer): Total number of line search/trust region function evaluations.
    *   `njetot` (integer): Total number of explicit Jacobian evaluations.
    *   `nfetot` (integer): Total number of function evaluations (includes Jacobian evaluations).
    *   `trmcod` (integer): Termination code indicating the reason for stopping.
*   **Characteristics**:
    *   **Strengths**: Extremely versatile, robust, and efficient for a wide range of non-linear systems, including ill-conditioned problems. It combines advanced optimization strategies and offers fine-grained control over the solution process, making it suitable for complex research and industrial applications. Provides extensive diagnostic outputs to monitor convergence and identify issues.
    *   **Limitations**: Its high number of parameters and intertwining algorithms make it complex to configure and understand without thorough study. Requires careful tuning for optimal performance on specific problems.
*   **Real-world Applications**:
    *   **Chemical Process Simulation**: Solving large systems of coupled material and energy balance equations in complex chemical reactors, distillation columns, or entire process flowsheets.
    *   **Structural Analysis**: Analyzing large-scale non-linear structural behavior (e.g., plasticity, buckling, large deformations).
    *   **Computational Fluid Dynamics (CFD)**: Solving systems arising from discretization of non-linear partial differential equations (PDEs).
    *   **Optimization**: A powerful backend solver for general non-linear optimization problems that can be reformulated as finding roots of gradients.
    *   **Advanced Control Systems**: Designing and analyzing complex non-linear control systems.
*   **References**: R.S. Bain (Copyright 1991). (Also mentions David M. Gay for F90 release).
*   **Associated Files**: `tnnes.pas` (driver program), `unnes.pas` (Part 1/3), `unnes1.pas` (Part 2/3), `unnes2.pas` (Part 3/3), `fcn.pas` (test function definitions), `nneshelp.txt` (detailed parameter explanations), `utils.pas` (utility functions).

### Nonlinear System of Two Variables

*   **Description**: This program is specifically designed to solve a system of two non-linear equations `f(x,y) = 0` and `g(x,y) = 0`. It employs a numerical approximation of Newton's method for two variables, using finite differences to estimate the partial derivatives.
*   **Mathematical Principle**: At each iteration, the program calculates the partial derivatives of `f` and `g` with respect to `x` and `y` numerically using a small step size `h` (e.g., `a = df/dx`, `b = dg/dx`, `c = df/dy`, `d = dg/dy`). These form a 2x2 Jacobian matrix. The linear system `J * delta_vec = -F_vec` (where `delta_vec = [dx, dy]^T` and `F_vec = [f(x,y), g(x,y)]^T`) is then solved using Cramer's rule to find the update steps `p` (`dx`) and `q` (`dy`). The solution `(x,y)` is updated `x = x - p`, `y = y - q` until convergence.
*   **Key Parameters**:
    *   `x0`, `y0` (REAL_AR): The starting point (initial guesses).
    *   `prec` (REAL_AR): The desired precision (convergence criterion for `ABS(p) + ABS(q)`).
    *   `maxiter` (integer): Maximal number of iterations.
*   **Outputs**:
    *   `x`, `y` (REAL_AR): The approximate solution `(x,y)`.
    *   `iter` (integer): The number of iterations performed.
    *   `error` (integer): Error code.
        *   `0`: OK (successful convergence).
        *   `1`: Error in evaluating `f(x,y)` or `g(x,y)` (e.g., `Ln(x)` for `x <= 0`).
        *   `2`: Singular system (Jacobian determinant `t` is too small).
*   **Characteristics**:
    *   **Strengths**: Provides a general approach for solving 2x2 non-linear systems. It does not require the user to provide analytical partial derivatives, as they are computed numerically. Conceptually simpler than general N-dimensional solvers.
    *   **Limitations**: Strictly limited to systems of two variables. The choice of the finite difference step size `h` can impact accuracy and stability. May be sensitive to initial guesses, and can suffer from slow convergence or failure if the Jacobian becomes singular or ill-conditioned.
*   **Real-world Applications**:
    *   **Interdisciplinary Modeling**: Any model that simplifies to two coupled non-linear equations, such as predator-prey models in biology, or simple economic models with two interacting variables.
    *   **Electrical Circuit Analysis**: Solving for node voltages or currents in a two-node non-linear circuit.
    *   **Geometry**: Finding the intersection points of two implicitly defined curves in 2D.
    *   **Basic Engineering Problems**: Where two parameters are related by non-linear equations (e.g., finding the intersection of two characteristic curves of a pump).
*   **References**: M. Ducamp and A. Reverchon, "Mathematiques en Turbo-Pascal part 1", Editions EYROLLES, Paris, 1991.
*   **Associated Files**: `nlinsyst.pas`

## Real-World Applications and Problem Solving

The numerical methods provided in this library are foundational tools with widespread applicability across various scientific and engineering disciplines. Understanding how to frame a problem as a root-finding task and selecting the appropriate algorithm is key to effective problem-solving.

Here are examples of how these algorithms can be used to solve real-world problems, along with considerations for their application:

*   **Engineering Design:**
    *   **Mechanical Engineering**:
        *   **Problem**: Finding equilibrium points of complex mechanical linkages or optimizing spring-damper systems.
        *   **Application**: Formulate the force balance equations as a system of non-linear equations `F(x)=0`. **NNES** or **Brown's Method** could be used.
        *   **Problem**: Calculating the natural frequencies of vibrating structures.
        *   **Application**: These often involve finding roots of characteristic polynomials. **Bairstow's Method** (for all roots) or **Bernouilli's Method** (for dominant frequencies) would be suitable.
    *   **Electrical Engineering**:
        *   **Problem**: Determining operating points of non-linear circuits (e.g., transistor biasing).
        *   **Application**: Set up Kirchhoff's laws as a system of non-linear equations. **Brown's Method** or **NNES** could be used. For simple 2-variable circuits, `nlinsyst.pas` is an option.
        *   **Problem**: Designing filters and control systems where frequency responses are non-linear.
        *   **Application**: Finding roots of transfer functions or characteristic equations. **Bairstow's Method** for complex roots, or **Mueller's Method (1D)** for complex roots of a single function.
    *   **Chemical Engineering**:
        *   **Problem**: Solving reaction-rate equations, calculating phase equilibria, or optimizing reactor design parameters.
        *   **Application**: These often lead to complex non-linear systems or implicit equations. **NNES** (for complex systems), **Levenberg-Marquardt** (for parameter fitting), or **Modified False Position/Lagrange** (for implicit single-variable equations) are applicable.
    *   **Civil Engineering**:
        *   **Problem**: Analyzing the stability of structures under various loads, optimizing bridge designs, or modeling fluid flow in pipes and channels.
        *   **Application**: Stability analysis can involve roots of polynomials (**Root4**, **Bairstow's**). Fluid flow equations (e.g., Colebrook equation) are often implicit and can be solved by **Regula Falsi** or **Brent's Method**.
*   **Physics and Applied Mathematics:**
    *   **Problem**: Solving transcendental equations for energy levels in quantum mechanics or orbital parameters in astrophysics.
    *   **Application**: These are typically single-variable non-linear equations. **Brent's Method**, **Zeroin Method**, or **Newton's Method** (if derivative available) are excellent choices.
    *   **Problem**: Solving non-linear Navier-Stokes equations or finding solutions to boundary layer problems in fluid dynamics.
    *   **Application**: Often involves large, complex systems of non-linear equations. **NNES** is designed for such high-dimensional problems.
    *   **Optimization**: Many optimization problems can be transformed into root-finding problems by setting the gradient of the objective function to zero.
    *   **Application**: If the gradient is a single function, any single-variable method applies. If it's a vector, system solvers like **NNES** or **Brown's Method** are used.
*   **Finance and Economics:**
    *   **Problem**: Valuing complex financial derivatives using non-linear Black-Scholes equations or similar models.
    *   **Application**: Often requires solving implicit equations (e.g., for implied volatility). **Brent's Method** or **Zeroin Method** are reliable.
    *   **Problem**: Finding market equilibrium points where supply and demand functions (often non-linear) intersect.
    *   **Application**: If a single market, **Bisection** or **Regula Falsi**. For multi-market equilibrium, **NNES** or **Brown's Method**.
*   **Computer Graphics and Robotics:**
    *   **Problem**: Determining the joint angles of a robot arm (inverse kinematics).
    *   **Application**: Formulate as a system of non-linear equations. **Brown's Method** or **NNES** can solve this. For 2-DOF, `mueller2.pas` or `nlinsyst.pas` are specific options.
    *   **Problem**: Calculating intersection points of rays with complex geometric surfaces defined by implicit functions (ray tracing).
    *   **Application**: This is often a single-variable root-finding problem along the ray. **Mueller's Method (1D)** or **Brent's Method** could be used.

### General Problem-Solving Workflow with `jpmMath` Library:

1.  **Formulate the Problem**: Clearly define the function(s) `f(x) = 0` or system `F(x) = 0` whose roots you need to find.
2.  **Analyze Function Properties**:
    *   Is `f(x)` continuous? Differentiable?
    *   Are there multiple roots expected? Are they real or complex?
    *   Is the function well-behaved (e.g., no sharp turns, flat regions, or asymptotes near roots)?
    *   Is the system well-conditioned (e.g., small changes in inputs lead to small changes in outputs)?
3.  **Initial Guesses**:
    *   **Bracketing Methods (Bisection, Regula Falsi, Brent, Zeroin, Pegasus, Lagrange)**: Requires an initial interval `[a,b]` where `f(a)` and `f(b)` have opposite signs. Ensure the root is truly bracketed.
    *   **Open Methods (Newton, Secant, Mueller)**: A closer initial guess leads to significantly faster and more reliable convergence. Use insights from `roottest.pas` if dealing with polynomials.
4.  **Convergence Criteria**: Choose appropriate values for tolerance (`E`, `tol`, `eps`, `prec`) and maximum iterations (`M`, `N`, `maxit`, `maxiter`). A very small tolerance might lead to excessive computation or convergence issues due to machine precision limitations.
5.  **Method Selection**:
    *   **Robustness First**: For initial exploration or when reliability is paramount, consider **Brent's Method** or **Zeroin Method** for single variables. For polynomials, **Bairstow's Method** is robust. For systems, **Brown's Method** or **NNES** are good starting points.
    *   **Speed When Possible**: If analytical derivatives are available and the function is well-behaved, **Newton's Method** is often the fastest. **Secant** and **Mueller** are good derivative-free alternatives. **Aitken Acceleration** can boost fixed-point iterations.
    *   **Specialized Needs**:
        *   For finding the largest root of a polynomial: **Bernouilli's Method**.
        *   For low-degree polynomials (2, 3, 4): `root4.pas` offers direct solutions.
        *   For multi-variable problems that can be reduced to 2D: `mueller2.pas` or `nlinsyst.pas`.
        *   For non-linear least squares fitting: **Levenberg-Marquardt**.
    *   **Complex or Multiple Roots**: **Mueller's Method (1D)**, **Allroot**, or **Bairstow's Method** are appropriate. Use `nextroot.pas` for finding additional roots after some are known.
6.  **Error Handling and Diagnostics**: Pay close attention to output warnings (e.g., "divide by zero", "process not convergent", "Jacobian singular", "overflow"). These messages are crucial for identifying numerical instability, ill-conditioned problems, or an inappropriate choice of method or initial guess. The **NNES** module, in particular, offers extensive diagnostic outputs that can help pinpoint issues.
7.  **Scaling (for Systems)**: For systems of equations, especially with **NNES**, adaptively scaling variables and functions is crucial to improve numerical stability and convergence. This is particularly important for problems where variables or function values have widely varying magnitudes.

By carefully considering these aspects, users can effectively leverage this library to tackle a wide array of challenging mathematical problems in their respective fields.

## Project Structure

The `jpmMath` library is organized into several Free Pascal files, typically with a main program (`.pas`) demonstrating a specific method and accompanying units (`.pas`) or explanatory text files (`.txt`).

*   `_info.txt`: Provides a high-level overview and descriptions of the programs and their corresponding `.txt` explanation files. This serves as an index to the library's contents.
*   `*.txt` files: These files contain detailed explanations of the numerical methods. They delve into the theoretical background, mathematical principles, advantages, limitations, and often include example results and practical considerations for each algorithm. They are crucial for a deep understanding of the algorithms implemented.
*   `*.pas` files: These contain the Free Pascal source code for the numerical methods.
    *   **Main Program Files**: Often start with `Demo_`, `Test_`, or simply the method name (e.g., `aitken.pas`, `zbrent.pas`, `test_unnes.pas`). These files typically include the `uses WinCrt` unit for console I/O, which is common in legacy Pascal environments. They serve as executable examples demonstrating the usage of a particular algorithm.
    *   **Unit Files**: (e.g., `brown.pas`, `fcn.pas`, `fpegasus.pas`, `fzeroin.pas`, `lm.pas`, `unnes.pas`, `unnes1.pas`, `unnes2.pas`, `utils.pas`, `fonction.pas`, `basis.pas`, `basis_r.pas`, `type_def.pas`, `wincrt1.pas`, `strings.pas`) These provide modular implementations of methods, functions, or utility routines. They are designed to be reusable across different programs and solvers, promoting code organization and maintainability. For instance, the NNES solver is split across `unnes.pas`, `unnes1.pas`, and `unnes2.pas` for modularity.

This modular structure allows for both easy demonstration of individual algorithms and flexible integration of complex solvers into larger applications.

## Utility Modules

These modules provide supporting functions and data structures used across various numerical methods in the library.

### Utils (`utils.pas`)

*   **Description**: This unit contains a collection of general-purpose utility functions that are frequently used by other modules in the JPM Math Library. These utilities aim to reduce code duplication and improve consistency across the codebase.
*   **Key Functions**:
    *   `Dot_Product(n, a, b)`: Calculates the dot product of two vectors `a` and `b` of size `n`.
    *   `IMAX(a, b)`, `IMIN(a, b)`: Returns the maximum or minimum of two integers.
    *   `Line0`, `Line1`, `Msg(text)`: Functions for formatted output to the console or file, often used for debugging, logging, and structuring output reports (e.g., adding borders, messages).
    *   `log(x)`: Calculates the logarithm base 10 of `x`. Includes a check for `x > 1E-12` to prevent `Ln(0)` errors.
    *   `MAX(a, b)`, `MIN(a, b)`: Returns the maximum or minimum of two real numbers.
    *   `mmulv(N, A, B, C)`: Performs matrix-vector multiplication: `C = A * B`, where `A` is an `N x N` matrix, and `B`, `C` are `N x 1` vectors.
    *   `mmulvt(N, A, B, C)`: Performs multiplication of the transpose of a matrix by a vector: `C = A^T * B`, where `A` is an `N x N` matrix, and `B`, `C` are `N x 1` vectors.
    *   `Power(x, n)`: Calculates `x` raised to the integer power `n`.
    *   `Power1(y, x)`: Calculates `y` raised to the real power `x` (`y^x`). Includes a check for `x < 0`.
    *   `Sign(a, b)`: Returns `a` with the sign of `b`.
    *   `PrintVec(n, Name, V)`: (For debug only) Prints a vector `V` with a given `Name`.
*   **Purpose**: To provide essential mathematical and programming utilities, promoting modularity and reusability within the library.

### Function Definitions (`fonction.pas`, `fcn.pas`)

These units provide standardized test functions and systems of equations for various numerical methods in the library, facilitating testing, validation, and demonstration.

*   **`fonction.pas`**:
    *   **Description**: This unit defines six real functions (`fct`) that serve as test cases for the Pegasus method (`fpegasus.pas`). These functions cover a range of behaviors (e.g., exponential, polynomial, trigonometric) for which roots can be sought.
    *   **Purpose**: To provide a convenient and standardized set of single-variable non-linear equations for testing and comparing root-finding algorithms. Each function can be selected by setting a `NumFunc` variable.
*   **`fcn.pas`**:
    *   **Description**: This unit defines example non-linear systems of equations (`fcn1`) of various sizes (`n=2`, `n=4`, `n=6`). These systems are specifically designed as test cases for multi-variable non-linear equation solvers like NNES (`unnes.pas`, etc.) and Levenberg-Marquardt (`lm.pas`). It also includes a placeholder for an analytical Jacobian (`jacob`) that users can implement if `JACTYP = 0` is chosen in NNES.
    *   **Purpose**: To provide standardized, ready-to-use multi-variable non-linear problems for testing and demonstrating the capabilities of systems solvers, ensuring consistent evaluation of their performance.
