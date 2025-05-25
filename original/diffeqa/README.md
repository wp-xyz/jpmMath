
# jpmMath/diffeqa - Differential Equations Library

This directory contains a collection of Pascal units and programs for solving various types of differential equations, including Ordinary Differential Equations (ODEs) as Initial Value Problems (IVPs) and Boundary Value Problems (BVPs), as well as Partial Differential Equations (PDEs). The implementations cover a range of numerical methods, from classical approaches like Runge-Kutta and Adams-Bashforth to more advanced techniques like Gear's method for stiff systems and the Bulirsch-Stoer extrapolation method.

The library emphasizes numerical accuracy and efficiency, often including explicit error control mechanisms and adaptive step-size adjustments.

## Implemented Numerical Methods

The `diffeqa` library provides implementations for a diverse set of numerical methods. Below is an overview of the key functionalities:

### Initial Value Problems (IVPs)

Initial Value Problems are typically defined by a differential equation and initial conditions at a single point. The library offers various methods to integrate these problems over an interval.

*   **Euler-Romberg Method**:
    *   **Description**: An iterative method that repeatedly applies Euler's method with halved step sizes and uses linear extrapolation to improve approximation accuracy. It aims to maintain consistent precision throughout the computation.
    *   **Files**: `eulromb.pas`, `eulromb.txt` (explanation file).
    *   **Reference**: Based on "Méthode de calcul numérique- Tome 2 - Programmes en Basic et en Pascal By Claude Nowakowski, Edition du P.S.I., 1984" [BIBLI 04].

*   **Adams-Bashforth Method**:
    *   **Description**: An explicit multi-step method for numerical integration of ordinary differential equations. It uses function evaluations from previous steps to compute the current step, improving speed compared to single-step methods like Runge-Kutta for a given accuracy.
    *   **Files**: `adambash.pas`, `adambash.txt` (explanation file).
    *   **Reference**: Based on "Méthode de calcul numérique- Tome 2 - Programmes en Basic et en Pascal By Claude Nowakowski, Edition du P.S.I., 1984" [BIBLI 04].

*   **Adams-Moulton Prediction-Correction Method**:
    *   **Description**: An implicit multi-step method that uses a predictor formula (Adams-Bashforth) to estimate the solution and then a corrector formula (Adams-Moulton) to refine it. This predictor-corrector approach offers higher accuracy than explicit methods.
    *   **Files**: `adammoul.pas`, `ab_mou.pas` (unit for Adams-Bashforth-Moulton with automatic step size control), `adambash.txt` (explanation file).
    *   **Reference**: `adammoul.pas` from "Méthode de calcul numérique- Tome 2 - Programmes en Basic et en Pascal By Claude Nowakowski, Edition du P.S.I., 1984" [BIBLI 04]. `ab_mou.pas` based on "Numerical Algorithms with C, By Gisela Engeln-Muellges and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].

*   **Runge-Kutta Methods**:
    *   **Description**: A family of widely used single-step methods for numerical integration. The library includes implementations for various orders and applications.
        *   **Order 4**: A popular choice offering a good balance between speed and accuracy (error is typically O(h^5) per step). It is used for single first-order equations, higher-order equations (by converting to a system of first-order equations), and systems of first-order equations.
        *   **Automatic Step Size Control**: Advanced Runge-Kutta implementations feature adaptive step-size control to maintain a desired precision and optimize computational effort.
        *   **Runge-Kutta-Fehlberg (RKF45)**: An embedded method that computes two approximations of different orders (4th and 5th) to estimate the local truncation error and adapt the step size accordingly.
    *   **Files**:
        *   Core units: `eqdif1.pas`, `equdif.pas`, `urkf45.pas`.
        *   Explanations: `rkutta.txt`.
        *   Test programs: `teqdif1.pas`, `teqdif1a.pas` (RK order 1 with graph), `teqdifn.pas`, `teqdifn1.pas` (RK order N with graph), `teqdifp.pas`, `teqdifp1.pas` (RK for P variables with graph), `tequdif.pas` (RK with time step control), `trk4.pas` (RK for 2-variable systems), `trkf45.pas` (RKF45 test).
    *   **Reference**: `eqdif1.pas` and `rkutta.txt` from "Analyse en Turbo Pascal versions 5.5 et 6.0" By Marc DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991 [BIBLI 03]. `equdif.pas` by Jean-Pierre Dumont. `urkf45.pas` based on H A Watts and L F Shampine, Sandia Laboratories.

*   **Bulirsch-Stoer-Gragg Extrapolation Method**:
    *   **Description**: A highly accurate and efficient method for non-stiff IVPs. It uses Richardson extrapolation on a sequence of numerical approximations (often from the midpoint method) to achieve very high orders of accuracy. It features robust step-size control.
    *   **Files**: `bulirsch.pas`, `tbulirsc.pas` (test program).
    *   **Reference**: Based on "Numerical Algorithms with C, By Gisela Engeln-Muellges and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].

*   **Stiff Systems Solvers**:
    *   **Description**: Stiff differential equations are numerically challenging because explicit methods require very small step sizes for stability, even if the solution is slowly varying. Implicit methods are generally used for such problems.
        *   **Implicit Gear Method (Order 4)**: A multi-step implicit method well-suited for stiff systems. It involves solving a system of linear equations at each step, often utilizing Jacobian matrices.
        *   **Rosenbrock Method (Order 3 or 4)**: A single-step implicit method that is often more efficient for stiff problems than Gear's method when the Jacobian is easily available or can be approximated well.
    *   **Files**:
        *   Gear: `gear.pas`, `mgear.pas` (test program).
        *   Rosenbrock: `tros4.pas`.
    *   **Reference**: Gear based on "Numerical Algorithms with C, By Gisela Engeln-Muellges and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11]. Rosenbrock from Numath Library By Tuan Dang Trong.

*   **Automatic Step Size Control (AWP)**:
    *   **Description**: A generic procedure (`awp`) that integrates a system of first-order ODEs with automatic step-size control using various embedded Runge-Kutta-type formulas. It supports different methods (Runge-Kutta 2nd/3rd order, England 4th/5th order, Prince-Dormand 4th/5th order) and includes stiffness detection.
    *   **Files**: `uawp.pas`, `tawp.pas` (test program).
    *   **Reference**: Based on "Numerical Algorithms with C, By Gisela Engeln-Muellges and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].

### Boundary Value Problems (BVPs)

Boundary Value Problems involve differential equations where conditions are specified at multiple points (boundaries) of the independent variable. The shooting method is a common approach for solving such problems by transforming them into IVPs.

*   **Shooting Method**:
    *   **Description**: This method converts a BVP into an IVP by guessing the unknown initial conditions, integrating the IVP, and then adjusting the guessed initial conditions until the boundary conditions at the other end are satisfied. This often involves solving a non-linear system (e.g., using Newton's method).
        *   **First Order Systems**: Solves systems of first-order DEs by determining an approximation for the initial values, relying on various IVP solvers (Runge-Kutta, Adams-Bashforth-Moulton, Bulirsch-Stoer).
        *   **Second Order DE**: Addresses second-order DEs by converting them to a system of first-order equations and using the shooting method.
    *   **Files**: `urwp.pas` (unit for first order systems), `m_rwp.pas` (driver program for first order systems), `limits.pas` (for second order DEs using RK), `tshoot.pas` (test for second order DEs).
    *   **Reference**: `urwp.pas` based on "Numerical Algorithms with C, By Gisela Engeln-Muellges and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11]. `limits.pas` and `tshoot.pas` based on "Méthode de calcul numérique- Tome 2 - Programmes en Basic et en Pascal By Claude Nowakowski, Edition du P.S.I., 1984" [BIBLI 04].

### Partial Differential Equations (PDEs)

The library includes examples for solving Laplace's equation, a fundamental elliptic PDE, using the relaxation method.

*   **Laplace Equation by Relaxation Method**:
    *   **Description**: Solves the Laplace equation (∇²T = 0) on a 2D domain using an iterative finite-difference relaxation technique. It approximates the solution at each grid point as the average of its neighbors, iteratively updating until convergence. Different boundary conditions (Dirichlet and Neumann) are demonstrated.
    *   **Files**: `laplace.pas` (Example 1: square plate, Dirichlet), `laplace1.pas` (Example 2: rectangular plate with a hole, Dirichlet), `laplace2.pas` (Example 3: square plate with mixed Dirichlet/Neumann conditions).
    *   **Reference**: Based on "Méthode de calcul numérique- Tome 2 - Programmes en Basic et en Pascal By Claude Nowakowski, Edition du P.S.I., 1984" [BIBLI 04].

## How to Use the Library (General Approach)

The numerical methods in this library are typically implemented as Pascal `UNIT`s, which provide reusable procedures and functions, and `PROGRAM`s, which serve as executable examples demonstrating the usage of these units.

To utilize these methods:

1.  **Define the Differential Equation**: For IVPs and BVPs, you need to implement a Pascal function or procedure that evaluates the right-hand side `F(x, y)` (or `F(x, y1, y2, ..., yn)` for systems) of your differential equation(s). This function will be passed as a parameter to the solver routine. Refer to the `dgl` procedure in `t_dgls.pas` for examples of how differential equations are defined for system solvers.
2.  **Set Initial/Boundary Conditions**: Provide the necessary initial values (e.g., `y(x0) = y0`) for IVPs or boundary conditions (e.g., `y(x0) = a`, `y(x1) = b`) for BVPs.
3.  **Specify Integration Parameters**: Define parameters such as the integration interval (`x0`, `xend`), desired step size (`h`), and error tolerances (`epsabs`, `epsrel`). For iterative methods, maximum iterations might also be required.
4.  **Call the Solver Routine**: Pass your defined function and parameters.
5.  **Process Results**: The solver will typically return the solution values in an array or vector, often alongside statistics like the number of function evaluations or the final step size. Test programs usually print these results to the console or a file.

For PDEs solved by relaxation methods, you typically define the domain, grid size, and boundary conditions directly within the program, and the iteration proceeds to fill the grid with the solution.

## Real-World Problem Solving with `diffeqa`

This library provides a foundational set of tools for modeling and solving various real-world phenomena described by differential equations.

*   **Physics and Engineering**:
    *   **Classical Mechanics**: Model the motion of objects under forces (e.g., projectile motion, pendulum dynamics) using second-order ODEs. `Equadiffn` or `limits.pas` could be applied.
    *   **Electrical Circuits**: Analyze RLC circuits described by systems of first-order ODEs. `Equadiffp` or `uawp.pas` could be used. Stiff circuits might require `gear.pas` or `tros4.pas`.
    *   **Heat Transfer and Diffusion**: Solve steady-state heat distribution problems in materials using Laplace's equation (`laplace.pas`, `laplace1.pas`, `laplace2.pas`). Transient problems would involve other PDE types (not directly in this section but similar numerical principles apply).
    *   **Fluid Dynamics**: Model simplified fluid flows or pressure distributions.
*   **Biology and Chemistry**:
    *   **Population Dynamics**: Simulate population growth models (e.g., logistic growth) or predator-prey systems using coupled first-order ODEs. `Equadiffp` or `uawp.pas` are suitable.
    *   **Chemical Kinetics**: Describe reaction rates and concentrations of chemical species over time.
*   **Finance**:
    *   **Option Pricing**: The Black-Scholes equation, while typically a PDE, can sometimes be simplified or approximated using ODE techniques for specific boundary conditions.
*   **General Modeling**:
    *   **System Dynamics**: Model the behavior of complex systems where the rate of change of one variable depends on other variables within the system. Adaptive step-size control methods (`uawp.pas`, `ab_mou.pas`, `bulirsch.pas`) are particularly useful here for efficiency and accuracy.
    *   **Control Systems**: Analyze the response of control systems.

### Solving Strategies and Considerations:

*   **Choice of Method**:
    *   For **non-stiff problems** requiring high accuracy, `bulirsch.pas` is often a good choice.
    *   For general-purpose IVP solving with automatic step-size control, `uawp.pas` provides flexible options.
    *   For problems known to be **stiff**, `gear.pas` or `tros4.pas` are recommended.
    *   For simpler problems or for teaching purposes, `eulromb.pas` or basic `eqdif1.pas` methods are effective.
*   **Error Control**: Pay close attention to `epsabs` (absolute error tolerance) and `epsrel` (relative error tolerance). A good strategy is to set `epsrel` based on the desired number of significant digits and `epsabs` to a small value relevant to the scale of `y` if `y` might pass through zero.
*   **Initial Step Size**: While adaptive methods adjust the step size, a reasonable initial guess for `h` (e.g., based on the scale of the problem or known dynamics) can help convergence.
*   **System Transformation**: Remember that higher-order ODEs can always be converted into a system of first-order ODEs. This is a common practice before applying most numerical solvers in this library.
*   **Boundary Conditions**: For BVPs, the shooting method often requires an iterative process (like Newton's method) to find the correct initial conditions. The `urwp.pas` unit encapsulates this complexity.

By combining the various methods provided, users can tackle a wide array of differential equation problems encountered in scientific and engineering domains.

## File Descriptions

This section provides a brief description of each Pascal source file within the `diffeqa` directory.

*   `_info.txt`: Provides a general description and a list of programs contained in the directory.
*   `adambash.txt`: Explains the theoretical background of Adams-Bashforth methods.
*   `adambash.pas`: Program demonstrating the Adams-Bashforth method for a first-order ODE.
*   `adammoul.pas`: Program demonstrating the Adams-Moulton Prediction-Correction method for a first-order ODE.
*   `ab_mou.pas`: Unit implementing the Adams-Bashforth-Moulton predictor-corrector method for first-order ODE systems with automatic step size control.
*   `bulirsch.pas`: Unit implementing the Bulirsch-Stoer-Gragg extrapolation method for first-order ODE systems.
*   `eulromb.txt`: Explains the theoretical background of the Euler-Romberg method.
*   `eulromb.pas`: Program demonstrating the Euler-Romberg method for a first-order ODE.
*   `eqdif1.pas`: Unit containing Runge-Kutta procedures (`Equadiff1`, `Equadiffn`, `Equadiffp`) for solving first-order ODEs (single, N-order, and systems of P variables).
*   `equdif.pas`: Unit providing `odeint1`, a Runge-Kutta integration procedure with time step control for first-order ODE systems.
*   `fgauss.pas`: Unit implementing the Gauss algorithm for solving linear equations, including LU decomposition. It is a dependency for several stiff system solvers.
*   `gear.pas`: Unit implementing the implicit Gear method of 4th order for solving stiff first-order ODE systems.
*   `laplace.pas`: Program solving Laplace's equation on a square plate with given boundary conditions using the relaxation method (Example 1).
*   `laplace1.pas`: Program solving Laplace's equation on a rectangular plate with a hole using the relaxation method (Example 2).
*   `laplace2.pas`: Program solving Laplace's equation on a square plate with mixed Dirichlet and Neumann boundary conditions using the relaxation method (Example 3).
*   `limits.pas`: Program solving a two-point boundary value problem for a second-order DE using a Runge-Kutta integration approach.
*   `m_rwp.pas`: Test program demonstrating the shooting method (`rwp` unit) for first-order boundary value problems, using various IVP solvers.
*   `mgear.pas`: Test program for the `gear.pas` unit, demonstrating the implicit Gear method for stiff systems.
*   `rkutta.txt`: Explains the theoretical background of the Runge-Kutta method.
*   `stormer.txt`: Explains the theoretical background of Stormer's method for second-order ODEs.
*   `stormer.pas`: Program demonstrating Stormer's method for a second-order ODE.
*   `t_dgls.pas`: Unit containing example definitions of right-hand side functions (`dgl`) for various systems of first-order differential equations, used by several test programs.
*   `tabmou.pas`: Test program for the `ab_mou.pas` unit (Adams-Bashforth-Moulton method with automatic step size control).
*   `tawp.pas`: Test program for the `uawp.pas` unit (automatic step size control with various embedded Runge-Kutta formulas).
*   `tbulirsc.pas`: Test program for the `bulirsch.pas` unit (Bulirsch-Stoer-Gragg method).
*   `teqdif1.pas`: Basic test program for `Equadiff1` (Runge-Kutta order 1).
*   `teqdif1a.pas`: Test program for `Equadiff1` with a graphical output option.
*   `teqdifc.pas`: Test program for a general prediction-correction method.
*   `teqdifn.pas`: Basic test program for `Equadiffn` (Runge-Kutta order N).
*   `teqdifn1.pas`: Test program for `Equadiffn` with a graphical output option.
*   `teqdifp.pas`: Basic test program for `Equadiffp` (Runge-Kutta for P variables).
*   `teqdifp1.pas`: Test program for `Equadiffp` with a graphical output option.
*   `tequdif.pas`: Test program for the `equdif.pas` unit (`odeint1` procedure, RK with time step control).
*   `trk4.pas`: Test program for solving a 2-variable first-order ODE system using Runge-Kutta.
*   `trkf45.pas`: Test program for the `urkf45.pas` unit (Runge-Kutta-Fehlberg method).
*   `tros4.pas`: Test program for the Rosenbrock method for stiff systems.
*   `tshoot.pas`: Test program demonstrating the shooting method for a second-order boundary value problem.
*   `uawp.pas`: Unit providing the `awp` procedure for automatic step size control, using different embedded Runge-Kutta formulas.
*   `urkf45.pas`: Unit implementing the Runge-Kutta-Fehlberg method (RKF45) for systems of ODEs.
*   `urwp.pas`: Unit implementing the shooting method (`rwp` procedure) for first-order boundary value problems, which relies on other IVP solvers.
*   `utils.pas`: (Renamed to `Utils1` in some uses statements) Contains general utility functions such as vector/matrix operations, input/output helpers, and basic mathematical functions.
*   `Type_def.pas`: Contains common type definitions (`REAL_AR`, `Table`, `RV`, etc.) used across multiple units, standardizing data structures for numerical operations.
*   `WinCrtMy.pas`: Provides Windows console-specific input/output and screen manipulation routines, often used by test programs for basic console interaction.
*   `WinProcs.pas`: Contains Windows API procedure declarations, typically for low-level system interactions, potentially used for graphical output or environment setup in older Pascal compilers.
*   `Graph_2D.pas`: Offers 2D graphics functionalities for plotting solutions, used by some test programs to visualize results.
