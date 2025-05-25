# jpmMath

## Introduction

**jpmMath** is a collection of numerical algorithms and mathematical routines originally developed by Jean-Pierre Moreau, now preserved and modernized for Free Pascal. Jean-Pierre Moreau’s extensive library of math code (written in Turbo Pascal and other languages) was at risk of disappearing when his personal website went offline in 2023. This repository aims to **prevent the loss** of those valuable algorithms by collecting, updating, and maintaining them in a modern Object Pascal framework. The project not only **archives the legacy code** but also **refactors** it for current platforms, making the routines easier to use, more modular, and compatible with Free Pascal/Lazarus.

## Repository Structure

The repository is organized into two main sections:

-   **`original/`** – An archive of Jean-Pierre Moreau’s original source code and documentation. This folder contains the original Turbo Pascal units and programs, organized by topic as they were on Moreau’s site. For example, you will find subfolders like `functions1`, `functions2`, `geometry`, `graphic`, `linearprog`, `lstsqr`, etc., containing the original `.pas` files, example programs (often prefixed with `t...` for “test”), and any accompanying notes (`.txt`) or PDFs describing the algorithms.
    
-   **`fpc/`** – The modernized Free Pascal version of the library. Code here is updated to work with **Free Pascal** (FPC) and Lazarus. The structure under `fpc/` mirrors the categories of the original code, but with each routine or algorithm in its own module for clarity and reuse. Notable subdirectories include:
    
    -   **`fpc/function1/`** – Contains ports of routines from the original “functions1” collection. Each subfolder under `function1` corresponds to a specific algorithm or routine, with its source code and demos. _For example:_ the `akima/` subfolder implements the **Akima interpolation** algorithm.
        
    -   **`fpc/stat/`** – Contains statistical routines. For instance, there are modules for probability distributions (e.g. normal distribution in `stat/normal/` and Student’s t-distribution in `stat/student/`), each with their own test/demo programs.
        
    -   _(As the project evolves, other categories from the original codebase – such as geometry, linear programming, special functions, etc. – will be added under `fpc/` with a similar structure.)_
        

Each algorithm’s folder under `fpc/functions1/` (and other categories) typically includes: the Pascal **unit** implementing the routine, an **explanation document** (if available from Moreau’s originals), and one or more **demo programs** showing how to use the routine. Demo programs might be provided in two forms – a simple console application (CLI demo) and a Lazarus project with a GUI – to illustrate usage.

## Modernizing the Legacy Code

The original code dates back to Turbo Pascal era and was written in a monolithic style. In jpmMath, we are **modernizing** these routines to improve their integration and maintainability. Key refactoring and modernization steps include:

-   **Eliminating global variables:** Original programs often used global state to pass data to routines. We refactor these into functions/procedures with parameters and return values. This makes each routine **reentrant** and easier to reuse in different contexts without side effects.
    
-   **Replacing static arrays with dynamic arrays:** Legacy Pascal code commonly used fixed-size or static arrays (often with 1-based indexing). We update these to use **dynamic arrays** (`array of ...`) or generics, allocated at runtime with `SetLength`. This allows the routines to handle variable-sized data sets and avoid arbitrary limits compiled into the code.
    
-   **Using 0-based indexing:** Along with dynamic arrays, we adopt 0-based indexing (the default in modern Object Pascal for dynamic arrays). This means adjusting array loops and algorithms that were originally written with 1-based indices. By doing this, the code aligns with Free Pascal conventions and avoids off-by-one errors when mixing old and new code. (Where appropriate, we also update algorithm logic that assumed 1-based indexing, to ensure the math still works correctly with 0-based arrays.)
    
-   **Encapsulation and modular design:** Each routine is turned into a **self-contained unit** or class. For example, a root-finding algorithm that was a standalone program is now a unit (e.g. `unit AkimaInterpolation`) providing a procedure or function (say `InterpolateAkima(...)`) that can be called from any application. Internal helper procedures are made **local** to the unit or turned into methods, and shared utilities are moved to a common `utils.pas` or similar module if needed. This modular approach improves **reusability** – you can include just the units you need in your project.
    
-   **Code style and safety:** We update code to comply with Free Pascal’s mode (usually ObjFPC or Delphi mode), use descriptive naming where possible, add comments, and remove deprecated or unsafe constructs. Memory management is updated (e.g., no manual `GetMem` with unbounded size; instead use dynamic allocation or managed types). These changes make the code more maintainable and safer (less prone to memory errors).
    

The modernization is done gradually – each algorithm is reviewed and refactored while ensuring that its **original functionality** (and accuracy of results) remains intact. We include the original sources and documentation (_info.txt or PDFs_) to verify that the refactored code produces the expected outputs.

## License ##
The original author does not specify any conditions how his code can be used. After inquiry about this (https://www.lazarusforum.de/viewtopic.php?p=150255#p150255), his son responded in a mail:

```
Hello,

My father, now very old, who created this website to capitalize and share his passion for
developing numerical computing algorithms, would be very happy if his work could benefit other people.

You can recover and reuse his work with great pleasure.

Kind regards,
François-Xavier Moreau
```

## Usage Guide

Using the jpmMath routines typically involves adding the appropriate unit to your uses clause and calling the provided functions. This section describes how to compile and run a module using the example of the **Akima interpolation** module:

-   **Compiling the demo program:** To see the Akima interpolation in action, you can compile the provided console demo. Navigate to `fpc/function1/akima/cli_demo/` and compile `akima.pas` using Free Pascal:
    

-   `fpc akima.pas` 
    
    This will produce an executable (on Windows, `akima.exe`) that runs a demonstration of the Akima spline interpolation. The program uses sample data points, calls the interpolation routine, and prints the interpolated results so you can verify the output.
    
-   **Running the GUI demo:** If you prefer a graphical demonstration, open the Lazarus project file `fpc/function1/akima/gui_demo/Akima_Demo.lpi` in the Lazarus IDE. This will load a small GUI application that showcases the Akima interpolation (for example, plotting the spline through a set of points). Compile and run it from Lazarus – a form will appear, demonstrating how the interpolation smooths the data. _(Note: you may need to have the Lazarus Component Library (LCL) installed, which comes with Lazarus by default, to run GUI demos.)_
    
-   **Using the Akima routine in your code:** The core interpolation code is in the unit `fpc/function1/akima/akima.pas`. To use it in your own project, add this unit to your program’s uses clause. The Akima unit provides a routine (or set of routines) to perform the interpolation. Typically, you would pass arrays of known data points (X and Y values) to an initialization function that computes the Akima coefficients, then call an evaluation function to get interpolated Y for any X. Consult the comments in `akima.pas` or the `akima.txt` explanation file for the exact function names and expected parameters. An example usage might look like:
    

-   `uses Akima;  // make sure the compiler can find akima.pas in its search path
    
    var
      X, Y: array of Double;
      interpValue: Double;
    begin
      // Suppose X and Y are filled with known data points...
      AkimaPrepare(X, Y);                   // (Pseudo-code) Compute interpolation coefficients for dataset
      interpValue := AkimaInterpolate(5.0); // Get interpolated Y value at X = 5.0
    end.` 
    
    The actual function names may differ, but the idea is that after preparation you can interpolate any value. Refer to the module’s documentation for details.
    
-   **Compilation notes:** All modernized code should compile with **Free Pascal (FPC)** 3.x or later. We recommend using the **Lazarus IDE** for an easy start, especially when running the GUI demos. Each demo comes with a `.lpi` project file for Lazarus. If compiling via command-line FPC, ensure you pass the `-MObjFPC` or `-MDelphi` switch if needed (depending on the mode used in the unit’s source). The code is mostly platform-independent and uses standard FPC units. When using any unit, make sure to include its source path in your project or compile with the `-Fu` flag pointing to the appropriate directory (e.g., `-Fu/path/to/jpmMath/fpc/function1/akima`).
    

## Implemented Routines and Algorithms

The jpmMath library covers a **broad range of mathematical algorithms**, from interpolation and integration to statistics and optimization. Below is a list of some of the key routines currently available (or in the process of being ported), grouped by topic, along with a brief description of each:

-   **Interpolation and Curve Fitting:**
    
    -   _Akima spline interpolation_ – A method of interpolating data with a piecewise smooth curve that minimizes oscillations (useful for smooth curve fitting without overshoot). **(Ported to FPC)**
        
    -   _Lagrange polynomial interpolation_ – Constructing a polynomial that passes through a given set of points. (Present in original code, to be modernized)
        
    -   _Chebyshev interpolation and approximation_ – Using Chebyshev polynomials to approximate functions (reduces Runge’s phenomenon).
        
    -   _Rational interpolation_ – e.g., algorithms for Thiele’s continued fraction interpolation (`tratint` in original).
        
    -   _Least-squares polynomial fitting_ – Routines for polynomial regression/approximation of data (see the `lstsqr` folder for linear least squares and non-linear regression examples).
        
-   **Optimization and Root Finding:**
    
    -   _Newton-Raphson method_ – Iterative root-finding algorithm for solving f(x)=0 using function derivatives.
        
    -   _Brent’s method_ – A hybrid root-finding algorithm (combining bisection, secant, and inverse quadratic interpolation) for finding roots bracketed in an interval with fast convergence.
        
    -   _Golden section search_ – One-dimensional optimization (finding the minimum of a unimodal function) using the golden ratio to reduce interval size.
        
    -   _Brent’s minimization_ – (Related to above) Brent’s algorithm for finding the minimum of a function without derivatives, often used in conjunction with golden section and parabolic interpolation.
        
    -   _Powell’s method_ – Multidimensional optimization algorithm that performs a directional search to find a local minimum of a function without needing derivatives.
        
    -   _Nelder-Mead simplex (“Amoeba”) method_ – Another derivative-free optimization that uses a simplex of points to converge to a minimum (great for complicated landscapes in several variables).
        
    -   _Steepest descent method_ – An optimization algorithm that uses the gradient (if available) to descend towards a minimum (Moreau’s code includes variants and demonstrations of this).
        
    -   _Function minimization utilities_ – e.g., routines to bracket minima (`mnbrak`) and test programs for the above algorithms.
        
-   **Numerical Integration (Quadrature):**
    
    -   _Romberg integration_ – Adaptive Richardson extrapolation scheme to improve trapezoidal rule approximations of definite integrals (the code integrates a function to a desired accuracy).
        
    -   _Simpson’s rule (adaptive Simpson)_ – Classical method dividing intervals into subintervals for better polynomial approximation of the integrand.
        
    -   _Gaussian quadrature_ – High-accuracy integration using Gauss-Legendre nodes and weights. The code includes routines to perform Gauss–Kronrod or Gauss–Legendre integration over an interval. (For example, modules named `kubgauss.pas` and `kubnec.pas` in the original likely implement Gauss and Newton–Cotes formulas for cubature/integration.)
        
    -   _Newton–Cotes formulas_ – Integration using evenly spaced points (trapezoidal, Simpson, Boole’s rule, etc.).
        
    -   _Clenshaw-Curtis quadrature_ – An integration technique using Chebyshev polynomials (there is a `clencurt.pas` in the original code).
        
    -   The library often provides **test programs** (like `tromberg.pas`, `tsimpson.pas`, `tgauss.pas`) that compare the performance and accuracy of these integration methods on example integrals.
        
-   **Differential Equations (ODE) Solvers:**
    
    -   _Runge-Kutta 4th order (RK4)_ – A fixed-step ODE solver for initial value problems.
        
    -   _Runge-Kutta-Fehlberg 45 (RKF45)_ – An adaptive step-size ODE solver (4th-5th order) for better accuracy control.
        
    -   _Adams-Bashforth / Adams-Moulton methods_ – Multi-step ODE solvers (the presence of units like `uawp.pas` suggests an Adams predictor-corrector implementation).
        
    -   _Shooting method for boundary-value problems_ – The code (e.g., `tshoot.pas`) demonstrates solving ODE boundary conditions by shooting and adjusting initial slopes.
        
    -   These solvers are accompanied by example programs (e.g., solving classic equations or systems) and utility units for common ODE tasks. For instance, `urkf45.pas` and `urwp.pas` are units that implement the integration logic, and `teqdif*.pas` are test programs for differential equations.
        
-   **Linear Algebra and Linear Programming:**
    
    -   _Linear system solvers:_ While not explicitly a separate module in the listing, some routines (like `sysmat.bas` mentioned in documentation) suggest solutions for linear systems or matrix operations might be included (possibly in supporting code for least squares or other analyses). This may be an area for future extension in FPC.
        
    -   _Linear Programming (Simplex method):_ The `linearprog` folder contains code for solving linear programming problems. This includes:
        
        -   **Simplex algorithm** – The classical method for optimizing a linear objective subject to linear constraints (see `simplex.pas` and its test `tsimplex.pas`).
            
        -   **Dantzig’s algorithm** – Related to simplex, possibly specialized or two-phase simplex approach (`dantzig.pas`).
            
        -   **Transportation problem** – A module for solving transportation optimization problems (`transpor.pas` and `tpert.pas` for perturbation?).
            
        -   These routines demonstrate operations research methods and can be used to solve resource allocation or cost minimization problems.
            
-   **Statistical Functions:**
    
    -   _Probability distributions:_ The repository now includes routines for statistical distributions. For example,
        
        -   **Normal (Gaussian) distribution** – Code to compute the Gaussian distribution function (and possibly its inverse or related error function). In `fpc/stat/normal/` the unit `unormal.pas` (modernized as `jpmnormal.pas`) provides the standard normal distribution calculations, and `normal_test.pas` demonstrates usage (e.g., evaluating the cumulative distribution or generating random normals).
            
        -   **Student’s t-distribution** – A routine to evaluate the Student’s _t_ distribution, useful for statistical analyses when sample sizes are small. The `fpc/stat/student/` module includes `student_test.pas` to show how to compute _t_ distribution values for given degrees of freedom.
            
        -   _(Planned)_ Other distributions such as Chi-square, Poisson, etc., could be added in the future to broaden the statistical toolkit. The groundwork is laid for a **statistics module** where each distribution or test will have its own unit and demo.
            
    -   _Statistical tests and utilities:_ Code for computing moments, basic stats (mean, variance), or performing curve fitting (chi-square goodness of fit, etc.) is present in the original library and will be progressively updated.
        
-   **Special Functions:**  
    The collection features implementations of various special mathematical functions, many of which are not readily available in standard libraries. These include:
    
    -   **Beta and Gamma functions:** e.g., routines for the Beta function (`mbeta.pas`) and possibly Gamma or related functions.
        
    -   **Bessel functions:** Documentation outside this repo notes Moreau had code for Bessel functions (based on work by Amos and others). It’s likely that Bessel or Airy functions appear in the collection (for instance, `mairya.pas` and `mairyzo.pas` suggest Airy function Ai and Ai’ computations).
        
    -   **Hypergeometric and Legendre functions:** There are modules like `mchgm.pas` which likely computes confluent hypergeometric functions, and `legendre.pas` for Legendre polynomials (and `eval_leg.pas` to evaluate them). Associated are `laguerre.pas` (Laguerre polynomials) and `hermite.pas` (Hermite polynomials) for orthogonal polynomial calculations.
        
    -   **Elliptic integrals:** The file `cliptic.pas` (with `cliptic.txt` explanation) points to algorithms for complete elliptic integrals or related elliptic functions.
        
    -   **Miscellaneous**: Bernoulli numbers (`mbernoa.pas`), Euler numbers, sine/cosine integrals (`sinint.pas` computes Si(x), Ci(x) integrals), and more. Each of these special functions comes with references or example usages in the original texts. As they are ported, they will be useful for scientific computing tasks requiring higher mathematics beyond basic math libraries.
        
-   **Chaos, Fractals, and Graphics (Demos):**  
    In addition to the core numerical routines, Jean-Pierre Moreau’s code included a variety of educational/demo programs illustrating mathematical concepts:
    
    -   _Fractals:_ e.g., generating the Mandelbrot set (`mandel.pas` / `mandbrot.pas`) or Julia sets (`julia.pas`).
        
    -   _Chaos and Dynamical Systems:_ The **Lorenz attractor** (`lorentz.pas`), **Henon map** (`henon.pas`), **Logistic map (Verhulst)** (`verhulst.pas`), and **Rossler attractor** (`roessler.pas`) are present, showing how to simulate these systems.
        
    -   _Games and puzzles:_ There’s a `knight.pas` (likely the Knight’s tour problem) and `queens.pas` (N-Queens problem), illustrating backtracking algorithms.
        
    -   _Geometry calculations:_ Programs like `triangle.pas` (triangle geometry solver), `arcircle.pas` (circle through three points), `conical.pas` (conic sections), and `surface.pas` (surface of revolution calculations) provide routines for geometric computations.
        
    -   _Graphics utilities:_ Low-level drawing routines (e.g., `graph_2.pas`, `graph_3d.pas`) for plotting 2D/3D graphs in text mode or simple graphics mode.
        
    
    These are mostly found under `original/graphic` or `original/geometry`. They serve as interesting examples or visual demonstrations rather than reusable library components. Modernizing these is a lower priority, but they remain part of the repository’s preserved legacy. (They can often still be compiled with FPC in a console or using the Graph unit, if one is interested in retro graphics.)
    

As of now, **Akima interpolation, normal distribution, and Student’s t-distribution** are fully ported to the new framework (with more on the way). Many other algorithms listed above are available in the `original` directory and will be incrementally modernized. The aim is to eventually have every routine from Moreau’s collection usable in Free Pascal with clean, updated code.

## Future Plans and Contributing

This project is ongoing and welcomes contributions from the community. Future goals include:

-   **Porting remaining routines:** We plan to methodically refactor and test each algorithm from the original codebase. There’s a rich assortment of algorithms (as listed above) still to be modernized. If you need a particular routine or just want to help, feel free to port one of the remaining routines. For each routine, create a new folder under the appropriate `fpc/` category (for example, `fpc/function1/<routine_name>/`) and implement it following the established patterns (no globals, use dynamic arrays, etc.). Include a small demo if possible to verify its correctness.
    
-   **Improving documentation:** While many routines have an associated text or PDF explanation from the original author, we welcome improved documentation in English. This could include elaborating on the mathematical background of the algorithm, describing the function signatures and expected inputs/outputs, and writing usage examples. Enhancing the README (or creating Wiki pages) for individual modules will help users understand and trust the routines.
    
-   **Writing tests:** We aim to build a suite of unit tests for these numerical algorithms. Contributors can help by adding FPC **unit tests** for the routines (for example, comparing the output of an integration routine against a known high-precision result, or checking that root-finders actually find roots of a test function). Robust tests will ensure that the modernization process hasn’t introduced errors and that future changes don’t break functionality.
    
-   **Performance tuning:** Some legacy code might benefit from optimizations (using modern constructs or algorithms improvements). As long as the algorithm’s integrity is maintained, performance improvements (especially for large data or many iterations) are appreciated. Free Pascal’s compiler is quite powerful, and we can also consider using multi-threading or other enhancements for heavy computations, where appropriate.
    
-   **Expanding the library:** Beyond porting the legacy code, jpmMath could grow with new routines that fit its spirit. If you have a numerical algorithm implemented in Pascal that isn’t easily available elsewhere, you could add it here for others to use. We would review it for quality and consistency with the project’s coding style.
    

**How to contribute:** If you’d like to contribute, you can fork the repository, make your changes (ensure they compile with FPC and include any necessary demo or test), and open a Pull Request. It’s a good idea to open an issue or discussion first if you plan a major change, so we can coordinate efforts and ensure no duplicate work. We also welcome **bug reports** or **feature requests** – if you encounter an issue with a routine or need a particular algorithm, please open an issue on GitHub.

By contributing to jpmMath, you’ll be helping to preserve and enhance a valuable set of numerical methods for Pascal programmers. This project not only saves historical code from obscurity but also updates it for modern use. We appreciate any help in keeping these mathematical tools alive and improving them for future users.

**Enjoy exploring jpmMath!** Whether you are studying numerical analysis, building engineering software, or just curious about algorithms, we hope this library proves useful. With your help, we will continue to expand and refine this collection, bridging the past and present of Pascal numerical computing.
