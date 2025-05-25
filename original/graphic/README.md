
# jpmMath Library: Graphic

This repository contains the source code for the graphic programs and supporting units developed as part of the jpmMath. These programs primarily focus on mathematical visualization, simulations, and fractal generation, leveraging Pascal for graphical output.

## Core Units

The following Pascal units provide foundational functionalities for the graphical programs in this directory.

### `CrtGr2D.pas`
*   **Description**: This unit provides a set of procedures for 2D graphics with manual scaling and clipping capabilities. It allows for drawing lines, circles, polygons, axes, grids, and frames within a defined window.
*   **Key Features**:
    *   `Fenetre`: Defines the limits in physical (real) coordinates.
    *   `Cloture`: Defines the limits in window pixels for drawing and clipping.
    *   `MoveXY`, `LineXY`: Draw lines between points in physical coordinates.
    *   `Bordure`: Draws a frame around the drawing zone.
    *   `Axes`, `Gradue`, `Grille`: Functions to draw and graduate coordinate axes and a grid.
    *   `Cercle`, `Circle1`: Procedures to draw circles. `Circle1` uses a faster algorithm.
*   **Role in Project**: `CrtGr2D.pas` is a fundamental unit for almost all 2D graphical applications in this collection, providing the core drawing primitives and coordinate system management.

### `Equa_dif.pas`
*   **Description**: This unit implements numerical integration procedures, specifically a fourth-order Runge-Kutta method (`RK4`) with adaptive step size control (`rkqc`) for solving first-order ordinary differential equation systems.
*   **Key Features**:
    *   `odeint1`: The main integration procedure, managing the overall integration process from a start time to an end time with a specified precision.
    *   `RK4`: Implements the classical fourth-order Runge-Kutta algorithm for a single step.
    *   `rkqc`: Provides adaptive step size control, adjusting the time step to meet a required precision.
    *   `Derivs`: A placeholder procedure that must be defined by the calling program to specify the differential equations (`dy/dt = F(y,t)`).
*   **Role in Project**: Essential for simulating dynamic systems and trajectories, as demonstrated in programs like `planets.pas` and `rocketv2.pas`.

### `Graph_2.pas`
*   **Description**: This unit provides basic 2D geometry procedures for vector and point manipulations.
*   **Key Features**:
    *   `Point2D`, `Vecteur2D`, `Transform2D`: Custom data types for 2D points, vectors, and transformation matrices.
    *   `Norme2D`: Calculates the magnitude (length) of a 2D vector.
    *   `Dist2D`: Calculates the distance between two 2D points.
    *   `GetVecteur2D`: Creates a 2D vector from two points.
*   **Role in Project**: Used by `apollo.pas` for geometric calculations related to circles and by `rocketv2.pas` through `Graph_2D` (likely referring to this unit).

### `Graph_3D.pas`
*   **Description**: This unit provides procedures for drawing 3D curves and surfaces, including support for perspective and parallel projections.
*   **Key Features**:
    *   `InitProj`: Initializes auxiliary variables for projection calculations based on observer position (`rho`, `theta`, `phi`).
    *   `Project`: Transforms 3D physical coordinates (x, y, z) into 2D projected coordinates on the screen.
    *   `MoveXYZ`, `DrawXYZ`: Draw lines between points in 3D physical coordinates, automatically applying the defined projection.
    *   `Border`: Draws a border around the 3D drawing zone.
    *   `Axes`: Draws and labels 3D coordinate axes.
*   **Role in Project**: Enables the visualization of three-dimensional mathematical concepts, as seen in `surfaces.pas` and `surfpara.pas`.

### `LargArra.pas`
*   **Description**: This unit provides a mechanism to handle dynamic arrays larger than 64 KB, overcoming memory limitations in older Pascal environments. It manages memory in segments and provides an `Index` method to access elements.
*   **Key Features**:
    *   `LgArray` Object: Manages the dynamic allocation and access of large arrays.
    *   `Init`, `Done`: Constructor and destructor for memory management.
    *   `Index`: Returns a pointer to an element at a given logical index.
*   **Role in Project**: Specifically used by `rocketv2.pas` to store extensive numerical results of the trajectory simulation over time.

### `Maths_2D.pas`
*   **Description**: This unit contains general 2D mathematical utility procedures.
*   **Key Features**:
    *   `real_arSwap`: Swaps the values of two real-arithmetic variables.
*   **Role in Project**: Provides basic mathematical helper functions, primarily used by `Graph_2.pas`.

### `Type_def.pas`
*   **Description**: This simple unit defines common data types used across various units and programs within the project, such as `real_ar` (a real number type with arbitrary precision, likely `Double` or `Extended` in Free Pascal) and `pi`.
*   **Role in Project**: Ensures consistent type definitions throughout the codebase.

### Supporting Windows Interface Units

The following units are part of the Windows-specific environment (likely Turbo Pascal for Windows or Free Pascal with WinCRT) and are essential for enabling the graphical output and user interaction of the programs. While they might include custom adaptations, their primary role is to bridge the Pascal application with the underlying Windows API for display, input, and printing.

*   `Savecrt.pas`: Handles saving and loading monochrome `.IMG` picture files, primarily used by `CrtGr2D.pas` for screenshot functionality and by `wvisu.pas` for viewing.
*   `WinCrtMy.pas`, `WinCrt1.pas`: Custom or adapted `WinCrt` units that provide console-like input/output and integrate with graphical procedures under Windows. They enable window creation and access to the device context (`HDC`).
*   `WinPrint.pas`: Provides functionalities for printing graphical output to a printer.
*   `WObjects.pas`: Likely part of a Windows Objects library, providing object-oriented interfaces for Windows programming elements.
*   `WinProcs.pas`: Contains Windows API procedure declarations.
*   `WinTypes.pas`: Defines Windows data types.
*   `Time.pas`: Provides timing functionalities, used to measure program execution duration.
*   `Strings.pas`: Standard string manipulation routines.
*   `WinDos.pas`: Provides DOS-like functions for Windows environments.
*   `BWCC.pas`, `StdDlgs.pas`: Likely related to Borland Windows Custom Controls and standard dialog boxes (e.g., file open/save dialogs).

## Graphical Programs

This section details the individual graphical programs available in this folder, their purpose, usage, and potential real-world applications.

### `apollo.pas`: Apollonius Circles
*   **Description**: This program dynamically generates and displays Apollonius circles. It randomly defines three initial circles tangent to each other and then recursively draws all internal tangent circles until their radius becomes too small.
*   **Functionality**:
    *   Random generation of initial tangent circles.
    *   Recursive drawing of Apollonius circles.
    *   Demonstrates clipping using `CrtGr2D.pas`.
    *   Displays the average drawing time per step.
*   **Usage**: The program runs automatically, displaying circles in steps. Press the space bar to stop after a step.
*   **Real-world Problem Solving Ideas**:
    *   **Packing Problems**: Visualizing optimal packing arrangements in 2D space, useful in logistics, material science, or manufacturing.
    *   **Fluid Dynamics**: Simulating particle distribution in confined spaces or flow patterns.
    *   **Computational Geometry**: Demonstrating fundamental concepts of circle packing and geometric recursion, relevant to CAD/CAM or graphics algorithms.

### `billard.pas`: Elliptical Billiard Simulation
*   **Description**: Simulates the trajectory of a ball on an elliptical billiard table, visualizing the envelope of trajectories for different initial conditions (throwing the ball between foci, between focus and edge, or through a focus).
*   **Functionality**:
    *   Graphical representation of an elliptical billiard table.
    *   User input for initial trajectory (line equation `Y = MX + N`).
    *   Calculation and visualization of ball rebounds (up to 100 rebounds).
    *   Highlights the foci of the ellipse.
*   **Usage**: Input `M` and `N` for the line equation. The program then draws the ellipse, foci, and the ball's trajectory.
*   **Real-world Problem Solving Ideas**:
    *   **Optics and Acoustics**: Illustrating reflection properties of elliptical shapes, applicable to design of reflectors, whispering galleries, or acoustic focusing systems.
    *   **Physics Simulations**: Demonstrating conservation of angular momentum or energy in non-trivial geometries.
    *   **Architectural Design**: Exploring sound or light propagation within elliptically shaped rooms or structures.

### `bolygone.pas`: Bolygone Figures
*   **Description**: This program draws "bolygones," which are mathematical figures created by connecting points on a circle. It's an envelope of chords joining an angle `i` to `n*i` on a circle.
*   **Functionality**:
    *   User input for the "order" `n` (typically 2 to 30) and "step" in degrees (typically 1 to 4).
    *   Generates and draws complex, intricate patterns based on simple geometric rules.
*   **Usage**: Enter the desired order `n` and step value. The program will draw the corresponding bolygone.
*   **Real-world Problem Solving Ideas**:
    *   **Art and Design**: Generating intricate patterns for textile design, architectural ornamentation, or digital art.
    *   **Geometric Modeling**: Exploring mathematical curves and their properties in design software.
    *   **Education**: Visualizing trigonometric functions and their graphical outcomes in a dynamic way.

### `chaos.pas`: Feigenbaum Diagram
*   **Description**: Visualizes the Feigenbaum bifurcation diagram, which illustrates the behavior of the logistic map (`x_{n+1} = (1+r)x_n - rx_n^2`) as the growth rate `r` varies. It demonstrates the transition from stable attractors to chaotic behavior.
*   **Functionality**:
    *   Plots the long-term behavior (attractors) of the logistic map.
    *   Shows the period-doubling bifurcations leading to chaos.
    *   Highlights Feigenbaum's constant.
*   **Usage**: The program automatically draws the Feigenbaum diagram. No user input required during drawing.
*   **Real-world Problem Solving Ideas**:
    *   **Population Dynamics**: Modeling and understanding chaotic behavior in biological populations.
    *   **Economic Models**: Analyzing nonlinear systems in economics, such as market fluctuations or investment growth.
    *   **Climate Science**: Studying complex, non-linear atmospheric or oceanic processes that exhibit chaotic characteristics.

### `cyclo.pas`: Cycloidal Curves
*   **Description**: Draws various cycloidal curves, including Epicycloids and Hypocycloids, as well as their extended (Epitrochoid) and shortened (Hypotrochoid) versions. These curves are generated by a point on a circle rolling along another fixed circle.
*   **Functionality**:
    *   User input for the radii of the fixed (`R1`) and mobile (`R2`) circles, and a parameter `L` for extension/shortening.
    *   Option to draw the generator circles.
    *   Calculates and displays the period of the curve.
*   **Usage**: Input `R1`, `R2`, and `L`. Choose whether to draw the generator circles. The program then plots the resulting cycloidal curve.
*   **Real-world Problem Solving Ideas**:
    *   **Mechanical Engineering**: Designing gears, cams, or other mechanical components where specific rolling motion trajectories are required.
    *   **Robotics**: Planning complex trajectories for robotic arms or mobile robots.
    *   **Typography and Graphic Design**: Creating unique and aesthetically pleasing curve forms for fonts, logos, or illustrations.

### `fractal.pas`: Fractal Curves
*   **Description**: Generates various fractal curves, including Von Koch's Snowflakes and H Fractals, using an iterative function system approach where a "basis" figure is transformed by a "generatrix" figure.
*   **Functionality**:
    *   Offers multiple fractal examples (Von Koch, H-fractal, and others).
    *   User input for the fractal "order" to control complexity.
    *   Displays the base figure, the generatrix, and the final fractal.
*   **Usage**: Select an example number and the fractal order. The program will draw the components and the resulting fractal.
*   **Real-world Problem Solving Ideas**:
    *   **Computer Graphics**: Generating realistic natural landscapes (coastlines, mountains), foliage, or textures.
    *   **Antenna Design**: Creating fractal antennas that are compact and multi-band.
    *   **Data Compression**: Utilizing fractal compression techniques for images.
    *   **Medical Imaging**: Analyzing fractal dimensions in biological structures (e.g., blood vessels, brain folds).

### `hanoi.pas`: Hanoi Towers
*   **Description**: Simulates the classic "Towers of Hanoi" puzzle, demonstrating the disk shifts required to move `N` disks from one peg to another following the rules, using a recursive algorithm.
*   **Functionality**:
    *   User input for the number of disks.
    *   Animates the movement of disks between pegs.
    *   Displays the steps taken to solve the puzzle.
*   **Usage**: Enter the number of disks. The program then graphically animates the solution.
*   **Real-world Problem Solving Ideas**:
    *   **Algorithm Visualization**: A clear demonstration of recursion, a fundamental concept in computer science.
    *   **Logistics and Planning**: Illustrating sequential dependencies and optimization in multi-step processes.
    *   **Cognitive Science**: Studying problem-solving strategies and cognitive load.

### `henon.pas`: Henon's Attractors
*   **Description**: Explores the Hénon map, a 2D discrete-time dynamical system that produces a chaotic attractor. It demonstrates how simple non-linear equations can lead to complex, seemingly random, but structured patterns.
*   **Functionality**:
    *   Plots the iterated points of the Hénon map.
    *   Allows user input for the angle parameter `a`.
    *   Illustrates the concept of strange attractors.
*   **Usage**: Input an angle in radians (e.g., 1.5732). The program then generates and plots the Henon attractor.
*   **Real-world Problem Solving Ideas**:
    *   **Cryptography**: Generating pseudorandom sequences for secure communication.
    *   **Financial Modeling**: Exploring chaotic behavior in stock market data or financial systems.
    *   **Ecological Modeling**: Understanding complex interactions in ecosystems.

### `julia.pas`: Julia Set
*   **Description**: Generates and visualizes the Julia set, a fractal defined in the complex plane based on an iterative function `z^2 + c`, where `c` is a fixed complex constant and `z` is a variable. The shape of the Julia set critically depends on the choice of `c`.
*   **Functionality**:
    *   User input for the complex plane intervals (X1, X2, Y1, Y2), maximum iterations, and the complex constant `c` (p, q).
    *   Renders intricate fractal patterns, showcasing different "shapes" (e.g., dragon, jewel, cactus).
*   **Usage**: Provide the plotting boundaries, number of iterations, and the complex constant `c`. The program will render the Julia set.
*   **Real-world Problem Solving Ideas**:
    *   **Image Processing**: Generating unique textures and patterns for computer graphics, special effects, or generative art.
    *   **Complex Systems Analysis**: Visualizing the basins of attraction in non-linear systems.
    *   **Mathematical Research**: Studying the properties of complex dynamics and fractals.

### `lorentz.pas`: Lorentz Attractor
*   **Description**: Visualizes the Lorentz attractor, a set of chaotic solutions for a simplified model of atmospheric convection. It's a classic example of a "strange attractor" and demonstrates sensitive dependence on initial conditions (the "butterfly effect").
*   **Functionality**:
    *   Simulates a 3D dynamic system defined by three differential equations.
    *   Plots the trajectory of the system in 2D projection (x-z plane).
    *   Demonstrates chaotic behavior and the formation of a strange attractor.
*   **Usage**: The program automatically runs and displays the Lorentz attractor. Press any key to stop.
*   **Real-world Problem Solving Ideas**:
    *   **Meteorology and Climatology**: Understanding the inherent unpredictability and chaotic nature of weather and climate systems.
    *   **Control Systems**: Designing robust control systems that account for chaotic inputs or system responses.
    *   **Theoretical Physics**: Exploring non-linear dynamics and the emergence of chaotic behavior in physical systems.

### `macfunc.pas`: Maclaurin's Series Approximation
*   **Description**: This program demonstrates the approximation of functions (specifically `sin(x)` and `exp(x)`) using Maclaurin's series expansion. It visually shows how the approximation improves with an increasing number of terms.
*   **Functionality**:
    *   User selects the function to approximate (`sin(x)` or `exp(x)`).
    *   User specifies the number of terms in the series.
    *   Plots the true function alongside its Maclaurin series approximation.
*   **Usage**: Choose the function and the number of terms. The program will incrementally draw the approximation, pausing after each term's contribution.
*   **Real-world Problem Solving Ideas**:
    *   **Numerical Methods**: Illustrating fundamental concepts of series expansion and numerical approximation, crucial for scientific computing and engineering.
    *   **Signal Processing**: Approximating complex waveforms using simpler series for analysis or synthesis.
    *   **Algorithm Optimization**: Understanding how computational complexity changes with the number of terms in an approximation.

### `mandbrot.pas`: Mandelbrot Set (Zoom)
*   **Description**: This program allows for an interactive "zoom" into specific regions of the Mandelbrot set, showcasing its infinite complexity and self-similarity at different scales. It's an extension of `mandel.pas`.
*   **Functionality**:
    *   Renders the Mandelbrot set within user-defined complex plane intervals.
    *   Supports different color mappings based on the number of iterations for divergence.
    *   Allows interactive exploration by redefining zoom regions.
*   **Usage**: Input the desired X and Y intervals for the complex plane and the number of iterations. The program will then render the zoomed view.
*   **Real-world Problem Solving Ideas**:
    *   **Generative Art**: Creating visually stunning and complex images for digital art, design, and media.
    *   **Image Compression**: As with other fractals, exploring fractal compression techniques.
    *   **Research in Dynamics**: Studying scaling properties and universal constants in complex systems.

### `mandel.pas`: Mandelbrot Set (Complete)
*   **Description**: Generates and displays the complete Mandelbrot set, a iconic fractal in the complex plane defined by the iteration `z_{n+1} = z_n^2 + c`, where `c` is a fixed complex constant and `z_0 = 0`. Points `c` for which the sequence remains bounded form the set.
*   **Functionality**:
    *   Renders the entire Mandelbrot set within a default complex plane region.
    *   Uses color to indicate divergence speed, visually distinguishing the set from its surroundings.
    *   Employs symmetry to optimize drawing.
*   **Usage**: The program automatically draws the Mandelbrot set. No initial user input is required beyond starting the program.
*   **Real-world Problem Solving Ideas**:
    *   **Benchmarking**: Can be used as a simple benchmark for floating-point performance in graphical rendering.
    *   **Algorithmic Optimization**: Exploring techniques for fast fractal rendering (e.g., exploiting symmetry, escape time algorithms).
    *   **Mathematical Education**: A powerful visual tool for teaching complex numbers, iteration, and the concept of fractals.

### `mira.pas`: Mira's Aesthetic Chaos
*   **Description**: Visualizes the Mira iterative model, a 2D discrete-time system that produces surprisingly beautiful and complex chaotic patterns depending on its parameters.
*   **Functionality**:
    *   Plots the iterations of the Mira map (`x_{n+1} = b.y_n + F(x_n)`, `y_{n+1} = -x_n + F(x_{n+1})`).
    *   User input for parameters `a`, `b`, initial points `x0`, `y0`, number of iterations, and view limits.
*   **Usage**: Input the parameters `a`, `b`, initial coordinates `x0`, `y0`, number of iterations, and the desired view limits. The program then generates the pattern.
*   **Real-world Problem Solving Ideas**:
    *   **Art Generation**: Creating abstract and aesthetically pleasing visual designs for various applications, from screen savers to digital art.
    *   **Pattern Recognition**: Studying the formation of complex patterns from simple rules.
    *   **Chaos Theory**: Further exploration of different chaotic systems and their unique attractors.

### `planets.pas`: Eight Planets Attraction Problem
*   **Description**: Simulates the gravitational interactions and motions of up to eight celestial bodies (planets) around a central star. It uses a Runge-Kutta numerical integration method with adaptive time stepping to maintain precision.
*   **Functionality**:
    *   Calculates 2D trajectories of multiple bodies under mutual gravitational attraction.
    *   Employs the `Equa_dif.pas` unit for robust numerical integration.
    *   Visualizes the orbital paths.
*   **Usage**: The program automatically sets up and simulates an 8-planet system. The user can stop the simulation by pressing a key.
*   **Real-world Problem Solving Ideas**:
    *   **Astrodynamics**: Simulating orbital mechanics for satellite trajectories, space mission planning, or understanding asteroid movements.
    *   **Computational Physics**: Demonstrating N-body simulations and numerical stability in classical mechanics problems.
    *   **Educational Tool**: Visualizing Kepler's laws and gravitational concepts in a dynamic way.

### `queens.pas`: Eight Queens Problem
*   **Description**: Solves the classic "Eight Queens Problem," which involves placing eight chess queens on an 8x8 chessboard such that no two queens threaten each other. The program uses a backtracking algorithm to find and display all possible solutions.
*   **Functionality**:
    *   Allows selection of board size (4 to 8 queens).
    *   Finds and lists all valid solutions using a backtracking method.
    *   Graphically displays chosen solutions on a chessboard.
*   **Usage**: Input the size of the chessboard. The program will calculate all solutions and then allow the user to select and view specific solutions.
*   **Real-world Problem Solving Ideas**:
    *   **Algorithm Design**: A prime example for teaching and understanding backtracking algorithms, which are crucial for solving constraint satisfaction problems.
    *   **Resource Allocation**: Solving allocation problems where resources have conflicting constraints.
    *   **Combinatorial Optimization**: Demonstrating how to systematically search for solutions in a large problem space.

### `rk.pas`: Differential Equations System Visualization
*   **Description**: Visualizes the integral curves (solutions) of a system of first-order ordinary differential equations of the form `y' = f(x,y,z)` and `z' = g(x,y,z)`. It uses the Runge-Kutta method to draw trajectories through various starting points in the phase plane.
*   **Functionality**:
    *   Plots integral curves for user-defined differential equations (example provided: `y' = x^2*y^2 - 1`, `z' = x^2+y^2 - 4`).
    *   Allows setting the physical window and drawing step.
    *   Handles "switching" direction to draw the curve in both positive and negative directions.
*   **Usage**: Input the desired physical window limits and the drawing step. The program then plots the integral curves.
*   **Real-world Problem Solving Ideas**:
    *   **Dynamical Systems Analysis**: Visualizing phase portraits and understanding the behavior of dynamic systems in engineering, biology, or economics.
    *   **Control Systems**: Analyzing system stability and response by observing integral curves.
    *   **Chemical Reaction Kinetics**: Modeling and visualizing the concentration changes of reactants over time.

### `rocketv2.pas`: V2 Rocket Trajectory Simulation
*   **Description**: Simulates the trajectory of a V2 rocket, accounting for gravitational forces, atmospheric drag (ARDC model), and Earth's rotation (Coriolis effect). It uses robust numerical integration and provides detailed output and graphical analysis of the trajectory.
*   **Functionality**:
    *   Models 3D rocket flight from launch to impact.
    *   Includes atmospheric density model, thrust calculations, and Coriolis correction.
    *   Stores and displays key trajectory parameters (altitude, speed, dynamic pressure, distance to target).
    *   Provides graphical views of vertical trajectory, horizontal trajectory, and a zoomed-in impact zone.
*   **Usage**: The simulation runs automatically with predefined parameters. After computation, a menu allows choosing different graphical views of the results.
*   **Real-world Problem Solving Ideas**:
    *   **Aerospace Engineering**: Simulating ballistic missile trajectories, rocket launches, or re-entry dynamics.
    *   **Defense Applications**: Trajectory prediction for ballistics, targeting, and interception.
    *   **Atmospheric Science**: Studying the effects of atmospheric conditions on flight.

### `roessler.pas`: Roessler Attractor
*   **Description**: Visualizes the Roessler attractor, another example of a strange attractor derived from a simplified system of three differential equations. It exhibits complex, non-periodic behavior similar to the Lorenz attractor but with a simpler topological structure.
*   **Functionality**:
    *   Simulates a 3D dynamic system defined by three differential equations.
    *   Plots the trajectory of the system in a 2D projection (x-y+z plane).
    *   Demonstrates chaotic behavior and the formation of a strange attractor.
*   **Usage**: The program automatically runs and displays the Roessler attractor. Press any key to stop.
*   **Real-world Problem Solving Ideas**:
    *   **Circuit Design**: Modeling non-linear electronic circuits that exhibit chaotic behavior.
    *   **Neuroscience**: Investigating simplified models of neuronal activity and chaotic firing patterns.
    *   **Complex Systems**: Exploring minimal models that exhibit chaos in various scientific disciplines.

### `surfaces.pas`: 3D Surfaces (Z=F(X,Y))
*   **Description**: Draws 3D surfaces defined by explicit equations of the form `Z = F(X,Y)`. It includes options for real perspective, ordinary parallel, dimetric parallel, and isometric parallel projections, with hidden section removal.
*   **Functionality**:
    *   User inputs for X and Y intervals, number of lines and points for rendering.
    *   Choice of projection type (perspective or various parallel projections).
    *   Removes hidden lines/sections for clearer visualization.
*   **Usage**: Input the intervals for X and Y, resolution parameters, view type, and projection type with corresponding angles/distances. The program then renders the surface.
*   **Real-world Problem Solving Ideas**:
    *   **Data Visualization**: Representing 3D data sets (e.g., topographical maps, scientific data from experiments, statistical distributions) as surfaces.
    *   **Engineering Design**: Visualizing component geometries, stress distributions, or fluid flow fields.
    *   **Architecture**: Designing and visualizing complex building envelopes or structural forms.

### `surfpara.pas`: 3D Parametric Surfaces (X=F(U,V), Y=G(U,V), Z=H(U,V))
*   **Description**: Draws 3D surfaces defined by parametric equations `x = f(u,v)`, `y = g(u,v)`, and `z = h(u,v)`. It provides examples for common shapes like ellipsoids, spheres, tori, and hyperboloids, with various projection options.
*   **Functionality**:
    *   User selects the surface type (e.g., Ellipsoid, Sphere, Torus, Hyperboloid).
    *   User inputs for projection kind (perspective/parallel) and viewing angles/distances.
    *   User defines the parameter intervals (U, V) and step sizes.
*   **Usage**: Select a surface type, projection parameters, and `U`, `V` ranges. The program then draws the parametric surface.
*   **Real-world Problem Solving Ideas**:
    *   **Computer-Aided Design (CAD)**: Modeling complex shapes for product design, automotive, or aerospace industries.
    *   **Reverse Engineering**: Reconstructing 3D models from point cloud data.
    *   **Art and Sculpture**: Designing intricate 3D forms for digital or physical fabrication.

### `verhulst.pas`: Verhulst Diagram
*   **Description**: Visualizes the Verhulst diagram, demonstrating the growth of a population under limited resources using the logistic map. It shows how population dynamics can stabilize, oscillate, or become chaotic depending on the growth rate.
*   **Functionality**:
    *   Plots the population `x` over time `t` based on the logistic map `x_{n+1} = (1+r)x_n - rx_n^2`.
    *   User input for the growth rate `r`.
    *   Illustrates the behavior of population models.
*   **Usage**: Input a value for the growth rate `R` (between 0 and 3). The program then plots the Verhulst diagram.
*   **Real-world Problem Solving Ideas**:
    *   **Ecology**: Modeling population growth and resource management in biological systems.
    *   **Epidemiology**: Simulating the spread of diseases.
    *   **Resource Management**: Understanding sustainable resource utilization and potential for over-exploitation.

### `wvisu.pas`: Monochrome Picture File Viewer
*   **Description**: A utility program designed to visualize monochrome `.IMG` picture files, which are generated by other graphical programs in this collection. It provides basic file operations (open, print) and screen clearing.
*   **Functionality**:
    *   Loads and displays `.IMG` files.
    *   Allows printing of the loaded image to a printer.
    *   Provides options to clear the screen or open new files.
*   **Usage**: Launch the program, then use the "File/Open" menu to select an `.IMG` file. The image will be displayed. Right-click the mouse to clear the screen.
*   **Real-world Problem Solving Ideas**:
    *   **Legacy Data Viewing**: Useful for viewing historical graphic data from older applications.
    *   **Batch Processing**: Could be integrated into a system to review or process `.IMG` outputs from simulations.
    *   **Simple Image Debugging**: A quick way to check graphical output of other programs during development.

