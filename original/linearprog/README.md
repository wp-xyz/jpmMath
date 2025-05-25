
# jpmMath - Programming Library

This folder contains a collection of Pascal programs demonstrating various operational research and linear programming methods from the `jpmMath` library. These programs serve as educational examples and provide practical implementations for solving decision-making and optimization problems in diverse fields.

## Folder Structure

The `linearprog` directory includes the following files:

*   `_info.txt`: General information and a program index.
*   `appoint.pas`: Pascal source code for the Appointment Method.
*   `appoint.txt`: Detailed explanation and an example problem for the Appointment Method.
*   `dantzig.pas`: Pascal source code for Dantzig's Model (optimal path finding).
*   `simplex.pas`: Pascal source code for a basic Simplex Method.
*   `simplex.txt`: Detailed explanation and an example problem for the basic Simplex Method.
*   `tanneal.pas`: Pascal source code for the Traveling Salesman Problem using Simulated Annealing.
*   `tpert.pas`: Pascal source code for the Time P.E.R.T. Model (project scheduling).
*   `transpor.pas`: Pascal source code for the Transport Model.
*   `transpor.txt`: Detailed explanation and an example problem for the Transport Model.
*   `tsimplex.pas`: Pascal source code for a more general Simplex Method handling various constraint types.

## Programs Overview

### 1. The Appointment Method

*   **Source File:** `appoint.pas`
*   **Explanation File:** `appoint.txt`

#### Purpose and Kinds of Problems

The Appointment method aims to find the optimal assignment configuration among a set of possible assignments. The primary objective is generally to minimize a total cost or maximize earnings/satisfaction. This method is suitable for scenarios where distinct entities need to be assigned to distinct tasks or locations with varying efficiencies or costs for each pairing. It can effectively handle unbalanced problems (where the number of entities and assignments differ) by introducing fictitious entries with zero values.

#### Model and Algorithm

The method utilizes a matrix representing appointment costs or quantities. It employs successive iterations based on the **Hungarian Algorithm** to find the optimal assignment configuration. For maximization problems, the input matrix must be converted into a "regret" matrix (e.g., by subtracting each satisfaction value from a value greater than the maximum satisfaction). The algorithm robustly handles basic encoding issues and aims for an optimal solution by systematically reducing matrix values.

#### Usage

The `appoint.pas` program prompts the user for the number of jobs/applicants, and subsequently for each applicant, the costs or regrets for each job. The program then outputs the optimal appointments (e.g., "APPLICANT #X => JOB #Y") and the total satisfaction.

**Example Input Flow (from `appoint.pas` sample run comments):**

```
NUMBER OF JOBS ? 4

INPUT APPOINTMENT COSTS/REGRETS OF APPLICANTS:

APPLICANT #1:
   JOB #1 ? 10
   JOB #2 ? 20
   JOB #3 ? 80
   JOB #4 ? 60

APPLICANT #2:
   JOB #1 ? 10
   JOB #2 ? 30
   JOB #3 ? 70
   JOB #4 ? 20

APPLICANT #3:
   JOB #1 ? 60
   JOB #2 ? 30
   JOB #3 ? 80
   JOB #4 ? 20

APPLICANT #4:
   JOB #1 ? 50
   JOB #2 ? 60
   JOB #3 ? 80
   JOB #4 ? 40
```

#### Real-World Problem-Solving Ideas

*   **Human Resources:** Assigning employees to projects or tasks to optimize efficiency or minimize labor costs based on skill sets and preferences.
*   **Logistics:** Allocating delivery drivers to routes to minimize total travel time or fuel consumption.
*   **Manufacturing:** Assigning production machines to specific jobs to minimize setup time or production cost.
*   **Event Planning:** Matching performers to stages or time slots to maximize audience satisfaction or minimize logistical conflicts.

#### Key Variables (from `appoint.pas` comments)

*   `NP`: Number of jobs/applicants.
*   `C(NP,NP)`: Appointment cost matrix.
*   `MP(NP,NP)`: Job/applicant matrix (used for marking assignments).
*   `IF1`: Flag to indicate if an optimal appointment is found.
*   `IZ`: Flag to mark a zero in the matrix.
*   `XMIN`: Stores the minimum value during row/column reduction.
*   `XCASE`: Used in appointment of a zero in relation to minimums.
*   `NZ`: Number of zeroes found.
*   `ZI, ZJ`: Coordinates of the current case being examined.
*   `M, N`: Coordinates of a marked zero.
*   `A, B`: Auxiliary variables for changing appointments.

### 2. Dantzig's Model

*   **Source File:** `dantzig.pas`

#### Purpose and Kinds of Problems

Dantzig's Model, as implemented in `dantzig.pas`, is designed for finding either the minimum or maximum path value through a network (graph) composed of "marks" (nodes) and "arcs" (edges) with associated values (e.g., distances, costs, or benefits). It is highly suitable for problems involving network optimization where the objective is to identify the most efficient, shortest, or most beneficial route.

#### Model and Algorithm

The program implements an algorithm to traverse a graph, starting from an initial mark and exploring paths to find the optimal one (minimum or maximum value). The core logic is found within the `DANTZIG1` procedure. The algorithm iteratively updates the optimal path value and the sequence of marks (nodes) that form this path.

#### Usage

The `dantzig.pas` program first prompts the user to specify whether they want to find the minimum or maximum path. Following this, it requests the total number of marks (nodes) in the network. For each mark, the user must input the number of outgoing arcs, along with the arrival mark and the value associated with that specific arc. The program's output will then display the determined optimal path and its calculated total value.

**Example Input Flow (from `dantzig.pas` sample run comments):**

```
DO YOU WANT THE MINIMUM PATH (Y/N) ? N

NUMBER OF MARKS ? 9

NUMBER OF ARCS FROM MARK #0 ? 3
ARRIVAL MARK, VALUE ? 1 5
ARRIVAL MARK, VALUE ? 2 11
ARRIVAL MARK, VALUE ? 3 10

NUMBER OF ARCS FROM MARK #1 ? 2
ARRIVAL MARK, VALUE ? 2 4
ARRIVAL MARK, VALUE ? 4 6

NUMBER OF ARCS FROM MARK #2 ? 2
ARRIVAL MARK, VALUE ? 4 5
ARRIVAL MARK, VALUE ? 6 9

NUMBER OF ARCS FROM MARK #3 ? 1
ARRIVAL MARK, VALUE ? 5 2

NUMBER OF ARCS FROM MARK #4 ? 1
ARRIVAL MARK, VALUE ? 7 5

NUMBER OF ARCS FROM MARK #5 ? 2
ARRIVAL MARK, VALUE ? 6 4
ARRIVAL MARK, VALUE ? 8 15

NUMBER OF ARCS FROM MARK #6 ? 1
ARRIVAL MARK, VALUE ? 7 6

NUMBER OF ARCS FROM MARK #7 ? 1
ARRIVAL MARK, VALUE ? 8 8
```

#### Real-World Problem-Solving Ideas

*   **Route Optimization:** Finding the shortest travel path between two locations in a road network for navigation systems or logistics planning.
*   **Critical Path Analysis:** Identifying the longest sequence of activities in a project, which determines the minimum project completion time (critical path).
*   **Network Flow Optimization:** Determining optimal data routing in telecommunication networks to minimize latency or maximize throughput.
*   **Resource Distribution:** Finding the most efficient distribution routes for goods or services from a source to multiple destinations.

#### Key Variables (from `dantzig.pas` comments)

*   `R`: 'Y' for minimum path, 'N' for maximum path.
*   `NR`: Total number of marks (nodes) in the network.
*   `V(NR,NR)`: Matrix storing path values (adjacency matrix of the graph).
*   `XL(NR)`: Array used for weighting each path.
*   `IE(NR*2)`: Set of adopted marks (nodes) in the optimal path.
*   `XM1, XM2`: Used to seek minimum/maximum path values.
*   `IASSO`: Arrival mark associated with the set `IE`.
*   `IO, ID`: Auxiliary variables for iterating through marks.
*   `IR`: Number of arcs coming from a specific mark.

### 3. The Simplex Method (Basic)

*   **Source File:** `simplex.pas`
*   **Explanation File:** `simplex.txt`

#### Purpose and Kinds of Problems

The Simplex Method is a foundational algorithm in linear programming used to optimize (maximize or minimize) a linear economic function subject to a set of linear constraints. This basic implementation focuses on problems where constraints are primarily upper limits (`<=`) and independent variables are non-negative.

Typical applications include:
*   Determining the optimal distribution of product components.
*   Dispatching machine working times efficiently.
*   Optimizing investment portfolios.
*   Planning production quantities to maximize output or minimize cost.
*   Maximizing benefits by optimizing selling prices.
*   Minimizing production costs.

#### Model and Algorithm

The `simplex.pas` program implements the core **Simplex Method**. It iteratively moves from one basic feasible solution to another, improving the objective function value at each step until an optimal solution is reached. Minimization problems are handled by transforming the objective function (negating its coefficients) and solving it as a maximization problem. The program identifies the optimal values for the independent variables that satisfy all given constraints.

#### Usage

The `simplex.pas` program prompts the user to choose whether to maximize or minimize the economic function. It then requests the number of variables in the economic function and the number of constraints. Users input the coefficients for the economic function and its right-hand side, followed by the coefficients and right-hand side for each constraint. The program outputs the optimal values for each variable and the final optimal value of the economic function.

**Example Input Flow (from `simplex.pas` sample run comments):**

```
MAXIMIZE (Y/N) ? Y

NUMBER OF VARIABLES OF ECONOMIC FUNCTION ? 3

NUMBER OF CONSTRAINTS ? 3

INPUT COEFFICIENTS OF ECONOMIC FUNCTION:
      #1 ? 15
      #2 ? 17
      #3 ? 20
      Right hand side ? 0

CONSTRAINT #1:
      #1 ? 0
      #2 ? 1
      #3 ? -1
      Right hand side ? 2

CONSTRAINT #2:
      #1 ? 3
      #2 ? 3
      #3 ? 5
      Right hand side ? 15

CONSTRAINT #3:
      #1 ? 3
      #2 ? 2
      #3 ? 1
      Right hand side ? 8
```

#### Real-World Problem-Solving Ideas

*   **Production Planning:** Maximizing profit by determining optimal production quantities for various products, given limitations on raw materials, labor, and machine hours.
*   **Resource Allocation:** Allocating limited resources (e.g., budget, personnel, time) across different activities or departments to achieve a specific objective (e.g., maximizing output, minimizing cost).
*   **Diet Formulation:** Creating a diet plan that meets nutritional requirements at the minimum cost by selecting appropriate quantities of different food items.
*   **Financial Portfolio Optimization:** Maximizing returns on investments while adhering to budget and risk tolerance constraints.

#### Key Variables (from `simplex.pas` comments)

*   `R`: 'Y' for maximize, 'N' for minimize.
*   `NV`: Number of variables of the economic function.
*   `NC`: Number of constraints.
*   `TS`: Simplex tableau, a 2D array storing coefficients.
*   `R1`: Multiplier (1 for maximize, -1 for minimize).
*   `R2`: Auxiliary variable for input.
*   `NOPTIMAL`: Flag indicating if iterations should continue (0 for optimal, 1 to continue).
*   `XMAX`: Stores the greater coefficient of the economic function during pivot selection.
*   `RAP`: Stores the smallest ratio (for pivot row selection).
*   `V`: Auxiliary variable.
*   `P1, P2`: Line and column indices of the pivot element.
*   `XERR`: Flag indicating if no solution exists.

### 4. The Simplex Method (General)

*   **Source File:** `tsimplex.pas`

#### Purpose and Kinds of Problems

The `tsimplex.pas` program provides a more comprehensive and robust implementation of the Simplex Method. It is designed to solve linear programming problems involving a mix of constraint types: less than or equal to (`<=`), greater than or equal to (`>=`), and exact equalities (`=`). This versatility makes it suitable for a broader range of real-world optimization challenges.

#### Model and Algorithm

This implementation employs a **Two-Phase Simplex** approach. In Phase I, it finds an initial basic feasible solution, often by introducing artificial variables for `>=` and `=` constraints. In Phase II, it proceeds to optimize the objective function, similar to the basic Simplex method, to find the global optimum. The program uses `simp1`, `simp2`, and `simp3` subroutines to manage the pivot operations and tableau updates efficiently.

#### Usage

The program first requires the user to input the number of variables in the economic function and the count for each type of constraint (`<=`, `>=`, `=`). Subsequently, the coefficients for the economic function and its constant term are entered. For each constraint, the coefficients for its variables and its constant term are then provided. The output indicates the maximum value of the economic function and the optimal values for each variable. It also reports if no feasible solution exists.

**Example Input Flow (from `tsimplex.pas` sample run comments):**

```
Number of variables in E.F.: 4
Number of <= inequalities..: 2
Number of >= inequalities..: 1
Number of = equalities.....: 1

Input Economic Function:
Coefficient # 1: 1
Coefficient # 2: 1
Coefficient # 3: 3
Coefficient # 4: -0.5
Constant term..: 0

Input constraint # 1:
Coefficient # 1: 1
Coefficient # 2: 0
Coefficient # 3: 2
Coefficient # 4: 0
Constant term..: 740

Input constraint # 2:
Coefficient # 1: 0
Coefficient # 2: 2
Coefficient # 3: 0
Coefficient # 4: -7
Constant term..: 0

Input constraint # 3:
Coefficient # 1: 0
Coefficient # 2: 1
Coefficient # 3: -1
Coefficient # 4: 2
Constant term..: 0.5

Input constraint # 4:
Coefficient # 1: 1
Coefficient # 2: 1
Coefficient # 3: 1
Coefficient # 4: 1
Constant term..: 9
```

#### Real-World Problem-Solving Ideas

*   **Comprehensive Resource Optimization:** Solving complex allocation problems where resources may be exactly required, have a minimum threshold, or have a maximum capacity.
*   **Blending Problems:** Determining the optimal mix of ingredients to produce a product with specific characteristics at minimum cost, accounting for various ingredient availabilities and quality constraints.
*   **Financial Portfolio Optimization:** Maximizing returns on investments considering diverse types of financial regulations, risk tolerance, and specific asset allocation targets.
*   **Production Scheduling:** Creating detailed production schedules that meet minimum output requirements, stay within maximum capacities, and adhere to specific processing time equalities.

#### Key Variables (from `tsimplex.pas` comments and context)

*   `M, N`: Total number of constraints (M) and variables (N).
*   `M1, M2, M3`: Number of `<=`, `>=`, and `=` constraints, respectively.
*   `A`: Simplex tableau matrix (`MAT` type).
*   `IPOSV, IZROV`: Integer vectors tracking indices of basic (left-hand) and non-basic (right-hand) variables.
*   `ICASE`: Output flag for solution status (0 = optimal, -1 = no feasible solution, 1 = unbounded solution).
*   `R`: Auxiliary variable for input.
*   `simp1, simp2, simp3`: Internal procedures for pivot selection and tableau transformation.

### 5. Time P.E.R.T. Model

*   **Source File:** `tpert.pas`

#### Purpose and Kinds of Problems

The PERT (Program Evaluation and Review Technique) model is a project management tool used for scheduling, organizing, and coordinating tasks within a project. It is particularly valuable for projects with uncertain activity durations. This program helps in determining the critical path, identifying potential bottlenecks, forecasting project completion times, and understanding activity flexibilities.

#### Model and Algorithm

The program calculates several key metrics for each activity:
*   **Probable Duration:** Calculated using a weighted average of optimistic, estimated, and pessimistic durations (formula: `(Optimistic + 4*Estimated + Pessimistic) / 6`).
*   **Soonest Begin Date (SBD):** The earliest time an activity can start without delaying subsequent activities.
*   **Latest Begin Date (LBD):** The latest time an activity can start without delaying the overall project completion.
*   **Soonest End Date (SED):** The earliest time an activity can finish.
*   **Latest End Date (LED):** The latest time an activity can finish without delaying the overall project.
*   **Total Margin (Slack):** The amount of time an activity can be delayed without affecting the overall project duration (`LBD - SBD`).
*   **Critical Activities (CA):** Activities with zero total margin, meaning any delay in these activities will directly delay the entire project.

The program also generates a simplified text-based **Gantt's Diagram** to visually represent the project schedule, showing probable duration (with `*`) and total margin (with `+`) for each activity.

#### Usage

The `tpert.pas` program first requests the total number of activities in the project and the total count of all preceding activity relationships. For each activity, the user inputs its estimated, optimistic, and pessimistic durations. Subsequently, for each activity, the user specifies its direct preceding activities. The program then provides a detailed breakdown of calculated dates and margins for each activity, culminating in a text-based Gantt's diagram.

**Example Input Flow (from `tpert.pas` sample run comments):**

```
NUMBER OF ACTIVITIES ? 9

NUMBER OF PRECEDING ACTIVITIES ? 11

ACTIVITY #1:
  ESTIMATED DURATION   ? 3
  OPTIMISTIC DURATION  ? 2
  PESSIMISTIC DURATION ? 4
  NUMBER OF PRECEDING ACTIVITIES ? 0

ACTIVITY #2:
  ESTIMATED DURATION   ? 7
  OPTIMISTIC DURATION  ? 4
  PESSIMISTIC DURATION ? 10
  NUMBER OF PRECEDING ACTIVITIES ? 1
  PREC. ACTIVITY #1 ? 1

ACTIVITY #3:
  ESTIMATED DURATION   ? 16
  OPTIMISTIC DURATION  ? 8
  PESSIMISTIC DURATION ? 24
  NUMBER OF PRECEDING ACTIVITIES ? 1
  PREC. ACTIVITY #1 ? 2

ACTIVITY #4:
  ESTIMATED DURATION   ? 4
  OPTIMISTIC DURATION  ? 2
  PESSIMISTIC DURATION ? 6
  NUMBER OF PRECEDING ACTIVITIES ? 1
  PREC. ACTIVITY #1 ? 2

ACTIVITY #5:
  ESTIMATED DURATION   ? 2
  OPTIMISTIC DURATION  ? 1
  PESSIMISTIC DURATION ? 3
  NUMBER OF PRECEDING ACTIVITIES ? 1
  PREC. ACTIVITY #1 ? 3

ACTIVITY #6:
  ESTIMATED DURATION   ? 24
  OPTIMISTIC DURATION  ? 16
  PESSIMISTIC DURATION ? 32
  NUMBER OF PRECEDING ACTIVITIES ? 2
  PREC. ACTIVITY #1 ? 3
  PREC. ACTIVITY #2 ? 4

ACTIVITY #7:
  ESTIMATED DURATION   ? 3
  OPTIMISTIC DURATION  ? 2
  PESSIMISTIC DURATION ? 4
  NUMBER OF PRECEDING ACTIVITIES ? 1
  PREC. ACTIVITY #1 ? 6

ACTIVITY #8:
  ESTIMATED DURATION   ? 16
  OPTIMISTIC DURATION  ? 8
  PESSIMISTIC DURATION ? 24
  NUMBER OF PRECEDING ACTIVITIES ? 2
  PREC. ACTIVITY #1 ? 5
  PREC. ACTIVITY #2 ? 7

ACTIVITY #9:
  ESTIMATED DURATION   ? 0
  OPTIMISTIC DURATION  ? 0
  PESSIMISTIC DURATION ? 0
  NUMBER OF PRECEDING ACTIVITIES ? 1
  PREC. ACTIVITY #1 ? 8
```

#### Real-World Problem-Solving Ideas

*   **Construction Project Management:** Scheduling construction phases (e.g., foundation, framing, roofing) and identifying critical dependencies to ensure timely completion.
*   **Software Development Life Cycle:** Planning software releases by managing dependencies between development, testing, and deployment tasks to estimate project duration and resource needs.
*   **Event Planning:** Organizing large-scale events by scheduling tasks like venue booking, vendor coordination, and marketing, to identify the critical path for a successful event.
*   **Research and Development:** Managing complex R&D projects by breaking them into smaller tasks, estimating durations, and tracking progress against a critical path to achieve research goals efficiently.

#### Key Variables (from `tpert.pas` comments)

*   `N`: Total number of activities.
*   `NA`: Total number of preceding activities (sum of all `IANTE`).
*   `D(N)`: Probable duration for each activity.
*   `PRED(N)`: Index of preceding activities (pointer to `ANTE` array).
*   `ANTE(NA)`: Array storing preceding activities for each activity.
*   `MARQ(2,N)`: Activity level and status during calculations.
*   `TOT(N)`: Soonest date for each activity.
*   `TARD(N)`: Latest date for each activity.
*   `DE`: Estimated duration.
*   `DOP`: Optimistic duration.
*   `DP`: Pessimistic duration.
*   `IANTE`: Number of preceding activities for a specific activity.
*   `K`: Number of activity levels / current level.
*   `R`: Used to seek a level (auxiliary variable).
*   `FTA`: Latest ending date.
*   `FOT`: Soonest ending date.
*   `LO`: Maximum width of Gantt's diagram.
*   `MAXT`: Probable maximum time (for Gantt scale).
*   `X`: Normalization coefficient for Gantt's diagram.
*   `Y`: Number of units for a sign (`*` or `+`) in Gantt's diagram.

### 6. The Transport Model

*   **Source File:** `transpor.pas`
*   **Explanation File:** `transpor.txt`

#### Purpose and Kinds of Problems

The Transport Model is a specialized form of linear programming designed to determine the most cost-effective dispatching of resources from multiple starting points (sources) to various destination points. Its primary objective is to minimize the total transportation cost. This model is ideal for situations where products or resources need to be moved across a network with known supply, demand, and unit transportation costs.

Common applications include:
*   Dispatching products from different factories to various customer locations at minimum cost.
*   Optimizing the distribution of liquid products through a pipeline network.
*   Optimal dispatching of electric power from generating stations to consumers.

The model explicitly handles scenarios where total supply does not equal total demand. If total demand exceeds supply, fictitious sources can be introduced. If total supply exceeds demand, a fictitious destination (representing excess stock) can be added, both with zero transportation costs, to balance the problem.

#### Model and Algorithm

The program first establishes an initial basic feasible solution using the **North-West Corner Rule**. Following this, it applies the **Stepping Stone Algorithm** to iteratively improve the solution. The Stepping Stone method identifies "unused" routes (paths) in the transportation matrix and evaluates their potential for cost reduction. By systematically moving quantities along these paths, the algorithm aims to reduce the total transport cost until an optimal distribution is achieved.

#### Usage

The `transpor.pas` program prompts the user to input the number of sources and destinations. Then, it requests the available quantities for each source and the required quantities for each destination. Finally, the user inputs the unitary transport costs for each source-destination pair. For routes where transport is impossible or undesirable, a very high value (e.g., 9999 as in the example) can be entered as the unitary cost. The program's output displays the optimal quantities to be transported from each source to each destination and the calculated total minimum transport cost.

**Example Input Flow (from `transpor.pas` sample run comments):**

```
NUMBER OF SOURCES ? 3

NUMBER OF DESTINATIONS ? 4

INPUT THE AVAILABLE SOURCE QUANTITIES:
 SOURCE #1 ? 5000
 SOURCE #2 ? 10000
 SOURCE #3 ? 6000

INPUT THE REQUIRED DESTINATION QUANTITIES:
 DESTINATION #1 ? 2000
 DESTINATION #2 ? 11000
 DESTINATION #3 ? 4000
 DESTINATION #4 ? 4000

INPUT TRANSPORT COSTS MATRIX:
 FROM SOURCE #1 TO DESTINATION #1 ? 200
 FROM SOURCE #1 TO DESTINATION #2 ? 700
 FROM SOURCE #1 TO DESTINATION #3 ? 800
 FROM SOURCE #1 TO DESTINATION #4 ? 9999
 FROM SOURCE #2 TO DESTINATION #1 ? 400
 FROM SOURCE #2 TO DESTINATION #2 ? 400
 FROM SOURCE #2 TO DESTINATION #3 ? 500
 FROM SOURCE #2 TO DESTINATION #4 ? 500
 FROM SOURCE #3 TO DESTINATION #1 ? 400
 FROM SOURCE #3 TO DESTINATION #2 ? 500
 FROM SOURCE #3 TO DESTINATION #3 ? 600
 FROM SOURCE #3 TO DESTINATION #4 ? 600
```

#### Real-World Problem-Solving Ideas

*   **Supply Chain Optimization:** Minimizing the cost of shipping goods from manufacturing plants to warehouses or from warehouses to retail stores.
*   **Logistics Planning:** Determining the most cost-effective way to move raw materials from suppliers to production sites and finished products to distribution centers.
*   **Disaster Relief:** Optimizing the distribution of aid supplies from various collection points to affected areas to minimize delivery time or cost.
*   **Waste Management:** Planning the transportation of waste from collection points to processing or disposal facilities to reduce operational expenses and environmental impact.

#### Key Variables (from `transpor.pas` comments)

*   `ISOU`: Number of sources.
*   `IDES`: Number of destinations.
*   `XM(ISOU,IDES)`: Matrix representing the resources distribution (quantities transported).
*   `C(ISOU,IDES)`: Unitary transport costs matrix.
*   `R(ISOU,IDES)`: Path search matrix (auxiliary for Stepping Stone).
*   `P(4,ISOU+IDES)`: Table of found paths.
*   `DAO(ISOU)`: Available source quantities.
*   `RAD(IDES)`: Required destination quantities.
*   `IX, IY`: Auxiliary variables for current position in transport path.
*   `CC`: Current path cost.
*   `XINFCC`: Minimum of path costs during optimization.
*   `NR`: Number of known marks.
*   `IF1`: Flag for marks still to eliminate.
*   `ID`: Displacement length.
*   `TT`: Transport total quantity.
*   `LT`: Total path length.
*   `CT`: Total transport cost.
*   `IOPTIMAL`: Flag indicating if optimality is found (1 for optimal).
*   `QT`: Transport quantity for one path.

### 7. The Traveling Salesman Problem (Simulated Annealing)

*   **Source File:** `tanneal.pas`

#### Purpose and Kinds of Problems

This program addresses the classic Traveling Salesman Problem (TSP): finding the shortest possible route that visits a given set of cities exactly once and returns to the origin city. It utilizes the **Simulated Annealing Method**, a metaheuristic inspired by the annealing process in metallurgy, to find a near-optimal solution. This method is particularly useful for complex optimization problems where finding the absolute global optimum is computationally infeasible within a reasonable timeframe.

#### Model and Algorithm

The problem is modeled as finding the shortest Hamiltonian cycle in a complete graph where cities are nodes and distances are edge weights. The program can calculate distances in two ways:
1.  **"As the crow flies" (Euclidean distance):** Calculates straight-line distances between cities based on their geographical coordinates provided in `villes.dat`.
2.  **"By road":** Uses pre-defined road distances between cities, hardcoded within the `InitDist1` and `InitDist2` procedures (specifically for a set of French towns).

The core **Simulated Annealing** algorithm works by:
*   Starting with an initial random itinerary (a permutation of cities).
*   Repeatedly proposing small changes to the current itinerary (e.g., swapping two cities, reversing a segment).
*   Accepting better itineraries.
*   Accepting worse itineraries with a certain probability that depends on a "temperature" parameter (`TP`). This temperature gradually decreases over time, allowing the algorithm to escape local optima early on and converge towards good solutions later.
The process continues for a specified number of iterations, gradually "cooling down" the "temperature" parameter. Due to its probabilistic nature, running the program multiple times might yield slightly different but often near-optimal results.

#### Usage

The `tanneal.pas` program first prompts the user to choose between "as the crow flies" (option 1) or "by road" (option 2) for distance calculation. Subsequently, it asks for the starting number of iterations, the maximum number of iterations, and the increment value for iterations. It reads town data (number of towns, names, and coordinates) from a file named `villes.dat` and writes the results (the shortest itinerary found and its total distance) to `villes.lst`.

**Example `villes.dat` format:**
```
38
AMIENS           525 115
ANGOULEME        365 585
... (each line: TOWN_NAME X_COORDINATE Y_COORDINATE)
```

**Example Input Flow (from `tanneal.pas` sample run comments):**

```
Distances between towns (km):
    1 : as the crow flies
    2 : by road
Your choice (1 or 2): 2

Starting number of iterations: 10
Maximum number of iterations : 2000
Increment value of iterations: 10
```

#### Real-World Problem-Solving Ideas

*   **Logistics and Delivery Routing:** Optimizing delivery routes for parcel services, postal services, or food delivery to minimize travel time, fuel consumption, and operational costs.
*   **Vehicle Routing Problems (VRP):** While TSP is a specific case, it forms the fundamental building block for more complex VRPs involving multiple vehicles, capacities, and time windows.
*   **Circuit Board Drilling:** Optimizing the path of a drill head to make holes in a circuit board, reducing manufacturing time and increasing efficiency.
*   **DNA Sequencing:** Ordering DNA fragments to reconstruct the original DNA strand, where minimizing "travel" between fragments corresponds to efficient sequence assembly.
*   **Tour Planning:** Creating efficient travel itineraries for tourists visiting multiple attractions in a region.

#### Key Variables (from `tanneal.pas` comments)

*   `NV`: Number of towns to visit.
*   `CH`: Table of town names.
*   `Vs`: Stores the shortest itinerary found.
*   `D`: Matrix of distances in km (for both options).
*   `ICO`: Geographical coordinates of towns (X, Y).
*   `IT1, LL`: Indices for steps of resolution/iterations.
*   `choice`: User's choice for distance calculation (1 for air, 2 for road).
*   `diter`: Increment value of iterations.
*   `itermin`: Stores the number of iterations at which the minimum distance was found.
*   `nbreiter`: Starting number of iterations.
*   `nmax`: Maximum number of iterations.
*   `DDMIN`: Overall minimum distance found across all runs.
*   `DMAX`: Maximum distance found in a single annealing run.
*   `DMIN`: Minimum distance found in a single annealing run.
*   `DT`: Total distance of the current itinerary.
*   `TP`: "Temperature" parameter in Simulated Annealing, which decreases over time.
*   `X, Y`: Auxiliary variables for coordinate differences (used in "as the crow flies" calculation).
