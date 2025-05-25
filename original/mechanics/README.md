
# jpmMath Mechanics

This repository contains a comprehensive collection of Pascal programs designed to solve various problems in mechanical engineering, encompassing structural analysis, vibrations, and composite materials. The programs are developed with a focus on clarity and demonstrate the application of fundamental mathematical and physical principles.

Each program is either self-contained or relies on a small set of auxiliary units provided within this library, making it suitable for educational purposes, study, modification, and extension. This documentation aims to provide a detailed description of each program's functionality, its theoretical underpinnings, mathematical principles, input requirements, output interpretation, and potential applications in real-world engineering scenarios.

---

## Programs Overview

The `mechanics` folder contains programs addressing the following areas:

*   **Single Degree of Freedom (1-DOF) Systems**:
    *   Response to sinusoidal input force with viscous or structural damping (`1dof01.pas`).
    *   Response to periodic input force using Fourier series (`1dof02.pas`).
*   **Beam Analysis**:
    *   Calculation of static beam deflection under various support and loading conditions (`beam.pas`).
    *   Determination of resonance frequencies, modal mass, stiffness, and deformation modes for bending beams (`beam1.pas`).
*   **Runge-Kutta Solvers & Analogous Systems**:
    *   LRC circuit or damped mass-spring problem (`circuit.pas`).
    *   General differential equations of order 1 (`eqdifp.pas` unit, used by `pendulum.pas`).
    *   Angular motion of an elementary mass pendulum (`pendulum.pas`).
*   **Composite Materials Analysis**:
    *   Stresses in a unidirectional composite layer (`compos01.pas`).
    *   Stiffness matrix of laminated composite materials (`compos03.pas`).
    *   Deformations and stresses in laminated composite materials under imposed efforts (`compos04.pas`).
*   **Multi-Degree of Freedom (N-DOF) Systems**:
    *   Frequencies and modes of undamped mass-spring systems via Transfer Method (`modes.pas`).
    *   Response to sinusoidal force using Transfer Matrices Method (`ndof01.pas`).
    *   Frequencies, eigenmodes, and modal properties of undamped mass-spring systems (`ndof02.pas`).
    *   Response to sinusoidal force using a Direct Method (`ndof03.pas`).
    *   Step-by-step solution of dynamic systems using the Wilson-Theta Method (`ndof04.pas`).
*   **Specialized Simulations**:
    *   Bouncing ball trajectory with damping (`rebounds.pas`).
*   **Finite Element Method (FEM) Demonstration**:
    *   EF3D - a small finite elements program for 3D structures (`ef3d/ef3d.pas`).
*   **Information File**:
    *   Program descriptions overview (`_info.txt`).

---

## Detailed Program Documentation

### 1. `1dof01.pas`: Response of a 1-DOF Mass-Spring System with Sinusoidal Input Force

**Program Description:**
This program calculates the steady-state displacement amplitude and phase response of a single degree-of-freedom (1-DOF) mass-spring system when subjected to a sinusoidal external force. It provides options to model the damping as either viscous or structural.

**Core Functionality:**
The program determines the steady-state amplitude and phase lag of the system's displacement as a function of excitation frequency. It prompts the user for system parameters (Mass `M`, Stiffness `K`, Force `F0`) and damping type (Viscous `C` or Structural `E`). It then performs a frequency sweep within a specified range and step.

**Key Variables:**
*   `M`: Mass
*   `K`: Stiffness
*   `F0`: Amplitude of the sinusoidal input force
*   `C`: Viscous Damping Coefficient (if `B='V'`)
*   `E`: Structural Damping Factor (if `B='S'`)
*   `F`: Excitation Frequency
*   `W`: Angular excitation frequency (`2 * PI * F`)
*   `M1`: Displacement Amplitude
*   `T3`: Phase (in degrees)

**Mathematical Principles:**
The program solves the steady-state response of a 1-DOF system to harmonic excitation, considering the relevant damping model (viscous or structural). The underlying equations of motion are:
*   For viscous damping: `M X" + C X' + K X = F0 sin(Wt)`
*   For structural damping: `M X" + K(1+iE) X = F0 sin(Wt)` (where `E` is the structural damping factor applied to stiffness, `i` is the imaginary unit)
It calculates the amplitude of displacement and the phase difference between the force and displacement for a range of frequencies using complex impedance methods.

**Sample Usage (from source comments):**
```
 Mass M..... = 50
 Stiffness K = 1e6
 Force F0... = 10000

 (V)iscous or (S)tructural Damping: S

 Structural damping factor: 0.05

 Frequency Scanning

 Resonance Frequency of undamped system = 22.508 Hz.

 Starting frequency = 20
 Ending frequency...= 25
 Frequency step.....= 0.25

 Frequency   Displacement   Phase (deg)
 --------------------------------------
  20.00      0.046234       13.366013
  20.25      0.050756       14.701446
  20.50      0.056293       16.347699
  20.75      0.063206       18.423074
  21.00      0.072037       21.111461
  21.25      0.083609       24.711299
  21.50      0.099181       29.729354
  21.75      0.120525       37.058260
  22.00      0.149218       48.252814
  22.25      0.181993       65.500492
  22.50      0.199980       89.194985
  22.75      0.183564      -66.609176
  23.00      0.149839      -48.520651
  23.25      0.119585      -36.721514
  23.50      0.097048      -29.028193
  23.75      0.080680      -23.790753
  24.00      0.068578      -20.053143
  24.25      0.059388      -17.273961
  24.50      0.052222      -15.136029
  24.75      0.046502      -13.444959
  25.00      0.041843      -12.076311
```

**Real-World Problem Solving Ideas:**
*   **Vibration Isolator Design**: Evaluate the effectiveness of a single-stage vibration isolation system (e.g., engine mounts, sensitive equipment mounts) by predicting its response amplitude and phase at various operating frequencies to minimize force transmission or displacement.
*   **Resonance Avoidance in Simple Structures**: Determine the dynamic response of a simplified component (e.g., a cantilevered beam, a panel) to harmonic forces, aiding in design modifications to shift resonance away from operating frequencies or to manage amplitude at resonance.
*   **Damping System Design**: Evaluate the effectiveness of different damping mechanisms (viscous vs. structural) in reducing vibration amplitudes. For instance, comparing the displacement reduction achieved by adding a viscous damper versus a material with structural damping.
*   **Accelerometer Calibration**: Understand frequency response characteristics of vibration sensors.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Mécanique des vibrations linéaires" By M. Lalanne, P. Berthier, J. Der Hagopian, Masson, Paris 1980 [BIBLI 16].

### 2. `1dof02.pas`: Response of a 1-DOF Mass-Spring System with Viscous Damping to a Periodic Input Force

**Program Description:**
This program computes the time-domain displacement response of a 1-DOF mass-spring system with viscous damping when subjected to a general periodic input force, represented by its Fourier series components.

**Core Functionality:**
The program takes as input the system's mass (`M`), stiffness (`K`), viscous damping coefficient (`C`), and the Fourier series coefficients (`A0`, `A(P)`, `B(P)`) of the periodic excitation force. It then calculates the displacement response `X(t)` over a specified time range and step.

**Key Variables:**
*   `M`: Mass
*   `K`: Stiffness
*   `C`: Viscous Damping Coefficient
*   `F8`: Fundamental excitation frequency (Hz)
*   `W`: Angular excitation frequency (`2 * PI * F8`)
*   `N`: Number of harmonic terms in the Fourier series (P from 1 to N)
*   `A0, A[P], B[P]`: Fourier series coefficients
*   `T`: Current time
*   `T1, T2`: Starting and ending time for scanning
*   `T3`: Time step
*   `X[J]`: Displacement at time `T`

**Mathematical Principles:**
The program solves the equation of motion `M X" + C X' + K X = F(t)`. The periodic force `F(t)` is represented as a Fourier series: `F(t) = A0/2 + Sum [A(P) COS(PWt) + B(P) SIN(PWt)]`. The solution is obtained by superimposing the steady-state responses to each harmonic component of the force. For each harmonic, the program calculates the complex steady-state response and then combines them in the time domain to get `X(t)`.

**Sample Usage (from source comments):**
```
 Mass M = 50
 Stiffness K = 1e6
 Viscous damping coefficient: 0.1

                N
 F(t) = A0/2 + Sum [A(P) COS(PWt) + B(P) SIN(PWt)]
               P=1

 N = 1

 Excitation frequency = 20

 A0 = 0

 A(1)= 0

 B(1)= 10000

 Time Scanning

 Excitation Period = 0.050

 Starting time = 0
 Ending time...= 0.05
 Time step.....= 0.0025

   Time     Displacement
 ------------------------
  0.0000    -0.000003
  0.0025     0.014682
  0.0050     0.027930
  0.0075     0.038444
  0.0100     0.045195
  0.0125     0.047521
  0.0150     0.045196
  0.0175     0.038447
  0.0200     0.027935
  0.0225     0.014688
  0.0250     0.000003
  0.0275    -0.014682
  0.0300    -0.027930
  0.0325    -0.038444
  0.0350    -0.045195
  0.0375    -0.047521
  0.0400    -0.045196
  0.0425    -0.038447
  0.0450    -0.027935
  0.0475    -0.014688
  0.0500    -0.000003
```

**Real-World Problem Solving Ideas:**
*   **Machine Component Vibrations**: Analyze the response of a machine part to operational forces generated by rotating or reciprocating elements, where the force profile can be complex but periodic (e.g., from an unbalanced shaft, gear meshing, or internal combustion engine forces).
*   **Structural Response to Wind Loads**: Model the dynamic response of a building to simplified periodic wind loads, which often have dominant harmonic components.
*   **Acoustic-Structure Interaction (Simplified)**: Predict the structural vibration induced by periodic sound waves.
*   **Control System Disturbance Analysis**: Evaluate the effectiveness of a control system in mitigating the effects of periodic disturbances on a mechanical system.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Mécanique des vibrations linéaires" By M. Lalanne, P. Berthier, J. Der Hagopian, Masson, Paris 1980 [BIBLI 16].

### 3. `beam.pas`: Calculate Beam Deflection for Four Support/Loading Systems

**Program Description:**
This program calculates the maximum static deflection of a beam for four standard support and loading configurations. These include simply supported beams with concentrated or uniformly distributed loads, and cantilever beams with concentrated or uniformly distributed loads.

**Core Functionality:**
The user inputs the beam's material elasticity (`elasticity`), moment of inertia (`moment_of_inertia`), and length (`length`). Then, a support/loading system is chosen from four predefined cases:
1.  Supported at each end, concentrated load.
2.  Supported at each end, uniformly distributed load.
3.  Supported at one end (cantilever), concentrated load at the free end.
4.  Supported at one end (cantilever), uniformly distributed load.
The program then calculates and displays the maximum deflection.

**Key Variables:**
*   `elasticity`: Young's modulus of the beam material (lb/in^2).
*   `moment_of_inertia`: Area moment of inertia of the beam's cross-section (in^4).
*   `length`: Beam length (converted from ft to in internally by multiplying by 12.0).
*   `load`: Concentrated force (lb) or distributed weight (lb/ft, converted to lb/in internally by the problem context).
*   `systemID`: Integer representing the chosen support/loading system (1-4).
*   `deflection`: Calculated maximum deflection (in).

**Mathematical Principles:**
The program applies direct analytical formulas derived from the Euler-Bernoulli beam theory. These formulas relate the maximum deflection to the beam's material properties (Modulus of Elasticity `E`), cross-sectional geometry (Moment of Inertia `I`), length `L`, and applied load `P` (concentrated) or `w` (distributed). Specific constants are used for each loading/support case.
*   Supported at each end, concentrated load at center: `deflection = -load * L^3 / (48.0 * E * I)`
*   Supported at each end, uniformly distributed load: `deflection = -5.0 * load * L^4 / (384.0 * E * I)`
*   Cantilever, concentrated load at free end: `deflection = -load * L^3 / (3.0 * E * I)`
*   Cantilever, distributed load: `deflection = -load * L^4 / (8.0 * E * I)`

**Sample Usage (from source comments):**
```
 Give elasticity (lb/in^2) and moment of inertia (in^4): 30e6 797
 Give the beam length in ft: 20
 Choose one of these support/loading systems:
 1 - supported at each end, concentrated load
 2 - supported at each end, uniformly distributed load
 3 - supported at one end, concentrated load at free end
 4 - supported at one end, distributed load
 Input your choice (1 to 4): 1
 Give the concentrated force: 50000

 Deflection =  -0.6023

 More? (y/n): y

 Give elasticity (lb/in^2) and moment of inertia (in^4): 30e6 797
 Give the beam length in ft: 20
 Choose one of these support/loading systems:
 1 - supported at each end, concentrated load
 2 - supported at each end, uniformly distributed load
 3 - supported at one end, concentrated load at free end
 4 - supported at one end, distributed load
 Input your choice (1 to 4): 3
 Give the concentrated force: 10000

 Deflection =  -1.9272

 More? (y/n): n
```

**Real-World Problem Solving Ideas:**
*   **Preliminary Structural Sizing**: Rapidly estimate the necessary cross-sectional dimensions of beams in buildings, bridges, or machine frames to meet serviceability deflection limits (e.g., maximum permissible sag).
*   **Material and Section Optimization**: Compare different materials (varying `E`) or beam cross-sections (varying `I`) to achieve desired deflection performance or minimize material usage.
*   **Load Capacity Assessment**: Determine the maximum allowable load for a beam given its geometric and material properties and a maximum permissible deflection.
*   **Educational Demonstrations**: Illustrate the direct application of fundamental solid mechanics principles for beam design.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Problem Solving with Fortran 90" By David R. Brooks, Springer-Verlag New York, 1997.

### 4. `beam1.pas`: Resonance Frequencies, Modal Mass and Stiffness, Deformation Modes, and Maximum Strain of a Bending Beam

**Program Description:**
This program determines the natural (resonance) frequencies, associated modal masses, modal stiffnesses, and the corresponding deformation mode shapes (eigenvectors) of a uniform bending beam. It supports six different end conditions and also calculates the slope and maximum strain distribution along the beam for each mode.

**Core Functionality:**
The user specifies the beam's limit conditions (Fixed-Free, Fixed-Supported, Fixed-Fixed, Free-Free, Free-Supported, Supported-Supported) and material/geometric properties (Young's Modulus `E`, Volumic Mass `R`, Beam Width `B`, Beam Thickness `H`, Beam Length `L`). The program then computes the first five natural frequencies. Optionally, it can display the deformation modes (deviation, slope, strain) for a specified number of points along the beam and the corresponding modal mass and stiffness for each mode.

**Key Variables:**
*   `M`: Integer representing the chosen limit condition (1-6).
*   `A1`: Array storing `Xn^2` values for different modes and boundary conditions.
*   `E`: Young's Modulus (Pa).
*   `R`: Volumic Mass (kg/m^3).
*   `B`: Beam Width (m).
*   `H`: Beam Thickness (m).
*   `L`: Beam Length (m).
*   `J`: Moment of Inertia of the cross-section (`B*H^3/12.0`).
*   `S`: Cross-sectional area (`B*H`).
*   `F`: Resonance Frequency (Hz).
*   `W`: Pulsation (rad/s).
*   `M1`: Modal Mass.
*   `K1`: Modal Stiffness.
*   `X`: Array of positions along the beam for detailed output.
*   `Y`: Deviation.
*   `P`: Slope.
*   `M1` (reused as local `M1`): Strain (`sigma`).

**Mathematical Principles:**
The program solves the eigenvalue problem for beam vibration. It applies the Euler-Bernoulli beam equation for free vibration: `EI d^4v/dx^4 + rho*S d^2v/dt^2 = 0`.
*   Natural frequencies `wn` are calculated using predefined `Xn^2` values for various boundary conditions: `wn = Xn^2 /L^2 * sqrt(EI/(rho*S))`.
*   The deformation mode shape `v(x)` is determined by the general solution `v(x) = C sin(beta.x) + D cos(beta.x) + E sh(beta.x) + F ch(beta.x)`, where coefficients C, D, E, F are determined by the boundary conditions and a normalization constant.
*   Modal mass and stiffness are derived from the mode shapes and system properties.
*   Maximum strain `sigma` is computed from the second spatial derivative of the deformation `d^2v(x)/dx^2` and the beam's Young's modulus `E`: `sigma = E * H/2 * d^2v(x)/dx^2`.

**Sample Usage (from source comments):**
```
 Fixed-Fixed Beam, M=3
 Young Modulus: 2E11
 Volumic Mass.: 7800
 Beam Width....: 0.4
 Beam Thickness: 0.5
 Beam Length...: 10


 Frequency (Hz) =   26.025373
 Frequency (Hz) =   71.739944
 Frequency (Hz) =  140.638973
 Frequency (Hz) =  232.483361
 Frequency (Hz) =  347.290041


 Do you want the modes, modal masses & Stiffnesses (y/n): y

 How many modes (Maximum 5): 5

 How many points for deviation, slope & max. strain: 11

 Do you want automatic divisions (y/n): y


        MODE #1
        -------

 OMEGA= 163.5222  FREQUENCY=  26.0254

 Modal Mass= 5.40002743637149E+0006  Modal Stiffness= 1.44394161046433E+0011

     X      Deviation   Slope    Strain (x10^6)
 ------------------------------------------------
    0.0000    0.0000    0.0000      22771.74
    1.0000    0.1925    0.3498      12232.17
    2.0000    0.6304    0.4927       2225.41
    3.0000    1.1155    0.4494      -6194.01
    4.0000    1.4814    0.2635     -11846.94
    5.0000    1.6164    0.0000     -13841.17
    6.0000    1.4814   -0.2635     -11846.94
    7.0000    1.1155   -0.4494      -6194.01
    8.0000    0.6304   -0.4927       2225.41
    9.0000    0.1925   -0.3498      12232.17
   10.0000   -0.0000   -0.0000      22771.74

        MODE #2
        -------

 OMEGA= 450.7554  FREQUENCY=  71.7399

 Modal Mass= 1.64384471757779E+0009  Modal Stiffness= 3.33997017195899E+0014

     X      Deviation   Slope    Strain (x10^6)
 ------------------------------------------------
    0.0000    0.0000    0.0000      61624.92
    1.0000    0.4554    0.7516      14015.32
    2.0000    1.2058    0.6202     -24480.45
    3.0000    1.5043   -0.0778     -40796.43
    4.0000    1.0338   -0.8261     -29766.62
    5.0000    0.0000   -1.1411         -0.00
    6.0000   -1.0338   -0.8261      29766.62
    7.0000   -1.5043   -0.0778      40796.43
    8.0000   -1.2058    0.6202      24480.46
    9.0000   -0.4554    0.7516     -14015.31
   10.0000   -0.0000    0.0000     -61624.92

        MODE #3
        -------

 OMEGA= 883.6607  FREQUENCY= 140.6390

 Modal Mass= 6.30221918268558E+0011  Modal Stiffness= 4.92112747787204E+0017

     X      Deviation   Slope    Strain (x10^6)
 ------------------------------------------------
    0.0000    0.0000    0.0000     120907.45
    1.0000    0.7701    1.1128      -6282.87
    2.0000    1.5079    0.1215     -77726.96
    3.0000    0.8687   -1.2982     -47991.93
    4.0000   -0.6284   -1.3976      39638.73
    5.0000   -1.4060   -0.0000      85988.24
    6.0000   -0.6284    1.3976      39638.73
    7.0000    0.8687    1.2982     -47991.92
    8.0000    1.5079   -0.1215     -77726.96
    9.0000    0.7701   -1.1128      -6282.87
   10.0000   -0.0000   -0.0000     120907.44

        MODE #4
        -------

 OMEGA=1460.7360  FREQUENCY= 232.4834

 Modal Mass= 2.62456494665287E+0014  Modal Stiffness= 5.60016486006314E+0020

     X      Deviation   Slope    Strain (x10^6)
 ------------------------------------------------
    0.0000    0.0000    0.0000     199859.15
    1.0000    1.0745    1.2736     -58760.62
    2.0000    1.3192   -0.9913     -120007.61
    3.0000   -0.4227   -1.9219      45103.94
    4.0000   -1.3935    0.3075     139911.06
    5.0000   -0.0000    1.9969          0.02
    6.0000    1.3935    0.3075     -139911.05
    7.0000    0.4227   -1.9219     -45103.98
    8.0000   -1.3192   -0.9913     120007.59
    9.0000   -1.0745    1.2736      58760.66
   10.0000   -0.0000    0.0000     -199859.09

        MODE #5
        -------

 OMEGA=2182.0877  FREQUENCY= 347.2900

 Modal Mass= 1.14990481414791E+0017  Modal Stiffness= 5.47527941977509E+0023

     X      Deviation   Slope    Strain (x10^6)
 ------------------------------------------------
    0.0000    0.0000    0.0000     298555.55
    1.0000    1.3218    1.1293     -144271.18
    2.0000    0.6736   -2.2318     -91130.35
    3.0000   -1.3394   -0.7648     201616.08
    4.0000   -0.2202    2.4118      33178.41
    5.0000    1.4146    0.0000     -211057.80
    6.0000   -0.2202   -2.4118      33178.38
    7.0000   -1.3394    0.7648     201616.09
    8.0000    0.6736    2.2318     -91130.32
    9.0000    1.3218   -1.1293     -144271.20
   10.0000   -0.0000   -0.0000     298555.50
```

**Real-World Problem Solving Ideas:**
*   **Vibration Control in Mechanical Systems**: Essential for designing structures and components to avoid resonance with operational or environmental excitation frequencies (e.g., rotating machinery, vehicle chassis, aircraft wings).
*   **Fatigue Life Prediction**: Identifying regions of maximum strain in vibrating components helps engineers predict potential fatigue failure locations and improve design durability.
*   **Structural Health Monitoring**: Understanding the natural frequencies and mode shapes of a structure provides a baseline for detecting damage or degradation over time (e.g., changes in mode shapes or frequencies can indicate structural integrity issues).
*   **Acoustic Design**: Used in the design of musical instruments, loudspeakers, and sound-absorbing panels where controlling vibrational modes is crucial for acoustic performance.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Mécanique des vibrations linéaires" By M. Lalanne, P. Berthier, J. Der Hagopian, Masson, Paris 1980 [BIBLI 16].

### 5. `circuit.pas`: Runge-Kutta Method for LRC Circuit or Damped Mass-Spring Problem

**Program Description:**
This program utilizes the fourth-order Runge-Kutta numerical integration method to solve the governing differential equations for either an LRC (Inductor-Resistor-Capacitor) electrical circuit or its analogous damped mass-spring mechanical system. It can simulate scenarios with or without an oscillating (time-dependent) forcing term.

**Core Functionality:**
The program allows the user to specify whether there's an oscillating term in the governing equation.
*   **Oscillating Term (Ld^2q/dt^2 + Rdq/dt + q/C = V):** Solves for charge `q` and current `i` (or displacement `x` and velocity `v`) over time. Inputs correspond to `V` (Force `F`), `L` (Mass `M`), `R` (Damping `D`), and `C` (inverse of Stiffness `1/K`). It also provides an analytic solution for comparison. A condition `4L/C - R^2 > 0` must be met for oscillatory behavior.
*   **No Oscillating Term (Ldi/dt + Ri = V):** Solves for current `i` (or velocity `v`) over time for a first-order system. Inputs correspond to `V` (Force `F`), `L` (Mass `M`), and `R` (Damping `D`). It also provides an analytic solution for comparison.

**Key Variables:**
*   `V`: Voltage (or Force `F`).
*   `L`: Inductance (or Mass `M`).
*   `R`: Resistance (or Damping Constant `D`).
*   `C`: Capacitance (or inverse of Spring Constant `1/K`).
*   `q`: Charge (or Displacement `x`).
*   `i`: Current (or Velocity `v`).
*   `t`: Current time.
*   `t_final`: End time for simulation.
*   `n`: Number of points for output.
*   `dt`: Time step (`t_final / n`).

**Mathematical Principles:**
*   **LRC Circuit Equation**: `L d^2q/dt^2 + R dq/dt + q/C = V(t)`.
*   **Mass-Spring System Equation**: `M X" + D X' + K X = F(t)`.
*   **Analogies**: The program leverages direct analogies between electrical and mechanical systems:
    *   Voltage `V` is analogous to Force `F`.
    *   Charge `q` is analogous to Displacement `x`.
    *   Current `i` is analogous to Velocity `v`.
    *   Inductance `L` is analogous to Mass `M`.
    *   Resistance `R` is analogous to Damping constant `D`.
    *   Reciprocal of Capacitance `1/C` is analogous to Spring constant `K`.
*   **Runge-Kutta Method (Order 4)**: The core of the program is the fourth-order Runge-Kutta method, which numerically approximates the solution to initial value problems of ordinary differential equations. It transforms the second-order ODE into a system of two first-order ODEs for `x` and `v` (or `q` and `i`) and applies the Runge-Kutta method for numerical integration over time.

**Sample Usage (from source comments):**
```
 (Solve mass-spring system with:
  F=75 N, M=50 kg, K=100 N (C=0.01), D=0.05)

 Specify [o]scillating or [n]o oscillating term...
 o
 Give V or F, L or M, R or D, C or 1/K (4L/C-R^2 > 0):
 75 50 0.05 0.01
 one period at t=   4.442879    s
 Give t_final and number of points:
 10 10

     time        q or d    i or speed  analytic q or d
   1.000000      0.7504      1.4993      0.6331
   2.000000      2.2485     -0.0017      1.4628
   3.000000      0.7483     -2.9955      1.0890
   4.000000     -2.2440      0.0055      0.1436
   5.000000      0.7555      5.9850      0.2223
   6.000000      6.7320     -0.0150      1.1913
   7.000000      0.7351    -11.9580      1.4147
   8.000000    -11.2021      0.0378      0.5163
   9.000000      0.7878     23.8922      0.0131
  10.000000     24.6302     -0.0915      0.7537

 (For a better accuracy, increase number of points).
```

**Real-World Problem Solving Ideas:**
*   **Electrical Engineering Design**: Simulate the transient response of RLC filters, oscillators, or power electronics circuits to various inputs (e.g., step voltage, pulse).
*   **Vehicle Suspension Modeling**: Analyze the ride dynamics of a vehicle's suspension system (approximated as a mass-spring-damper) when encountering bumps or changes in road conditions.
*   **Biomedical System Modeling**: Model simplified biological systems, such as the response of a limb to external forces or the dynamics of a blood flow system.
*   **Control System Development**: Evaluate the performance of control algorithms designed to damp oscillations or achieve desired transient responses in mechanical or electrical systems.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Problem Solving with Fortran 90" By David R. Brooks, Springer-Verlag New York, 1997.

### 6. `compos01.pas`: Stresses in a Unidirectional Composite Layer

**Program Description:**
This program calculates the stress components (`sxx`, `syy`, `sxy`) within a single, homogeneous unidirectional layer of a composite material. The stresses are provided in two coordinate systems: the global (x,y) plane and the material's principal axes (L,W), where 'L' is aligned with the fiber direction and 'W' is transverse to it. The analysis assumes a plane stress state.

**Core Functionality:**
The program takes as input the material's engineering constants (Young's Moduli `EL`, `ET`, Poisson's Ratio `NULT`, Shear Modulus `GLT`) and the applied strains (`EXX`, `EYY`, `GXY`), along with the fiber orientation angle (`TH`). It then computes:
1.  The "reduced stiffness matrix" in the main material axes (L,W).
2.  The "reduced stiffness matrix" transformed to the global (x,y) axes.
3.  The stresses (`sxx`, `syy`, `sxy`) in the global (x,y) axes.
4.  The stresses transformed back to the main material axes (L,W).

**Key Variables:**
*   `EL`, `ET`: Young's Moduli in longitudinal (fiber) and transverse directions (GPa).
*   `NULT`: Major Poisson's ratio.
*   `GLT`: Shear Modulus in the L-T plane (GPa).
*   `TH`: Angle of fibers in x-direction (degrees, converted to radians for calculations).
*   `EXX`, `EYY`, `GXY`: Normal and shear strains in the global (x,y) axes.
*   `Q`, `Q1`: Reduced stiffness matrices.
*   `S`, `S1`: Stress vectors.

**Mathematical Principles:**
The program applies the constitutive relations for orthotropic materials under plane stress.
1.  **Reduced Stiffness Matrix (`Q`)**: The material's stiffness matrix `Q` in its principal material axes (L,W) is calculated from the engineering constants (`EL`, `ET`, `NULT`, `GLT`).
2.  **Transformation Matrix**: A transformation matrix is derived based on the fiber angle `theta` to rotate stresses and strains between the global (x,y) and local (L,W) coordinate systems.
3.  **Transformed Reduced Stiffness Matrix (`Q_bar` or `Q1`)**: The `Q` matrix is transformed to the global (x,y) coordinate system to obtain `Q_bar`. This matrix relates strains in the global system (`exx`, `eyy`, `gxy`) to stresses in the global system (`sxx`, `syy`, `sxy`): `{sigma_xy} = [Q_bar] {epsilon_xy}`.
4.  **Stress Calculation**: Stresses in global axes `S` are computed by multiplying the transformed stiffness matrix `Q1` by the strain vector `E`. Stresses in main axes `S1` are then transformed back using another transformation matrix `Q`.

**Sample Usage (from source comments):**
```
   STRESSES IN A UNIDIRECTIONAL COMPOSITE LAYER

      Angle of fibers:  30.000 °
      eps. xx        :  10.000 mm
      eps. yy        :  -5.000 mm
      gam. xy        :  20.000 mm

   REDUCED STIFFNESS MATRIX (MAIN AXES):

     41.0509     3.2841     0.0000
      3.2841    10.2627     0.0000
      0.0000     0.0000     4.5000

   REDUCED STIFFNESS MATRIX (X, Y AXES):

     28.3391     8.2989     9.5611
      8.2989    12.9450     3.7706
      9.5611     3.7706     9.5148

   Stresses in axes x, y:

    433.1189 MPa
     93.6746 MPa
    267.0540 MPa

   Stresses in main axes:

    579.5334 MPa
   -52.7399 MPa
   -13.4567 MPa
```

**Real-World Problem Solving Ideas:**
*   **Micromechanics of Composites**: Understand how macroscopic loads translate to stresses at the fiber and matrix level, which is essential for predicting material failure.
*   **Ply-by-Ply Failure Analysis**: Before analyzing a full laminate, this program helps to understand the stress capacity of individual layers under various strain states.
*   **Material Characterization**: Validate experimental strain measurements against theoretical stress predictions for new composite materials.
*   **Design of Hybrid Composites**: Aid in selecting optimal fiber orientations and materials for individual plies in a multi-layered composite to achieve desired performance and avoid local failure.

**Dependencies**: `WinCrt`.

**Reference:**
*   "MATERIAUX COMPOSITES - Comportement mécanique et analyse des structures" By J.-M. Berthelot - MASSON 1996.

### 7. `compos03.pas`: Stiffness Matrix of a Laminated Composite Material

**Program Description:**
This program calculates the **A, B, and D matrices** which form the global stiffness matrix of a laminated composite material. These matrices represent the extensional, coupling, and bending stiffnesses, respectively, linking the resultant forces and moments (Nx, Ny, Nxy, Mx, My, Mxy) to the mid-plane strains and curvatures (Eps xx, Eps yy, gam xy, Kx, Ky, Kxy). The laminate is composed of multiple unidirectional composite layers, each with a specific angle and thickness.

**Core Functionality:**
The user inputs the number of layers (`NCOUCHES`), and for each layer, its fiber angle (`TH`) and thickness (`EP`). The program also requires the basic material engineering constants (`EL`, `ET`, `NULT`, `GLT`) for the unidirectional layers.
The program then calculates and displays:
1.  The "basic reduced stiffness matrix" (`Q0`) for a 0-degree ply.
2.  The transformed reduced stiffness matrix (`Q`) for each individual layer, considering its fiber angle.
3.  The `A`, `B`, and `D` matrices, which form the 6x6 laminate stiffness matrix `R`.
    *   `A` matrix (extensional stiffness): Relates in-plane forces to in-plane strains.
    *   `B` matrix (coupling stiffness): Couples in-plane forces to curvatures, and moments to in-plane strains.
    *   `D` matrix (bending-torsion stiffness): Relates moments to curvatures.

**Key Variables:**
*   `EL`, `ET`, `NULT`, `GLT`: Material engineering constants (GPa).
*   `NCOUCHES`: Number of layers.
*   `TH[I]`: Fiber angle of layer `I` (degrees, converted to radians).
*   `EP[I]`: Thickness of layer `I` (m).
*   `Q0`: Basic reduced stiffness matrix (for theta=0).
*   `Q`: Array of transformed reduced stiffness matrices for each layer.
*   `A`, `B`, `D`: Sub-matrices of the laminate stiffness matrix.
*   `R`: The full 6x6 laminate stiffness matrix.

**Mathematical Principles:**
The program strictly adheres to the Classical Laminate Plate Theory (CLPT).
1.  **Lamina Stiffness (`Q_bar` or `Q`)**: For each layer, the transformed reduced stiffness matrix (`Q_bar`) in the global (x,y) coordinate system is calculated, considering its fiber orientation.
2.  **Laminate Stiffness Matrices (A, B, D)**: These matrices are formed by integrating the `Q_bar` matrices through the laminate thickness.
    *   `A_ij = SUM (Q_bar_ij)_k * (h_k - h_(k-1))` for k layers. This matrix relates in-plane forces to in-plane strains.
    *   `B_ij = SUM (Q_bar_ij)_k * 0.5 * (h_k^2 - h_(k-1)^2)`. This matrix represents the coupling between in-plane forces and bending curvatures, and vice versa.
    *   `D_ij = SUM (Q_bar_ij)_k * (1/3) * (h_k^3 - h_(k-1)^3)`. This matrix relates moments to curvatures.
3.  The program also assembles these into a single 6x6 `R` matrix: `[[A, B], [B, D]]`.

**Sample Usage (from source comments):**
```
 CALCULATE THE STIFFNESS MATRIX OF A LAMINATE

 Number of layers: 4
  Layer 1: angle  -30.0 deg. - thickness    1.0 mm
  Layer 2: angle   15.0 deg. - thickness    1.5 mm
  Layer 3: angle  -15.0 deg. - thickness    1.5 mm
  Layer 4: angle   30.0 deg. - thickness    1.0 mm

 Material Parameters:
    El =   38.0    Et =    9.0
  NUlt =   0.32    Glt=    3.6


 BASIC REDUCED STIFFNESS MATRIX IN GPa:
   38.9445     2.9516     0.0000
    2.9516     9.2237     0.0000
    0.0000     0.0000     3.6000

 Layer #1
 REDUCED STIFFNESS MATRIX IN GPa:
   26.2896     8.1763    -9.4512
    8.1763    11.4292    -3.4183
   -9.4512    -3.4183     8.8247

 Layer #2
 REDUCED STIFFNESS MATRIX IN GPa:
   35.2120     4.6931     6.7316
    4.6931     9.4731     0.6986
    6.7316     0.6986     5.3416

 Layer #3
 REDUCED STIFFNESS MATRIX IN GPa:
   35.2120     4.6931    -6.7316
    4.6931     9.4731    -0.6986
   -6.7316    -0.6986     5.3416

 Layer #4
 REDUCED STIFFNESS MATRIX IN GPa:
   26.2896     8.1763     9.4512
    8.1763    11.4292     3.4183
    9.4512     3.4183     8.8247

 A MATRIX IN 1E6 N/M:
  158.2153    30.4320     0.0000
   30.4320    51.2776     0.0000
    0.0000     0.0000    33.6741

 B MATRIX IN 1E3 N:
    0.0000     0.0000    22.6588
    0.0000     0.0000    12.1012
   22.6588    12.1012     0.0000

 D MATRIX IN N.DM:
  293.9255    77.3325     0.0000
   77.3325   114.6529     0.0000
    0.0000     0.0000    84.0869

 STIFFNESS R MATRIX:
 1.5822e+08  3.0432e+07  0.00000000  0.00000000  0.00000000  2.2659e+04
 3.0432e+07  5.1278e+07  0.00000000  0.00000000  0.00000000  1.2101e+04
 0.00000000  0.00000000  3.3674e+07  2.2659e+04  1.2101e+04  0.00000000
 0.00000000  0.00000000  2.2659e+04  293.925540  77.3325243  0.00000000
 0.00000000  0.00000000  1.2101e+04  77.3325243  114.652882  0.00000000
 2.2659e+04  1.2101e+04  0.00000000  0.00000000  0.00000000  84.0868611
```

**Real-World Problem Solving Ideas:**
*   **Composite Structural Design**: Crucial for designing aerospace, automotive, or marine composite components (e.g., fuselage panels, car hoods, boat hulls) by allowing engineers to predict their overall stiffness and response to external loads.
*   **Finite Element Analysis (FEA) Pre-processing**: The calculated A, B, D matrices are often required inputs for commercial FEA software when modeling composite shell or plate elements. This program can serve as a validation or pre-computation tool.
*   **Anisotropic Material Characterization**: Aids in understanding and quantifying the highly anisotropic (direction-dependent) behavior of composite laminates.
*   **Laminate Optimization**: By varying ply angles and stacking sequences, engineers can use this program to optimize laminates for specific stiffness requirements (e.g., high bending strength, minimal hygrothermal warping, or specific coupling effects).

**Dependencies**: `WinCrt`.

**Reference:**
*   "MATERIAUX COMPOSITES - Comportement mécanique et analyse des structures" By J.-M. Berthelot - MASSON 1996.

### 8. `compos04.pas` - Deformations and Stresses in Laminated Composite Material

**Program Description:**
This program determines the resulting deformations (mid-plane strains and curvatures) and the detailed ply-by-ply stress distributions (in both global and local fiber-aligned axes) for a laminated composite material. The analysis is performed based on known resultant imposed efforts (forces per unit length and moments per unit length). The output is saved to `compos04.lst`.

**Core Functionality:**
The program takes as input the imposed resultant forces (`Nx`, `Ny`, `Nxy`) and moments (`Mx`, `My`, `Mxy`), along with the laminate's stacking sequence (number of layers, angles, and thicknesses) and the basic material properties.
It first calculates the full 6x6 laminate stiffness matrix (similar to `compos03.pas`). Then, it inverts this matrix to obtain the laminate compliance matrix. By multiplying the compliance matrix with the imposed effort vector, it calculates the laminate's mid-plane strains and curvatures (`E0` vector). Finally, for each layer and at different `z` positions (top and bottom surfaces of the layer), it computes:
1.  Deformations in the global (x,y) axes.
2.  Deformations in the main material axes (L,W).
3.  Stresses in the global (x,y) axes.
4.  Stresses in the main material axes (L,W).
All results are output to a file named `compos04.lst`.

**Key Variables:**
*   `NX`, `NY`, `NXY`: Imposed resultant forces per unit length (N).
*   `MX`, `MY`, `MXY`: Imposed resultant moments per unit length (N.m).
*   `EL`, `ET`, `NULT`, `GLT`: Material engineering constants (GPa).
*   `NCOUCHES`: Number of layers.
*   `TH[I]`: Fiber angle of layer `I` (degrees, converted to radians).
*   `EP[I]`: Thickness of layer `I` (m).
*   `R`: The full 6x6 laminate stiffness matrix.
*   `DET`: Determinant of `R` matrix.
*   `N`: Input effort vector (`NX`, `NY`, `NXY`, `MX`, `MY`, `MXY`).
*   `E0`: Vector of mid-plane strains and curvatures.
*   `EDEB`, `EFIN`: Deformation/Stress vectors at top and bottom surfaces of a layer.
*   `sigma L`, `sigma T`, `sigma LT`: Stresses in main material axes (MPa).
*   `fp`: File pointer to `compos04.lst`.

**Mathematical Principles:**
This program also relies on Classical Laminate Theory (CLT). It extends the calculation of the laminate stiffness matrix (ABD matrix) by inverting it to obtain the compliance matrix. The compliance matrix directly relates the resultant forces and moments to the laminate's overall strains and curvatures. Once these global deformations are known, the individual strain and stress components within each ply are calculated using through-thickness strain variation (`epsilon(z) = epsilon0 + z * kappa`) and ply constitutive relations, accounting for ply orientation. The `MATINV` procedure implements the Gauss-Jordan method for matrix inversion.

**Sample Usage (from source comments):**
```
   The output file Compos04.lst contains (extract):

  CALCULATE DEFORMATIONS AND STRESSES
       IN A LAMINATED MATERIAL

  Imposed resulting efforts:
     NX =   1.0000e+06  NY =   5.0000e+05  NXY =   2.5000e+05
     MX =   0.00000000  MY =   0.00000000  MXY =   0.00000000

  Number of layers: 4

   Layer 1: angle   15.0 deg. - thickness    1.5 mm

   Layer 2: angle  -30.0 deg. - thickness    1.0 mm

   Layer 3: angle  -15.0 deg. - thickness    1.5 mm

   Layer 4: angle   30.0 deg. - thickness    1.0 mm


  Material Parameters:
     El =   38.0    Et =    9.0
   NUlt =   0.32    Glt=    3.6


  BASIC REDUCED STIFFNESS MATRIX IN GPa:
    38.9445     2.9516     0.0000
     2.9516     9.2237     0.0000
     0.0000     0.0000     3.6000

  Layer #:1

  REDUCED STIFFNESS MATRIX IN GPa:
    35.2120     4.6931     6.7316
     4.6931     9.4731     0.6986
     6.7316     0.6986     5.3416

  Layer #:2

  REDUCED STIFFNESS MATRIX IN GPa:
    26.2896     8.1763    -9.4512
     8.1763    11.4292    -3.4183
    -9.4512    -3.4183     8.8247

  Layer #:3

  REDUCED STIFFNESS MATRIX IN GPa:
    35.2120     4.6931    -6.7316
     4.6931     9.4731    -0.6986
    -6.7316    -0.6986     5.3416
  --/--
  Stresses in each layer in main axes:

  Layer #1
   For z =     -2.50000 mm:
    sigma L :    282.73562  MPa
    sigma T :     78.07995  MPa
    sigma LT:     45.15936  MPa
   For z =     -1.00000 mm:
    sigma L :    288.58230  MPa
    sigma T :     70.53304  MPa
    sigma LT:     34.70107  MPa

  Layer #2
   For z =     -1.00000 mm:
    sigma L :     86.43422  MPa
    sigma T :    105.75930  MPa
    sigma LT:      5.73649  MPa
   For z =      0.00000 mm:
    sigma L :    111.92641  MPa
    sigma T :     96.96499  MPa
    sigma LT:      8.38896  MPa

  Layer #3
   For z =      0.00000 mm:
    sigma L :    151.46586  MPa
    sigma T :     90.07486  MPa
    sigma LT:     21.12949  MPa
   For z =      1.50000 mm:
    sigma L :    192.64459  MPa
    sigma T :     76.37100  MPa
    sigma LT:     19.34600  MPa

  Layer #4
   For z =      1.50000 mm:
    sigma L :    333.21160  MPa
    sigma T :     51.87584  MPa
    sigma LT:      8.77294  MPa
   For z =      2.50000 mm:
    sigma L :    317.90586  MPa
    sigma T :     50.19097  MPa
    sigma LT:      1.40860  MPa
```

**Real-World Problem Solving Ideas:**
*   **Composite Failure Analysis**: The most critical application is identifying which plies are most highly stressed and in which material direction, allowing engineers to predict potential failure (e.g., fiber breakage, matrix cracking) by comparing ply stresses to material strength limits.
*   **Delamination Risk Assessment**: High interlaminar stresses (especially shear stresses) predicted at ply interfaces can indicate a risk of delamination, a common and critical failure mode in composite laminates.
*   **Design Validation and Refinement**: Verifying that a composite laminate design can safely withstand expected operational loads and refining stacking sequences for improved performance or weight reduction.
*   **Post-Processing of FEA Results**: Can be used to verify detailed stress results obtained from more complex Finite Element Analyses.

**Dependencies**: `WinCrt`.

**Reference:**
*   "MATERIAUX COMPOSITES - Comportement mécanique et analyse des structures" By J.-M. Berthelot - MASSON 1996.

### 9. `ef3d/ef3d.pas`: EF3D Finite Elements Demonstration Program

**Program Description:**
EF3D (Finite Elements 3D) is a demonstrative program that calculates the eigenfrequencies and eigenmodes (natural frequencies and mode shapes) of three-dimensional structures. It supports various structural elements including 3D beams, axial bars, and torsion bars. The program incorporates features like imposed degrees of freedom (fixed supports), the addition of local (point) masses, and discrete spring elements. It utilizes advanced linear algebra techniques for solving the underlying generalized eigenvalue problem.

**Core Functionality:**
The program allows for the definition of structural models composed of various 3D elements: beam, bar, and torsion bar elements, along with the ability to add local masses and springs. It takes nodal coordinates, element properties (material, geometry), and boundary conditions (fixed/free degrees of freedom) as input.
It then assembles the global stiffness (`K`) and mass (`M`) matrices for the entire structure. Using these matrices, it solves the generalized eigenvalue problem `[K]X = \lambda[M]X` to determine the natural frequencies (`omega2` or `lambda`) and corresponding mode shapes (`Phi`).
The program offers several post-processing options:
*   Displaying general model parameters.
*   Presenting a table of modes with frequencies and generalized masses.
*   Displaying the numerical shape of modes (nodal displacements/rotations).
*   Graphical representation of a selected mode's shape (requires a graphical environment).
*   Performing an orthogonality test on the computed mode shapes.
*   Saving and loading model data and results.

**Key Procedures/Units:**
*   `Open_Input_File`: Reads model data from a `.DAT` file (nodes, coordinates, boundary conditions, elements, local masses, springs).
*   `Points_Nodaux`: Processes nodal data.
*   `Elements`: Processes element data (type, material properties, connectivity).
*   `Partitionner`: Organizes degrees of freedom into free and prescribed.
*   `Assemble_K` / `Assemble_M`: Assembles global stiffness and mass matrices from element contributions, applying coordinate transformations.
*   `Add_Local_Masses` / `Add_Springs`: Adds point masses and springs to specific degrees of freedom.
*   `Frequences_et_Masses_Generalisees`: Calculates natural frequencies and generalized masses from eigenvalues and eigenvectors.
*   `Table_of_Modes`: Displays a summary table of modes, frequencies, and generalized masses.
*   `Shape_of_Modes`: Displays numerical values of mode shapes for each node and DOF.
*   `Draw_Mode`: Graphically displays a selected mode shape for a specific DOF.
*   `Test_Orthogonality`: Verifies the orthogonality of mode shapes.

**Mathematical Principles:**
The program implements the Finite Element Method (FEM).
1.  **Discretization**: The continuous structure is discretized into finite elements (beams, bars).
2.  **Element Matrix Formulation**: For each element, local stiffness (`[k_e]`) and mass (`[m_e]`) matrices are formulated based on the element's geometry, material properties, and assumed displacement functions.
3.  **Coordinate Transformation**: Local element matrices are transformed from local to global coordinates using rotation matrices (`Change_Base`).
4.  **Assembly**: Individual element matrices are assembled into global stiffness (`[K]`) and mass (`[M]`) matrices for the entire structure.
5.  **Boundary Conditions**: Imposed degrees of freedom are handled during assembly, reducing the system size.
6.  **Generalized Eigenvalue Problem**: The free vibration problem is expressed as `([K] - w^2 [M]) {phi} = {0}`. The program solves this using the `Reduc1` procedure from `ef3d/algeblin.pas` to transform it into a standard symmetric eigenvalue problem `[P]Z = w^2Z`. This transformation involves Cholesky decomposition of the mass matrix (`M`).
7.  **Eigenvalue/Eigenvector Solution**: `Tred2` and `Tql2` (from `ef3d/algeblin.pas`) are used. `Tred2` performs Householder reduction to tridiagonalize the matrix, and `Tql2` then finds eigenvalues (`omega2` - natural frequencies squared) and eigenvectors (`Phi` - mode shapes) of the transformed system.
8.  **Back-Transformation**: `Rebaka` (from `ef3d/algeblin.pas`) transforms the eigenvectors back to the original problem's coordinate system.
9.  **Normalization**: Eigenvectors are normalized to unity based on the largest translation or rotation component to simplify interpretation.
10. **Generalized Mass/Stiffness**: Calculated from the mode shapes and mass/stiffness matrices.

**Sample Usage (from source comments):**
The program's main menu provides options for interacting with a loaded problem.
```
EF3D FINITE ELEMENTS PROGRAM
   Release TPW 3.0, May 2006
(C)Copyright J.-P. Dumont 1987,1991
             J.-P. Moreau 2006

Is it a new problem (y/n)? y
Name of data file ? mybeam.dat (user input, example data file not provided)

(Program reads data, assembles matrices, solves eigenvalue problem)

Assembling Ok.
Do you want to save the results (y/n)? y

MENU OF RESULTS
  General Parameters: 1
  Table of modes    : 2
  Shape of modes    : 3
  Draw a mode       : 4
  Orthogonality Test: 5
  Exit              : 6
  Other Problem     : 7

Option (1 to 7): 2

Mode  Frequency  Generalized Mass
----  ---------  ----------------
   1  26.025     5.4000e+06
   2  71.740     1.6438e+09
   3  140.639    6.3022e+11
...

Option (1 to 7): 3
Input number of desired mode (99 to exit)
Mode number? 1

Shape of Mode 1
------------
Mode # 1   Eigenvalue=   1.63522238466632E+0002
           Frequency=         26.02537330 Hz
           Generalized Mass= 5.40002743637149E+0006
Eigenvector:
------------
# Node u         v         w     theta x   theta y   theta z
  1 0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
  ... (nodal displacements and rotations)
```

**Real-World Problem Solving Ideas:**
*   **Structural Dynamics**: Analyzing the dynamic response of complex 3D frame structures, such as building frames, offshore platforms, or machine structures, to understand their inherent vibration characteristics.
*   **Aerospace Structures**: Predicting the natural frequencies of aircraft components (e.g., fuselage, wing spars) to avoid resonance with engine vibrations or aerodynamic loads.
*   **Robotics and Automation**: Designing and optimizing robotic manipulators and other automated machinery by ensuring their natural frequencies are well-separated from operating frequencies, minimizing vibrations and improving precision.
*   **Fatigue Analysis**: Identifying critical vibration modes and high-stress regions that could lead to fatigue failure in dynamically loaded components.
*   **Educational Tool**: Serves as an excellent example for students and researchers to understand the practical implementation of FEM algorithms, particularly for dynamic analysis and eigenvalue problems.

**Dependencies**: `WinCrt1`, `WinTypes`, `WinProcs`, `Strings`, `Errors`, `ErrorHan`, `Declara1`, `AlgebLin`, `Utilit`, `Graph_2D`. (Note: Some of these units are located in the `ef3d` subfolder).

### 10. `eqdifp.pas`: Unit for Solving Differential Equations (Runge-Kutta)

**Unit Description:**
This is a Pascal unit providing procedures for solving systems of first-order ordinary differential equations (ODEs) using the fourth-order Runge-Kutta method. It is designed to be integrated into other programs that require numerical solutions to differential equations.

**Core Functionality:**
The unit defines a `Equadiffp` procedure that takes a function pointer (`fp`) representing the system of differential equations (`yi' = f(y1, y2, ..., yn)`), initial conditions (`Yi`), and simulation parameters (time range `xi`, `xf`, number of points `m`, number of variables `p`, and `finesse` for intermediate calculation steps). It then calculates and stores the solutions for the first two variables (`t1`, `t2`) over the specified time range.
It also includes a `Display` procedure to print the results to the console and optionally save them to external `.asc` files (`theta.asc`, `theta1.asc`), which can be used for plotting.

**Key Procedures/Functions:**
*   `fp(k:INTEGER;X:real_ar;Y:Table):real_ar`: A function type (`Fp1`) that defines the right-hand side of the system of differential equations. `k` refers to the index of the equation, `X` is the independent variable (e.g., time), and `Y` is an array of dependent variables (`Table`).
*   `Equadiffp`: Implements the fourth-order Runge-Kutta algorithm to numerically integrate the system of ODEs.
*   `Display`: Outputs the calculated results (time, y1, y2) to the console and files.

**Key Variables (within `eqdifp.pas` procedures):**
*   `h`: Step size.
*   `x`: Independent variable.
*   `t1, t2`: Output dynamic arrays (`RV` type) for the first two dependent variables.
*   `xi, xf`: Initial and final values of the independent variable.
*   `Yi`: Initial conditions array (`Table` type).
*   `m`: Number of desired output points.
*   `p`: Number of independent variables (equations).
*   `fi`: Finesse (number of intermediate steps).
*   `ta, tb, tc, td`: Runge-Kutta intermediate coefficient arrays.
*   `y, z`: Temporary arrays for dependent variables.
*   `fp1, fp2`: File pointers for `theta.asc` and `theta1.asc`.

**Mathematical Principles:**
The unit employs the classical fourth-order Runge-Kutta method. This is a widely used iterative numerical technique for approximating the solution of initial value problems for ordinary differential equations. It is known for its balance of accuracy and computational efficiency.

**Real-World Problem Solving Ideas:**
*   **Any System Modeled by ODEs**: This unit is a fundamental building block for a wide range of engineering and scientific problems. Examples include:
    *   **Chemical Reaction Kinetics**: Simulate the concentration changes of reactants and products over time.
    *   **Population Dynamics**: Model changes in population sizes based on birth/death rates.
    *   **Thermodynamics**: Analyze heat transfer or transient temperature changes in systems.
    *   **Pharmacokinetics**: Model drug concentration in the body over time.
    *   **Financial Modeling**: Simulate changes in financial variables based on differential equations.
*   **Custom Dynamic System Simulations**: Engineers can write custom `fp` functions to represent their specific physical systems (e.g., an electrical circuit, a mechanical vibrator, a control system) and use `Equadiffp` to analyze their time-domain behavior.

**Dependencies**: `WinCrt`, `Type_def`. (Note: `Type_def` is provided in `ef3d/type_def.pas` but is also accessible if the full library structure is maintained).

**Reference:**
*   "Analyse en Turbo Pascal versions 5.5 et 6.0" By Marc DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991 [BIBLI 03].

### 11. `modes.pas`: Frequencies and Modes of Mass-Spring Systems (Transfer Method)

**Program Description:**
This program identifies the natural (resonance) frequencies and corresponding mode shapes (eigenvectors) of undamped multi-degree-of-freedom (MDOF) mass-spring systems. It achieves this by employing the Transfer Matrix Method. The program supports various boundary conditions at the system's ends.

**Core Functionality:**
The user defines a 1D system by specifying the number and type of elements (masses or springs) in sequence. Masses and spring stiffnesses are provided. The program supports four common boundary conditions: Fixed-Fixed, Fixed-Free, Free-Fixed, and Free-Free.
It can operate in two primary modes:
1.  **Frequency Sweep**: Computes the determinant of the overall transfer matrix (which is related to the characteristic equation) over a specified frequency range. Resonance frequencies correspond to frequencies where this determinant approaches zero or changes sign.
2.  **Mode Calculation**: For a given resonance frequency, it calculates the corresponding mode shape (displacement and force at each node).

**Key Variables:**
*   `N`: Number of elements.
*   `C[I]`: Type of element `I` (1: Spring, 2: Mass).
*   `M1[I]`: Mass of element `I`.
*   `K1[I]`: Stiffness of element `I`.
*   `C8`: Integer code for system boundary conditions (1-4).
*   `F`: Current frequency (Hz).
*   `W`: Pulsation (`2 * PI * F`).
*   `W2`: `W * W`.
*   `T`, `T1`, `A`: 2x2 matrices used in the transfer method.
*   `X`: Determinant value (varies depending on boundary conditions).
*   `V1`, `V2`: 2-element vectors representing state variables (Force and Displacement or vice-versa) at each node for mode calculation.

**Mathematical Principles:**
The program utilizes the Transfer Matrix Method.
1.  **Element Transfer Matrix**: For each element (mass or spring), a 2x2 transfer matrix `[A]` is defined that relates the state vector `[Force, Displacement]^T` at one end of the element to the state vector at the other end.
    *   For a mass `M`: `A = [[1, -M*W^2], [0, 1]]` (relates displacement and force before and after the mass)
    *   For a spring `K`: `A = [[1, 0], [1/K, 1]]` (relates displacement and force before and after the spring)
2.  **Global Transfer Matrix**: These individual transfer matrices are then multiplied sequentially to obtain a global transfer matrix `[T]` for the entire system: `[T] = [A_N] * ... * [A_2] * [A_1]`.
3.  **Characteristic Equation**: The natural frequencies are found by setting a specific element of the global transfer matrix (determined by the boundary conditions) to zero, forming the characteristic equation (`X` in the output).
    *   Fixed-Fixed (Code=1): `T[2,1]` (Displacement-Force element of `T` matrix).
    *   Fixed-Free (Code=2): `T[1,1]` (Force-Force element of `T` matrix).
    *   Free-Fixed (Code=3): `T[2,2]` (Displacement-Displacement element of `T` matrix).
    *   Free-Free (Code=4): `T[1,2]` (Force-Displacement element of `T` matrix).
4.  **Mode Shape Calculation**: Once a natural frequency is identified, the corresponding mode shape (the relative displacements of the masses) is determined by propagating an initial state vector (e.g., assuming unit displacement or unit force at the first node, depending on boundary conditions) through the system using the transfer matrices at that specific natural frequency.

**Sample Usage (from source comments):**
This example models a 2-DOF system: `k - 3m - k - m - k` with fixed-fixed ends.
```
 Frequencies and Modes of Mass-spring Systems
 By Transfer Method.

 Kind of elements
 ----------------
 C(I)=1: Spring
 C(I)=2: Mass

 Number of Elements: 5

 Kind of element 1 C(1) = 1
 Kind of element 2 C(2) = 2
 Kind of element 3 C(3) = 1
 Kind of element 4 C(4) = 2
 Kind of element 5 C(1) = 1

 Fixed-Fixed System, code=1.
 Fixed-Free System,  code=2.
 Free-Fixed System,  code=3.
 Free-Free System,   code=4.

 Code = 1

 Mass #1 = 3
 Mass #2 = 1

 Spring #1 = 1
 Spring #2 = 1
 Spring #3 = 1

 Frequency Sweep: (S)WEEP
 Calculate Mode.: (M)ODE
 Exit Program...: (E)XIT

 Answer = S

 Frequency Sweep
 ---------------
 Starting Frequency: 0
 Ending Frequency..: 1
 Frequency Step....: 0.1

 Frequency= 0.00000000000000E+0000 Hz  Determinant= 3.00000000000000E+0000
 Frequency= 1.00000000000000E-0001 Hz  Determinant= 3.09290228614616E-0001
 Frequency= 2.00000000000000E-0001 Hz  Determinant=-2.15207544198299E+0000
 ... (indicates a resonance near 0.1 Hz)

 Note: A minimum value of matrix determinant (here T[2,1]) occurs near a
       resonance frequency, here f=0.1 hz. Refine sweep to improve solution.

 Answer = S (Refined Sweep)
 Starting Frequency: 0.10
 Ending Frequency..: 0.12
 Frequency Step....: 0.0005
 ... (further refinement shows a minimum near 0.107 Hz)

 Answer = M (for Mode)
 Resonance Frequency: 0.107

 Node #1 Force= 1.00000000000000E+0000 Displacement= 0.00000000000000E+0000
 Node #2 Force= 1.00000000000000E+0000 Displacement= 1.00000000000000E+0000
 Node #3 Force=-3.55965209456865E-0001 Displacement= 1.00000000000000E+0000
 Node #4 Force=-3.55965209456865E-0001 Displacement= 6.44034790543135E-0001
 Node #5 Force=-6.47061466008975E-0001 Displacement= 6.44034790543135E-0001
 Node #6 Force=-6.47061466008975E-0001 Displacement=-3.02667546583979E-0003
```

**Real-World Problem Solving Ideas:**
*   **Shaft Dynamics**: Analyze torsional or lateral vibrations in multi-disk shafts, gear train systems, or power transmission systems to identify critical speeds and prevent resonance.
*   **Train Dynamics**: Model simplified train car systems connected by couplers (springs) and having distributed masses, to study their longitudinal vibrations and stability.
*   **Conveyor Systems**: Understand the dynamic characteristics of conveyor belts and rollers, which can influence material transport efficiency and component longevity.
*   **Conceptual Design of MDOF Systems**: Provides an analytical tool for initial conceptual design of systems where complex FEA models are not yet necessary.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Mécanique des vibrations linéaires" By M. Lalanne, P. Berthier, J. Der Hagopian, Masson, Paris 1980 [BIBLI 16].

### 12. `ndof01.pas`: Response of N-DOF Mass-Spring System with Damping (Transfer Matrices Method)

**Program Description:**
This program calculates the frequency-dependent response (displacement amplitude and phase) of a multi-degree-of-freedom (N-DOF) mass-spring system, including various types of damping, subjected to a sinusoidal force. It employs the Transfer Matrix Method to analyze the system.

**Core Functionality:**
The user defines the system by specifying the number and type of elements (Spring, Mass, Viscous Damper, Spring+Viscous Damper in parallel, Spring with Structural Damping, Sinusoidal Force). Material properties (masses, spring stiffnesses, damping coefficients/factors) and sinusoidal force characteristics (amplitude, phase) are provided. The program supports four common boundary conditions (Fixed-Fixed, Fixed-Free, Free-Fixed, Free-Free).
The primary output is the frequency response of a specific node (displacement amplitude and phase) over a defined frequency sweep.

**Key Variables:**
*   `N`: Number of elements.
*   `C[I]`: Type of element `I` (1-6).
*   `M1[I]`: Mass.
*   `K1[I]`: Spring stiffness.
*   `C1[I]`: Viscous Damping Coefficient.
*   `K2[I]`, `C2[I]`: Spring stiffness and viscous damping coefficient for parallel element.
*   `K3[I]`, `E3[I]`: Spring stiffness and structural damping factor.
*   `F1[I]`, `F2[I]`: Real and Imaginary components of the sinusoidal force.
*   `N1`: Node number for which response is desired.
*   `C8`: Integer code for system boundary conditions (1-4).
*   `F`: Current frequency (Hz).
*   `W`: Pulsation (`2 * PI * F`).
*   `W2`: `W * W`.
*   `T`, `T1`, `A`, `T2`: 5x5 matrices used in the complex transfer method.
*   `V1`, `V2`: 5-element vectors used for state propagation.
*   `M2`: Displacement amplitude at node `N1`.
*   `T3`: Phase (in degrees) at node `N1`.

**Mathematical Principles:**
This program extends the Transfer Matrix Method (as seen in `modes.pas`) to include damping and external forces.
1.  **Generalized Transfer Matrix**: For each element (spring, mass, damper, or force application point), a 5x5 complex transfer matrix `[A]` is constructed. This matrix relates the complex state vector (typically including displacement, force, complex conjugate displacement, complex conjugate force, and external force components) at one end of the element to the state vector at the other end.
2.  **Global System Matrix**: These individual transfer matrices are multiplied in sequence (`[T] = [A_N] * ... * [A_1]`) to obtain a global transfer matrix for the entire system.
3.  **Boundary Conditions and Excitation**: The program effectively solves a system of complex linear equations derived from the global transfer matrix, applying the boundary conditions (e.g., fixed ends mean zero displacement, free ends mean zero force) and the external sinusoidal force components. This often involves converting the complex system into a real-valued system of twice the original size (e.g., 2N x 2N) before solving.
4.  **Response Calculation**: After solving for the system's overall response, the program extracts the displacement amplitude and phase for the specified node (`N1`) by analyzing the real and imaginary components of its displacement.

**Sample Usage (from source comments):**
This example models a 3-DOF system: `k - 2m - 2k - m - k - 3m` with a force on the first mass, and optionally structural damping in the last spring.
```
 Kind of elements
 ----------------
  1: Spring
  2: Mass
  3: Viscous Damper
  4: Spring + Viscous Damper in parallel
  5: Spring with structural Damping
  6: Sinusoidal Force

 Number of elements: 7
 Kind of element 1: 1
 Kind of element 2: 2
 Kind of element 3: 1
 Kind of element 4: 2
 Kind of element 5: 1
 Kind of element 6: 2
 Kind of element 7: 6

 Mass #1 = 2
 Mass #2 = 1
 Mass #3 = 3

 Spring #1 = 1
 Spring #2 = 2
 Spring #3 = 1

 Excitation Force #1 F.COS(PHI) = 1000
 Excitation Force #1 F.SIN(PHI) = 0

 For which node number do you want the response: 3

 Fixed-Fixed System, code=1.
 Fixed-Free System,  code=2.
 Free-Fixed System,  code=3.
 Free-Free System,   code=4.

 Code = 2

 Frequency Sweep
 ---------------
 Starting Frequency: 0.30
 Ending Frequency..: 0.32
 Frequency Step....: 0.001

 Freq=  0.300  Displacement=     96.1  Phase=  0.0 Deg.
 Freq=  0.301  Displacement=    101.2  Phase=  0.0 Deg.
 Freq=  0.302  Displacement=    107.1  Phase=  0.0 Deg.
 Freq=  0.303  Displacement=    114.0  Phase=  0.0 Deg.
 Freq=  0.304  Displacement=    122.2  Phase=  0.0 Deg.
 Freq=  0.305  Displacement=    132.1  Phase=  0.0 Deg.
 Freq=  0.306  Displacement=    144.1  Phase=  0.0 Deg.
 Freq=  0.307  Displacement=    159.3  Phase=  0.0 Deg.
 Freq=  0.308  Displacement=    178.8  Phase=  0.0 Deg.
 Freq=  0.309  Displacement=    204.7  Phase=  0.0 Deg.
 Freq=  0.310  Displacement=    240.9  Phase=  0.0 Deg.
 Freq=  0.311  Displacement=    294.8  Phase=  0.0 Deg.
 Freq=  0.312  Displacement=    383.9  Phase=  0.0 Deg.
 Freq=  0.313  Displacement=    558.2  Phase=  0.0 Deg.
 Freq=  0.314  Displacement=   1051.8  Phase=  0.0 Deg.
 Freq=  0.315  Displacement=  12214.5  Phase=  0.0 Deg. <-- 3rd resonance
 Freq=  0.316  Displacement=   1226.3  Phase=180.0 Deg.
 Freq=  0.317  Displacement=    574.1  Phase=180.0 Deg.
 Freq=  0.318  Displacement=    370.7  Phase=180.0 Deg.
 Freq=  0.319  Displacement=    271.5  Phase=180.0 Deg.
 Freq=  0.320  Displacement=    212.8  Phase=180.0 Deg.

 Example #2 (with structural damping):
  Kind of element 5: 5 (Spring with structural Damping)
  Spring with structural Damping #1  K = 1
  Spring with structural Damping #1  E = 0.05

 Freq=  0.300  Displacement=     95.7  Phase=  3.5 Deg.
 ...
 Freq=  0.315  Displacement=    644.2  Phase= 84.1 Deg. <-- 3rd resonance
 ...
```

**Real-World Problem Solving Ideas:**
*   **Dynamic Analysis of Machine Foundations**: Predict the steady-state vibration of heavy machinery foundations, which can be modeled as N-DOF systems, under harmonic excitation from rotating or reciprocating components.
*   **Bridge Response to Traffic**: Analyze the dynamic amplification of bridge structures under repetitive vehicle loads, where bridges can be simplified to N-DOF systems.
*   **Vibration Absorber/Damper Effectiveness**: Evaluate how different damping strategies (viscous, structural, tuned dampers) affect the overall dynamic response of complex systems at various frequencies.
*   **Aerospace Component Vibration**: Predict the forced response of aircraft or spacecraft components to engine vibrations or aerodynamic forces, particularly important for fatigue life.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Mécanique des vibrations linéaires" By M. Lalanne, P. Berthier, J. Der Hagopian, Masson, Paris 1980 [BIBLI 16].

### 13. `ndof02.pas`: Frequencies, Eigenmodes, Modal Masses & Stiffnesses of Undamped MDOF System

**Program Description:**
This program calculates the natural frequencies (eigenvalues), corresponding eigenmodes (eigenvectors), and associated modal masses and modal stiffnesses for an undamped multi-degree-of-freedom (MDOF) mass-spring system. The system's dynamics are described by the generalized eigenvalue problem `[M] . [X]" + [K] . [X] = [0]`. It uses an iterative power method to solve for these properties.

**Core Functionality:**
The user inputs the order of the system (`N`), the symmetric mass matrix `[M]`, and the symmetric stiffness matrix `[K]`. It then asks for the number of eigenvectors to be calculated, the precision for eigenvalues, and the maximum number of iterations.
The program uses an iterative method to find the eigenvalues (natural frequencies squared) and corresponding eigenvectors (mode shapes). It also calculates and displays the modal mass matrix and modal stiffness matrix.

**Key Variables:**
*   `N`: Order of the system (number of degrees of freedom).
*   `M`, `K`: Mass and Stiffness matrices (using pointers `pMAT`).
*   `A`, `D`, `K1`, `S`, `S1`, `T`, `T1`, `V`: Matrices used in the iterative solution.
*   `L`: Eigenvalues (natural frequencies squared).
*   `X`: Current eigenvector iteration vector.
*   `F`, `F1`: Temporary vectors.
*   `W2`: `w^2` during iterations.
*   `P`: Precision for eigenvalue convergence (`P` squared internally).
*   `N1`: Number of eigenvectors asked for.
*   `N3`: Maximum number of iterations.
*   `DET`: Determinant (from `MATINV`).

**Mathematical Principles:**
The program solves the generalized eigenvalue problem `[K] {phi} = w^2 [M] {phi}` which arises from the undamped free vibration equation `[M] . [X]" + [K] . [X] = [0]`.
1.  **Transformation**: The problem is transformed by adding a shift `A0` and inverting `[K1] = [K + A0*M]` to get `[K1]^-1 [M] X = (1/(w^2 + A0)) X`. This converts it to a standard eigenvalue problem `[D] X = \lambda_shifted X` where `D = [K1]^-1 [M]`.
2.  **Iterative Power Method**: The program then uses an iterative power method (specifically, inverse iteration) on the matrix `[D]`. When applied to `[D]`, the power method converges to the eigenvector corresponding to the dominant eigenvalue of `[D]`, which corresponds to the smallest `(w^2 + A0)` value, thus the smallest `w^2` (lowest natural frequency). The process is repeated for subsequent modes by deflation (removing the contribution of already found modes).
3.  **Matrix Operations**: `MATINV` procedure (Gauss-Jordan method with full pivoting) is implemented directly within the program to compute matrix inverses. `MATMUL` performs matrix multiplication.
4.  **Modal Properties**: After obtaining eigenvalues `L = w^2` and eigenvectors `V = phi`, the modal mass matrix `[M_modal] = {phi}^T [M] {phi}` and modal stiffness matrix `[K_modal] = {phi}^T [K] {phi}` are calculated. For an undamped system, these matrices are ideally diagonal.

**Sample Usage (from source comments):**
This example models a 3-DOF system: `k - 2m - 2k - m - k - 3m`.
```
 Order of system: 3

 M(1,1) = 2
 M(1,2) = 0
 M(1,3) = 0
 M(2,2) = 1
 M(2,3) = 0
 M(3,3) = 3

 K(1,1) = 3
 K(1,2) = -2
 K(1,3) = 0
 K(2,2) = 3
 K(2,3) = -1
 K(3,3) = 1

 Number of eigenvectors asked for: 3

 What precision for eigenvalues: 1e-5

 Maximum number of iterations: 100


 Eigenvector #1, Convergence after 10 iterations.
 Eigenvector #2, Convergence after 9 iterations.
 Eigenvector #3, Convergence after 2 iterations.


 Eigenvalues:

 L(1) =  1.05173321566663E-0001
 L(2) =  8.08612063101068E-0001
 L(3) =  3.91945398307909E+0000

 Pulsations:

 W(1) =  3.24304365629980E-0001
 W(2) =  8.99228593351584E-0001
 W(3) =  1.97976109242481E+0000

 Frequencies:

 F(1) =  5.16146428562926E-0002
 F(2) =  1.43116675601476E-0001
 F(3) =  3.15088764000420E-0001

 Eigenvectors:

 E.V.  1
              1.00000000000000E+0000
              1.39482667843334E+0000
              2.03778199949639E+0000

 E.V.  2
              1.00000000000000E+0000
              0.69138346664719E+0000
             -0.48490130871242E+0000

 E.V.  3
              1.00000000000000E+0000
             -2.41959287715014E+0000
              2.24903464725100E+0000

 Modal Mass Matrix:

  16.40321  -0.00001   0.00000
  -0.00001   3.18340  -0.00003
   0.00000  -0.00003   8.00617

 Modal Stiffness Matrix:

   1.72517   0.00000  -0.00000
   0.00000   2.57413   0.00001
  -0.00000   0.00001  31.38059
```

**Real-World Problem Solving Ideas:**
*   **Structural Dynamics Analysis**: Fundamental for understanding the free vibration characteristics of complex structures (e.g., buildings, bridges, aircraft) modeled as N-DOF systems. This is a crucial step in vibration analysis and dynamic design.
*   **Seismic Analysis**: Predict how a structure will deform and oscillate when subjected to earthquake ground motion. Natural frequencies dictate which ground motion frequencies will excite the structure most.
*   **Rotating Machinery Dynamics**: Analyze critical speeds of rotating shafts where resonance might occur due to imbalance or external excitations.
*   **System Health Monitoring**: Changes in natural frequencies over time can indicate structural degradation or damage. This program provides the foundation for such analysis.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Mécanique des vibrations linéaires" By M. Lalanne, P. Berthier, J. Der Hagopian, Masson, Paris 1980 [BIBLI 16].

### 14. `ndof03.pas`: Response of N-DOF Mass-Spring System with Damping (Direct Method)

**Program Description:**
This program determines the frequency response (displacement amplitude and phase) of a multi-degree of freedom (N-DOF) mass-spring system with damping subjected to a sinusoidal force, using a direct method. It supports both viscous and structural damping models.

**Core Functionality:**
The user inputs the system's order (`N`), the mass matrix `[M]`, stiffness matrix `[K]` (or `K1`), and damping matrix `[C]` (or `K2`, if structural damping). The complex amplitude of the excitation force vector (real and imaginary parts `F1`, `F2`, combined into `B` vector) is also provided.
The program then performs a frequency sweep, calculating the displacement amplitude and phase for a specific degree of freedom over a specified frequency range (linear or logarithmic step).

**Key Variables:**
*   `N`: Number of degrees of freedom.
*   `M`: Mass matrix.
*   `K1`: Stiffness matrix (or real part of complex stiffness for structural damping).
*   `K2`: Damping matrix (or imaginary part of complex stiffness for structural damping).
*   `C9`: Damping type (1: Viscous, 2: Structural).
*   `A`: Augmented `2N x 2N` dynamic stiffness matrix (or impedance matrix).
*   `B`: Augmented `2N`-element force vector (combining real and imaginary parts of applied forces).
*   `X`: Augmented `2N`-element displacement vector (combining real and imaginary parts of displacements).
*   `N3`: Degree of freedom for which displacement and phase are output.
*   `F`: Current frequency (Hz).
*   `W`: Pulsation (`2 * PI * F`).
*   `M1`: Displacement amplitude at `N3`.
*   `T1`: Phase (in degrees) at `N3`.
*   `F1`, `F2`: Starting and ending frequencies for sweep.
*   `F3`: Frequency step value (Hz for LIN, factor for LOG).

**Mathematical Principles:**
The program solves the frequency-domain steady-state response of a damped N-DOF system directly. For a sinusoidal excitation `F(t) = F_0 * exp(i\omega t)` and a response `X(t) = X_0 * exp(i\omega t)`, the equation of motion `[M]X" + [C]X' + [K]X = F(t)` transforms into a system of linear algebraic equations involving complex numbers:
*   For viscous damping: `([-M\omega^2 + iC\omega + K])X_0 = F_0`
*   For structural damping: `([-M\omega^2 + K_1 + iK_2])X_0 = F_0` (where `K_2` is the structural damping matrix).
This complex system is converted into a larger real-valued system of twice the original size (`2N x 2N`) by separating real and imaginary parts. The program constructs this augmented dynamic stiffness matrix `[A]` and the augmented force vector `[B]`. It then solves the linear system `[A]X = B` for the real and imaginary parts of the displacements. The `MATINV` procedure (Gauss-Jordan elimination) is used for matrix inversion. The amplitude and phase of the desired displacement are then calculated from its real and imaginary components.

**Sample Usage (from source comments):**
This example models a 3-DOF system with viscous damping: `k - 2m - 2k - m - k - 3m`, with a viscous damper on the last mass.
```
 How many degrees of freedom (d.o.f.): 3

 Input Mass Matrix [M]:
  M(1,1) = 2
  M(1,2) = 0
  M(1,3) = 0
  M(2,2) = 1
  M(2,3) = 0
  M(3,3) = 3

 Viscous Damping...: VIS
 Structural Damping: STR

 Your choice: VIS

 Input Stiffness Matrix [K]:
  K(1,1) = 3
  K(1,2) = -2
  K(1,3) = 0
  K(2,2) = 3
  K(2,3) = -1
  K(3,3) = 1

 Input Damping Matrix [C]:
  C(1,1) = 0
  C(1,2) = 0
  C(1,3) = 0
  C(2,2) = 0
  C(2,3) = 0
  C(3,3) = 0.5

 Input Excitation Vector:
  F1(1) = 1000
  F1(2) = 0
  F1(3) = 0

  F2(1) = 0
  F2(2) = 0
  F2(3) = 0

 Number of d.o.f. to calculate: 3


 Frequency Sweep
 ---------------
 Starting Frequency: 0.31
 Ending Frequency..: 0.32

 Linear frequency step.....: LIN
 Logarithmic frequency step: LOG

 Your choice: LIN


 Frequency Step....: 0.0005


 Frequency (Hz)   Displacement   Phase (deg.)
 --------------------------------------------
     0.3100          240.66        0.0
     0.3105          264.73        0.0
     0.3110          294.72        0.0
     0.3115          333.08        0.0
     0.3120          383.85        0.0
     0.3125          454.17        0.0
     0.3130          557.90        0.0
     0.3135          725.80        0.0
     0.3140         1041.65        0.0
     0.3145         1824.42        0.0
     0.3150         4347.56        0.0  <-- 3rd resonance
     0.3155         2248.74       -0.0
     0.3160         1151.46       -0.0
     0.3165          757.64       -0.0
     0.3170          560.61       -0.0
     0.3175          443.06       -0.0
     0.3180          365.13       -0.0
     0.3185          309.74       -0.0
     0.3190          268.36       -0.0
     0.3195          236.29       -0.0
     0.3200          210.72       -0.0
```

**Real-World Problem Solving Ideas:**
*   **Forced Vibration Response**: Determine the steady-state response of mechanical systems (e.g., turbine blades, vehicle bodies, machine components) subjected to continuous harmonic loads, taking into account the effects of damping.
*   **Frequency Response Function (FRF) Generation**: Generate points on a Frequency Response Function plot (amplitude vs. frequency, phase vs. frequency) for a specific degree of freedom, which is crucial for understanding system behavior and for experimental modal analysis.
*   **Resonance Characterization**: Observe how damping influences the amplitude at resonance and the phase shift, which is vital for designing effective vibration control strategies.
*   **Isolation System Design**: Evaluate the effectiveness of vibration isolation systems by predicting the transmitted force or displacement across a range of frequencies.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Mécanique des vibrations linéaires" By M. Lalanne, P. Berthier, J. Der Hagopian, Masson, Paris 1980 [BIBLI 16].

### 15. `ndof04.pas`: Step-by-Step Solution of MDOF System (Wilson-Theta Method)

**Program Description:**
This program provides a time-domain, step-by-step numerical solution for the dynamic response of a multi-degree-of-freedom (N-DOF) mass-spring system, including viscous damping, subjected to general time-varying external forces. It utilizes the "Wilson-Theta" direct integration method, which is known for its unconditional stability. The output is saved to `ndof04.txt`.

**Core Functionality:**
The user inputs the system's order (`N`), the symmetric mass matrix `[M]`, damping matrix `[C]`, and stiffness matrix `[K]`. Initial conditions for displacement (`X0`), velocity (`V0`), and acceleration (`G0`) are required. The program also needs the simulation time range (`T0`, `T3`), time increment (`D`), a parameter `Theta` (referred to as `T2` in code, commonly 1.4 for unconditional stability), and the characteristics of the external dynamic force (`Fmax`, `freq`, and `N1` for the force application node).
The program numerically integrates the equations of motion over time, outputting the displacement, velocity, and acceleration for each degree of freedom at each time step to an external file (`ndof4.txt`). It also tracks the maximum displacement for node #1.

**Key Variables:**
*   `N`: Number of degrees of freedom.
*   `M`, `C`, `K`: Mass, Damping, and Stiffness matrices.
*   `X0`, `V0`, `G0`: Initial displacement, velocity, and acceleration vectors at `T0`.
*   `X1`, `V1`, `G1`: Displacement, velocity, and acceleration vectors at `T1`.
*   `T0`, `T3`: Start and end times.
*   `D`: Time increment (`\Delta t`).
*   `T2` (Theta): Wilson-Theta integration parameter (typically 1.4).
*   `Fmax`, `freq`: Amplitude and frequency of the sinusoidal forcing function.
*   `N1`: Node number where the force is applied.
*   `F`: Current force value.
*   `F0`, `F1`, `F2`: Force vectors at different time instances (`t`, `t+dt`, `t+theta*dt`).
*   `K1`: Effective stiffness matrix.
*   `X2`: Calculated incremental displacement vector.
*   `XMAX`: Maximum displacement observed at node #1.
*   `fp`: File pointer to `ndof4.txt`.

**Mathematical Principles:**
The program implements the Wilson-Theta method, an implicit, unconditionally stable, direct integration method for solving transient dynamic problems. This method assumes a linear variation of acceleration over an extended time interval (`\Theta \Delta t`), leading to a set of linear algebraic equations that must be solved at each time step.
The solution involves:
1.  **Effective Stiffness Matrix**: Formation of an effective stiffness matrix `[K_eff] = [K] + a_0[M] + a_1[C]`, where `a_0 = 6/((\Theta \Delta t)^2)` and `a_1 = 3/(\Theta \Delta t)` are integration constants depending on `\Delta t` and `Theta`.
2.  **Effective Force Vector**: Calculation of an effective force vector `[F_eff]` which incorporates the applied external force at `t + \Theta \Delta t` and terms from the previous time step's displacements, velocities, and accelerations.
3.  **Linear System Solution**: Solving `[K_eff] X_{new} = F_{eff}` for the unknown displacement vector `X_{new}` (or for the incremental displacement). The `MATINV` procedure (Gauss-Jordan elimination) is used for matrix inversion.
4.  **Update**: Updating the velocities and accelerations at the current time step using the newly calculated displacements.

**Sample Usage (from source comments):**
This example models a 3-DOF system with viscous damping: `k - 2m - 2k - m - k - 3m`, with a viscous damper on the last mass.
```
 How many degrees of freedom (d.o.f.): 3

 Mass Matrix:
  M(1,1) = 2
  M(1,2) = 0
  M(1,3) = 0
  M(2,2) = 1
  M(2,3) = 0
  M(3,3) = 3

 Damping Matrix:
  C(1,1) = 0
  C(1,2) = 0
  C(1,3) = 0
  C(2,2) = 0
  C(2,3) = 0
  C(3,3) = 0.5

 Stiffness Matrix:
  K(1,1) = 3
  K(1,2) = -2
  K(1,3) = 0
  K(2,2) = 3
  K(2,3) = -1
  K(3,3) = 1

 Starting time (sec.) = 0
 Ending time......... = 60
 Time increment...... = 0.1

 Starting motion X(1) = 0
 Starting motion X(2) = 0
 Starting motion X(3) = 0

 Starting speed V(1) = 0
 Starting speed V(2) = 0
 Starting speed V(3) = 0

 Starting acceleration G(1) = 0
 Starting acceleration G(2) = 0
 Starting acceleration G(3) = 0

 Theta = 1.4

 Number of force componant: 1

 Force (maximum) = 1000
 Force frequency = 0.315

 Output file ndof04.txt contains:

 Time= 0.00000000000000E+0000 Force= 0.00000000000000E+0000
 Node  Displacement  Speed  Acceleration
 ---------------------------------------
 1  0.00000000000000E+0000  0.00000000000000E+0000  0.00000000000000E+0000
 2  0.00000000000000E+0000  0.00000000000000E+0000  0.00000000000000E+0000
 3  0.00000000000000E+0000  0.00000000000000E+0000  0.00000000000000E+0000
 Time= 1.00000000000000E-0001 Force= 1.96630694615420E+0002
 1  1.63063348145476E-0001  4.89190044436429E+0000  9.78380088872859E+0001
 2  1.05501179771710E-0003  3.16503539315131E-0002  6.33007078630263E-0001
 3  1.13432172045768E-0006  3.40296516137303E-0005  6.80593032274606E-0004
 ... (further time steps)
 Time= 6.00000000000006E+0001 Force=-5.87785252291540E+0002
 1 -9.89652416589468E+0002 -4.43026677849652E+0003  3.68939572939423E+0003
 2  2.68713644244187E+0003  9.63603447643417E+0003 -1.02316256295891E+0004
 3 -2.73627842586534E+0002 -7.89759695744624E+0002  1.10877730102250E+0003

 Maximum Displacement (node #1) =  2.43218400402500E+0003
```

**Real-World Problem Solving Ideas:**
*   **Transient Dynamic Analysis**: Simulate the time-history response of structures to transient loads like impact, sudden application/removal of forces, or blast loads. This is crucial for understanding how structures behave under non-harmonic, short-duration events.
*   **Earthquake Engineering**: Analyze the response of buildings and bridges to recorded earthquake ground motion time histories. The Wilson-Theta method is a robust choice for such analyses.
*   **Machine Startup/Shutdown**: Simulate the dynamic behavior of machinery during start-up or shut-down phases where forces might be non-periodic or rapidly changing.
*   **Control System Response**: Predict the precise time-domain behavior of a control system under various input signals.

**Dependencies**: `WinCrt`.

**Reference:**
*   "Mécanique des vibrations linéaires" By M. Lalanne, P. Berthier, J. Der Hagopian, Masson, Paris 1980 [BIBLI 16].

### 16. `pendulum.pas`: Mass Pendulum Angular Motion

**Program Description:**
This program calculates the angular motion (angle `theta` and angular velocity `theta'`) of an elementary mass pendulum over time. It uses the fourth-order Runge-Kutta method by leveraging the `eqdifp.pas` unit.

**Core Functionality:**
The program sets up and solves a system of two first-order differential equations that describe the motion of a damped mass pendulum:
*   `y1' = y2` (where `y1` is `theta` and `y2` is `theta'`)
*   `y2' = -K * sin(y1) - K' * y2` (representing the forces/moments acting on the pendulum, with `K` relating to gravity and pendulum length, and `K'` to fluid resistance).
The user provides initial conditions (`theta` and `theta'` at `t0`), the total simulation time range (`xi`, `xf`), the number of output points (`ndata`), and a `finesse` parameter (number of intermediate steps for Runge-Kutta). Results are displayed on the console and optionally saved to `theta.asc` and `theta1.asc` files for plotting.

**Key Variables:**
*   `xi`, `xf`: Start and end times for the simulation.
*   `yi[0]`: Initial angle `theta` at `t0` (radians).
*   `yi[1]`: Initial angular velocity `theta'` at `t0` (radians/second).
*   `ndata`: Number of output points.
*   `fi` (finesse): Number of intermediate Runge-Kutta steps per output interval.
*   `v1`: Output vector for angle `theta` (using `RV` type).
*   `v2`: Output vector for angular velocity `theta'` (using `RV` type).
*   `p`: Number of variables (fixed at 2).
*   `K`: Constant (`pi*pi/4`) for gravitational restoring force.
*   `K'`: Constant (`0.1`) for viscous damping.

**Mathematical Principles:**
The program utilizes the `Equadiffp` procedure from the `eqdifp.pas` unit, which implements the fourth-order Runge-Kutta method. The second-order non-linear pendulum equation: `theta" = -K sin(theta) - K' theta'` is converted into two coupled first-order equations, which are then numerically integrated.

**Sample Usage (from source comments):**
```
    DIFFERENTIAL EQUATION WITH P VARIABLE OF ORDER 1
        of type yi' = f(y1,y2,...,yn), i=1..n
              (Case of a Mass Pendulum)

  Begin value time    : 0
  End value time      : 12
  theta value at t0   : 0.25
  theta' value at t0  : 0
  Number of points    : 25
  Finesse             : 10

    Time         Theta        Theta'
--------------------------------------
  0.000000     0.250000     0.000000
  0.500000     0.178615    -0.268789
  1.000000     0.009250    -0.372792
  ... (further time steps)
 12.000000     0.136627     0.011313

(Another example for large initial angle)
  Begin value time    : 0
  End value time      : 40
  theta value at t0   : 3.14
  theta' value at t0  : 0
  Number of points    : 41
  Finesse             : 10

    Time     Theta     Theta'
--------------------------------------
  0.000000  3.140000  0.000000
  1.000000  3.137678 -0.005478
  ... (further time steps)
 40.000000  0.261876 -0.412985
```

**Real-World Problem Solving Ideas:**
*   **Simple Harmonic Motion Analysis**: Study the behavior of pendulums, which are classic examples of oscillatory systems. This helps understand concepts like period, amplitude, and damping.
*   **Damping Effects**: Investigate how different damping coefficients (`K'`) affect the decay of oscillations in a pendulum, applicable to understanding shock absorbers or vibration isolators.
*   **Non-Linear Dynamics**: Explore the behavior of a pendulum at large initial angles where the `sin(theta)` term becomes significant, demonstrating non-linear effects.
*   **Robotics and Control**: Model the swing dynamics of a single-joint robotic arm or a control system for an inverted pendulum.

**Dependencies**: `WinCrt`, `Type_def`, `Eqdifp`.

**Reference:**
*   "Analyse en Turbo Pascal versions 5.5 et 6.0" By Marc DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991 [BIBLI 03].

### 17. `rebounds.pas`: The Bouncing Ball Simulation

**Program Description:**
This program simulates the trajectory of an elastic ball bouncing on a flat ground, dropped with an initial horizontal speed from a given height. It calculates and visualizes successive parabolic rebound trajectories, accounting for energy loss due to damping.

**Core Functionality:**
The user provides the initial horizontal speed (`vx`), a damping coefficient (`amort`) (representing energy loss during impact), and the number of desired rebounds (`nrebounds`). The program starts with a ball falling from a fixed height (`Height = 12` units).
It then calculates and draws each parabolic section of the ball's trajectory. For each rebound, the impact velocity is scaled by the damping coefficient, and new parabolic coefficients are calculated to define the subsequent trajectory.

**Key Variables:**
*   `Height`: Initial drop height (fixed at 12 units).
*   `G`: Gravitational acceleration (fixed at 9.81 units).
*   `vx`: Initial horizontal speed.
*   `amort`: Damping coefficient (0 to 1).
*   `nrebounds`: Number of simulated rebounds.
*   `v0`: Speed at ground level (initial for each rebound).
*   `vy`: Vertical speed component.
*   `angle`: Angle of velocity vector with the horizontal (initial for each rebound).
*   `x1`, `x2`: Start and end horizontal positions for a parabolic segment.
*   `xs`, `ys`: Coordinates of the parabola's summit.
*   `a`, `b`, `c`: Coefficients of the parabolic equation `y = a*x^2 + b*x + c`.
*   `bounds`: Counter for completed rebounds.

**Mathematical Principles:**
The simulation combines principles of projectile motion and energy dissipation.
1.  **Projectile Motion**: Each segment of the trajectory between bounces is modeled as a parabolic path under constant gravitational acceleration, assuming no air resistance. The initial vertical velocity is derived from the drop height: `Vy = sqrt(2GH)`.
2.  **Impact Dynamics**: Upon impact with the ground, the total velocity `V0` is reduced by multiplying it with the `amort` damping coefficient (`v0:=v0*amort`), simulating inelastic collisions. The angle of rebound (from the horizontal) is assumed to be conserved, simplifying the model.
3.  **Parabolic Equation**: The trajectory of each bounce is represented by `y = a*x^2 + b*x + c`. The coefficients `a`, `b`, `c` of the parabolic trajectory are determined by solving a system of linear equations using Cramer's method. The three points used are:
    *   The start point of the bounce `(x1, 0)`.
    *   The end point of the bounce `(x2, 0)`.
    *   The apex (highest point) of the parabola `(xs, ys)`, where `ys = (V0^2 * sin^2(theta)) / (2G)` and `xs = (x1 + x2) / 2`.

**Sample Usage (from source comments):**
```
  The bouncing ball
  -----------------

  Initial horizontal speed: 5
  Damping coefficient.....: 0.9
  Number of rebounds......: 5

  (Graphical output showing the trajectory)
  Height= 12  Speed= 5.0  Total length= 10.39  Damping= 0.90
```

**Real-World Problem Solving Ideas:**
*   **Sports Science**: Analyze the physics of various sports (e.g., ball games like basketball, tennis, golf) by modeling the bounce characteristics, which can inform equipment design or player technique.
*   **Robotics/Automation**: Inform the design of robotic systems that interact with surfaces through bouncing or impacting actions.
*   **Physics Education**: A clear visual demonstration of projectile motion, energy conservation (or dissipation), and the effects of damping.
*   **Material Testing**: Qualitatively assess the elasticity of materials by observing the rebound characteristics of a ball dropped onto them.

**Dependencies**: `WinCrtMy` (renamed from `ef3d/wincrt1.pas`), `WinProcs`, `Strings`, `Type_def` (from `ef3d/type_def.pas`), `CrtGr2D` (renamed from `ef3d/graph_2d.pas`).

**Reference:**
*   "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0" By R. Dony - MASSON 1990 page 113 [BIBLI 12].

### 18. `_info.txt`: Program Descriptions Overview

**File Description:**
This is a descriptive text file providing a brief overview of the programs available in the `mechanics` folder. It serves as a quick reference for the functionalities offered by the library.

**Content:**
A bulleted list of the different calculation types and problems that the Pascal programs in this directory are designed to solve. It is not an executable program.

**Purpose:**
Provides a high-level summary for users to quickly identify relevant programs for their needs.

---

## Units and Libraries

Several Pascal units are included to support the main programs, providing common functionalities and abstracting complex logic.

### 1. `ef3d/algeblin.pas` - Linear Algebra Procedures

**Unit Description:**
This unit provides a collection of essential linear algebra procedures used primarily by the EF3D finite element program. It includes routines for matrix manipulation, decomposition, and solving linear systems and eigenvalue problems.

**Key Procedures:**
*   `Transpose_Matrice_Carree(n:Integer; Var A:pMatc)`: Transposes a square matrix `A` of size `n x n`.
*   `Matvec(n:Integer; Var K: pMatc; Var Ve,Vs: pVect; option:integer)`: Calculates the product of a square matrix `K` (of size `n x n`) by a column vector `Ve`, storing the result in `Vs`. `option` specifies matrix storage (0: full, 1: upper triangle, -1: lower triangle).
*   `Produit_Scalaire(n:Integer; Var U,V: pVect; Var result : arg_reel)`: Computes the scalar product of two vectors `U` and `V` of size `n`.
*   `Multiplie_Matrix_3x3(Var A,B,C: Matrix_3x3)`: Multiplies two 3x3 matrices `A` and `B`, storing the result in `C`.
*   `Transpose_Matrix_3x3(Var A,ATr:Matrix_3x3)`: Transposes a 3x3 matrix `A`, storing the result in `ATr`.
*   `Choldet1(n:integer; Var a:pMatc; Var p:pVect; Var d1:arg_reel; Var d2:Integer; Var Erreur:Boolean)`: Performs Cholesky decomposition of a symmetric positive definite matrix `a` (of size `n x n`) and computes its determinant (`d1` and `d2` for scientific notation). Sets `Erreur` if the matrix is not positive definite.
*   `Cholsol1(n,r : integer; Var a : pMatc; Var b,x: Matrice; Var p : pVect)`: Solves a linear system `L L^T x = b` (or `A x = b` where `A = L L^T`) using the Cholesky decomposition (`L` is stored in `a`, `p` is the diagonal inverse factor). `r` is the number of right-hand side vectors.
*   `Tred2(n:Integer; tol:arg_reel; Var a,z:pMatc; Var d,e:pVect)`: Reduces a symmetric matrix `a` (of size `n x n`) to a tridiagonal form using Householder transformations. The orthogonal transformation matrix is stored in `z`. The diagonal elements of the tridiagonal matrix are in `d`, and sub-diagonal elements in `e`.
*   `Tql2(n:Integer; macheps:arg_reel; Var d,e:pVect; Var z :pMatc; Var Erreur: Boolean)`: Calculates all eigenvalues (`d`) and eigenvectors (`z`) of a symmetric tridiagonal matrix. `macheps` is machine epsilon. Sets `Erreur` if convergence is not achieved.
*   `Reduc1(n:Integer; var a,b:pMatc; var dl:pVect; var erreur: Boolean)`: Reduces the generalized eigenvalue problem `[A]x = \lambda[B]x` to the symmetric standard form `[P]z = \lambda z`, where `[A]` and `[B]` are symmetric, and `[B]` is positive definite. `[P]` is stored in `a`, and the Cholesky factor of `[B]` (matrix `L`) is stored in `b` (lower triangle) and its diagonal in `dl`.
*   `Rebaka(n,m1,m2:Integer; b:pMatc; dl:pVect; Var z: pMatc)`: Constructs the eigenvectors of the original generalized problem from the eigenvectors of the symmetric standard problem. `m1` and `m2` define the range of eigenvectors. `b` and `dl` are the Cholesky factors from `Reduc1`.

**Relevance to Problem Solving:**
These procedures are fundamental to numerical methods in engineering and scientific computing. They underpin various structural analysis, vibration, and heat transfer simulations. For example, `Choldet1` and `Cholsol1` are essential for solving systems of equations arising from static structural analysis. `Tred2`, `Tql2`, `Reduc1`, and `Rebaka` are critical for solving eigenvalue problems in structural dynamics to find natural frequencies and mode shapes.

### 2. `ef3d/declara1.pas` - General Constants and Types Declarations

**Unit Description:**
This unit defines global constants and custom data types used across various programs within the `ef3d` subfolder. It centralizes common definitions to ensure consistency and modularity.

**Key Declarations:**
*   **Constants**:
    *   `Macheps`: Machine epsilon (1E-15).
    *   `MaxSize`: General maximum size for arrays (4096).
    *   `Size`: A smaller common size (1024).
    *   `npoin`: Max nodal points (512).
    *   `maxc`: Max. number of unknowns/columns (89, due to Turbo Pascal 64KB segment limit).
    *   `maxndp`: Max. number of nodal points (89).
    *   `maxddl`: Max. number of degrees of freedom (534, which is 6 * `maxc`).
*   **Types**:
    *   `arg_reel`: `Double` for high precision floating-point numbers.
    *   `Matrice`, `Matrice_carree`, `Matrice_n_mPlus1`, `Matrice_nxr`, `Matrix_3x3`, `Block12x12`: Various matrix types for different dimensions.
    *   `ddl_status`: `Array[1..6,1..maxndp] of Byte` for DOF status (0=free, 1=prescribed).
    *   `connexion`: `Array[1..2,1..maxc] of Byte` for element connectivity.
    *   `Classe`: `Array[1..Maxndp] of Integer` for arbitrary node numbering.
    *   `Vecteur_ligne`, `pVect`, `Vecteur_colonne`, `Vecteur_Position`, `Vecteur_Reel`, `Vecteur_Ent`, `Vecteur_Byte_long`, `Vecteur_Byte_court`: Various vector types.
    *   String types: `Mot_Cle`, `Descriptif`, `Enregistrement`, `Nom_Fichier`, `Suffixe`.
    *   `Gen_Data_Record`: A record type for storing general problem data, enabling saving/loading of entire model definitions for EF3D.

**Relevance to Problem Solving:**
This unit acts as a common data dictionary, making the code more readable, maintainable, and less prone to type inconsistencies. It allows complex data structures, such as those needed for finite element models, to be defined and reused across different program modules, ensuring data integrity and simplifying development.

### 3. `ef3d/errorhan.pas` and `ef3d/errors.pas` - Error Handling Units

**Unit Descriptions:**
*   `errors.pas`: Defines a structured way to categorize and describe common runtime and I/O errors.
*   `errorhan.pas`: Provides a generic error handler procedure that displays detailed information about errors caught during program execution.

**Key Components:**
*   `Errors` unit: Declares `Error_Type` record (Error ID, Error Message) and `Run_Errors`, `IO_Errors` constant arrays listing specific error codes and their descriptive messages (e.g., 'Floating point overflow', 'File does not exist').
*   `ErrorHan` unit: Contains the `ErrorHandler` procedure. This procedure is set as the program's exit procedure, meaning it is automatically called upon an unhandled runtime error. It interprets the error number and displays a user-friendly message, including the error type, ID, and memory address where the error occurred, along with advice to report the error.

**Relevance to Problem Solving:**
Robust error handling is crucial for any computational software. These units enhance the usability of the programs by providing informative error messages instead of cryptic system errors. This helps users and developers diagnose issues more effectively, whether it's a file not found, a mathematical singularity, or memory exhaustion, without requiring deep knowledge of the underlying Pascal runtime.

### 4. `ef3d/graph_2d.pas` - 2D Graphing Utilities

**Unit Description:**
This unit provides a set of procedures for drawing 2D curves (`y=f(x)`) within a specified window, with options for linear or logarithmic scales, automatic or manual axis adjustments, and basic drawing primitives. It's designed to facilitate graphical output for numerical results.

**Key Procedures:**
*   `Fenetre(P: HDC; num: word)`: Defines and opens one of 11 predefined partial windows (e.g., quarters, halves, full screen) within the current program window (`P` is the Device Context).
*   `PleinEcran`: Resets the drawing area to the full current window.
*   `Echelle(VAR ech :real_ar)`: Internal procedure to adjust scale factors (e.g., to 10, 5, 2.5, 2, 1.25, 1) for axis graduations.
*   `TracerLesAxes(P:HDC)`: Draws X and Y axes.
*   `GraduerLesAxes(P:HDC)`: Adds graduations and labels to the axes (including grid lines).
*   `InitFenetre(P,fntr:integer;xmn,xmx,ymn,ymx:real_ar)`: Initializes a virtual drawing window with physical coordinates (`xmn, xmx, ymn, ymx`), sets scales (automatic or manual), and draws axes/graduations. `fntr` selects the window type.
*   `MoveXY(P:HDC; x,y: real_ar)`: Moves the drawing pen to a specified physical coordinate `(x,y)`.
*   `LineXY(P:HDC; x,y: real_ar)`: Draws a line from the current pen position to `(x,y)`.
*   `TextXY(P:HDC;x,y:real_ar;text:PChar)`: Writes a text string `text` at physical coordinate `(x,y)`.
*   `CroixXY(P:HDC;x,y:real_ar)`: Draws a cross symbol at physical coordinate `(x,y)`.
*   `Legendes (P:HDC; titre,titrex,titrey:Pchar)`: Adds a main title to the graph (`titre`) and labels for X (`titrex`) and Y (`titrey`) axes.
*   `MinMax(n:integer;Y:RV; VAR ymin,ymax:real_ar)`: Finds minimum and maximum values in a dynamic `Real_Vector` `Y` of `n` elements.
*   `CourbeXY(P:HDC;n,fntr: integer; Y: RV; xn,xm: real_ar)`: Draws a 2D curve from `n` tabulated values in `Y` over the x-range `[xn, xm]`, including initialization of the drawing window (calls `InitFenetre`).
*   `TracerXY(P:HDC;n: integer; Y: RV; xn,xm: real_ar)`: Draws an additional 2D curve on an already initialized window, ideal for superimposing multiple curves.
*   `Circle(P:HDC;xc,yc,r: real_ar; trait: boolean)`: Draws a circle with center `(xc,yc)` and radius `r`, as a solid or dotted line.
*   `WinCrtInit(Nom:PChar)`: Initializes the underlying `WinCrt1` window with a title `Nom`.
*   `SortieGraphique`: Provides options for saving the generated graph to a file (`.IMG`), printing, or continuing/exiting.

**Relevance to Problem Solving:**
Graphical representation of results is vital for interpreting complex numerical data. This unit allows users to visualize time-history responses, frequency response functions, mode shapes, or any `y=f(x)` relationship computed by the mechanics programs. This aids in quickly understanding trends, identifying critical points (e.g., resonance peaks), and verifying computational accuracy.

### 5. `ef3d/savecrt.pas` - Screen Save/Load Utilities

**Unit Description:**
This unit provides procedures for saving and loading black-and-white screenshots of the `WinCrt1` window to and from disk files.

**Key Procedures:**
*   `WCrttoFile( P : HDC; name : string)`: Captures the current `WinCrt1` window content (pixel by pixel, converting to B&W based on pixel color) and saves it to a specified file with `.IMG` extension.
*   `WLoadCrt( P : HDC; name : string)`: Loads a previously saved `.IMG` file and displays it on the `WinCrt1` window, drawing blue pixels where black pixels were saved.

**Relevance to Problem Solving:**
This utility is useful for documenting results or creating a visual archive of simulations without needing a dedicated screenshot tool. Users can save graphical outputs (e.g., mode shapes, pendulum trajectories, bouncing ball path) for later review or inclusion in reports.

### 6. `ef3d/type_def.pas` - Generic Type Definitions for Graphing

**Unit Description:**
This unit defines generic types and global variables primarily used by the `graph_2d.pas` unit and other programs that utilize 2D graphing functionalities.

**Key Declarations:**
*   **Constants**:
    *   `Macheps`: Machine epsilon (1E-12, a general floating-point precision constant).
    *   `Size`: Default array size for vectors (2048, can be modified by user).
*   **Types**:
    *   `real_ar`: `Double` (for high precision real numbers).
    *   `complex`: Array for complex numbers (`Array[1..2] of REAL`).
    *   `Real_Vector`, `Integer_Vector`, `Complex_Vector`: Dynamic array types for vectors.
    *   `RV`, `CV`: Pointers to `Real_Vector` and `Complex_Vector` for dynamic memory allocation.
    *   `Descriptif`, `FileName`, `Title`: String types for labels and filenames.
*   **Variables**:
    *   `MaxX`, `MaxY`: Screen resolution (e.g., 800x600 for SVGA by default, used by `graph_2d`).
    *   `Log_X`, `Log_Y`, `Ech_auto`: Booleans to control logarithmic scaling for X/Y axes and automatic axis adjustment.
    *   `X_mini`, `X_maxi`, `Y_mini`, `Y_maxi`: Manual axis limits for plotting.
    *   `CrtDC`: Device context handle for the `WinCrt1` window, used for direct GDI drawing.

**Relevance to Problem Solving:**
Similar to `declara1.pas`, this unit provides standardized types and global configuration options for graphical output across multiple programs. It ensures consistency in data representation for plotting and allows for easy configuration of graph appearance. This promotes modularity and reusability across the library's graphical programs.

### 7. `ef3d/utilit.pas` - Utility Procedures for I/O

**Unit Description:**
This unit contains utility procedures for reading and writing real numbers and integers from/to the console and text files, along with basic file handling.

**Key Procedures/Functions:**
*   `evalue_log(l : ar_reel) : INTEGER`, `power10(n : INTEGER) : ar_reel`: Internal helper functions for `disp_real` to handle scientific notation.
*   `Read_integer(min:INTEGER;max:longint):INTEGER`: Reads an integer from input, performing error control to ensure it's within a specified range and is a valid number.
*   `Read_filename(VAR fich:TEXT)`: Prompts the user for a filename, assigns it to a `TEXT` file variable `fich`, and attempts to open it for reading. Includes error checking for file existence.
*   `disp_real(l2 : ar_reel)`: Displays a `double` (real number) to the console using a controlled 10-character wide format, automatically switching to scientific notation for very large or very small numbers.
*   `f_disp_real(VAR f:TEXT; l2 : ar_reel)`: Same as `disp_real` but writes the formatted real number to a specified text file `f`.
*   `WriteHead(VAR f:TEXT; s:STRING)`: Writes a formatted header (separator lines and a string `s`) to a text file `f`.
*   `WriteEnd(VAR f:TEXT)`: Writes a separator line to a text file `f`, typically used as a footer.

**Relevance to Problem Solving:**
These utilities simplify common input/output tasks, making the programs more robust and user-friendly by handling common errors like invalid numerical input or non-existent files. The formatted real number display ensures consistent and readable numerical output, which is important for presenting quantitative results in a clear and professional manner.

### 8. `ef3d/wincrt1.pas` - Modified Windows CRT Interface

**Unit Description:**
This is a modified version of Borland's `WinCrt` unit, designed to provide a console-like text and graphics window environment within Windows applications. It's essential for the user interface and basic graphical capabilities of the programs.

**Key Functionality:**
*   **Window Management**: Creates and manages a standard Windows console window with a customizable title, initial position, and size.
*   **Text Input/Output**: Provides basic text input/output functions, including `WriteBuf` (writing a buffer of characters), `WriteChar` (writing a single character), `ReadKey` (reading a single character, blocking until input is available), and `ReadBuf` (reading a buffer of characters, handling backspace).
*   **Cursor Control**: Manages cursor visibility, position (`GotoXY`, `WhereX`, `WhereY`), and ensures the cursor is visible during input (`TrackCursor`).
*   **Screen Manipulation**: Offers screen clearing capabilities (`ClrScr`, `ClrEol`) and basic scrolling.
*   **Keyboard Handling**: Processes keyboard messages, including special keys (e.g., arrow keys, Ctrl-C for break).
*   **Graphics Context**: Provides a device context handle (`CrtDC`) which allows other units (like `graph_2d.pas`) to perform direct GDI (Graphics Device Interface) drawing operations within the console window.
*   **Event Loop**: Contains an internal message loop to process Windows messages, ensuring the window remains responsive.

**Relevance to Problem Solving:**
This unit forms the backbone of the user interaction for many programs in the `mechanics` folder. It provides a simple, consistent environment for displaying prompts, reading user input, and presenting textual results. Its integration with graphical capabilities allows for mixed text and graphical output within the same application window, which is crucial for programs like EF3D and `rebounds.pas` that need to display both numerical tables and visual plots. This unit effectively abstracts away the complexities of direct Windows API programming for basic console and graphical I/O.
