
# jpmMath - Series

## Contents

The `Series` module provides implementations for the following mathematical programs and subroutines:

*   **Asymptotic Error Function (`ASYMERF`)**: Calculates `erf(x)` and `erfc(x)` using an asymptotic series, particularly effective for large `x`.
*   **Chebyshev Economization (`CHEBECON`)**: Reduces the degree of a polynomial approximation while preserving accuracy, leveraging Chebyshev polynomials.
*   **Chebyshev Series Coefficients (`CHEBYSER`)**: Generates the coefficients for Chebyshev polynomials of the first kind.
*   **Chi-Square Distribution Functions (`CHI_SQ`, `CHI_SQR`)**: Computes both the Probability Density Function (PDF) and Cumulative Distribution Function (CDF) for the Chi-Square statistic.
*   **Logarithm of Factorial X (`LN_FACTX`)**: Provides a stable and accurate approximation for `ln(x!)`, crucial for large arguments in statistical calculations.
*   **Complex Series Evaluation (`CMPLXSER`)**: Extends polynomial evaluation to complex arguments, crucial for many engineering and physics applications.
*   **Horner's Rule (`HORNER`)**: Implements an algorithm for shifting the expansion point of a polynomial, significantly improving numerical stability.
*   **Inverse Normal Distribution (`INVNORM`)**: Approximates the inverse of the complementary error function, widely used to find Z-scores for given tail probabilities.
*   **Polynomial Inversion (`RECIPRO`)**: Calculates the coefficients of a polynomial that approximates the reciprocal of another polynomial (`1/P(x)`).
*   **Series Reversion (`REVERSE`)**: Computes the coefficients of the inverse power series `X(Y)` given `Y(X)`.

---

## Detailed Program/Subroutine Documentation

Each section below provides a detailed description of a program or subroutine, covering its purpose, underlying mathematical principles, algorithmic details, input/output specifications, usage notes with sample data (where available in the source code), and potential real-world applications.

### 1. ASYMERF (Asymptotic Error Function)

*   **Program File**: `asymerf.pas`
*   **Explanation File**: `asymerf.txt`

**Description:**
The `ASYMERF` program calculates the value of the error function, `erf(x)`, and its complement, `erfc(x) = 1 - erf(x)`, using an asymptotic series expansion. This method is particularly effective for larger values of the argument `x` (typically `x > 3`), where standard Taylor series approximations might converge slowly or numerical tables might lack sufficient detail. The core `Asym_erf` subroutine determines the level of maximum accuracy for the approximation.

**Mathematical Background:**
The error function `erf(x)` is defined as:
```
              x   -t^2 dt
erf(x) = 2/sqrt(pi) ∫ e
              0
```
For large positive values of `x`, the `Asym_erf` procedure implements an asymptotic expansion derived by Peirce. The expansion for `erf(x)` is:
```
erf(x) = 1 - [1/(2x^2)] + [1*3/(2x^2)^2] - [1*3*5/(2x^2)^3] + ...
```
As with most asymptotic series, this expression eventually diverges. The procedure terminates the summation at the smallest term (where the magnitude of the next term begins to increase). The magnitude of this first neglected term provides an estimate of the error bound. For `x <= 3`, a rapidly convergent Taylor series approximation could be used, though this program specifically focuses on the asymptotic approach for larger `x`.

**Algorithm (`Asym_erf` procedure in `asymerf.pas`):**
The `Asym_erf` procedure computes `erfc(x)` by iteratively summing terms of its asymptotic series.
1.  Initializes `y = 1`, `n = 1`, and the first term `c2 = 1 / (2 * x * x)`.
2.  Enters a loop where it subtracts `c2` from `y`, increments `n` by 2, and updates `c1` to `c2`.
3.  The next term `c2` is calculated as `-c1 * n / (2 * x * x)`.
4.  The loop continues as long as `ABS(c2) <= ABS(c1)`, which is the condition for convergence (or rather, for summing terms up to the smallest one before divergence begins).
5.  Once the loop terminates, `n` is adjusted (`(n + 1) DIV 2`).
6.  An intermediate factor `e` is calculated as `EXP(-x * x) / (x * SQRT(PI))`.
7.  `y1` (approximating `erfc(x)`) is computed as `y * e`.
8.  `y` (approximating `erf(x)`) is then `1.0 - y1`.
9.  The error estimate `e` is updated to `e * c2`.

**Inputs (Global Variables):**
*   `x` (`DOUBLE`): The argument for which to evaluate `erf(x)`. It must be `x > 0`.

**Outputs (Global Variables):**
*   `y` (`DOUBLE`): The calculated `erf(x)`.
*   `y1` (`DOUBLE`): The calculated `erfc(x)`.
*   `e` (`DOUBLE`): An estimate of the error, approximately equal to the first term neglected in the series summation.
*   `n` (`INTEGER`): The number of terms evaluated in the series.

**Usage Notes/Sample:**
The main program prompts for the value of `X`. The output provides `ERF(X)`, the error estimate, and the number of terms evaluated.

```
SAMPLE RUN:
Input X ? 3
ERF(X)= .9999779 with error estimate=-0.00000000
Number of terms evaluated was 10

Input X ? 4
ERF(X)= 1.0000000 with error estimate= 0.0000000
Number of terms evaluated was 17
```

**Real-world Applications:**
*   **Heat Conduction and Diffusion**: Fundamental for solving transient heat transfer problems (e.g., cooling of semi-infinite solids, heat penetration into a material) and modeling ionic diffusion processes in materials science or chemical engineering.
*   **Probability and Statistics**: Used in calculations related to the normal distribution's cumulative density function, particularly for probabilities in the tails of the distribution. It's related to the Z-score and p-value calculations for large deviations.
*   **Fluid Dynamics**: Analyzing boundary layer flows or dispersion of pollutants in fluids.

**Limitations:**
*   The `Asym_erf` procedure is specifically designed for `x > 0`.
*   No explicit error checks are performed on the input `x` to ensure it is positive.

### 2. CHEBECON & CHEBYSER (Chebyshev Economization and Series Coefficients)

This pair of programs provides robust tools for polynomial approximation. `CHEBYSER` is a foundational utility for generating Chebyshev polynomial coefficients, which are then utilized by `CHEBECON` to perform polynomial economization. This technique reduces the degree of a polynomial approximation while maintaining a controlled level of accuracy, significantly optimizing computational performance.

#### 2.1. CHEBYSER (Chebyshev Series Coefficients Evaluation)

*   **Program File**: `chebyser.pas`

**Description:**
The `CHEBYSER` program contains the `Cheby_Ser` procedure, which calculates the coefficients for Chebyshev polynomials of the first kind (`T_n(x)`) when expressed as a sum of powers of `x`. Chebyshev polynomials possess unique properties that make them ideal for approximating functions and for polynomial economization, as they distribute the approximation error uniformly across a given interval.

**Mathematical Background:**
Chebyshev polynomials of the first kind, `T_n(x)`, are defined by the recurrence relation:
*   `T_0(x) = 1`
*   `T_1(x) = x`
*   `T_n(x) = 2x * T_{n-1}(x) - T_{n-2}(x)` for `n >= 2`.

When expanded, `T_n(x)` can be written as a polynomial in `x`:
`T_n(x) = Σ_{j=0}^{n} B_{n,j} * x^j`
The `Cheby_Ser` procedure computes these `B_{n,j}` coefficients for a given `n`.

**Algorithm (`Cheby_Ser` procedure in `chebyser.pas` and `chebecon.pas`):**
The `Cheby_Ser` procedure implements the recurrence relation to determine the coefficients `B[i,j]`, where `i` is the degree of the Chebyshev polynomial `T_i(x)` and `j` is the power of `x` (i.e., `x^j`).
1.  **Base Cases**: Initializes coefficients for `T_0(x)` (`B[0,0] = 1`) and `T_1(x)` (`B[1,0] = 0`, `B[1,1] = 1`).
2.  **Recursive Calculation**: For `i` from 2 up to `n`, it calculates `B[i,j]` using the recurrence relation: `B[i, j] := 2 * B[i - 1, j - 1] - B[i - 2, j]`. The constant term is also determined as `B[i,0] := -B[i-2,0]`.

**Inputs (Global Variable):**
*   `n` (`INTEGER`): The desired degree of the Chebyshev polynomial for which to generate coefficients.

**Outputs (Global Variable):**
*   `B` (`2D Array of DOUBLE`): Stores the computed coefficients. `B[i, j]` represents the `j`-th coefficient of the `i`-th degree Chebyshev polynomial.

**Usage Notes/Sample:**
The main program in `chebyser.pas` iterates from `n = 2` to `SIZE` (defined as 10), calling `Cheby_Ser` for each `n` and then printing the calculated coefficients for `T_n(x)`.

```
SAMPLE RUN (excerpt):
Chebyshev polynomial coefficients for degree 2
  A(0) = -1
  A(1) = 0
  A(2) = 2

Chebyshev polynomial coefficients for degree 3
  A(0) = 0
  A(1) = -3
  A(2) = 0
  A(3) = 4
```

**Real-world Applications:**
*   **Function Approximation**: A cornerstone for approximating functions with polynomials that exhibit uniform error bounds.
*   **Numerical Integration and Differentiation**: Chebyshev methods offer highly accurate techniques for these operations.
*   **Spectral Methods**: Used as basis functions in solving differential equations and in analyzing oscillatory data.
*   **Filter Design**: In electrical engineering, Chebyshev filters are designed based on these polynomials due to their optimal approximation properties.

**Limitations:**
*   The `B` array size is fixed (`SIZE = 10` in `chebyser.pas`, `SIZE = 40` in `chebecon.pas`), limiting the maximum degree of Chebyshev polynomial whose coefficients can be generated.

#### 2.2. CHEBECON (Chebyshev Economization)

*   **Program File**: `chebecon.pas`

**Description:**
The `CHEBECON` program demonstrates Chebyshev economization, a technique to reduce the degree of a polynomial approximation while maintaining a controlled level of accuracy over a specified interval. This is achieved by converting the original polynomial's coefficients into Chebyshev series coefficients, truncating the series by discarding higher-order terms, and then converting the resulting lower-degree Chebyshev series back into a power series. This method is highly effective for optimizing polynomial evaluations and is crucial for balancing accuracy with computational efficiency.

**Mathematical Background:**
Chebyshev economization leverages the properties of Chebyshev polynomials to minimize the maximum error of polynomial approximations over an interval. The process involves:
1.  **Scaling**: The input polynomial `P(x)` defined over an interval `[-x0, x0]` is transformed to operate on the standard Chebyshev interval `[-1, 1]` via a change of variable `x' = x/x0`.
2.  **Conversion to Chebyshev Series**: The scaled power series `P(x') = Σ C'_i (x')^i` is converted into a Chebyshev series `P(x') = Σ A_i T_i(x')`.
3.  **Truncation**: The Chebyshev series is truncated by discarding terms with coefficients `A_i` for `i > m1`, where `m1` is the desired lower degree. The error introduced by truncation is bounded and well-distributed due to the properties of Chebyshev polynomials.
4.  **Conversion back to Power Series**: The truncated Chebyshev series `P_economized(x') = Σ_{i=0}^{m1} A_i T_i(x')` is converted back into a standard power series `P_economized(x') = Σ C''_i (x')^i`.
5.  **Un-scaling**: The resulting polynomial is un-scaled to operate on the original interval `[-x0, x0]`.

**Algorithm (`Cheby_Econ` procedure in `chebecon.pas`):**
1.  **Scaling Input Coefficients**: Multiplies `C[i]` by `x0^i` to scale them for the `[-1, 1]` interval.
2.  **Conversion to Chebyshev Coefficients (`A` array)**: Iterates `n` from `m` (original degree) down to 0. For each `n`:
    *   Calls `Cheby_Ser(n)` to get the `x`-power coefficients for `T_n(x)` in `B[n,j]`.
    *   Calculates `A[n] = C[n] / B[n,n]`.
    *   Subtracts `A[n] * T_n(x)` from the current polynomial `C(x)` (by subtracting `A[n] * B[n,l]` from `C[l]`) to effectively remove the `T_n(x)` component.
3.  **Truncation and Reconversion to Power Series (`C` array)**: Iterates `i` from 0 to `m1` (desired degree):
    *   For each `i`, sums up `A[i] * B[i,j]` into `C[j]`, effectively reconstructing the economized polynomial.
4.  **Un-scaling Output Coefficients**: Divides `C[i]` by `x0^i` to convert them back to the original `[-x0, x0]` range.

**Inputs (Global Variables):**
*   `m` (`INTEGER`): The initial degree of the input polynomial.
*   `m1` (`INTEGER`): The desired degree of the economized polynomial (`m1 < m`).
*   `x0` (`DOUBLE`): The maximum range of `x` for which the polynomial is defined (e.g., if valid for `[-2, 2]`, `x0` is `2`), used for scaling.
*   `C` (`Array[0..SIZE] of DOUBLE`): The coefficients `C[0]` to `C[m]` of the input polynomial. (Note: This array is modified in-place to store the economized polynomial coefficients).

**Outputs (Global Variables):**
*   `A` (`Array[0..SIZE] of DOUBLE`): The Chebyshev series coefficients `A[0]` to `A[m]`.
*   `C` (`Array[0..SIZE] of DOUBLE`): The coefficients `C[0]` to `C[m1]` of the economized polynomial.

**Usage Notes/Sample:**
The program prompts for the initial and desired degrees, the range `x0`, and then the coefficients of the input polynomial. It then outputs the intermediate Chebyshev coefficients and the final economized polynomial coefficients.

```
SAMPLE RUN (excerpt):
What is the degree of the input polynomial: ? 15
What is the degree of the desired economized polynomial: ? 9
What is the range of the input polynomial: ? 1.57
Input the coefficients:
C( 0) = ? 0
C( 1) = ? 1
C( 2) = ? 0
C( 3) = ? -0.166666666
... (other coefficients for sin(x) Taylor series)

The Chebyshev series coefficients are:
A( 0) =  0.0000000000
A( 1) =  1.1334708982
...
A(15) =  0.0000000000

The economized polynomial coefficients are:
C( 0) =  0.0000000000
C( 1) =  0.9999999767
...
C( 9) =  0.0000025906
```

**Real-world Applications:**
*   **Numerical Libraries**: Optimizing polynomial approximations for standard mathematical functions (e.g., `sin(x)`, `cos(x)`, `exp(x)`) in scientific computing libraries, leading to faster execution and reduced floating-point errors.
*   **Embedded Systems**: Implementing mathematical functions on resource-constrained devices where computational speed, memory usage, and power consumption are critical.
*   **Computer Graphics and Games**: Efficiently evaluating curves and surfaces defined by polynomials, such as Bezier curves, for faster and smoother rendering without perceptible loss of detail.
*   **Data Analysis**: Reducing model complexity while maintaining predictive power, e.g., in regression models or surrogate modeling.

**Limitations:**
*   The `C` and `A` arrays are hardcoded to `SIZE = 40`, limiting the maximum degree of polynomial that can be processed.
*   The input `C` array is overwritten with the economized coefficients.

### 3. CHI-SQUARE DISTRIBUTION FUNCTIONS: PDF and CDF

This subsection provides robust utilities for working with the Chi-Square probability density function (PDF) and cumulative distribution function (CDF). These functions are fundamental in statistical hypothesis testing. The module includes an accurate auxiliary subroutine for the logarithm of the factorial function, which is crucial for calculating the Gamma function terms required in the Chi-Square formulas, especially for large arguments.

#### 3.1. LN_FACTX (Logarithm of Factorial X)

*   **Program File**: `logn!.pas`
*   **Dependency for**: `chi_sq.pas`, `chi_sqr.pas`

**Description:**
The `LN_FACTX` subroutine calculates the natural logarithm of the factorial of a number `x` (i.e., `ln(x!)`). This function is invaluable in scenarios where `x!` itself would exceed the maximum representable number in standard floating-point types, particularly for large `x`. By operating in the logarithmic domain, calculations involving very large factorials become feasible, which are common in advanced statistical and combinatorial problems.

**Mathematical Background:**
The Gamma function, `Gamma(z+1) = z!`, extends the factorial concept to real and complex numbers. For large `x`, `ln(x!) = ln(Gamma(x+1))` can be accurately approximated using an asymptotic expansion that is a refined version of Stirling's approximation. The formula implemented is:
```
ln(x!) = (x + 0.5) * ln(x) - x * (1 - 1/(12x^2) + 1/(360x^4) - 1/(1260x^6) + 1/(1680x^8)) + 0.918938533205
```
The constant `0.918938533205` is approximately `1/2 * ln(2*pi)`. This expansion provides high accuracy for larger values of `x`.

**Algorithm (`LN_FACTX` procedure in `logn!.pas`, `chi_sq.pas`, `chi_sqr.pas`):**
The procedure directly implements the truncated asymptotic expansion for `ln(x!)`. It calculates `x1 = 1 / (x * x)` and then uses powers of `x1` in the series sum to compute the terms within the parenthesis.

**Inputs (Global Variable):**
*   `x` (`DOUBLE`): The number for which to calculate `ln(x!)`. `x` should be `> 0`. For optimal accuracy of this specific series, `x >= 3` is recommended, and accuracy improves significantly for `x > 10`.

**Outputs (Global Variable):**
*   `y` (`DOUBLE`): The calculated value of `ln(x!)`.

**Usage Notes/Sample:**
The `logn!.pas` program demonstrates its accuracy by calculating `LN(X!)` for `X` from 1 to 15 and also `EXP(LN(X!))` to show the original factorial value.

```
SAMPLE RUN (excerpt):
   X       LN(X!)      EXP(LN(X!))
----------------------------------
    1    -0.000307                1
    2     0.693146                2
    3     1.791759                6
    4     3.178054               24
    5     4.787492              120
...
   15    27.899271    1307674368000
```

**Real-world Applications:**
*   **Probability and Statistics**: Essential for computing probabilities involving large combinations and permutations, such as in population genetics, statistical mechanics, and Bayesian inference (e.g., for likelihood functions of complex models). It is a core component for computing various statistical distributions, including the Chi-Square, Gamma, and Beta distributions.
*   **Thermodynamics and Statistical Physics**: Used in entropy calculations for large systems, where the number of microstates involves factorials that are astronomically large.
*   **Data Science and Bioinformatics**: When dealing with very large sample sizes or combinatorial problems in areas like machine learning, network analysis, and genomic sequence analysis.

**Limitations:**
*   The approximation's accuracy improves with larger `x`. For very small `x` (e.g., `x < 3`), direct calculation or other specific approximations might be more appropriate or accurate, although the routine produces results.
*   No explicit error handling for invalid inputs like `x <= 0`.

#### 3.2. Chi_Square (Chi-Square Probability Density Function - PDF)

*   **Program File**: `chi_sqr.pas`
*   **Explanation File**: `chi_sq.txt`
*   **Dependency**: `LN_FACTX`

**Description:**
The `Chi_Square` subroutine calculates the value of the Chi-Square probability density function (PDF), `p(x)`, for a given number of degrees of freedom `M` and a specific Chi-Square value `x`. This function is fundamental in statistical analyses for describing the probability distribution of the sum of squared standard normal variates, especially in goodness-of-fit tests.

**Mathematical Background:**
The Chi-Square probability density function `p(x)` is defined as:
```
           e^(-x/2) * x^(M/2 - 1)
p(x) = -----------------------------   (Equation 2.3.4 from chi_sq.txt)
       (2^(M/2)) * Gamma(M/2)
```
where `x = chi^2` and `M` represents the number of degrees of freedom.
To handle potential numerical overflow when computing `Gamma(M/2)` for large `M` (e.g., `Gamma(100) = 99!`), the calculation is performed in the logarithmic domain. The natural logarithm of `p(x)` is:
```
ln[p(x)] = -x/2 + (M/2 - 1)ln(x) - (M/2)ln(2) - ln[Gamma(M/2)]   (Equation 2.3.6)
```
The `ln[Gamma(M/2)]` term is approximated using the `LN_FACTX` subroutine (since `Gamma(z+1) = z!`). Finally, `p(x)` is obtained by exponentiating the result: `p(x) = exp{ln[p(x)]}` (Equation 2.3.7).

**Algorithm (`Chi_Square` procedure in `chi_sqr.pas`):**
1.  Saves the input `x` to `m1`.
2.  Sets `x` to `m/2 - 1` and calls the `LN_FACTX` subroutine to compute `ln((m/2 - 1)!)`, which approximates `ln[Gamma(m/2)]`. The result is stored in the global `y` variable.
3.  Restores `x` from `m1`.
4.  Calculates `c` using the logarithmic form of the PDF equation: `c := -x / 2 + (m / 2 - 1) * LN(x) - (m / 2) * LN(2) - y`.
5.  The final result `y` (the Chi-Square PDF) is then `EXP(c)`.

**Inputs (Global Variables):**
*   `m` (`INTEGER`): The number of degrees of freedom. Must be a non-negative integer.
*   `x` (`DOUBLE`): The value of the chi-square statistic. Must be positive.

**Outputs (Global Variable):**
*   `y` (`DOUBLE`): The calculated value of the chi-square probability density function, `p(x)`.

**Usage Notes/Sample:**
The `chi_sqr.pas` program prompts for the number of degrees of freedom (`m`), a range (`X1`, `X2`), and a table step size (`X3`). It then prints a table of `X` values and their corresponding Chi-Square PDF values.

```
SAMPLE RUN (excerpt):
How many degrees of freedom: 100
What is the range (X1,X2):
  X1: 50
  X2: 150
What is the table step size: 5

   X       Chi-Square PDF
------------------------
   50        0.00000
   55        0.00003
   60        0.00018
...
  100        0.02816
...
  150        0.00017
```

**Real-world Applications:**
*   **Goodness-of-Fit Tests**: Directly provides the density for observed chi-square values, aiding in determining the probability of such observations under a null hypothesis, e.g., in tests to see if observed data fits a theoretical distribution (like normal, Poisson, etc.).
*   **Statistical Modeling**: Used as a building block for more complex statistical models, particularly those that involve distributions of variances or sums of squared normal variables.
*   **Educational Tools**: Demonstrates the computation of a fundamental statistical distribution, useful for teaching and understanding statistical theory.

**Limitations:**
*   No explicit error checks on the inputs (e.g., `m` as a non-negative integer, `x` as a positive real). The calling program is responsible for validating these to prevent mathematical domain errors.
*   Accuracy is dependent on the `LN_FACTX` approximation, which is very good even for `M` as low as 4.

#### 3.3. Chi_Square_Cumul (Chi-Square Cumulative Distribution Function - CDF)

*   **Program File**: `chi_sq.pas`
*   **Explanation File**: `chi_sq.txt`
*   **Dependencies**: `Chi_Square`, `LN_FACTX`

**Description:**
The `Chi_Square_Cumul` subroutine approximates the Chi-Square cumulative distribution function (CDF), `P(x)`. This function represents the probability that a Chi-Square random variable with `m` degrees of freedom will take a value less than or equal to `x`. This is vital for statistical hypothesis testing, where `P(x)` helps determine the significance of observed data.

**Mathematical Background:**
The Chi-Square CDF `P(x)` is approximated using the following asymptotic series expansion (Equation 2.3.5 from `chi_sq.txt`):
```
             2x           inf                x^k
P(x) = --- * p(x) * [1 + Sum --------------------------- ]
             M            k=1 (M + 2)(M + 4) ... (M + 2k)
```
where `p(x)` is the chi-square probability density function (calculated by `Chi_Square`). The sum iterates until the contribution of the current term (`x2` in the code) falls below a specified convergence factor or error bound (`ee`). Since the series summed in this equation contains terms of the same sign, little round-off error is expected in the final result.

**Algorithm (`Chi_Square_Cumul` procedure in `chi_sq.pas`):**
1.  Initializes `y1 = 1` (the first term of the series `[1 + Sum ...]`) and `x2 = x / (m + 2)`. `m2` is set to `m + 2`.
2.  Enters a loop where it adds `x2` to `y1`.
3.  If `x2` is less than the error bound `ee`, it exits the loop.
4.  Otherwise, it increments `m2` by 2 and calculates the next term `x2` using `x2 := x2 * (x / m2)`. This form helps avoid potential overflow for large `x` or `m`.
5.  After the summation loop, it calls the `Chi_Square` subroutine to obtain `p(x)`. The result is stored in the global `y`.
6.  Finally, it calculates the CDF `y` as `(y1 * y * 2) * (x / m)`.

**Inputs (Global Variables):**
*   `m` (`INTEGER`): The number of degrees of freedom.
*   `x` (`DOUBLE`): The value of the chi-square statistic for which to calculate the CDF.
*   `ee` (`DOUBLE`): The summation truncation error bound, controlling the accuracy of the series summation.

**Outputs (Global Variable):**
*   `y` (`DOUBLE`): The calculated value of the chi-square cumulative distribution function, `P(x)`.

**Usage Notes/Sample:**
The `chi_sq.pas` program demonstrates its use by prompting for degrees of freedom (`m`), a range (`X1`, `X2`), a step size (`X3`), and a summation truncation error bound (`ee`). It then generates a table of `X` values and their corresponding Chi-Square CDF values.

```
SAMPLE RUN (excerpt):
How many degrees of freedom: 100
What is the range (X1,X2):
  X1: 50
  X2: 150
What is the table step size: 5
Summation truncation error bound: 1e-6

   X       Chi-Square CDF
------------------------
   50       0.00001
   55       0.00007
   60       0.00052
...
  100       0.51881
...
  150       0.99910
```

**Real-world Applications:**
*   **Statistical Hypothesis Testing**: Directly used to find p-values in Chi-Square tests for goodness-of-fit (e.g., testing if observed data fits a theoretical distribution) and independence (e.g., analyzing contingency tables). This is fundamental for making statistical decisions.
*   **Experimental Design**: Determining critical values for statistical tests to establish thresholds for significance in experimental outcomes.
*   **Medical Research**: Analyzing categorical data from clinical trials to assess relationships between variables (e.g., treatment efficacy across different patient groups).
*   **Quality Control**: Assessing if defect rates conform to a specific distribution or if changes in a process have a statistically significant impact.

**Limitations:**
*   Relies on the `Chi_Square` subroutine, which in turn relies on `LN_FACTX`. It inherits their accuracy characteristics.
*   No explicit error checks on the inputs `m` or `x`.

### 4. CMPLXSER (Complex Series Evaluation)

*   **Program File**: `cmplxser.pas`
*   **Explanation File**: `cmplxser.txt`

**Description:**
The `CMPLXSER` program is designed to evaluate a polynomial `P(Z)` with real coefficients `A(i)` for a complex argument `Z = x + iy`. This extends the utility of standard polynomial evaluation routines into the complex plane, which is critical in various engineering, physics, and applied mathematics domains where complex numbers are inherent to problem formulations.

**Mathematical Background:**
A polynomial `P(Z)` of degree `m` with real coefficients `A[i]` is defined as:
```
P(Z) = A[0] + A[1]Z + A[2]Z^2 + ... + A[m]Z^m
```
where `Z` is a complex number (`Z = x + iy`). The core computational challenge is efficiently calculating `Z^n` for a complex number `Z`. This is accomplished by converting `Z` to its polar form `u * (cos(v) + i*sin(v))`, where `u` is the magnitude (`u = |Z| = sqrt(x^2 + y^2)`) and `v` is the phase angle (`v = arg(Z) = atan2(y, x)`).
Then, De Moivre's Theorem is applied to compute powers:
```
Z^n = u^n * (cos(nv) + i*sin(nv))
```
The result is then converted back to rectangular form `(real + i * imaginary)` for summation.

**Algorithm (`Complex_Series` procedure in `cmplxser.pas`):**
1.  Initializes the real (`z1`) and imaginary (`z2`) parts of the sum with the constant term `A[0]`.
2.  Saves the original real (`x`) and imaginary (`y`) parts of `Z` into `a1` and `a2` respectively.
3.  Iterates `n` from 1 to `m` (the degree of the polynomial):
    *   Restores the current `x` and `y` from `a1` and `a2` (since the internal power calculation modifies them).
    *   Calls the `RectPower` subroutine to calculate `(x + iy)^n`. The result overwrites `x` and `y` with its real and imaginary parts. `RectPower` itself orchestrates the following:
        *   `RectPol`: Converts `(x, y)` to polar `(u, v)`.
        *   `PolPower`: Calculates `u^n` and `n*v`, storing them in `u1` and `v1`. It also normalizes `v1` to be within `0` and `2*PI`.
        *   `PolRect`: Converts `(u1, v1)` back to rectangular `(x, y)`.
    *   Adds `A[n] * x` to `z1` (real part) and `A[n] * y` to `z2` (imaginary part) to form the partial sum.
4.  Restores the original `x` and `y` values before exiting the procedure.

**Internal Helper Subroutines:**
*   `ATAN(Numerateur, denominateur: DOUBLE; Var Phase: DOUBLE)`: Custom `atan2` function to return phase between `0` and `2*PI`.
*   `Power(x: DOUBLE; n: INTEGER): DOUBLE`: Calculates `x` to the power of `n` for real `x` and integer `n`.
*   `RectPol`: Converts rectangular coordinates `(x, y)` to polar coordinates `(u, v)`.
*   `PolRect`: Converts polar coordinates `(u, v)` to rectangular coordinates `(x, y)`.
*   `PolPower`: Calculates the `n`-th power of a complex number in polar form, resulting in `u1` and `v1`.
*   `RectPower`: Calculates the `n`-th power of a complex number in rectangular form by orchestrating calls to `RectPol`, `PolPower`, and `PolRect`.

**Inputs (Global Variables):**
*   `m` (`INTEGER`): The degree of the polynomial.
*   `A` (`Array[0..5] of DOUBLE`): The real coefficients `A[0]` to `A[m]` of the polynomial. (Note: The example `Coeff` procedure hardcodes coefficients for a polynomial of degree 5).
*   `x` (`DOUBLE`): The real part of the complex argument `Z`.
*   `y` (`DOUBLE`): The imaginary part of the complex argument `Z`.

**Outputs (Global Variables):**
*   `z1` (`DOUBLE`): The real component of the evaluated polynomial `P(Z)`.
*   `z2` (`DOUBLE`): The imaginary component of the evaluated polynomial `P(Z)`.

**Usage Notes/Sample:**
The example `cmplxser.pas` program's `Coeff` procedure defines the coefficients for `P(Z) = (1+Z)^5`. It prompts the user for the real and complex parts of the argument `Z`.

```
SAMPLE RUN:
Input the complex number as prompted:
  Real part    = 1
  Complex part = 0

Results are:
  Z1 =  32.0000
  Z2 =   0.0000
```
This output `Z1 = 32.0000` is correct for `P(1+0i) = (1+1)^5 = 2^5 = 32`.

**Real-world Applications:**
*   **Electrical Engineering**: Analyzing AC circuits, impedance calculations, frequency response of systems, and filter design where transfer functions or signals are complex-valued polynomials.
*   **Signal Processing**: Fundamental for Fourier analysis, digital filter design, and transformations in the frequency domain, where signals are often represented by complex exponentials.
*   **Quantum Mechanics**: Solving various equations (e.g., Schrödinger's equation) that involve complex wave functions and operators.
*   **Control Systems**: Stability analysis of linear systems using techniques like root locus or Nyquist plots, which operate in the complex s-plane.
*   **Applied Mathematics**: General complex analysis, conformal mapping, and evaluation of analytic functions.

**Limitations:**
*   The `A` array size is hardcoded to `Array[0..5]`, limiting the polynomial degree to 5 in the example.
*   Small non-zero values (`1e-30`) are used to guard against division by zero in phase calculations (`RectPol`), which can affect precision for very small `x` or `y` values.
*   The coefficients are assumed to be real.

### 5. HORNER (Horner's Rule for Polynomial Shifting)

*   **Program File**: `horner.pas`
*   **Explanation File**: `horner.txt`

**Description:**
The `HORNER` program demonstrates Horner's rule, a method for efficiently evaluating and restructuring polynomials. Specifically, it implements the `Horner_Shift` subroutine to change the expansion point of a polynomial from `x = 0` to a new point `x0`. This transformation converts a polynomial of the form `y(x) = a_0 + a_1x + a_2x^2 + ...` into `y(x) = b_0 + b_1(x-x_0) + b_2(x-x_0)^2 + ...`. This re-expression can significantly minimize round-off errors and improve numerical stability when evaluating the polynomial near `x_0`, especially if the polynomial's value approaches zero at that point.

**Mathematical Background:**
Horner's rule, also known as Horner's method or nested multiplication, provides an algorithm for polynomial evaluation that minimizes the number of multiplications and additions. When applied to shifting the expansion point, it leverages principles similar to synthetic division. For an Nth-degree polynomial `y(x) = Σ_{i=0}^{N} a_i x^i`, the new coefficients `b_j` for the expansion `y(x) = Σ_{j=0}^{N} b_j (x-x_0)^j` can be found iteratively. The algorithm effectively performs repeated division by `(x - x_0)`.

The procedure uses a 2D array `C` for intermediate calculations following the definition:
`b_n^(-1) = a_n` for `n = 0, ..., N` (where `N` is the degree of the polynomial, 4 in this implementation).
Then, for `m = 0, ..., N`:
`b_0^m = b_0^(m-1)`
`b_j^m = x_0 * b_{j-1}^m + b_j^(m-1)` for `j = 1, ..., N-m`.
Finally, `b_m = b_m^(N-m)`.

**Algorithm (`Horner_Shift` procedure in `horner.pas`):**
The `Horner_Shift` procedure provided is specifically tailored for quartic polynomials (degree 4). It uses a 2D array `C` for intermediate calculations.
1.  **Initialization**: The first column of `C` (`C[j,0]`) is initialized with the input coefficients `A[i]` in reverse order of their powers (i.e., `C[0,0] = A[4]`, `C[1,0] = A[3]`, ..., `C[4,0] = A[0]`).
2.  **Iterative Calculation**: It iterates `i` from 0 to 4 (representing steps of the synthetic division process):
    *   The first element in the next column `C[0,i+1]` is set equal to `C[0,i]`.
    *   For `j` from 1 up to `4-i`, it calculates `C[j,i+1] := x0 * C[j-1,i+1] + C[j,i]`. This is the core recurrence for the coefficients.
3.  **Extraction of Shifted Coefficients**: The final shifted coefficients `B[i]` are extracted from the diagonal elements of the `C` array (specifically `B[4-i] := C[i,4-i+1]`).

**Inputs (Global Variables):**
*   `A` (`Array[0..10] of DOUBLE`): The original coefficients `A[0]` to `A[4]` of the quartic polynomial. `A[i]` is the coefficient of `x^i`.
*   `x0` (`DOUBLE`): The new expansion point.

**Outputs (Global Variable):**
*   `B` (`Array[0..10] of DOUBLE`): The shifted coefficients `B[0]` to `B[4]` of the polynomial with respect to `(x - x0)`. `B[i]` is the coefficient of `(x-x0)^i`.

**Usage Notes/Sample:**
The program prompts for five coefficients `A(0)` to `A(4)` and the expansion point `x0`. It then outputs the shifted coefficients `B(i)`.

```
SAMPLE RUN:
Input the five coefficients:
   A( 0) =  1
   A( 1) =  4
   A( 2) =  6
   A( 3) =  4
   A( 4) =  1
What is the expansion point: 1

The shifted coefficients are:
   B( 0) =  16
   B( 1) =  32
   B( 2) =  24
   B( 3) =  8
   B( 4) =  1
```
This example shows `y(x) = 1 + 4x + 6x^2 + 4x^3 + x^4 = (x+1)^4`. Shifting around `x0 = 1` yields `y(x) = ((x-1)+2)^4`, whose expansion is `16 + 32(x-1) + 24(x-1)^2 + 8(x-1)^3 + 1(x-1)^4`, which matches the output coefficients.

**Real-world Applications:**
*   **Numerical Stability**: Crucial in scientific and engineering simulations where polynomial evaluations must maintain high accuracy, particularly when evaluating functions near their roots or critical points where direct evaluation could lead to catastrophic cancellation errors.
*   **Computer-Aided Design (CAD)**: Efficient evaluation of splines and Bezier curves, which are typically defined by polynomials, improving rendering speed and precision in design software.
*   **Compiler Optimization**: Compilers can use Horner's rule to optimize polynomial evaluations, leading to more efficient machine code and faster program execution.
*   **Root Finding Algorithms**: Forms the basis for algorithms like Newton-Raphson for finding polynomial roots, as the coefficients of the Taylor expansion (which is what shifting creates) are directly related to the derivatives of the polynomial.

**Limitations:**
*   The current implementation (`horner.pas`) of `Horner_Shift` is hardcoded for quartic polynomials (degree 4). Extending it to higher degrees would require adjusting array bounds and loop limits.
*   No input validation is performed to ensure the correct number of coefficients or valid numerical inputs.

### 6. INVNORM (Inverse Normal Distribution)

*   **Program File**: `invnorm.pas`
*   **Explanation File**: `invnorm.txt`

**Description:**
The `INVNORM` program provides the `Inverse_Normal` subroutine, which calculates an approximation to the inverse of the complementary error function. This is closely related to finding the inverse of the standard normal distribution's cumulative tail probability. Given a probability `y` (representing `P(Z > X)`), it returns the value `x` (often called the z-score) such that the area under the standard normal curve from `x` to infinity is equal to `y`. This is a fundamental operation in statistical analysis.

**Mathematical Background:**
The complementary error function `erfc(x)` is related to the standard normal distribution tail probability `Q(x)`:
```
Q(x) = (1/sqrt(2pi)) * Integral_from_x_to_inf (e^(-t^2/2) dt) = 1/2 * erfc(x/sqrt(2))
```
The program uses a rational polynomial approximation, specifically chosen for its accuracy within the range `0 < y <= 0.5`. The formula, adapted from Abramowitz and Stegun, is:
```
Xo = t - (c0 + c1*t + c2*t^2) / (1 + d1*t + d2*t^2 + d3*t^3)
```
where `t = sqrt(-LN(y * y))` (which simplifies to `t = sqrt(ln(1/y^2))`).
The constants `c0, c1, c2, d1, d2, d3` are predefined empirical coefficients that optimize the approximation. The error in this approximation is stated as `|E(Q)| < 0.0005` (the error is referenced to `Q`, not `Xo`).

**Algorithm (`Inverse_Normal` procedure in `invnorm.pas`):**
1.  **Coefficient Definition**: Defines the fixed empirical coefficients (`c0, c1, c2, d1, d2, d3`).
2.  **Input Validation/Edge Case**: Checks if `y <= 0`. If true, it sets `x` to a very large number (`1e13`), approximating positive infinity, as the tail probability approaches zero.
3.  **Transformation**: Calculates the intermediate variable `z` (denoted `t` in the mathematical background) using `z := SQRT(-LN(y * y))`.
4.  **Rational Polynomial Evaluation**: Computes the numerator and denominator of the rational polynomial based on `z`.
5.  **Final Calculation**: Calculates `x` using the inverse normal approximation formula `x := z - (c0 + c1 * z + c2 * z * z) / (1.0 + d1 * z + d2 * z * z + d3 * z * z * z)`.

**Inputs (Global Variable):**
*   `y` (`DOUBLE`): The tail probability `Q(x)` (i.e., `P(Z > x)`). The valid range for this approximation is `0 < y <= 0.5`.

**Outputs (Global Variable):**
*   `x` (`DOUBLE`): The corresponding value `X` (z-score) such that the probability of a standard normal variable being greater than `X` is `y`.

**Usage Notes/Sample:**
The main program generates a table of `P(Z>X)` values (from `0.5` down to approximately `0`) and their corresponding `X` values.

```
SAMPLE RUN (excerpt):
 P(Z>X)     X
----------------
 0.50    0.0000
 0.48    0.0500
 0.46    0.1002
...
 0.02    2.0542
-0.00    INF.
```

**Real-world Applications:**
*   **Statistical Analysis**: Finding critical values (e.g., z-scores) for hypothesis tests, constructing confidence intervals, and performing power analysis in various statistical studies.
*   **Financial Risk Management**: Essential for quantitative finance, particularly in option pricing models (like Black-Scholes), Value-at-Risk (VaR) calculations, and credit risk modeling, where probabilities of extreme events are critical.
*   **Quality Control**: Determining control limits for processes that follow a normal distribution to monitor and improve product quality.
*   **Machine Learning**: Used in various probabilistic models and algorithms that rely on normal distributions, such as Gaussian Mixture Models or Bayesian inference with Gaussian priors.

**Limitations:**
*   The approximation is valid only for `0 < y <= 0.5`. For `y > 0.5`, the symmetry of the normal distribution can be used (`x_inverse(y) = -x_inverse(1-y)`).
*   The stated accuracy is `|E(Q)| < 0.0005`, which might not be sufficient for all applications requiring extremely high precision.
*   The program defines `x := 1e13` for `y <= 0`, indicating a practical "infinity", but doesn't explicitly handle negative `y` as an invalid input in a strict error-checking sense.

### 7. RECIPRO (Polynomial Inversion)

*   **Program File**: `recipro.pas`
*   **Explanation File**: `recipro.txt`

**Description:**
The `RECIPRO` program provides the `Reciprocal` procedure, which computes the coefficients of a polynomial `Q(x)` that approximates the reciprocal of another polynomial `P(x)`. This means that `P(x) * Q(x)` should ideally approximate `1`. This routine is useful for finding inverse relationships or for implementing polynomial division in a series context.

**Mathematical Background:**
Given an input polynomial `P(x) = a_0 + a_1x + a_2x^2 + ... + a_n x^n` and an output reciprocal polynomial `Q(x) = b_0 + b_1x + b_2x^2 + ... + b_m x^m`, the condition `P(x)Q(x) = 1` implies that the coefficients of all powers of `x` (except `x^0`) in the product `P(x)Q(x)` must be zero. By equating coefficients of like powers of `x`, a system of linear equations for `b_i` is formed:
*   `x^0`: `a_0 * b_0 = 1`  => `b_0 = 1 / a_0` (provided `a_0 <> 0`)
*   `x^1`: `a_0 * b_1 + a_1 * b_0 = 0` => `b_1 = -(a_1 * b_0) / a_0`
*   `x^2`: `a_0 * b_2 + a_1 * b_1 + a_2 * b_0 = 0` => `b_2 = -(a_1 * b_1 + a_2 * b_0) / a_0`
*   ...
*   `x^k`: `a_0 * b_k + a_1 * b_{k-1} + ... + a_k * b_0 = 0`
    This general recurrence relation can be rearranged to solve for `b_k`:
    `b_k = -(1/a_0) * Σ_{j=1}^{k} a_j * b_{k-j}` (for `k > 0`)

The procedure iteratively solves for `b_k` starting from `b_0`.

**Algorithm (`Reciprocal` procedure in `recipro.pas`):**
1.  **Normalization**: The procedure first normalizes the coefficients of `P(x)` by dividing all `A[i]` coefficients by `A[0]` (stored in `l`). This effectively makes the new `A[0]` equal to 1, simplifying subsequent calculations. `B[0]` is initialized to 1 (corresponding to `1/A[0]` after normalization).
2.  **Array Clearing**: Initializes relevant parts of the `A` and `B` arrays to zero for degrees beyond `n` to ensure correct summation.
3.  **Coefficient Calculation (`B` array)**: Iterates `i` from 1 to `m` (the desired degree of `Q(x)`):
    *   Initializes `B[i]` to 0.
    *   Enters an inner loop that calculates `B[i] := B[i] - A[j] * B[i - j]` for `j` from 1 to `i`. This implements the summation `Σ_{j=1}^{k} a_j * b_{k-j}`.
    *   Since `A[0]` was normalized to 1, the division by `A[0]` in the formula `b_k = -(1/a_0) * ...` becomes simply `b_k = -( ... )`.
4.  **Un-normalization**: After calculating all `B[i]` coefficients in the normalized domain, it multiplies the `A[i]` coefficients back by `l` and divides the `B[i]` coefficients by `l` to restore their original scaling relative to `P(x)`.

**Inputs (Global Variables):**
*   `n` (`INTEGER`): The degree of the input polynomial `P(x)`.
*   `m` (`INTEGER`): The desired degree of the reciprocal polynomial `Q(x)`.
*   `A` (`Array[0..10] of DOUBLE`): The coefficients `A[0]` to `A[n]` of the input polynomial `P(x)`.

**Outputs (Global Variable):**
*   `B` (`Array[0..10] of DOUBLE`): The coefficients `B[0]` to `B[m]` of the inverted polynomial `Q(x)`.

**Usage Notes/Sample:**
The program prompts for the degree of the input polynomial (`n`), the desired degree of the inverted polynomial (`m`), and then the coefficients `A(0)` to `A(n)`.

```
SAMPLE RUN:
What is the degree of the input polynomial: 2
What is the degree of the inverted polynomial: 6
Input the polynomial coefficients:
   A( 0) = ? 1
   A( 1) = ? .1
   A( 2) = ? .01

The inverted polynomial coefficients are:
   B( 0) =  1.000000
   B( 1) = -0.100000
   B( 2) =  0.000000
   B( 3) =  0.001000
   B( 4) = -0.000100
   B( 5) =  0.000000
   B( 6) =  0.000001
```
For the given input `P(x) = 1 + 0.1x + 0.01x^2`, the output `Q(x)` coefficients are calculated.

**Real-world Applications:**
*   **Control Systems**: Finding the inverse of system transfer functions, which are often rational polynomials, to design controllers for desired system behavior (e.g., feedforward control).
*   **Digital Signal Processing**: Designing inverse filters for deconvolution (e.g., removing blurring from images or echoes from audio signals) where the inverse of a system's response is needed.
*   **Numerical Algebra**: Used in some iterative methods for solving linear equations or approximating matrix inverses in numerical linear algebra.
*   **Approximation Theory**: Constructing rational function approximations by finding the reciprocal of a given polynomial.
*   **Model Identification**: In situations where a system's output is related to its input by a polynomial, and the inverse relationship is desired to determine input from observed output.

**Limitations:**
*   **Critical Requirement**: `A[0]` (the constant term of `P(x)`) must not be zero. If `A[0] = 0`, a divide-by-zero error will occur, as the normalization step relies on it. This means `P(x)` must not have a root at `x = 0`.
*   The infinite series for `Q(x)` will not converge for `x` values that correspond to roots of `P(x)`. The truncated approximation to `Q(x)` behaves like a truncated Taylor series—its error tends to grow with increasing `x`. Therefore, it is wise to empirically check the range of validity of the calculated `Q(x)` before extensive use.
*   The arrays `A` and `B` are fixed to a maximum size of 10 (`Array[0..10] of DOUBLE`), limiting the maximum degree of the polynomials.

### 8. REVERSE (Series Reversion)

*   **Program File**: `reverse.pas`

**Description:**
The `REVERSE` program provides the `Reverse_coeff` procedure, which performs series reversion. Given a power series `Y = A[0] + A[1] * X + A[2] * X^2 + ...`, this routine computes the coefficients of the inverse power series `X = B[0] + B[1] * Y + B[2] * Y^2 + ...`. This is particularly useful when it is difficult or impossible to find an explicit inverse function analytically, but the original function's power series expansion is known.

**Mathematical Background:**
Series reversion involves expressing the independent variable `X` as a power series in `Y`. The coefficients `B[i]` of the reverted series are determined by substituting the `X` series into the original `Y` series and equating coefficients of powers of `Y`. This can be derived from methods such as the Lagrange Inversion Theorem. The `Reverse_coeff` procedure directly implements the algebraic formulas for `B[i]` derived from these equations.

A critical condition for the reversion to be well-behaved at `X=0` (or `Y=A[0]`) is that `A[1]` (the coefficient of the linear term `X` in the input polynomial) must be non-zero. If `A[1] = 0`, the inverse function's derivative at `Y=A[0]` is undefined or infinite, making a simple power series expansion at that point impossible.

The implemented formulas are explicit for the first few coefficients `B[i]` (up to degree 7). For example:
*   `B_1 = 1 / A_1`
*   `B_2 = -A_2 / A_1^3`
*   `B_3 = (2 * A_2^2 - A_1 * A_3) / A_1^5`
*   ...and so on for higher terms, becoming progressively more complex.

**Algorithm (`Reverse_coeff` procedure in `reverse.pas`):**
The procedure directly calculates the coefficients `B[1]` through `B[7]` using the explicit formulas mentioned above.
1.  **Check `A[1]`**: It first checks if `A[1]` is zero. If it is, a "Divide zero error" message is printed, and the procedure exits, as the reversion is not possible.
2.  **Calculate `B[1]`**: `B[1]` is set to `1.0 / A[1]`. Intermediate variables `aa` and `bb` are initialized based on `1/A[1]` for efficiency in subsequent calculations.
3.  **Calculate `B[2]` to `B[7]`**: Each `B[i]` (for `i` from 2 to 7) is calculated using its respective direct algebraic formula, which involves various `A[j]` coefficients and powers of `aa` and `bb`.
4.  **Calculate `B[0]`**: After `B[1]` to `B[7]` are found, `B[0]` is calculated as `B[0] = -sum(B[i] * A[0]^i for i=1 to 7)`. This accounts for the constant term `A[0]` in the original series to ensure `X=B(Y)` gives the correct `X` value when `Y` is given.
5.  **Small Value Handling**: Coefficients `B[i]` that are very close to zero (less than `1e-20` in absolute value) are set to `0.0` for cleaner output.

**Inputs (Global Variables):**
*   `n` (`INTEGER`): The degree of the input polynomial `Y(X)`.
*   `A` (`Array[0..10] of DOUBLE`): The coefficients `A[0]` to `A[n]` of the input polynomial `Y(X)`. `A[1]` must not be zero.

**Outputs (Global Variable):**
*   `B` (`Array[0..10] of DOUBLE`): The coefficients `B[0]` to `B[7]` of the reverted polynomial `X(Y)`.

**Usage Notes/Sample:**
The program prompts for the degree of the input polynomial (`n`) and its coefficients `A(0)` to `A(n)`. It then outputs the calculated reverted polynomial coefficients `B(i)`.

```
SAMPLE RUN:
What is the degree of the input polynomial: 3
Input the coefficients as prompted:
   A( 0) = ? 1
   A( 1) = ? 1
   A( 2) = ? 1
   A( 3) = ? 1

The reversed polynomial coefficients are:
   B( 0) =  9
   B( 1) =  1
   B( 2) =  0
   B( 3) = -1
   B( 4) =  0
   B( 5) =  3
   B( 6) =  0
   B( 7) = -12
```
For `Y = 1 + X + X^2 + X^3`, the expected reversion around `A[0]=1` would typically yield `X = (Y-1) - (Y-1)^2 + 2(Y-1)^3 - ...`. The `B[0]=9` in the sample is unusual for a standard series reversion starting at `A[0]`. The `B[0]` calculation in the code (`B[0] := B[0] - B[i] * aa; aa := aa * A[0]`) suggests it is calculating `B(0) = - (B_1 A_0 + B_2 A_0^2 + ... + B_7 A_0^7)`, which results in a constant term for the reversed series when evaluated at `Y=0`. This is a specific interpretation.

**Real-world Applications:**
*   **Solving Implicit Equations**: If an equation is given in the form `f(X) = Y` where `f` is a power series, series reversion can approximate `X` as a function of `Y`, which is invaluable when `f` is difficult to invert directly.
*   **Function Inversion**: Useful for inverting transcendental functions or special functions that have known power series expansions but are difficult to invert directly (e.g., in some branches of physics, engineering, or applied mathematics).
*   **Perturbation Theory**: In physics and engineering, sometimes a complex relationship can be expressed as a series expansion, and its inverse is needed for perturbation analysis or sensitivity studies.
*   **Numerical Methods**: Can be used as part of iterative schemes to find roots or fixed points of equations.

**Limitations:**
*   **Critical Requirement**: `A[1]` must be non-zero. This is a fundamental mathematical condition for series reversion at the expansion point.
*   The degree of reversion is hardcoded and limited to seven. This means only `B[0]` through `B[7]` are calculated, and extending it to higher degrees would require deriving and implementing more complex explicit formulas.
*   No input validation for the degree `n` or other coefficients is performed.
*   The `A` and `B` arrays are fixed to a maximum size of 10.

---

## General Usage Notes

*   **Data Types**: All programs primarily utilize standard Pascal data types (`DOUBLE` for floating-point numbers, `INTEGER` for degrees, counts, and indices).
*   **Input/Output**: Input values are typically read from the console (e.g., using `read` or `readln`), and results are printed to the console (e.g., using `Writeln`).
*   **Global Variables**: Many subroutines utilize global variables for inputs and outputs (e.g., `x`, `y`, `m`, `A`, `B`, `C`). Users should be aware of this scope when integrating these procedures into larger applications.
*   **Accuracy**: The numerical accuracy of each routine is discussed in its respective explanation file and within this documentation, often referencing the original sources. Users should consult these notes to understand the limitations and appropriate use cases.
*   **Error Handling**: Error handling for invalid inputs (e.g., negative degrees of freedom for Chi-Square, `A[0]=0` for `RECIPRO`, `A[1]=0` for `REVERSE`) is generally limited or noted as the responsibility of the calling program or user input.

## References

This module is based on the `BASIC Scientific Subroutines, Vol. II` by F.R. Ruckdeschel, published by BYTE/McGraw-Hill in 1981 [1]. The Pascal versions were created by J.-P. Moreau (www.jpmoreau.fr). Specific mathematical references are noted in the individual explanation files and summarized below:

*   [1] Ruckdeschel, F.R. *BASIC Scientific Subroutines, Vol. II*. BYTE/McGraw-Hill, 1981.
*   Abramowitz, M., and Stegun, I.A. *Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables*. National Bureau of Standards, Applied Mathematics Series 55, 1964. (Referenced for `INVNORM`).
*   Fike, C.T. *Computer Evaluation of Mathematical Functions*. Englewood Cliffs, NJ: Prentice-Hall, 1968. (Referenced for rational polynomials).
*   Henrici, P. *Computational analysis*. Wiley, 1982. (Referenced for `RECIPRO`).
*   Peirce, B.O. *A Short Table of Integrals*. Ginn and Company, 1957. (Referenced for `ASYMERF`).
*   *CRC Standard Mathematical Tables*, 24th Edition. (Referenced for `LN_FACTX` and `REVERSE`).
*   Texas Instruments SR-51 Owners Manual, 1974. (Referenced for `CHI_SQR`).
*   Hewlett-Packard statistics programs, 1974. (Referenced for `CHI_SQ`).
