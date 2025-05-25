
# jpmMath - Statistical Functions

## Table of Contents

1.  [About This Library](#about-this-library)
2.  [Special Mathematical Functions](#special-mathematical-functions)
    *   [`gamma.pas` - Gamma Function Demonstration](#gammapas---gamma-function-demonstration)
    *   [`ibeta.pas` - Incomplete Beta Function Ix(a,b)](#ibetapas---incomplete-beta-function-ixa-b)
3.  [Probability Distributions](#probability-distributions)
    *   [`distri.pas` - Statistical Distributions Calculator](#distripas---statistical-distributions-calculator)
    *   [`chi2.pas` - Chi-square and Inverse Chi-square Functions](#chi2pas---chi-square-and-inverse-chi-square-functions)
    *   [`fdistri.pas` - F Distribution Q(F|nu1,nu2)](#fdistripas---f-distribution-qfnu1nu2)
    *   [`normal.pas` - Normal and Inverse Normal Probability Functions](#normalpas---normal-and-inverse-normal-probability-functions)
    *   [`tnormal.pas` - Standardized Normal Law Probabilities](#tnormalpas---standardized-normal-law-probabilities)
    *   [`student.pas` - Student T-Probability Law](#studentpas---student-t-probability-law)
4.  [Descriptive Statistics & Moments](#descriptive-statistics--moments)
    *   [`fstat.pas` - Statistical Functions for One or Two Variables](#fstatpas---statistical-functions-for-one-or-two-variables)
    *   [`moment.pas` - Means and Moments of a Statistical Variable with Frequencies](#momentpas---means-and-moments-of-a-statistical-variable-with-frequencies)
    *   [`momnts.pas` - Mean and First Third Moments of a Set of Data (Incremental)](#momntspas---mean-and-first-third-moments-of-a-set-of-data-incremental)
    *   [`tmoment.pas` - Statistical Moments of a Distribution (Comprehensive)](#tmomentpas---statistical-moments-of-a-distribution-comprehensive)
5.  [Median Calculation](#median-calculation)
    *   [`tmdian.pas` - Median Value with Heapsort Method](#tmdianpas---median-value-with-heapsort-method)
    *   [`tmdian1.pas` - Median Value Without Sorting](#tmdian1pas---median-value-without-sorting)
6.  [Regression Analysis](#regression-analysis)
    *   [`simlp.pas` - Simple Linear Regression Minimizing Sum of Absolute Deviation (SAD)](#simlppas---simple-linear-regression-minimizing-sum-of-absolute-deviation-sad)
7.  [Probability Expression Evaluation](#probability-expression-evaluation)
    *   [`prob.pas` - Evaluate Probabilities by Parsing String Expression](#probpas---evaluate-probabilities-by-parsing-string-expression)
8.  [Statistical Unit (`Stats.pas`)](#statistical-unit-statspas)
    *   [`stats.pas` - Unit for Basic Statistical Functions](#statspas---unit-for-basic-statistical-functions)
    *   [`tstats.pas` - Driver Program to Test Unit STATS](#tstatspas---driver-program-to-test-unit-stats)
9.  [References](#references)

---

## 1. About This Library

The `stat` folder contains a suite of Pascal programs and a reusable unit for performing various statistical computations. These tools are designed to demonstrate fundamental statistical algorithms and their practical applications, offering a foundational toolkit for quantitative analysis in scientific and engineering domains.

The programs are characterized by:
*   **Modularity**: Each program typically addresses a specific statistical concept or calculation.
*   **Clarity**: The source code is well-commented, emphasizing the underlying mathematical and algorithmic principles.
*   **Utility**: Providing core functions for statistical distributions, moments, regression, and data description.

This library serves as a valuable resource for learning statistics, inspecting classical algorithms, and as a starting point for building more complex statistical applications.

## 2. Special Mathematical Functions

This section details programs that implement fundamental mathematical functions often utilized as building blocks in statistical computations.

### `gamma.pas` - Gamma Function Demonstration

*   **Description**: This program demonstrates the calculation of the Gamma function, $\Gamma(x)$, for a range of positive real numbers. The Gamma function is a mathematical function that extends the concept of the factorial function to real and complex numbers. For positive integers $n$, $\Gamma(n) = (n-1)!$.
*   **Core Functionality**:
    *   `Function Gamma(xx:double):double`: Computes $\Gamma(x)$ for $x > 0$ using an approximation based on the natural logarithm of the Gamma function, `EXP(tmp+LN(stp*ser))`, which provides robust calculations.
*   **Usage**: The program iteratively calculates and prints the value of `Gamma(X)` for `X` ranging from `0.5` to `5.0` in increments of `0.5`. No user input is required during execution.
*   **Sample Run**:
    ```
          X        Gamma(X)
     -------------------------
        0.5000      1.772454
        1.0000      1.000000
        1.5000      0.886227
        2.0000      1.000000
        2.5000      1.329340
        3.0000      2.000000
        3.5000      3.323351
        4.0000      6.000000
        4.5000     11.631728
        5.0000     24.000000
    ```
*   **Real-world Applications**:
    *   **Probability Theory**: The Gamma function is a cornerstone of many continuous probability distributions, including the Gamma distribution itself (modeling waiting times), the Chi-square distribution, Student's t-distribution, and the Beta distribution.
    *   **Mathematical Physics and Engineering**: It appears in diverse fields, including quantum mechanics, fluid dynamics, and various integral transforms, where developers might integrate this function into larger simulations or analytical tools.
    *   **Special Functions**: It is a fundamental building block for many other special functions in applied mathematics.
*   **Reference**: "Numerical Recipes, By W.H. Press, B.P. Flannery, S.A. Teukolsky and T. Vetterling, Cambridge University Press, 1986" [BIBLI 08].

### `ibeta.pas` - Incomplete Beta Function Ix(a,b)

*   **Description**: This program calculates the regularized incomplete Beta function, denoted as $I_x(a,b)$. This function is essential for calculating cumulative probabilities of various statistical distributions and is a generalization of the Beta function.
*   **Core Functionality**:
    *   `Function GAMMLN(XX:Double): Double`: Computes $\ln(\Gamma(x))$ for $x > 0$. This function is a prerequisite for `BETAI` and is identical to the one found in `fdistri.pas` and `gamma.pas`.
    *   `Function BETACF(A,B,X:Double): Double`: Calculates the continued fraction for the incomplete Beta function. This is an internal helper function used by `BETAI` for numerical stability and convergence.
    *   `Function BETAI(A,B,X:Double): Double`: The primary function that computes $I_x(a,b)$. It leverages the `GAMMLN` and `BETACF` functions and employs a symmetry transformation to ensure accurate calculation across the domain $0 \le x \le 1$.
*   **Usage**: The program computes `Ix(a,b)` for pre-defined values of `a`, `b`, and `x` (specifically `a=0.50`, `b=5.0`, `x=0.20`). No direct user input is requested during execution.
*   **Sample Run**:
    ```
     Incomplete Beta Function Ix(A,B)

     A = 5.00000000000000E-0001
     B = 5.00000000000000E+0000
     x = 2.00000000000000E-0001

     Y = 8.55072399728097E-0001
    ```
*   **Real-world Applications**:
    *   **Cumulative Distribution Functions (CDFs)**: Crucial for calculating CDFs for:
        *   **Beta Distribution**: Modeling probabilities of probabilities, especially in Bayesian statistics (e.g., inferring success rates).
        *   **F-distribution**: As seen in `fdistri.pas`, the incomplete Beta function is used to calculate F-distribution probabilities for ANOVA and regression.
        *   **Student's t-distribution**: Probabilities for the t-distribution can also be expressed in terms of the incomplete Beta function.
    *   **Statistical Hypothesis Testing**: Directly involved in computing p-values for statistical tests that utilize the F or t distributions.
    *   **Quality Control**: Assessing the reliability of systems or processes by modeling distributions of success/failure rates, where outcomes are bounded proportions.
*   **Reference**: "Numerical Recipes, By W.H. Press, B.P. Flannery, S.A. Teukolsky and T. Vetterling, Cambridge University Press, 1986" [BIBLI 08].

## 3. Probability Distributions

This section documents programs dedicated to calculating probabilities for various standard statistical distributions.

### `distri.pas` - Statistical Distributions Calculator

*   **Program Name**: `stat_distributions`
*   **Description**: This program provides a menu-driven interface to compute probabilities for six common statistical distributions: binomial, Poisson, normal (univariate and bivariate), chi-square, and Student's T. It serves as a unified calculator for a range of probability calculations.
*   **Core Functionality**:
    *   `Function PowerI(x:double; n:integer): double`: Calculates $x^n$ for integer $n$.
    *   `Function Power(y,x:DOUBLE):DOUBLE`: Calculates $y^x$ for real exponents.
    *   `Procedure Binomial(p, n, x, var fx,px,qx)`: Computes the probability mass function ($P(X=x)$), cumulative distribution function ($P(X \le x)$), and complementary cumulative distribution function ($P(X \ge x)$) for a binomial distribution.
    *   `Procedure Poisson(xm, x, var fx,px,qx)`: Computes the probability mass function, CDF, and complementary CDF for a Poisson distribution.
    *   `Procedure normal(xm, s, x, var fx,px,qx)`: Computes the probability density function ($f(x)$), CDF ($P(X \le x)$), and complementary CDF ($P(X \ge x)$) for a univariate normal distribution.
    *   `Procedure Normal2(xm,ym,sx,sy,ro,x,y, var fxy)`: Computes the probability density function for a bivariate normal distribution, given means (`xm`, `ym`), standard deviations (`sx`, `sy`), and correlation coefficient (`ro`).
    *   `Procedure Chisqr(nu, x, var fx,px,qx)`: Computes the probability density function, CDF, and complementary CDF for a chi-square distribution.
    *   `Procedure Student(nu, x, var bt,px,qx)`: Computes $P(-X \le \text{random variable} \le X)$, CDF, and complementary CDF for Student's T distribution.
*   **Usage**: The program first presents a tutorial and a menu of distributions. The user selects a distribution by entering a number (1-6) and then provides the necessary parameters and the value(s) of the random variable for the chosen distribution.
*   **Sample Run (Normal distribution)**:
    ```
           Statistical distributions

     Tutorial

     1. Input distribution number:

        1: binomial distribution
        2: Poisson distribution
        3: normal distribution
        4: normal distribution (2 variables)
        5: chi-square distribution
        6: Student T distribution

     2. Define the parameters of chosen distribution

     3. Input value of random variable


     Input distribution number (1 to 6): 3

       Normal distribution

       MU=mean
       S =standard deviation
       X =value of random variable

       MU = 2
       S  = 3
       X  = 5

       Probability of random variable  = X: .0806569
       Probability of random variable <= X: .8413447
       Probability of random variable >= X: .1586552
    ```
*   **Real-world Applications**:
    *   **Risk Assessment**: Quantifying the probability of events in various fields (e.g., number of defects in a manufacturing batch using Binomial, or number of customer calls per hour using Poisson).
    *   **Data Analysis**: Understanding the likelihood of observing specific data points given a known distribution (e.g., how likely is a student to score above a certain mark if grades are normally distributed?).
    *   **Statistical Modeling**: Serving as a foundational tool for building more complex statistical models, especially in finance (e.g., option pricing with Normal distribution, valuing correlated assets with Normal2) and quality control.
    *   **Hypothesis Testing**: Directly computing p-values for various statistical tests where the underlying data conforms to one of these distributions (e.g., t-tests for small samples, chi-square tests).
*   **Reference**: "Mathematiques et statistiques By H. Haut, PSI Editions, France, 1981" [BIBLI 13].

### `chi2.pas` - Chi-square and Inverse Chi-square Functions

*   **Program Name**: `Stats01`
*   **Description**: This program calculates the Chi-squared probability distribution function for a given chi-squared statistic (`X`) and degrees of freedom (`df`). It also provides the inverse function, finding the chi-squared value corresponding to a given probability and degrees of freedom.
*   **Core Functionality**:
    *   `FUNCTION poz (z:DOUBLE): Double`: Computes the probability of a standard normal z value. This is an auxiliary function used internally by `chi2`.
    *   `FUNCTION chi2 (X: DOUBLE; df: INTEGER): Double`: Calculates the cumulative probability for a given Chi-square value ($X$) and degrees of freedom ($df$). It employs Algorithm 299 from Collected Algorithms for the CACM.
    *   `FUNCTION invchi2 (P: DOUBLE; idf: INTEGER): Double`: Calculates the inverse Chi-square function, returning the Chi-square value ($X$) for a given cumulative probability ($P$) and degrees of freedom ($idf$). It uses an iterative bisection method.
*   **Usage**: The program prompts the user for the `Value of variable` (chi-squared statistic `X`) and `Degree of freedom` (`idf`). It then calculates and outputs the `Probability` and verifies the inverse function.
*   **Sample Run**:
    ```
      CHI2 Law:
      ========
    Value of variable: 100
    Degree of freedom: 100

    Probability =  4.81191684527957E-0001

    Verify Inverse Law:
      X =  1.00000000148592E+0002
    ```
*   **Real-world Applications**:
    *   **Hypothesis Testing**: This function is fundamental to chi-squared tests, such as:
        *   **Goodness-of-Fit Tests**: Determining if an observed frequency distribution significantly differs from an expected distribution (e.g., testing if dice rolls are fair, or if survey responses align with a theoretical proportion).
        *   **Tests of Independence**: Assessing if there is a statistically significant relationship between two categorical variables (e.g., is there an association between gender and preferred type of media, or vaccination status and disease incidence?).
    *   **Confidence Intervals**: Used to construct confidence intervals for population variance or standard deviation based on sample data, assuming the underlying population is normally distributed.
*   **Reference**: Adapted from JavaScript functions by John Walker, originally C implementations by Gary Perlman of Wang Institute, Tyngsboro. The reference URL is `http://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html`.

### `fdistri.pas` - F Distribution Q(F|nu1,nu2)

*   **Program Name**: `FDISTRI`
*   **Description**: This program calculates the F-distribution cumulative probability, denoted as $Q(F|\nu_1, \nu_2)$. This function is critical for evaluating the significance of F-statistics, commonly used in ANOVA (Analysis of Variance) and regression analysis.
*   **Core Functionality**:
    *   `FUNCTION gammln(xx:REAL): REAL`: Computes $\ln(\Gamma(x))$. (Shared utility function, also found in `gamma.pas` and `ibeta.pas`).
    *   `FUNCTION betacf(a,b,x:REAL):REAL`: Computes the continued fraction for the incomplete Beta function. (Shared utility function, also found in `ibeta.pas`).
    *   `FUNCTION betai(a, b, x: REAL): REAL`: Computes the incomplete Beta function $I_x(a,b)$. (Shared utility function, also found in `ibeta.pas`).
    *   The F-distribution probability is calculated using the relationship: $Q(F|\nu_1, \nu_2) = I_x(\nu_2/2, \nu_1/2, \nu_2/(\nu_2+\nu_1 F))$.
*   **Usage**: The program prompts the user for the `Ratio F of dispersion of 1st sample / dispersion of 2nd sample` (`F`), the `Degree of freedom nu1 of first sample` (`nu1`), and the `Degree of freedom nu2 of second sample` (`nu2`). It then calculates and outputs the cumulative probability `Q`.
*   **Sample Run**:
    ```
     Calculate F Distribution Function Q(F|nu1,nu2) for
     F, nu1 and nu2 given.

     Ratio F of dispersion of 1st sample / dispersion of 2nd sample: 1.25
     Degree of freedom nu1 of first  sample: 25
     Degree of freedom nu2 of second sample: 15


     Number of iterations: 6


     Q = 0.332373
    ```
*   **Real-world Applications**:
    *   **Analysis of Variance (ANOVA)**: Used to test for significant differences between the means of three or more groups (e.g., comparing the effectiveness of multiple teaching methods on student performance, or yields from different agricultural treatments).
    *   **Regression Analysis**: Evaluating the overall significance of a regression model by testing if the model explains a significant amount of variance in the dependent variable (e.g., is a set of predictors better than simply using the mean to predict an outcome?).
    *   **Comparing Variances**: Performing F-tests to compare the variances of two populations, which is important for assumptions in other statistical tests (e.g., t-tests).
*   **Reference**: "Numerical Recipes, By W.H. Press, B.P. Flannery, S.A. Teukolsky and T. Vetterling, Cambridge University Press, 1986, page 166" [BIBLI 08].

### `normal.pas` - Normal and Inverse Normal Probability Functions

*   **Program Name**: `Normal_Law`
*   **Description**: This program implements functions for the Standard Normal Probability Distribution, specifically the Probability Density Function (PDF) and its inverse (quantile function). It works with standardized normal values (mean=0, standard deviation=1).
*   **Core Functionality**:
    *   `FUNCTION PHI (u: DOUBLE): Double`: Computes the value of the standard normal PDF at a given point $u$. This describes the relative likelihood for a given random variable to take on a given value.
    *   `FUNCTION NORMAL (P: DOUBLE): Double`: Computes the inverse of the standard normal function. Given a probability `P`, it iteratively searches for the Z-score `X` such that `PHI(X)` (the probability density at `X`) equals `P`. While standard inverse normal functions typically find the Z-score for a given *cumulative* probability, this implementation appears to find the Z-score for a given *density* value `P`.
*   **Usage**: The program prompts the user for a `Value of variable` (a standardized normal `X` value or z-score). It then outputs the calculated `Probability` (density at `X`) and verifies the inverse function by finding the `X` value corresponding to the obtained probability.
*   **Sample Run**:
    ```
     Value of variable: 0.2

     Probability =  3.91042693975456E-0001

     Verify:
       X =  2.00000286102295E-0001
    ```
*   **Real-world Applications**:
    *   **Statistical Inference**: Calculating the likelihood of observing specific z-scores or finding the z-score corresponding to a particular probability density.
    *   **Quality Control**: Monitoring processes by calculating z-scores for measurements and interpreting their likelihood relative to a target distribution.
    *   **Research**: When combined with cumulative probability functions (like those in `tnormal.pas`), these functions are vital for hypothesis testing and constructing confidence intervals for normally distributed data.
*   **Reference**: This program includes functions commonly found in basic statistical libraries, focusing on direct implementation of probability and inverse probability for the standard normal distribution.

### `tnormal.pas` - Standardized Normal Law Probabilities

*   **Program Name**: `TNormal`
*   **Description**: This program calculates probabilities and probability densities for the standardized normal distribution (mean=0, standard deviation=1). It demonstrates three different algorithms for these calculations: `alnorm`, `NORMP`, and `NPROB`, each offering varying levels of accuracy or computational approaches.
*   **Core Functionality**:
    *   `Function alnorm(x:double; upper:boolean): Double`: Implements Algorithm AS66 from Applied Statistics. It evaluates the tail area of the standardized normal curve (from $x$ to infinity if `upper` is true, or from negative infinity to $x$ if `upper` is false).
    *   `Procedure NORMP(Z:Double; Var P, Q, PDF:Double)`: Implements Algorithm 5666 (from Hart et al.). Calculates the cumulative probability to the left of $Z$ ($P$), to the right of $Z$ ($Q$), and the probability density function (PDF).
    *   `Procedure NPROB(Z:Double; Var P,Q,PDF:Double)`: Implements Algorithm 39 (from Adams). Provides another method for calculating probabilities $P$, $Q$, and PDF for the standard normal distribution.
*   **Usage**: The program uses a hardcoded `X` value (e.g., `5.0`) representing a z-score for demonstration. It then calls each of the three functions/procedures (`alnorm`, `NORMP`, `NPROB`) to calculate probabilities and probability densities based on this `X` value, showing their respective outputs.
*   **Sample Run**:
    ```
      P1= 9.99999713348428E-0001
      P2= 2.86651571867398E-0007

      P1= 9.99999713348428E-0001
      P2= 2.86651571892333E-0007
      Q=  1.48671951473430E-0006

      P1= 9.99999713348428E-0001
      P2= 2.86651571867398E-0007
      Q=  1.48671951467306E-0006
    ```
*   **Real-world Applications**:
    *   **Statistical Inference**: Performing precise calculations for Z-tests (e.g., comparing a sample mean to a population mean when the population standard deviation is known).
    *   **Confidence Intervals**: Determining critical z-values for constructing confidence intervals for means or proportions, assuming a normal distribution.
    *   **Option Pricing**: The cumulative standard normal distribution function is a core component of the Black-Scholes model and other financial models for valuing options and other derivatives.
    *   **Quality Control**: Assessing the probability of a product falling outside specification limits, assuming a normal distribution of measurements.
    *   **Data Transformation**: Converting any normally distributed variable to a standard normal variable (z-score) allows for easy comparison across different scales and facilitates the use of standard normal tables or functions.
*   **References**:
    *   `alnorm`: "Algorithm AS66 Applied Statistics (1973) vol22 no.3".
    *   `NORMP`: "Hart, J.F. et al, 'Computer Approximations', Wiley 1968".
    *   `NPROB`: "Adams,A.G. Areas under the Normal Curve, ALGORITHM 39, COMPUTER J., VOL. 12, 197-8, 1969".

### `student.pas` - Student T-Probability Law

*   **Program Name**: `STUDENT`
*   **Description**: This program calculates probabilities associated with the Student's t-distribution. It provides two distinct functions for calculating tail areas: `PROBST` for the lower tail and `STUDNT` for the upper tail. The t-distribution is crucial for statistical inference when dealing with small sample sizes or when the population standard deviation is unknown.
*   **Core Functionality**:
    *   `FUNCTION Power(x:REAL;n:Integer):REAL`: Calculates $x^n$.
    *   `FUNCTION PROBST(T:REAL; IDF:Integer; Var IFAULT:Integer): REAL`: Implements Algorithm AS 3 from Applied Statistics. Calculates the lower tail probability of the Student's T-distribution for a given T-statistic and degrees of freedom ($IDF$).
    *   `FUNCTION STUDNT (T, DOFF:REAL; Var IFAULT:Integer): REAL`: Implements Algorithm AS 27 from Applied Statistics. Calculates the upper tail area under Student's T-distribution for a given T-statistic and degrees of freedom ($DOFF$).
*   **Usage**: The program uses hardcoded values for the t-statistic (`X=0.257`) and degrees of freedom (`NU=19`) for demonstration. It then calls both `PROBST` and `STUDNT` functions.
*   **Sample Run**:
    ```
      X= 0.2570000
      PROB1= 0.6000294
      ERROR=0
      X= 0.2570000
      PROB2= 0.3999705
      ERROR=0
      PROB1+PROB2= 0.9999998
    ```
*   **Real-world Applications**:
    *   **Hypothesis Testing (t-tests)**: Fundamentally used in t-tests for comparing means of two groups (independent or paired samples) or comparing a sample mean to a hypothesized population mean, particularly when sample sizes are small (`< 30`) and the population standard deviation is unknown. Examples include A/B testing in marketing with small sample sizes or clinical trials to compare drug efficacy.
    *   **Confidence Intervals**: Constructing confidence intervals for population means based on small samples, providing a range within which the true population mean is likely to fall.
    *   **Quality Control**: Assessing variations in small production batches where full population data is not available.
*   **Reference**: "JOURNAL OF APPLIED STATISTICS (1968) VOL.17, P.189, & VOL.19, NO.1".

## 4. Descriptive Statistics & Moments

This section includes programs focused on calculating various descriptive statistics that summarize key features of a dataset.

### `fstat.pas` - Statistical Functions for One or Two Variables

*   **Program Name**: `FStat`
*   **Description**: This program computes common statistical functions for a single set of data (e.g., $X_i$) or for a pair of related data sets ($X_i, Y_i$). It provides measures of central tendency, dispersion, and relationships between variables.
*   **Core Functionality**:
    *   `Procedure Stat_functions`: The main computational procedure that performs the statistical calculations. It distinguishes between single-variable and two-variable analysis based on user input.
        *   **For One Variable**: Calculates the mean of $X_i$, `(n-1)` standard deviation of $X_i$ (sample), `(n)` standard deviation of $X_i$ (population), and the standard error of the mean for both `(n-1)` and `(n)` denominators.
        *   **For Two Variables**: Calculates the means of $X_i$ and $Y_i$, `(n-1)` and `(n)` standard deviations for both $X_i$ and $Y_i$, standard deviations of their respective means, `(n-1)` and `(n)` covariance of $X,Y$, and the correlation coefficient.
*   **Usage**: The program first prompts the user to select the type of calculation: 1 for one variable, 2 for two variables. After selection, the user inputs the `Number of data` points (`n`), followed by the `n` values for `X(i)` (and `Y(i)` if type 2 is selected).
*   **Sample Run (Two variables)**:
    ```
     TUTORIAL

     1. Define type of calculus:"
        1: Statistical functions for a set X(i)
        2: Statistical functions for a set X(i), Y(i)

     2. Input number n of data

     3. Input successively the n values X(i) [and Y(i)]

     Type of calculus (1 or 2): 2

     Number of data: 5

       1  1 12
       2  2 9
       3  3 7
       4  4 15
       5  5 6


     Mean of X(i)....................: 3.00000000

     (n-1) standard deviation of X(i): 1.58113883
       (n) standard deviation of X(i): 1.41421356

     (n-1) standard dev. of X mean...: 0.70710678
       (n) standard dev. of X mean...: 0.63245553

     Mean of Y(i)....................: 9.80000000

     (n-1) standard deviation of Y(i): 3.70135111
       (n) standard deviation of Y(i): 3.31058908

     (n-1) standard dev. of Y mean...: 1.65529454
       (n) standard dev. of Y mean...: 1.48054045


     (n-1) covariance of X,Y.........: -1.50000000
       (n) covariance of X,Y.........: -1.20000000

     Correlation coefficient.........: -0.25630730
    ```
*   **Real-world Applications**:
    *   **Exploratory Data Analysis (EDA)**: Provides immediate insights into the characteristics of datasets, such as their central location (mean), spread (standard deviation), and linear relationship (covariance, correlation).
    *   **Process Monitoring**: In manufacturing, monitoring the mean and standard deviation of product dimensions to ensure quality control.
    *   **Financial Analysis**: Calculating the correlation coefficient between two assets to assess diversification benefits in a portfolio.
    *   **Research Studies**: Summarizing basic demographic or experimental data and identifying potential relationships between variables before more complex modeling.
*   **Reference**: "Mathematiques et statistiques by H. Haut, PSI Editions, France, 1981" [BIBLI 13].

### `moment.pas` - Means and Moments of a Statistical Variable with Frequencies

*   **Program Name**: `Moments`
*   **Description**: This program computes a comprehensive set of descriptive statistics for a statistical variable, including various means (arithmetic, geometric, harmonic), and the first four central moments (M1-M4). It also computes derived coefficients such as flatness (kurtosis) and asymmetry (skewness), accounting for data frequencies.
*   **Core Functionality**:
    *   `Function Power(y,x:DOUBLE):DOUBLE`: Auxiliary function to calculate $y^x$.
    *   `Procedure Calculate`: The main computational procedure that:
        *   Calculates necessary sums for weighted means and moments.
        *   Computes the arithmetic mean, geometric mean, harmonic mean (notes if undefined, e.g., when `X[i]=0`).
        *   Calculates the generalized mean $M(t)$ for a user-specified coefficient `t`.
        *   Determines the first four central moments ($M_1, M_2, M_3, M_4$).
        *   Computes the flatness coefficient (normalized $M_4$, related to kurtosis) and the coefficient of asymmetry (normalized $M_3$, related to skewness).
*   **Usage**: The program prompts the user for the `Number of data (n)` (count of distinct data points), then successively for `n` pairs of `X[i]` (data value) and `F[i]` (frequency of `X[i]`). Finally, it asks for `Calculate generalized mean for t` (a coefficient `t`).
*   **Sample Run**:
    ```
     TUTORIAL

     Means and moments of a statistical variable X[i]
     with frequency F[i]

     1. Input number n of data

     2. Input sucessively the n values X[i] and F[i]

     Number of data: 3

       1  1 6
       2  3 4
       3  5 8

     Calculate generalized mean for t= 2


     Arithmetic mean:  3.22222222

     Geometric mean :  2.61023904

     Harmonic mean  :  2.01492537

     Moments:

       M1=    3.22222222
       M2=    3.06172840
       M3=   -1.16323733
       M4=   12.56881570

     Flatness coefficient....:   1.34079084
     Coefficient of asymmetry:  -0.21712925

     Gen. mean M( 2): 3.66666667
    ```
*   **Real-world Applications**:
    *   **Financial Risk Analysis**: Evaluating the risk and return profiles of investments, where higher moments (skewness and kurtosis) provide insights into non-normal distribution characteristics (e.g., fat tails, asymmetry of returns).
    *   **Data Profiling**: Gaining a deeper understanding of the shape, symmetry, and peakedness of data distributions in various scientific and engineering datasets.
    *   **Signal Processing**: Characterizing the statistical properties of signals and noise.
    *   **Quality Control**: Monitoring process variations beyond just mean and standard deviation, helping to identify non-standard process behaviors.
    *   **Demographics**: Calculating weighted means (e.g., average income weighted by household size) or analyzing the distribution of wealth and its moments.
*   **Reference**: "Mathematiques et statistiques By H. Haut, PSI Editions, France, 1981" [BIBLI 13].

### `momnts.pas` - Mean and First Third Moments of a Set of Data (Incremental)

*   **Program Name**: `TMOMNTS`
*   **Description**: This program calculates the mean ($S1$) and the sums of powers of deviations ($S2, S3, S4$) from the mean for a set of data points `Y(i)`. This implementation is based on Algorithm AS 52, which is designed for numerical stability and incremental calculation of moments.
*   **Core Functionality**:
    *   `Procedure MOMNTS(X:REAL; K, N:Integer; Var S1, S2, S3, S4:REAL; Var IFAULT:Integer)`: This is the core incremental procedure. It adds a new value `X` to the calculation, updating the current mean ($S1$) and the sums of powers of deviations ($S2, S3, S4$).
        *   `S1`: Current mean.
        *   `S2`: Sum of $(Y_i - S1)^2$.
        *   `S3`: Sum of $(Y_i - S1)^3$.
        *   `S4`: Sum of $(Y_i - S1)^4$.
        *   `K`: Parameter controlling which moments are updated (1=mean only, 2=mean+S2, 3=mean+S2+S3, 4=mean+S2+S3+S4).
        *   `N`: The current count of observations (should be incremented for each call).
*   **Usage**: The program prompts the user for the `Number of data (ndata)`, then successively for `ndata` values of `Y[i]`. It then iteratively calls `MOMNTS` for each data point.
*   **Outputs**: The final values of $S1$ (mean), $S2$ (sum of squared deviations), $S3$ (sum of cubed deviations), and $S4$ (sum of fourth power deviations), along with an `Error code`.
*   **Sample Run**:
    ```
      Number of data: 5
      1: 12
      2: 9
      3: 7
      4: 15
      5: 6

      S1=  9.80000000000291E+0000
      S2=  5.48000000000466E+0001
      S3=  7.39200000003912E+0001
      S4=  1.02497600000165E+0003

      Error code: 0
    ```
*   **Real-world Applications**:
    *   **Online/Streaming Data Analysis**: This algorithm's incremental nature makes it suitable for calculating moments of large datasets or data streams without needing to store all data points simultaneously. This is crucial for applications where data arrives continuously (e.g., sensor data, financial tick data, network traffic analysis).
    *   **Numerical Stability**: The algorithm is designed to reduce the impact of catastrophic cancellation or precision loss that can occur when calculating moments from raw sums, especially with large numbers or very small variances.
    *   **Statistical Process Control**: Continuously monitoring the mean and variability of a process, as well as detecting shifts in distribution shape, by updating moments as new observations become available.
*   **Reference**: "Journal of Applied Statistics (1972) vol.21, page 226".

### `tmoment.pas` - Statistical Moments of a Distribution (Comprehensive)

*   **Program Name**: `TMOMENT`
*   **Description**: This program calculates various statistical moments and related measures for a given set of data. These include the average (mean), average deviation, standard deviation, variance, skewness, and kurtosis, providing a comprehensive description of the data's central tendency, dispersion, and shape.
*   **Core Functionality**:
    *   `Procedure Moment(data:pVEC; n:Integer; Var ave,adev,sdev,var0,skew,curt:Double)`: The primary procedure that computes:
        *   `ave`: The arithmetic mean.
        *   `adev`: The average absolute deviation from the mean.
        *   `var0`: The variance (calculated with `n-1` in the denominator for sample variance).
        *   `sdev`: The standard deviation (square root of variance).
        *   `skew`: The skewness (normalized third central moment, indicating asymmetry).
        *   `curt`: The kurtosis (normalized fourth central moment, typically reported as excess kurtosis by subtracting 3.0, indicating peakedness/tailedness relative to a normal distribution).
*   **Usage**: The program prompts the user for the `Number of data (ndata)` (must be at least 2), then successively for `ndata` values of `Y[i]`.
*   **Sample Run**:
    ```
     Number of data: 5
       1: 12
       2: 9
       3: 7
       4: 15
       5: 6

     Average ..........:   9.80000000000000E+0000
     Average  Deviation:   2.96000000000000E+0000
     Standard Deviation:   3.70135110466435E+0000
     Variance .........:   13.7000000000000E+0000
     Skewness .........:   0.29154869588874E+0000
     Kurtosis .........:  -1.90779903031595E-0000
    ```
*   **Real-world Applications**:
    *   **Data Exploration and Profiling**: Providing a quick and comprehensive summary of a dataset's characteristics, helping analysts understand the underlying distribution's shape and properties.
    *   **Risk Management**: In finance, skewness indicates the asymmetry of returns (positive for upside potential, negative for downside risk), while kurtosis quantifies the likelihood of extreme events (fat tails).
    *   **Quality Control and Process Improvement**: Monitoring process stability and identifying deviations from expected distributions. For example, a high kurtosis in manufacturing defect rates might indicate intermittent issues or abnormal operating conditions.
    *   **Statistical Modeling**: Informing the choice of appropriate statistical models by revealing properties like normality, skewness, or heavy tails in the data, which can affect model assumptions.
*   **Reference**: "Numerical Recipes, by W.H. Press, B.P. Flannery, S.A. Teukolsky and T. Vetterling, Cambridge University Press, 1986".

## 5. Median Calculation

This section details programs for calculating the median of a dataset, offering different algorithmic approaches.

### `tmdian.pas` - Median Value with Heapsort Method

*   **Program Name**: `TMDIAN`
*   **Description**: This program demonstrates how to calculate the median of a dataset by first sorting the data using the Heapsort algorithm. Heapsort is an efficient, in-place sorting algorithm suitable for large arrays (N log N complexity), making this method robust for finding the median even in substantial datasets where a sorted order is also desired.
*   **Core Functionality**:
    *   `PROCEDURE HPSORT(N:INTEGER; RA:pVEC)`: Implements the Heapsort algorithm to sort an array `RA` of length `N` in ascending order.
    *   `Procedure MDIAN(X:pVEC; N:Integer; Var XMED: REAL)`: This procedure first calls `HPSORT` to sort the input array `X`. Once sorted, it calculates the median: for an even number of elements, it's the average of the two middle elements; for an odd number, it's the single middle element. The input array `X` is modified (sorted).
*   **Usage**: The program generates a random table of numbers (e.g., `N=80` values between 0 and 1000) for demonstration purposes. No direct user input is required during execution. It displays the original table, the calculated median, and the sorted table.
*   **Sample Run**:
    ```
     Given table:

        1.0   63.8  278.9  406.2  546.8  657.7  638.4  324.6  745.5  852.3
      165.0  950.6  142.1  319.3  120.4  587.6  166.4  736.8  451.7  656.9
      605.7  312.7  565.0  614.3  326.3  660.0  933.0  494.3  349.6  559.1
      964.5  299.4  252.3  575.6  455.5   48.1  986.1  225.1  346.4   41.6
      283.1  288.0  999.4   44.4  815.1   20.3  452.0  699.7  460.0  584.8
      886.0  413.1  638.8  815.3   90.1  712.9    4.3  488.6  648.7  592.4
      172.8  456.4  991.6  237.6  981.1  854.4   94.6  631.3  678.4   16.2
      590.3   23.1  452.9  501.9  396.4  288.3  173.8  432.6  822.4  271.6

     Median value: 455.9250

     Sorted table (Heapsort method):

        1.0    4.3   16.2   20.3   23.1   41.6   44.4   48.1   63.8   90.1
       94.6  120.4  142.1  165.0  166.4  172.8  173.8  225.1  237.6  252.3
      271.6  278.9  283.1  288.0  288.3  299.4  312.7  319.3  324.6  326.3
      346.4  349.6  396.4  406.2  413.1  432.6  451.7  452.0  452.9  455.5
      456.4  460.0  488.6  494.3  501.9  546.8  559.1  565.0  575.6  584.8
      587.6  590.3  592.4  605.7  614.3  631.3  638.4  638.8  648.7  656.9
      657.7  660.0  678.4  699.7  712.9  736.8  745.5  815.1  815.3  822.4
      852.3  854.4  886.0  933.0  950.6  964.5  981.1  986.1  991.6  999.4
    ```
*   **Real-world Applications**:
    *   **Robust Central Tendency**: The median is a measure of central tendency that is less affected by outliers than the mean. This is particularly useful in fields like economics (e.g., median income, median housing prices) or environmental science where data distributions can be highly skewed.
    *   **Data Preparation**: Often, the first step in data analysis involves understanding the distribution and identifying the median as a representative value.
    *   **Efficient Sorting for Other Tasks**: While primarily for median calculation, the included `HPSORT` procedure provides a robust sorting utility that can be reused for other applications requiring efficient in-place sorting (e.g., finding percentiles, rank statistics).
*   **Reference**: "NUMERICAL RECIPES By W.H. Press, B.P. Flannery, S.A. Teukolsky and W.T. Vetterling, Cambridge University Press, 1986" [BIBLI 08].

### `tmdian1.pas` - Median Value Without Sorting

*   **Program Name**: `TMDIAN1`
*   **Description**: This program calculates the median value of an array of numbers using an iterative method that does *not* require sorting the entire array. This approach is beneficial for very large datasets where a full sort would be computationally expensive or memory-intensive, especially when only the median (or specific quantiles) is needed.
*   **Core Functionality**:
    *   `Function MAX(a,b:Double):Double`, `Function MIN(a,b:Double):Double`: Auxiliary functions for finding maximum and minimum.
    *   `Procedure MDIAN1(X:pVEC; N:Integer; Var XMED:Double)`: The primary procedure that computes the median. It uses an iterative approach based on successive guesses and adjustments (similar to a root-finding algorithm). This method maintains the original order of the input array.
*   **Usage**: The program reads `N` data points from an external file named `tmdian1.dat`. The file format is expected to first contain the number of data points ($N$), followed by $N$ floating-point numbers. It then displays the input table, calculates the median, and prints the computed median value.
*   **Sample Run (assuming `tmdian1.dat` contains the depicted data)**:
    ```
     N=80
     Given table:

     407.8 192.8 851.1 604.4 932.3 799.4 914.5 965.8 453.7 295.1
     154.5 977.4 410.2 916.2 934.7 504.8 823.2 225.2 456.6  49.0
     933.5 663.0 335.3 346.6 568.7 956.1 654.7 300.7 379.6 591.9
     992.9 689.6 644.7 305.4 148.2 257.2 664.6 612.1 713.0  99.7
      46.5 167.6 984.6 847.2  55.4  82.7 999.0  10.7 877.7 929.4
     398.1 972.8 874.1 755.1 472.1 122.8 671.4  35.5 128.8  76.8
     454.2 959.2 510.1 791.3 122.8 176.6 237.9 995.8 548.3 309.8
     162.6 996.5 750.0 250.6 577.7 761.1 101.9 797.1 539.0 723.5

     Median value: 558.5000
    ```
*   **Real-world Applications**:
    *   **Big Data Analytics**: For very large datasets that are too large to fit entirely into memory or where a full sort is impractical, this method provides an efficient way to find the median. This is common in fields like genomics, sensor networks, or large-scale simulations.
    *   **Stream Processing**: For real-time applications where data arrives continuously, this algorithm can be adapted to estimate the median without storing the entire history of observations, critical for performance and resource management.
    *   **Performance-Critical Applications**: In scenarios where computational time is critical and only the median (or specific quantiles) is needed, algorithms that avoid full sorting offer significant speed advantages.
*   **Reference**: "NUMERICAL RECIPES by W.H. Press, B.P. Flannery, S.A. Teukolsky and W.T. Vetterling, Cambridge University Press, 1986".

## 6. Regression Analysis

This section details a program for performing linear regression with a focus on robustness.

### `simlp.pas` - Simple Linear Regression Minimizing Sum of Absolute Deviation (SAD)

*   **Program Name**: `TSIMLP`
*   **Description**: This program performs a simple linear regression to fit a line of the form `Y = Alpha + Beta * X` to a set of data points. Unlike the common Ordinary Least Squares (OLS) method, this implementation minimizes the Sum of Absolute Deviations (SAD), also known as $L_1$ regression. This makes it a robust regression technique, less sensitive to outliers in the data.
*   **Core Functionality**:
    *   `Function ISign(a,b : Integer) : Integer`, `Function Sign(a,b : Double) : Double`: Auxiliary functions for determining the sign of a number.
    *   `Procedure SIMLP(N:Integer; X, Y:pVEC; Var SAD, ALPHA, BETA:Double; D: pVEC; Var ITER:Integer; INEXT:pIVEC; Var IFAULT:Integer)`: Implements Algorithm AS 132 from Applied Statistics. This iterative algorithm finds the regression coefficients (`ALPHA` and `BETA`) that minimize the sum of absolute differences between observed and predicted values. It employs a simplex-like approach.
*   **Usage**: The program uses hardcoded sample data (10 pairs of X and Y values) to demonstrate its functionality. It then calls the `SIMLP` procedure, which performs the iterative regression.
*   **Outputs**: The estimated intercept (`Alpha`), the estimated slope (`Beta`), the minimized Sum of Absolute Deviations (`SAD`), the number of iterations (`ITER`) the algorithm took to converge, and an `Error code` (0 for no error).
*   **Sample Run**:
    ```
     (Find Line Alpha + Beta X, minimizing SAD for next set
      of data):

     N = 10
     X(1)=1  Y(1)=1.15
     X(2)=2  Y(2)=2.22
     X(3)=3  Y(3)=2.94
     X(4)=4  Y(4)=3.85
     X(5)=5  Y(5)=5.10
     X(6)=6  Y(6)=5.99
     X(7)=7  Y(7)=7.25
     X(8)=8  Y(8)=8.04
     X(9)=9  Y(9)=9.003
     X(10)=10 Y(10)=9.999

       Alpha =  0.1667778
       Beta  =  0.9832222

       SAD =  0.8270000
       ITER = 3
       Error code: 0
    ```
*   **Real-world Applications**:
    *   **Robust Data Modeling**: Ideal for datasets that may contain outliers or measurement errors, where a standard OLS regression might be unduly influenced by extreme values. For instance, in financial modeling where market shocks or unusual events can occur, or in sensor data with occasional erroneous readings.
    *   **Econometrics**: In situations where the assumption of normally distributed errors (required for OLS) is violated, minimizing SAD can provide more reliable parameter estimates.
    *   **Signal Processing**: Filtering and fitting models to noisy signals where impulse noise or sparse errors are prevalent.
    *   **Medical and Biological Data**: Analyzing data where a few extreme values might represent unusual patient conditions or experimental anomalies, for which robust estimation is preferred.
*   **Reference**: "Journal of Applied Statistics (1978) vol.27 no.3".

## 7. Probability Expression Evaluation

This section documents a program designed to evaluate specific types of probability expressions.

### `prob.pas` - Evaluate Probabilities by Parsing String Expression

*   **Program Name**: `Prob`
*   **Description**: This program offers a unique functionality: it evaluates a mathematical expression provided as a string, specifically designed for combining binomial probability terms. It allows users to construct complex probability calculations from a simple input string.
*   **Core Functionality**:
    *   `Function Fact(x:real):Real`: Calculates the factorial of `x`.
    *   `Function C(n,k:integer): Integer`: Calculates combinations (`n` choose `k`), $C(n,k) = n! / (k! \cdot (n-k)!)$.
    *   `Function Power(x:real; n:integer): Real`: Calculates $x$ raised to the power of $n$.
    *   The main program loop parses the input string, identifying terms of the form `c(n,k,prob_a)`. This syntax is interpreted as $C(n,k) \cdot prob\_a^k \cdot (1.0 - prob\_a)^{(n-k)}$, which is the probability of exactly `k` successes in `n` trials with success probability `prob_a`. Terms can be combined using `+` or `-` operators.
*   **Supported Syntax**: Expressions like `'1-c(10,0,0.1)-c(10,1,0.1)'` are supported. (Note: The program indicates "no syntax checking" in its prompt, implying user responsibility for valid input format).
*   **Usage**: The program prompts the user to `Type expression to be evaluated`. An example is provided in the prompt: `1-c(10,0,0.2)-c(10,1,0.2)`.
*   **Outputs**: Prints the intermediate combination values $C(n,k)$ encountered during parsing and the final calculated `probability`.
*   **Sample Run**:
    ```
     Type expression to be evaluated (no syntax checking).
     Example: 1-c(10,0,0.2)-c(10,1,0.2)

     'c(10,0,.1)+c(10,1,0.1)+c(10,2,0.1)'

     C(10,0) = 1
     C(10,1) = 10
     C(10,2) = 45
     probability =    0.9298
    ```
*   **Real-world Applications**:
    *   **Custom Probabilistic Modeling**: Useful for scenarios where probabilities of multiple, independent binomial events need to be combined (e.g., calculating the probability of a system success given probabilities of individual component failures).
    *   **Risk Assessment**: Aggregating probabilities of different failure modes in engineering or financial systems, especially in discrete scenarios.
    *   **Educational Tool**: Demonstrating the application of binomial probability formula and combinations in a dynamic, interactive manner, allowing students to test various scenarios.
    *   **Quality Control**: Directly applicable to scenarios like calculating the probability of a certain range of defective items in a sample from a production line.
*   **Reference**: "Problem Solving with Fortran 90 By David R.Brooks, Springer-Verlag New York, 1997".

## 8. Statistical Unit (`Stats.pas`)

This section describes a reusable Free Pascal unit (`Stats.pas`) containing common statistical functions, along with its test driver program.

### `stats.pas` - Unit for Basic Statistical Functions

*   **Description**: `Stats.pas` is a Free Pascal unit (a module similar to a library) that provides a collection of reusable procedures for fundamental statistical operations. These include calculating descriptive statistics for normal distributions, performing linear regression, and generating arrays of normally distributed random numbers. This unit is designed to be integrated into larger Free Pascal projects.
*   **Dependencies**: The unit uses `WinCrt1` (likely a compatibility version of `WinCrt` for Free Pascal environments) for console I/O, though its primary function is to provide callable statistical routines.
*   **Core Functionality (Exported Procedures)**:
    *   **Types**: Defines `pVec` (pointer to a vector) and `VEC` (array of `Real`) for flexible data handling.
    *   `Function Dot_Product(n:Integer; x,y:pVec): Real`: Calculates the dot product of two real vectors. (Internal function).
    *   `Function Sum(n:Integer; x:pVec): Real`: Calculates the sum of elements in a real vector. (Internal function).
    *   `Procedure NormalStats(a:pVec; n:Integer; flag:Char; Var avg, std_dev: Real)`: Computes the `avg` (mean) and `std_dev` (standard deviation) for an array `a` of `n` real numbers. The `flag` parameter determines whether to calculate `'p'` (population) or `'s'` (sample) statistics. It handles cases of negative variance by setting `std_dev` to -1.0.
    *   `Procedure NormalArray(a:pVec; n:Integer)`: Fills an array `a` of `n` real numbers with pseudo-random numbers that follow a standard normal distribution (mean 0, standard deviation 1). It uses the Box-Muller transform, converting pairs of uniform random numbers into normal random numbers.
    *   `Procedure LinearReg(x,y:pVec; n:Integer; flag:Char; Var a,b,s_yx,r: Real)`: Performs a simple linear regression to fit a line $y = ax + b$ to `n` pairs of `x` and `y` data points. It calculates the intercept (`a`), slope (`b`), standard error of the estimate of $y$ on $x$ (`s_yx`), and the sample correlation coefficient (`r`). If `a` is initialized to `0.0` on input, the regression is forced to pass through the origin (0,0). Handles error conditions (e.g., zero standard deviation) by setting `r` to 0.0.
*   **Usage**: As a unit, it is included in other Pascal programs via the `Uses Stats;` directive. It exposes its procedures to be called directly from client code.
*   **Real-world Applications (for developers using this unit)**:
    *   **Custom Statistical Applications**: Developers can use these functions as building blocks for creating tailored statistical software or tools for specific domains (e.g., scientific simulations, data visualization tools, internal data analysis scripts).
    *   **Simulation and Modeling**: `NormalArray` is invaluable for generating synthetic data that follows a normal distribution, essential for Monte Carlo simulations, testing algorithms, or creating realistic input for models (e.g., financial market simulations, engineering system stress tests).
    *   **Data Preprocessing**: `NormalStats` can be used to quickly summarize and validate input data before more complex analyses or for standardizing data (e.g., z-score normalization for machine learning inputs).
    *   **Rapid Prototyping**: The unit provides pre-built, tested functions, accelerating the development of applications requiring basic statistical analyses.
*   **Reference**: "Problem Solving with Fortran 90 By David R.Brooks, Springer-Verlag New York, 1997".

### `tstats.pas` - Driver Program to Test Unit STATS

*   **Program Name**: `Tstats`
*   **Description**: This program acts as a comprehensive test driver for the `Stats.pas` unit. It demonstrates the proper usage and verifies the outputs of the `NormalStats`, `NormalArray`, and `LinearReg` procedures defined within the `Stats` unit. This program is useful for developers to understand how to integrate and use the `Stats` unit in their own applications.
*   **Dependencies**: This program explicitly `Uses WinCrt1, Stats;`, demonstrating its reliance on the `Stats` unit and `WinCrt1` for console output.
*   **Core Functionality**:
    *   Initializes the random number generator.
    *   **Tests Basic Statistics**: Generates `n` (e.g., 100) normally distributed random numbers using `NormalArray`. Calculates and prints both population and sample means and standard deviations using `NormalStats`.
    *   **Tests Linear Regression**: Generates new data (`x`, `y`) with a linear relationship plus added normal noise.
        *   Performs a full linear regression (`y=ax+b`) by initializing `a:=1.0` and calling `LinearReg`, then prints the regression coefficients, standard error, and correlation coefficient.
        *   Performs a linear regression forced through the origin (`y=bx`) by initializing `a:=0.0` and calling `LinearReg`, then prints the results.
*   **Usage**: The program runs automatically upon execution. No user input is required, as data is generated internally for testing purposes.
*   **Sample Run**:
    ```
     Population mean and standard deviation:   0.075133   0.879506

     Sample mean and standard deviation:   0.075133   0.883904

     For full regression:
                regression coefficients:   1.006200   1.986324
     standard error of estimate of y on x:   9.918043
                correlation coefficient:   0.985044

     For regression forced through (0,0):
                regression coefficients:   0.000000   2.006249
     standard error of estimate of y on x:   9.935045
                correlation coefficient:   0.984993
    ```
*   **Real-world Applications**:
    *   **Component Testing and Validation**: Provides a template for unit testing and verifying the correctness of statistical algorithms, which is crucial for ensuring the reliability of any software involving numerical computations.
    *   **Learning Resource**: Serves as an excellent practical example for developers to understand how to call and interpret the results from the `Stats` unit's procedures within a larger program.
    *   **Debugging and Benchmarking**: Can be adapted to debug specific issues or benchmark the performance and numerical stability of the statistical functions under different conditions.

## 9. References

This section consolidates the references cited across the various Pascal source files within the `stat` folder. These references provide the theoretical and algorithmic foundations for the implementations in this library.

*   **[BIBLI 08]**: Press, W.H., Flannery, B.P., Teukolsky, S.A., & Vetterling, W.T. (1986). *Numerical Recipes*. Cambridge University Press. (Cited in `fdistri.pas`, `gamma.pas`, `ibeta.pas`, `tmdian.pas`, `tmoment.pas`).
*   **[BIBLI 13]**: Haut, H. (1981). *Mathematiques et statistiques*. PSI Editions, France. (Cited in `distri.pas`, `fstat.pas`, `moment.pas`).
*   Adams, A.G. (1969). Areas under the normal curve, Algorithm 39. *Computer J.*, *12*, 197-198. (Cited in `tnormal.pas`).
*   Brooks, D.R. (1997). *Problem Solving with Fortran 90*. Springer-Verlag New York. (Cited in `prob.pas`, `tstats.pas`, `stats.pas`).
*   Hart, J.F. et al. (1968). *Computer Approximations*. Wiley. (Cited in `tnormal.pas`).
*   Hill, I.D., & Pike, M.C. (1967). Algorithm 299. *Collected Algorithms for the CACM*, p. 243. (Updated for rounding errors based on remark in ACM TOMS June 1985, page 185). (Cited in `chi2.pas`).
*   Ibbetson, D. (1963). Algorithm 209. *Collected Algorithms of the CACM*, p. 616. (Cited in `chi2.pas`).
*   Journal of Applied Statistics (1968) vol.17, p.189. (Cited in `student.pas`).
*   Journal of Applied Statistics (1972) vol.21, p.226. (Cited in `momnts.pas`).
*   Journal of Applied Statistics (1973) vol.22 no.3. (Cited in `tnormal.pas`).
*   Journal of Applied Statistics (1978) vol.27 no.3. (Cited in `simlp.pas`).
*   Walker, John (adapted from Gary Perlman). CHI2 and Inverse CHI2 Functions. Online resource: `http://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html`. (Cited in `chi2.pas`).
