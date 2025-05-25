
# jpmMath Library - arith - Basic Arithmetic and Number Theory

## Core Principles

The `arith` folder focuses on providing fundamental tools for tasks ranging from base conversions and high-precision arithmetic to equation solving and expression parsing, emphasizing direct mathematical implementation within the Pascal environment.

## Programs Overview

### 1. Numerical Base Conversion (`base.pas`)

*   **Description**: This program facilitates the conversion of numerical values from one base to another. It supports bases ranging from 2 (binary) to 36 (using digits 0-9 and letters A-Z for values 10-35). The conversion process involves an intermediate step where the number is first converted to base 10 (decimal) and then to the target base.

*   **Key Concepts**: Positional numeral systems, base 10 conversion, iterative division/multiplication for base conversion. The `Decodebase` function converts from an arbitrary base to base 10, while `Codebase` converts from base 10 to an arbitrary base.

*   **How to Use**:
    1.  Execute the `base.pas` program.
    2.  The program will prompt for the "Start base" (e.g., 5 for base 5, 10 for decimal, 16 for hexadecimal).
    3.  Next, it will ask for the "Arrival base".
    4.  You can then enter numbers in the "Start base" repeatedly. Enter a null string (just press Enter) to exit.

*   **Example Usage (from source code):**
    ```
    BASE CONVERSION

    Start  base (2 to 36): 5
    Arival base (2 to 36): 3

    Enter number in start base: 3421
    In base  3: 200000
    Enter number in start base: 3420
    In base  3: 122222


    BASE CONVERSION

    Start  base (2 to 36): 10
    Arival base (2 to 36): 16

    Enter number in start base: 65535
    In base 16: FFFF
    Enter number in start base: 100
    In base 16: 64


    BASE CONVERSION

    Start  base (2 to 36): 16
    Arival base (2 to 36): 2

    Enter number in start base: FF
    In base  2: 11111111
    Enter number in start base: 3A
    In base  2: 111010
    Enter number in start base: E2
    In base  2: 11100010
    ```

*   **Limitations**: The program uses `REAL` (floating-point) types for intermediate base 10 conversion (`MAXREAL = 1E40`), which limits the precision and magnitude of numbers that can be accurately converted, especially for very large numbers or conversions between bases that significantly expand the number of digits. For calculations requiring higher precision or larger numbers, `pcalc.pas` should be considered.

### 2. Arithmetic Operations in Any Base (`basisop.pas`)

*   **Description**: This program enables performing basic arithmetic operations (`+`, `-`, `*`, `/`, `^` for power) on numbers expressed in any base between 2 and 36. Similar to `base.pas`, it converts operands to base 10, performs the operation, and then converts the result back to the specified base.

*   **Key Concepts**: Extends base conversion to arithmetic operations. It re-implements `Decodebase` and `Codebase` functions locally, identical to those in `base.pas`.

*   **How to Use**:
    1.  Execute the `basisop.pas` program.
    2.  First, you will be prompted to "Enter base (2 to 36)".
    3.  Then, "Enter 1st number" in that base.
    4.  Subsequently, you can "Enter next number" and an "Operator (+-*/^)". The program will apply the operation to the current result and the new number, displaying the updated result. Enter a null string for the next number to exit.

*   **Example Usage (from source code):**
    ```
    ARITHMETIC OPERATIONS IN ANY BASE BETWEEN 2 AND 36

    Enter base (2 to 36): 25
    Enter 1st number : AM3G

    Enter next number: 28
    Operator (+-*/^) : +
    Result           : AM50
    Enter next number: J
    Operator (+-*/^) : *
    Result           : 86MD6
    Enter next number: 2B
    Operator (+-*/^) : /
    Result           : 39JM


    ARITHMETIC OPERATIONS IN ANY BASE BETWEEN 2 AND 36

    Enter base (2 to 36): 2
    Enter 1st number : 11001010

    Enter next number: 1100
    Operator (+-*/^) : -
    Result           : 10111110
    Enter next number: 1110
    Operator (+-*/^) : -
    Result           : 10110000
    Enter next number: 10101
    Operator (+-*/^) : +
    Result           : 11000101
    Enter next number: 10
    Operator (+-*/^) : /
    Result           : 1100010
    Enter next number: 11
    Operator (+-*/^) : *
    Result           : 100100110


    ARITHMETIC OPERATIONS IN ANY BASE BETWEEN 2 AND 36

    Enter base (2 to 36): 16
    Enter 1st number : FEAA

    Enter next number: FF
    Operator (+-*/^) : +
    Result           : FFA9
    Enter next number: 1799
    Operator (+-*/^) : -
    Result           : E810
    Enter next number: 2
    Operator (+-*/^) : /
    Result           : 7408
    ```

*   **Limitations**: Shares the same `REAL` precision limitations as `base.pas` (`MAXREAL = 1E40`), which means operations on very large numbers or those producing very long results may suffer from accuracy issues. As noted in the source, the number of useful digits (`udg`) in base `b` is `INT(40*log(2)/log(b))`, resulting in `udg = 12` for base 10, `10` for base 16, and `7` for base 36.

### 3. Combinatory Analysis (`combi.pas`)

*   **Description**: This program provides functions for common combinatory analysis operations: factorial (N!), combinations (C(n,p)), and permutations (A(n,p)). It handles large factorial results by presenting them in scientific notation (mantissa and exponent) for N >= 25, leveraging properties of logarithms.

*   **Key Functions**:
    *   `Factorial(n: INTEGER; VAR mantissa, exponent: REAL): REAL`: Calculates `n!`. For `n < 25`, it returns the direct factorial. For larger `n`, it uses Stirling's approximation (`(Ln(2*PI*n)/2 + n*Ln(n)-n)/Ln(10)`) to return the mantissa and exponent (base 10) of the result.
    *   `Cnp(n, p: INTEGER): REAL`: Calculates the number of combinations (C(n,p)). Returns 0 if `p > n`.
    *   `Anp(n, p: INTEGER): REAL`: Calculates the number of permutations (A(n,p)). Returns 0 if `p > n`.

*   **Constants**:
    *   `MAXREAL = 1E4000`: A very large real number used for overflow checking in combination and permutation calculations.

*   **How to Use**:
    1.  Execute the `combi.pas` program.
    2.  A menu will be displayed with options:
        *   1: Factorial n!
        *   2: Combination C n,p
        *   3: Permutation A n,p
        *   0: Quit
    3.  Enter your choice. For factorial, enter N. For combinations and permutations, enter N and P.

*   **Example Usage (from source code):**
    ```
    COMBINATORY ANALYSIS

      1:  Factorial n!
      2:  Combination C n,p
      3:  Permutation A n,p
      0:  Quit

      Your choice (0 to 3): 1

      N = 100
      N! = 9.3248476268  10^ 157

      Your choice (0 to 3): 2

      N = 7
      P = 3
      Cnp =           35

      Your choice (0 to 3): 3

      N = 10
      P = 6
      Anp =       151200
    ```

*   **Limitations**: While `Factorial` supports large numbers using mantissa/exponent representation, `Cnp` and `Anp` still rely on `REAL` calculations with `MAXREAL = 1E4000`. Very large combinations/permutations might still exceed this limit or lose precision.

### 4. Coincidence of Letters (`complet.pas`)

*   **Description**: This program identifies and lists common letters shared between two input words.

*   **Key Concepts**: String iteration, character comparison, basic boolean logic. It modifies the second word internally by replacing found common letters with a space to prevent re-reporting the same coincidence if a letter appears multiple times in the first word.

*   **How to Use**:
    1.  Execute the `complet.pas` program.
    2.  You will be prompted to "Input first word" and "Input second word".
    3.  The program will then print each letter found to be common to both words.

*   **Example Usage (from source code):**
    ```
    Input first word : Bonjour
    Input second word: Monsieur

    The letter o is common to the two words.
    The letter n is common to the two words.
    The letter u is common to the two words.
    The letter r is common to the two words.
    ```

*   **Limitations**: Words are limited to a maximum length of 26 characters (`Lmax = 26`). The internal modification of the second word means it's not suitable for scenarios where the original content of both strings must be preserved after the check.

### 5. Diophantine Equation Solver (`diophan.pas`)

*   **Description**: This program solves linear Diophantine equations of the form `Ax + By = C`, where A, B, C, x, and y are integer numbers. It finds general integer solutions (x0 + k\*q, y0 - k\*p) if they exist.

*   **Key Functions**:
    *   `GCD(a, b: REAL): REAL`: (Identical to `GCD` in `gcd.pas`) Computes the greatest common divisor of two integer numbers using the Euclidean algorithm.
    *   `Diophantian(a, b, c: REAL; VAR x0, y0, p, q: REAL): BOOLEAN`: Solves the equation. Returns `TRUE` if solutions exist (i.e., GCD(a,b) divides C), and sets `x0, y0, p, q`. Returns `FALSE` otherwise.

*   **How to Use**:
    1.  Execute the `diophan.pas` program.
    2.  Input the integer coefficients A, B, and C when prompted.
    3.  The program will output the general solutions for X and Y in terms of an arbitrary integer K, or indicate if no integer solutions exist.

*   **Example Usage (from source code):**
    ```
    SOLVING IN Z EQUATION AX + BY = C

       A = 3
       B = -2
       C = 7

    Solutions are:

      X=     1 +     2 *K
      Y=    -2 +     3 *K
    ```

*   **Limitations**: The calculations within `GCD` and `Diophantian` functions use `REAL` types, which might introduce floating-point inaccuracies when dealing with extremely large integer coefficients, although `abs(r)<1E-10` is used for comparison robustness.

### 6. Integer Factorization (`factors.pas`)

*   **Description**: This program factorizes a given integer number into its prime factors. It implements a trial division algorithm optimized by checking divisibility by 2, 3, and then numbers of the form 6i-1 and 6i+1.

*   **Key Concepts**: Prime factorization, trial division, basic optimization for prime number checks. The algorithm iteratively divides the number by its smallest prime factors until it becomes 1.

*   **How to Use**:
    1.  Execute the `factors.pas` program.
    2.  Enter the integer number you wish to factorize when prompted with `?`.
    3.  The program will print the prime factors separated by spaces.

*   **Example Usage (from source code):**
    ```
    ? 394616

    2 2 2 107 461
    ```

*   **Limitations**: The input number `n` is a `REAL` type, which severely limits the maximum integer value that can be accurately factorized due to floating-point precision constraints. The `eps` value (`1e-6`) used for equality checks with remainders also implies this limitation. It is not suitable for cryptographic-grade factorization of large numbers.

### 7. Operations on Fractions (`fraction.pas`)

*   **Description**: This program handles elementary arithmetic operations on fractions (addition, subtraction, multiplication, and division). It prompts for numerators and denominators and displays simplified results.

*   **Procedures**:
    *   `P300`: Performs subtraction of fractions.
    *   `P310`: Performs addition of fractions.
    *   `P330`: Performs division of fractions (inverts the divisor and multiplies).
    *   `P360`: Performs multiplication of fractions.
    *   `P390`: Simplifies the resulting fraction by dividing the numerator and denominator by their greatest common divisor.

*   **How to Use**:
    1.  Execute the `fraction.pas` program.
    2.  Enter the numerator and denominator of the first fraction (e.g., `3 4` for 3/4). If you enter `0` for numerator or denominator, the program exits.
    3.  Enter an operator:
        *   `+`, `-`, `*`, `/`: Perform the respective operation with the next fraction.
        *   `=`: Display the current result.
        *   `s`: Save the current result.
    4.  After an operator, enter the next fraction (numerator and denominator). If you enter `0` for the denominator, it will use the previously `s`aved result.
    5.  The program will continuously prompt for operations and new fractions until `0` is entered.

*   **Example Usage (from source code):**
    ```
    3 4
    *
    8 9
    +
    -4 5
    =
    Result is -2/15
    ```

*   **Limitations**: The numerators and denominators are stored as `INTEGER` types, which restricts the magnitude of the fractions that can be processed. Very large numerators or denominators will cause overflow.

### 8. GCD and SCM of Integer Numbers (`gcd.pas`)

*   **Description**: This program calculates the Greatest Common Divisor (GCD) and the Smallest Common Multiple (SCM) of multiple integer numbers provided by the user. The calculation is iterative, taking one number at a time and updating the GCD and SCM of the set.

*   **Key Function**:
    *   `GCD(a, b: REAL): REAL`: Computes the greatest common divisor of two integer numbers using the Euclidean algorithm. It handles absolute values and checks for potential overflow for very large inputs.

*   **How to Use**:
    1.  Execute the `gcd.pas` program.
    2.  Enter the "First number".
    3.  Continue entering "Next number" until you enter `0` to signal the end of input.
    4.  The program will then display the calculated GCD and SCM for all entered numbers.

*   **Example Usage (from source code):**
    ```
    GCD and SCM of integer numbers:

    First number: 9936
    Next number : 414
    Next number : 3174
    Next number : 0

    GCD =      138
    SCM =   228528
    ```

*   **Limitations**: The numbers `x`, `pg`, `pp` are `REAL` types, limiting the magnitude of integers for which GCD and SCM can be accurately calculated. Precision issues with `REAL` can affect the accuracy of the result for very large numbers.

### 9. Number of Days Between Two Dates (`ndays.pas`)

*   **Description**: This program computes the number of days between two given dates. It uses a simplified formula and a pre-calculated array for months.

*   **Procedure**:
    *   `Calcul`: Reads month, day, and year, then calculates `N`, the number of days from a fictitious start date (01/00/1901), based on a formula using an array `A` for month offsets and `365.25 * (Y - 1901)`.

*   **How to Use**:
    1.  Execute the `ndays.pas` program.
    2.  You will be prompted to "Enter first date (M D Y)" (Month Day Year, separated by spaces).
    3.  Then, "Enter second date (M D Y)".
    4.  The program will output the "Number of days" between these two dates.

*   **Example Usage (Conceptual, based on source description):**
    ```
    Compute number of days between two dates
        (from 01/01/1901 to 12/31/2099)

    Enter first date (M D Y): 1 1 2000
    Enter second date (M D Y): 1 1 2001

    Number of days:    366
    ```

*   **Limitations**:
    *   **Date Range**: Restricted to dates between January 1, 1901, and December 31, 2099. This avoids complexities related to "secular" leap years (e.g., 1900 not being a leap year, 2000 being one).
    *   **Precision**: Uses `REAL` for day calculations, which might introduce minor floating-point inaccuracies, though `Round(N)` is used to obtain integer days.

### 10. Recursive Descent Parser (`parser.pas`)

*   **Description**: This is a powerful recursive descent parser capable of evaluating mathematical expressions, handling variables with multiple characters, and supporting a set of built-in mathematical functions. It allows for variable assignment, complex expressions, and session management (saving/loading defined variables and formulas).

*   **Key Concepts**: Recursive descent parsing, lexical analysis (tokenizing), operator precedence, variable assignment, function calls, session management. Uses `WinCrt`, `WinTypes`, `WinProcs`, `WObjects`, `Strings` units for UI and string manipulation.

*   **Core Features**:
    *   **Variables:** Supports variables with multiple characters (e.g., `Masse`, `Area`). Variables can be assigned numerical values or formulas.
    *   **Arithmetic Operations:** Supports `+`, `-`, `*`, `/`, `%` (modulo), and `^` (power).
    *   **Built-in Functions:** Includes `arctan`, `cos`, `exp`, `log` (natural logarithm), `root` (square root), `sin`, `tan`.
    *   **Constants:** `pi` is a recognized constant.
    *   **Session Management:** Ability to save and load defined variables and their associated formulas to/from disk files.
    *   **Error Handling:** Provides basic error messages for syntax issues, division by zero, or invalid function arguments.

*   **Internal Working (High-Level):**
    The parser uses a recursive descent approach, breaking down expressions into smaller, self-similar components (terms, factors, primitives) and uses a set of mutually recursive procedures to evaluate them according to operator precedence.

*   **Tokenization (`get_token`):** This procedure scans the input string (`prog`) and extracts the next "token," which can be a delimiter (operator, parenthesis), a variable name, or a number. It determines the `tok_TYPE` (DELIMITER, VARIABLE, NUMBER).
*   **Operator Precedence (Level Functions):**
    *   `level2`: Handles addition and subtraction (`+`, `-`).
    *   `level3`: Handles multiplication, division, and modulo (`*`, `/`, `%`).
    *   `level4`: Handles exponentiation (`^`).
    *   `level5`: Handles unary minus (`-`) and built-in functions (e.g., `sin`, `cos`, `exp`).
    *   `level6`: Processes expressions enclosed in parentheses `()`.
    *   `primitive`: Extracts the actual numerical value of a number or a variable.
*   **Variable and Formula Management:** The parser maintains arrays (`vrbl`, `formule`, `flg`) to store variable names, their associated formulas (expressions), and flags indicating if a variable's value has been accessed. This allows for defining and recalling variables, enabling complex calculations by building upon previously defined expressions.
*   **Error Handling:** Includes basic error reporting for syntax errors, mismatched parentheses, and mathematical errors (e.g., division by zero, invalid logarithm/square root arguments).

*   **How to Use**:
    1.  Execute the `parser.pas` program.
    2.  A main menu will be displayed:
        *   `l`: List current variables and their formulas/values.
        *   `c`: Enter calculation mode. In this mode, you can type expressions like `A=10/4`, `F=A*X^2+B*X+C`, `sin(pi/2)`. An empty line exits calculation mode.
        *   `r`: Read a session (variables and formulas) from a disk file.
        *   `s`: Save the current session to a disk file.
        *   `x`: Suppress (delete) a defined variable.
        *   `q`: Quit the parser (prompts to save changes).

*   **Example Usage (from source code):**
    ```
        ********************************************
        *              PARSER                      *
        *                                          *
        *  List variables and formulas       :  l  *
        *  Make calculations                 :  c  *
        *  Read a session from disk          :  r  *
        *  Save a session to disk            :  s  *
        *  Suppress a variable               :  x  *
        *  --------------------------------------  *
        *  Quit parser                       :  q  *
        *                                          *
        ********************************************
            Your choice ( l, c, r, s, x or q ) : c

        *** CALCULATE ***

      A=10/4
      B=25.478
      C=-1
      X=0.72
      F=A*X^2+B*X+C
      F
        formula: A*X^2+B*X+C
        A=  2.50000000
        X=  0.72000000
        B= 25.47800000
        C= -1.00000000
        17.92476800
      exp(1)
        2.71828183
    ```

*   **Limitations**:
    *   Maximum of `NB_VAR = 40` variables.
    *   Mathematical operations are performed using `double` precision floating-point numbers, which may have precision limitations for certain types of computations.
    *   Error handling for invalid mathematical operations (e.g., log of non-positive number, sqrt of negative number) prints an error message and sets an internal `ok` flag to `false`.

### 11. Pascal's Triangle (`pastri.pas`)

*   **Description**: This program generates and displays Pascal's Triangle, a triangular array of the binomial coefficients. It can display up to a maximum of 15 lines.

*   **Key Concepts**: Binomial coefficients, combinatorial mathematics, iterative generation of triangle rows (each number is the sum of the two numbers directly above it).

*   **How to Use**:
    1.  Execute the `pastri.pas` program.
    2.  Enter the "How many lines" you wish to display. If the input exceeds 15, it will default to 15.
    3.  The program will then print the Pascal's Triangle.

*   **Example Usage (from source code):**
    ```
    How many lines: 10


    p      0    1    2    3    4    5    6    7    8    9
    n
    -----------------------------------------------------
    0 --   1
    1 --   1    1
    2 --   1    2    1
    3 --   1    3    3    1
    4 --   1    4    6    4    1
    5 --   1    5   10   10    5    1
    6 --   1    6   15   20   15    6    1
    7 --   1    7   21   35   35   21    7    1
    8 --   1    8   28   56   70   56   28    8    1
    9 --   1    9   36   84  126  126   84   36    9    1
    --------------------------------------------------------
    ```

*   **Limitations**: The maximum number of lines is hardcoded to 15 (`Nmax = 15`). The internal array `t` uses `Integer` type, which means numbers in the triangle that exceed the maximum value of an `Integer` will cause overflow, limiting the practical number of rows even further if `Nmax` were increased.

### 12. Prime Number Test (`prime.pas`)

*   **Description**: This program determines whether a given integer number is a prime number or not. It uses an optimized trial division method, checking for divisibility by 2, 3, and then numbers of the form 6k-1 and 6k+1 up to the square root of the number.

*   **Key Concepts**: Primality testing, trial division, optimized prime checking (skipping multiples of 2 and 3).

*   **How to Use**:
    1.  Execute the `prime.pas` program.
    2.  Enter the integer number you want to test when prompted with `?`.
    3.  The program will output either "Prime number." or "Not a prime number.".

*   **Example Usage (Conceptual):**
    ```
    ? 17

    Prime number.

    ? 18

    Not a prime number.
    ```

*   **Limitations**: The input number `n` is a `REAL` type. This significantly restricts the size of integers that can be reliably tested for primality due to floating-point precision. For cryptographically relevant prime numbers, much more advanced algorithms (e.g., Miller-Rabin) are required, which are beyond the scope of this implementation.

### 13. Table of Prime Numbers (`primes.pas`)

*   **Description**: This program generates and displays a table of prime numbers up to a predefined limit (N = 2000 in the current implementation). It implements a variation of the Sieve of Eratosthenes algorithm to efficiently find all primes within the specified range.

*   **Key Concepts**: Sieve of Eratosthenes, prime number generation, array manipulation for marking composite numbers.

*   **Algorithm (`Sieve of Eratosthenes` variant):**
    1.  Initialize an array `Raye` (meaning "crossed out") where `Raye[i]` is 0 if `i` is potentially prime, and 1 if `i` is composite.
    2.  Mark 1 as composite.
    3.  Iterate through numbers starting from `prem = 2`.
    4.  If `prem` is not marked as composite (`Raye[prem] = 0`), it's a prime. Mark all its multiples (`2*prem, 3*prem, ...`) as composite in the `Raye` array.
    5.  Repeat until `prem * prem` exceeds `N`.
    6.  Finally, collect all numbers `i` for which `Raye[i]` is 0 into the `Premier` array and display them.

*   **How to Use**:
    1.  Execute the `primes.pas` program.
    2.  The program will automatically generate primes up to its internal limit (2000) and display them in a formatted table.

*   **Example Usage (from source code - partial output):**
    ```
    Between 1 and 2000, there are 303 prime numbers:

        2    3    5    7   11   13   17   19   23   29   31   37   41   43   47
       53   59   61   67   71   73   79   83   89   97  101  103  107  109  113
       ... (continues up to 2000)
    ```

*   **Limitations**: The upper limit for prime number generation (`N = 2000`) and the sizes of the `Premier` and `Raye` arrays are hardcoded. To generate primes beyond 2000, these constants and array sizes would need to be modified and the program recompiled.

## Technical Details

This section elaborates on the core mechanisms and significant functions within the `arith` folder, providing insight into their implementation.

### PCalcFun Unit (`pcalc_fun.txt`)

The `PCalcFun` unit is crucial for handling "huge integers" in `pcalc.pas`. Instead of using standard integer types, it manipulates numbers as `PChar` (pointer to a character array, effectively a null-terminated string), allowing for arbitrary precision limited only by available memory.

Key functions in `PCalcFun`:

*   `SubChar(C1, C2: char; Var borrow: boolean): Char;` and `AddChar(C1, C2: char; Var carry: boolean): Char;`: These are assembler-level functions for single-digit arithmetic, handling borrows and carries efficiently.
*   `TrimLead0(P: PChar; Len: word);`: Removes leading zeros from a `PChar` number string.
*   `add(A, B: PChar; MaxLen: Word; Var Status: Word): PChar;`: Performs addition of two large numbers. It handles string manipulation to align digits and propagates carries. `addWTrail` is an internal helper for adding numbers with trailing zeros (used in multiplication).
*   `sub(A, B: PChar; MaxLen: Word; Var Status: Word): PChar;`: Performs subtraction. It first compares the numbers to ensure subtraction of the smaller from the larger, then adjusts the sign if necessary.
*   `prod(A, B: PChar; MaxLen: Word; Var Status: Word): PChar;`: Implements multiplication by repeatedly adding the multiplicand. It iterates through each digit of the multiplier, adding scaled copies of the multiplicand.
*   `divide(A, B, remainder: PChar; MaxLen: Word; Var Status: Word): PChar;`: Performs division using repeated subtraction. It iteratively subtracts scaled versions of the divisor from the remainder to build the quotient.
*   `fact(A: PChar; MaxLen: Word; Var Status: Word): PChar;`: Calculates the factorial of a large number by iteratively multiplying the number by decreasing integers until 1.
*   `Power(B, E: PChar; MaxLen: Word; Var Status: Word): PChar;`: Computes base `B` raised to exponent `E` for large numbers. It uses a "halving and squaring" algorithm for efficiency.
*   `AddComma(WW: PChar; MaxLen: Word): PChar;` and `StripComma(WW: PChar): PChar;`: Utility functions to add or remove comma separators for readability in large numbers.
*   `AllNums(P: PChar): Boolean;`: Checks if a `PChar` string contains only numeric digits.
*   `Compare(X, Y: PChar): ShortInt;`: Compares two large integer strings.

*   **Usage**: This unit is not a standalone program. It is designed to be used by other Free Pascal programs that require arbitrary-precision arithmetic capabilities. The `pcalc.pas` program is an example of how this unit is utilized. **Important Note on Compilation and Accessibility:** For the `PCalcFun` unit to be successfully compiled and used by programs like `pcalc.pas`, the `pcalc_fun.txt` file *must* first be renamed to `pcalc_fun.pas`. After renaming, this `.pas` file needs to be made accessible to the Free Pascal compiler. Common methods for ensuring accessibility include:
    *   **Placing it in the same directory** as the main program (`pcalc.pas`) that `Uses` it.
    *   **Placing it in one of Free Pascal's standard unit search paths** (e.g., in a `units` subdirectory within your Free Pascal installation, or a directory specified by the `fpc.cfg` configuration file).
    *   **Specifying the unit's directory during compilation** using the `-Fu` compiler option (e.g., `fpc -Fumy_units_path/ pcalc.pas`).
This ensures the compiler can locate, compile, and link the unit's interface and implementation when building your application.

*   **Limitations**: `MaxLen` parameter for functions indicates the maximum allowed length of the result string, typically up to 5000 characters as used in `pcalc.pas`. Operations can return status codes indicating overflow, low memory, or division by zero.

### PChar-Based Calculator Program (`pcalc.pas`)

*   **Description**: This program is a command-line calculator that leverages the `PCalcFun` unit to perform arithmetic operations on "huge" integers (up to 5000 digits). It supports addition, subtraction, multiplication, division, factorial, and exponentiation.

*   **Dependencies**:
    *   `PCalcFun` unit (`pcalc_fun.pas`)

*   **Key Concepts**: Demonstrates the practical application of the `PCalcFun` unit for large integer arithmetic. Command-line parsing for input expressions.

*   **How to Use**:
    1.  Execute the `pcalc.pas` program.
    2.  You will be prompted with `>`. Enter expressions in the format `n1 op n2` (e.g., `12345 + 67890`, `2 ^ 100`) or `n1 !` for factorial (e.g., `200 !`).
    3.  Enter a null string (just press Enter) or a two-character input to exit.

*   **Example Usage (from source code):**
    ```
    >2 ^ 100

    POWER: 1,267,650,600,228,229,401,496,703,205,376

    >200 !

    FACTORIAL: 30,414,093,201,713,378,043,612,608,166,064,768,844,377,
               641,568,960,512,000,000,000,000
    ```

*   **Limitations**: The maximum length of the integer result is constrained by the `Max` constant, typically 5000 characters. Input parsing expects specific formats (`n1 op n2` or `n1 !`).

### Base Conversion Logic (`base.pas` and `basisop.pas`)

Both `base.pas` and `basisop.pas` utilize the `Decodebase` and `Codebase` functions for converting numbers between an arbitrary base (2 to 36) and base 10 (real number representation).

*   `Decodebase(x: STRING; b: INTEGER; VAR y: REAL): BOOLEAN;`: Converts a number from `base b` (string `x`) to `base 10` (real `y`). It iterates through the input string from right to left, multiplying each digit by increasing powers of the base. It handles digits '0'-'9' and 'A'-'Z' (for values 10-35).
*   `Codebase(x: REAL; b: INTEGER; VAR y: STRING): BOOLEAN;`: Converts a number from `base 10` (real `x`) to `base b` (string `y`). It repeatedly takes the remainder of division by `b` to get the digits in the new base, then prepend them to the result string.

A limitation noted in `basisop.pas` is the precision of `REAL` numbers in Pascal, which has a mantissa of 40 digits in base 2. This translates to fewer useful digits in higher bases (e.g., 12 for base 10, 10 for base 16, 7 for base 36), limiting the magnitude of numbers that can be accurately processed without exceeding `MAXREAL`. For calculations requiring higher precision, the `pcalc.pas` program with `PCalcFun` should be used.

### Parser Design (`parser.pas`)

The `parser.pas` program implements a recursive descent parser. This parsing technique breaks down the input expression into smaller, self-similar components (terms, factors, primitives) and uses a set of mutually recursive procedures to evaluate them according to operator precedence.

*   **Tokenization (`get_token`):** This procedure scans the input string (`prog`) and extracts the next "token," which can be a delimiter (operator, parenthesis), a variable name, or a number. It determines the `tok_TYPE` (DELIMITER, VARIABLE, NUMBER).
*   **Operator Precedence (Level Functions):**
    *   `level2`: Handles addition and subtraction (`+`, `-`).
    *   `level3`: Handles multiplication, division, and modulo (`*`, `/`, `%`).
    *   `level4`: Handles exponentiation (`^`).
    *   `level5`: Handles unary minus (`-`) and built-in functions (e.g., `sin`, `cos`, `exp`).
    *   `level6`: Processes expressions enclosed in parentheses `()`.
    *   `primitive`: Extracts the actual numerical value of a number or a variable.
*   **Variable and Formula Management:** The parser maintains arrays (`vrbl`, `formule`, `flg`) to store variable names, their associated formulas (expressions), and flags indicating if a variable's value has been accessed. This allows for defining and recalling variables, enabling complex calculations by building upon previously defined expressions.
*   **Error Handling:** Includes basic error reporting for syntax errors, mismatched parentheses, and mathematical errors (e.g., division by zero, invalid logarithm/square root arguments).

## Real-World Problem Solving Ideas

The `jpmMath` library, particularly the `arith` folder, offers foundational mathematical tools that can be applied to a variety of real-world scenarios:

*   **Financial Modeling and Analytics:**
    *   **`pcalc.pas` (Huge Integer Calculator) & `PCalcFun`:** Essential for high-precision financial calculations where standard floating-point arithmetic might lead to rounding errors, especially for very large sums or complex interest calculations over extended periods. Examples include calculating compound interest on national debts or simulating long-term investment growth.
    *   **`fraction.pas` (Fraction Operations):** Useful for scenarios requiring exact fractional results, such as prorating expenses, calculating stock splits, or precise ratio analysis where decimal approximations are undesirable.
    *   **`ndays.pas` (Number of Days Between Dates):** For calculating interest accrual periods, bond yields, or duration of financial contracts.
*   **Cryptography and Security:**
    *   **`factors.pas` (Integer Factorization), `prime.pas` (Prime Number Check), `primes.pas` (Table of Prime Numbers):** Fundamental to understanding and implementing public-key cryptographic systems like RSA, which rely on the difficulty of factoring large numbers into their prime components. These tools can be used for educational purposes to demonstrate the underlying principles or for generating smaller primes for testing cryptographic algorithms.
    *   **`base.pas` (Numerical Base Conversion) & `basisop.pas` (Operations in Any Base):** While not directly a cryptographic algorithm, base conversion can be used in data encoding and obfuscation techniques. Operations in different bases might be relevant for specific, non-standard cryptographic transformations or coding challenges.
*   **Scientific and Engineering Applications:**
    *   **`combi.pas` (Combinatory Analysis):** For calculating probabilities in statistical analysis, determining the number of possible configurations in system design, or analyzing permutations in scheduling algorithms.
    *   **`parser.pas` (Recursive Descent Parser):** Can serve as the backend for a scientific calculator application, a simple scripting language for data analysis, or a formula evaluator in simulation software where users need to define custom equations.
    *   **`diophan.pas` (Diophantine Equation Solver):** Useful in optimization problems where resources must be allocated in integer quantities, such as distributing goods, scheduling tasks, or designing systems with discrete components. For example, finding combinations of quantities that satisfy specific cost or resource constraints.
    *   **`gcd.pas` (GCD and SCM Calculation):** Applicable in areas like gear design (finding optimal gear ratios), timing synchronization in electrical engineering, or even in musical theory for understanding harmonic intervals.
*   **Educational Tools and Simulations:**
    *   **`base.pas` and `basisop.pas`:** Excellent for teaching number systems and arithmetic across different bases, aiding students in understanding concepts like binary, hexadecimal, and other numerical representations.
    *   **`pastri.pas` (Pascal's Triangle):** A visual aid for teaching combinatorics, probability, and binomial expansion.
    *   **`complet.pas` (Letter Coincidences):** A simple tool for demonstrating string manipulation and set theory concepts in an approachable manner.
*   **Logistics and Resource Management:**
    *   **`gcd.pas` (GCD and SCM Calculation):** Can be used to optimize packing or cutting problems (e.g., finding the largest common size to cut materials without waste) or to determine efficient scheduling cycles for recurring events.
    *   **`diophan.pas` (Diophantine Equation Solver):** Useful for resource allocation problems where quantities must be integers, ensuring exact distributions without remainders.

## File Manifest

This section lists all the source code files provided in the `arith` folder, along with a brief description of their content and purpose.

*   `_info.txt`: Provides a high-level description of the programs available in the `arith` folder.
*   `base.pas`: Pascal program for converting numbers between different bases (2-36).
*   `basisop.pas`: Pascal program for performing arithmetic operations in arbitrary bases (2-36).
*   `combi.pas`: Pascal program for combinatory analysis (Factorial, Combination, Permutation).
*   `complet.pas`: Pascal program to find common letters between two words.
*   `diophan.pas`: Pascal program to solve Diophantine equations (Ax + By = C) and includes a GCD function.
*   `factors.pas`: Pascal program for factoring an integer number into its prime components.
*   `fraction.pas`: Pascal program for elementary operations on fractions.
*   `gcd.pas`: Pascal program to compute the Greatest Common Divisor (GCD) and Smallest Common Multiple (SCM) of several integers.
*   `ndays.pas`: Pascal program to compute the number of days between two dates.
*   `parser.pas`: Pascal program implementing a recursive descent parser for mathematical expressions with variables and functions.
*   `pastri.pas`: Pascal program to display Pascal's Triangle.
*   `pcalc.pas`: Pascal program for a PChar-based calculator handling huge integers.
*   `pcalc_fun.txt`: Pascal Unit (`PCalcFun`) containing functions for arithmetical operations on large integers represented as PChar strings. **This file should be named `pcalc_fun.pas` for compilation.**
*   `prime.pas`: Pascal program to test if an integer number is a prime number.
*   `primes.pas`: Pascal program to write a table of prime numbers from 1 to N.
