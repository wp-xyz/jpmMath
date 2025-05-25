
# jpmMath Library - matrices - Numerical Linear Algebra

The `matrices` sub-library offers a versatile collection of routines for various numerical linear algebra operations, including:
*   Solving linear systems of equations (dense, sparse, symmetric, banded, and special structures).
*   Calculating matrix properties (determinants, characteristic polynomials).
*   Performing matrix inversions.
*   Determining eigenvalues and eigenvectors using a variety of algorithms.

The implementations often prioritize numerical stability through techniques like pivoting, and some are designed for memory efficiency or optimized performance for specific matrix structures.

## Table of Contents

1.  [General Usage Guidelines](#1-general-usage-guidelines)
    *   [Matrix and Vector Storage Conventions](#11-matrix-and-vector-storage-conventions)
    *   [Memory Management](#12-memory-management)
    *   [Error Handling](#13-error-handling)
    *   [Input/Output Operations](#14-inputoutput-operations)
2.  [Core Concepts and Utility Modules](#2-core-concepts-and-utility-modules)
    *   [`basis.pas`](#21-basis.pas)
    *   [`ucomplex1.pas`](#22-ucomplex1.pas)
    *   [`Type_def.pas`](#23-type_def.pas)
3.  [Linear System Solvers](#3-linear-system-solvers)
    *   [General Dense Linear Systems](#31-general-dense-linear-systems)
        *   Gauss-Jordan Method (`sysmat.pas`)
        *   LU Decomposition (`lu.pas`)
        *   Direct Factorization (`tdlittl.pas`)
        *   Triangularization (Gauss Elimination) (`tlinear.pas`)
        *   Solving with Reduced Storage (`dple.pas`)
        *   Singular Value Decomposition (SVD) (`tsvbksb.pas`)
    *   [Iterative Linear System Solvers](#32-iterative-linear-system-solvers)
        *   Iterative Gauss-Seidel Method (`fseidel.pas`)
        *   Conjugate Gradient Method (`cgtst1.pas`)
        *   Conjugate Gradient for Sparse Systems (`tsparse.pas`)
    *   [Symmetric Linear Systems](#33-symmetric-linear-systems)
        *   Gauss Method (`syslin.pas`)
        *   Cholesky Decomposition (`fcholy.pas`)
        *   SYMSOL Method (`tsymsol.pas`)
    *   [Specialized Linear Systems](#34-specialized-linear-systems)
        *   Banded Matrices (with Pivots) (`fband.pas`)
        *   Banded Matrices (without Pivots) (`fbando.pas`)
        *   Non-Symmetric Banded Matrices (`nsbslv.pas`)
        *   Tridiagonal Systems (`tridiag.pas`)
        *   Vandermonde Systems (`tvander.pas`)
        *   Toeplitz Systems (`toeplitz.pas`)
4.  [Matrix Operations](#4-matrix-operations)
    *   [Matrix Inversion](#41-matrix-inversion)
        *   Inversion by LU Decomposition (`inv_lu.pas`)
        *   Inversion by Householder's Method (`househol.pas`)
        *   Inversion of a Symmetric Positive Definite Matrix by Cholesky Method (`choles.pas`)
    *   [Determinant Calculation](#42-determinant-calculation)
        *   Gauss Method (`deter.pas`)
        *   LU Decomposition Method (`deter1.pas`)
        *   Recursive Kramer's Rule (`deter2.pas`)
        *   `FindDet` Function (`tfinddet.pas`)
5.  [Eigenvalue Problems](#5-eigenvalue-problems)
    *   [Characteristic Polynomial Calculation](#51-characteristic-polynomial-calculation)
        *   Real Tridiagonal Matrix (`carpol.pas`)
        *   Complex Square Matrix (`carpol1.pas`)
        *   Real Square Matrix (`carpol2.pas`)
        *   Real Symmetric Square Matrix (Lanczos's Method) (`carpol3.pas`)
    *   [Power Methods for Eigenvalues](#52-power-methods-for-eigenvalues)
        *   Greatest Eigenvalue (`tpwm.pas`)
        *   Smallest Eigenvalue (`tpwimgt.pas`)
    *   [Jacobi's Method](#53-jacobis-method)
        *   Real Symmetric Matrices (`ujacobi.pas`)
        *   Hermitian Matrices (`tephj.pas`)
    *   [QL Method for Symmetric Matrices](#54-ql-method-for-symmetric-matrices)
        *   Real Symmetric Tridiagonal Matrices (`elprotd.pas`, `ttql2.pas`)
        *   Real Symmetric Matrices (Householder Reduction + QL) (`ttred2.pas`)
    *   [Rutishauser's Method and Inverse Iteration](#55-rutishauser-s-method-and-inverse-iteration) (`elpro.pas`)
    *   [QR Algorithm for Non-Symmetric Matrices](#56-qr-algorithm-for-non-symmetric-matrices) (`feigen0.pas`, `linpack.pas`, `test_hqr.pas`, `thqr.pas`)

---

## 1. General Usage Guidelines

This section provides essential guidance for effectively utilizing the `matrices` sub-library, focusing on common conventions and best practices for memory management, error handling, and data representation.

### 1.1. Matrix and Vector Storage Conventions

The library employs various conventions for storing matrices and vectors, primarily to optimize memory usage and access patterns in Pascal. Users must be aware of these conventions as they are specific to each unit and procedure.

*   **Static 2D Arrays (`Array[1..NMAX, 1..NMAX]` or `Array[0..SIZE, 0..SIZE]`):**
    *   Common in many older or simpler programs (e.g., `syslin.pas`, `deter.pas`, `tseidel.pas`).
    *   Matrices are typically indexed `A[row, col]`. The base index (0 or 1) and maximum size (`NMAX` or `SIZE`) are defined by constants within the respective unit.
    *   Memory is allocated on the stack or in the global data segment, which can be restrictive for very large matrices.

*   **Pointers to 1D Arrays for 2D Matrix Storage (`pVECT`):**
    *   A prevalent dynamic allocation strategy for larger matrices or where memory efficiency is critical (e.g., `lu.pas`, `basis.pas`, `fband.pas`).
    *   A conceptual matrix `A(N, N)` is linearized into a 1D array pointed to by `pVECT`.
    *   **Common Convention (1-indexed):** As seen in `lu.pas` (`NOTA2`), the element `A[i, j]` is accessed as `A^[i * (N + 1) + j]`. Here, `N+1` acts as the stride for moving to the next row, enabling 1-based indexing for both rows and columns.
    *   **Alternative Convention (0-indexed):** Some `basis.pas` procedures (e.g., `ReadMat`, `WriteMat`, `SetMat`) assume 0-indexed linear storage where `A[i,j]` is accessed as `A^[i * N + j]`.
    *   **Crucial Note:** Always verify the expected indexing convention (0- or 1-based) and the linear storage formula for the specific procedure being used. The maximum capacity for these `pVECT` types is often defined by `SIZE` or `NMAX` constants (e.g., `Array[0..SIZE] of REAL`), meaning your matrix dimensions should not exceed this predefined limit to avoid access violations.

*   **Pointers to Records/Arrays (`pMAT`, `pVEC`, `pMC`, `pIVECT`):**
    *   Many units define their own `TYPE` for `MATRIX` (e.g., `MATR` in `fseidel.pas`, `MATRIX` in `fcholy.pas`, `MATC` for complex in `tephj.pas`) or `VECTOR` (`VEC`, `VECTOR`), and then declare pointers to these types (`pMAT`, `pVEC`, `pMC`, `pIVECT`).
    *   These types often encapsulate standard 2D arrays (either 0-indexed or 1-indexed, `[0..SIZE, 0..SIZE]` or `[1..NMAX, 1..NMAX]`).
    *   **Complex Matrices:** Units dealing with complex numbers (e.g., `tephj.pas`) define `Complex` records and `pCMat` for matrices of `Complex` elements.
    *   **Always Consult Source:** For precise indexing and internal structure of these types, refer to the specific unit's source code declarations.

### 1.2. Memory Management

Most procedures in this library utilize dynamic memory allocation for matrices and vectors via pointers. Proper memory management is paramount to prevent memory leaks, ensure program stability, and optimize resource usage.

*   **Allocation (`New()`):** Before using any pointer-type variable (e.g., `pMAT`, `pVEC`), you **must** allocate memory for it using the `New()` procedure. This reserves the necessary space on the heap.
    ```pascal
    VAR
      myMatrixPtr: pMAT;
      myVectorPtr: pVEC;
    BEGIN
      New(myMatrixPtr); // Allocates memory for the matrix structure
      New(myVectorPtr); // Allocates memory for the vector structure
      // ... Proceed to fill data into myMatrixPtr^ and myVectorPtr^
    END;
    ```
    *   **Important:** The `New()` call allocates memory for the *structure* the pointer points to (e.g., `Array[0..SIZE,0..SIZE]`). Ensure your actual matrix/vector dimensions do not exceed the `SIZE` or `NMAX` constants defined in the respective unit, as this can lead to runtime errors (array out of bounds) or memory corruption.

*   **Deallocation (`Dispose()`):** After you have finished using a dynamically allocated variable, you **must** release its memory back to the system using the `Dispose()` procedure. Failing to do so causes memory leaks, which can degrade system performance over time, especially in long-running applications.
    ```pascal
    BEGIN
      // ... Your code using myMatrixPtr and myVectorPtr
      Dispose(myMatrixPtr);
      Dispose(myVectorPtr);
      // Pointers should not be used after Dispose() without re-allocating
    END;
    ```

*   **In-Place Operations:** Many procedures in this library are designed to operate "in-place," meaning they modify their input matrices and vectors directly. This is done to save memory and improve performance by avoiding redundant data copying.
    *   If you need to preserve the original data, make a copy of the matrix/vector before passing it to an in-place procedure. For example, `copy_vector` from `basis.pas` can be used for this purpose.
    *   Consult the documentation for each procedure to understand if it modifies its inputs.

### 1.3. Error Handling

Most numerical procedures are designed to return an integer `code` or `rc` (return code) parameter to indicate the success or failure of the operation. Robust error handling is crucial for reliable numerical software.

*   **Checking Return Codes:** A return code of `0` typically indicates successful completion. Non-zero values usually signify an error condition (e.g., singular matrix, non-convergence, invalid input parameters, division by zero).
    *   Always check the `code`/`rc` parameter immediately after a procedure call.
    *   Refer to the specific procedure's documentation (or source code comments) for the detailed meaning of each non-zero error code.
    ```pascal
    VAR
      resultCode: INTEGER;
      // ...
    BEGIN
      // ... Call a procedure that returns an error code
      LUDCMP(A, N, INDX, D, resultCode);
      IF resultCode <> 0 THEN
      BEGIN
        // Handle the error gracefully
        LogError('LU Decomposition failed with code: ' + IntToStr(resultCode));
        // Optionally, exit the program or attempt recovery
        Halt(1);
      END;
      // ... Continue with program logic only if successful
    END;
    ```

*   **Error Logging:** The `basis.pas` unit provides utility procedures for reporting errors:
    *   `LogError(text: STRING)`: Prints an error message to standard output and halts program execution.
    *   `LogError1(VAR fp: TEXT; text: STRING)`: Prints an error message to a specified text file (without halting by default, allowing for custom error recovery).

### 1.4. Input/Output Operations

Many of the example programs rely on reading input data from text files (often with `.dat` extension) and writing computational results to output text files (typically with `.lst` extension). The `basis.pas` unit provides convenient and formatted routines for these I/O operations.

*   **File Handling:**
    *   `Assign(fileVar, filename)`: Associates a file variable with a physical file name.
    *   `Reset(fileVar)`: Opens an existing file for reading.
    *   `Rewrite(fileVar)`: Creates a new file (or overwrites an existing one) for writing.
    *   `Close(fileVar)`: Closes an opened file, ensuring data is written and resources are released.
*   **Reading/Writing Data:**
    *   `ReadVec`, `WriteVec`: For single vectors.
    *   `ReadMat`, `WriteMat`: For matrices (which are internally stored as 1D vectors in many cases).
    *   Variations exist (e.g., `ReadVec1`, `WriteMat1`) to support both 0-based and 1-based indexing, and to specify the number of items per line (`ic`) for formatting.
*   **Formatted Numerical Output:**
    *   `f_aff_reel(VAR fp: TEXT; l2: REAL)`: This crucial procedure, typically from `basis.pas`, ensures that real numbers are displayed in a consistent 10-character format. It automatically switches to scientific notation (e.g., `1.23E+003`) for very large or very small numbers, preventing display issues and ensuring readability in numerical reports.
*   **Headers and Footers:**
    *   `WriteHead(VAR fp: TEXT; nom: string)`: Writes a standardized header to an output file.
    *   `WriteEnd(VAR fp: TEXT)`: Writes a separator line to an output file, indicating section breaks or end of results.

By adhering to these general guidelines, developers can effectively integrate and utilize the `jpmMath` matrices library components within their Free Pascal projects.

---

## 2. Core Concepts and Utility Modules

This section describes the foundational units that define common data structures and provide essential utility routines, forming the backbone for numerical operations across the `matrices` sub-library.

### 2.1. `basis.pas`

This unit (sometimes referred to as `basis_r.pas` in include statements) provides a collection of fundamental routines for vector and matrix manipulation. It is a highly utilized unit, serving as a common dependency for many other modules in the library.

**Key Procedures/Functions:**

*   `norm_max(VAR vektor: pVECT; n: INTEGER): REAL;`
    *   **Purpose**: Computes the maximum (infinity) norm of a real vector, defined as the maximum absolute value of its elements.
*   `copy_vector(VAR ziel: pVECT; VAR quelle: pVECT; n: INTEGER);`
    *   **Purpose**: Copies `n` elements from a source vector `quelle` to a destination vector `ziel`. This is essential for preserving original data before in-place operations.
*   `LogError(text:STRING);`
    *   **Purpose**: Prints a critical error message to standard output and halts program execution. Use for unrecoverable errors.
*   `LogError1(VAR fp:TEXT; text:STRING);`
    *   **Purpose**: Prints an error message to a specified text file. Useful for logging non-fatal errors or detailed error reports.
*   `ReadVec(VAR fp:TEXT; n, ic:INTEGER; VAR x:pVECT);`
    *   **Purpose**: Reads `n` real numbers into a vector `x` from a text file `fp`. Assumes 0-based indexing for `x`. `ic` specifies how many items are expected per line in the file.
*   `ReadVec1(VAR fp:TEXT; n, ic:INTEGER; VAR x:pVECT);`
    *   **Purpose**: Reads `n` real numbers into a vector `x` from a text file `fp`. Assumes 1-based indexing for `x`. `ic` specifies items per line.
*   `SetVec(n:INTEGER; VAR x:pVECT; val:REAL);`
    *   **Purpose**: Initializes all `n` elements of a vector `x` with a constant `val`.
*   `WriteVec(VAR fp:TEXT; n, ic:INTEGER; VAR x:pVECT);`
    *   **Purpose**: Writes `n` elements of vector `x` to a text file `fp`, using 0-based indexing. `ic` specifies how many items to write per line for formatting.
*   `WriteVec1(VAR fp:TEXT; n, ic:INTEGER; VAR x:pVECT);`
    *   **Purpose**: Writes `n` elements of vector `x` to a text file `fp`, using 1-based indexing. `ic` specifies items per line.
*   `ReadMat(VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);`
    *   **Purpose**: Reads an `m x n` matrix from a text file `fp`. The matrix is expected to be stored linearly (row-major) in the `pVECT` `a`, with 0-based indexing. `ic` specifies items per line.
*   `ReadMat1(VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);`
    *   **Purpose**: Reads an `m x n` matrix from a text file `fp`. The matrix is expected to be stored linearly (row-major) in the `pVECT` `a`, with 1-based indexing (`a^[i*(n+1)+j]`). `ic` specifies items per line.
*   `WriteMat(VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);`
    *   **Purpose**: Writes an `m x n` matrix to a text file `fp`. The matrix is stored linearly (row-major) in `pVECT` `a`, with 0-based indexing. `ic` specifies items per line.
*   `WriteMat1(VAR fp:TEXT; m, n, ic:INTEGER; VAR a:pVECT);`
    *   **Purpose**: Writes an `m x n` matrix to a text file `fp`. The matrix is stored linearly (row-major) in `pVECT` `a`, with 1-based indexing. `ic` specifies items per line.
*   `SetMat(m, n:INTEGER; VAR a:pVECT; val:REAL);`
    *   **Purpose**: Initializes all elements of an `m x n` matrix `a` with a constant `val`.
*   `WriteHead(VAR fp:TEXT; nom : string);`
    *   **Purpose**: Writes a formatted header with a given string `nom` to a text file `fp`.
*   `WriteEnd(VAR fp:TEXT);`
    *   **Purpose**: Writes a standardized separator line to a text file `fp`, typically indicating the end of a report section.
*   `min(a,b:INTEGER):INTEGER;`
    *   **Purpose**: Returns the smaller of two integers.
*   `max(a,b:INTEGER):INTEGER;`
    *   **Purpose**: Returns the larger of two integers.
*   `Swap(VAR a:REAL; VAR b:REAL);`
    *   **Purpose**: Exchanges the values of two real numbers.
*   `f_aff_reel(VAR fp:TEXT; l2 : REAL);`
    *   **Purpose**: Displays a real number `l2` to a text file `fp` in a standardized 10-character wide format. It automatically uses scientific notation (`E+nnn`) for numbers outside a common range (e.g., `[0.001, 10000]`), ensuring readability for varied magnitudes.

**Real-world Application Ideas:**
`basis.pas` serves as the foundational utility unit for any numerical computation involving matrices and vectors. Its comprehensive I/O and basic manipulation routines are essential for:
*   **Data Acquisition and Preparation**: Loading experimental data from text files into suitable matrix/vector structures.
*   **Reporting and Visualization**: Generating formatted output files for detailed analysis or input for plotting software.
*   **Code Debugging and Verification**: Quick inspection of intermediate computation results through `WriteMat`/`WriteVec` and `f_aff_reel`.
*   **Initializing Problem States**: Setting up initial conditions for iterative solvers or physical simulations.

### 2.2. `ucomplex1.pas`

This unit (`UComplex1.pas`, often aliased as `UComplex` in `USES` clauses) defines a `Complex` record type and provides a set of fundamental arithmetic operations for complex numbers. It is crucial for modules that deal with complex-valued matrices, such as those found in quantum mechanics or electrical engineering.

**Key Types and Procedures:**

*   **`Complex = Record r, i: Double End;`**
    *   **Purpose**: Defines a record structure to represent a complex number. `r` stores the real part and `i` stores the imaginary part, both as `Double` precision floating-point numbers.
*   `CADD(c1, c2: Complex; Var c3: Complex);`
    *   **Purpose**: Performs complex addition: `c3 = c1 + c2`.
    *   **Formula**: `c3.r = c1.r + c2.r`, `c3.i = c1.i + c2.i`.
*   `CDIF(c1, c2: Complex; Var c3: Complex);`
    *   **Purpose**: Performs complex subtraction: `c3 = c1 - c2`.
    *   **Formula**: `c3.r = c1.r - c2.r`, `c3.i = c1.i - c2.i`.
*   `CMUL(c1, c2: Complex; Var c3: Complex);`
    *   **Purpose**: Performs complex multiplication: `c3 = c1 * c2`.
    *   **Formula**: `c3.r = c1.r * c2.r - c1.i * c2.i`, `c3.i = c1.r * c2.i + c1.i * c2.r`.
*   `CPRO(alpha: Double; C: Complex; VAR c1: Complex);`
    *   **Purpose**: Multiplies a complex number `C` by a real scalar `alpha`: `c1 = alpha * C`.
    *   **Formula**: `c1.r = alpha * C.r`, `c1.i = alpha * C.i`.

**Real-world Application Ideas:**
`ucomplex1.pas` enables numerical algorithms to handle complex-valued data, which is essential in domains such as:
*   **Electrical Engineering (AC Circuits):** Representing impedances, voltages, and currents as complex numbers for frequency domain analysis.
*   **Quantum Mechanics:** Simulating quantum systems where wave functions and operators are inherently complex.
*   **Signal Processing:** Performing Discrete Fourier Transforms (DFT) and analyzing complex-valued signals (e.g., audio, radio frequencies).
*   **Control Systems:** Analyzing frequency response and stability of systems using complex transfer functions.

### 2.3. `Type_def.pas`

While the source code for `Type_def.pas` is not explicitly provided, its usage in modules like `feigen0.pas`, `linpack.pas`, and `test_hqr.pas` indicates its role in defining global data types.

**Presumed Key Types:**
*   **`REAL_AR = Double;`**
    *   **Purpose**: Standardizes the floating-point precision used across a set of modules to `Double` (64-bit floating-point). This ensures consistency and high precision for numerical stability in sensitive calculations like eigenvalue problems.
*   `Square_Matrix`, `Rect_Matrix`, `Real_Vector`, `Integer_Vector`:
    *   **Purpose**: Likely defines various matrix and vector types based on `REAL_AR` (and `INTEGER` for integer vectors), typically as fixed-size arrays (e.g., `ARRAY[0..Maxc, 0..Maxc] OF REAL_AR`). `Maxc`, `Nsol`, `Ndiag` are likely constants defining maximum dimensions for these statically allocated types.

**Real-world Application Ideas:**
The use of `Type_def.pas` promotes:
*   **Code Maintainability**: Centralizing type definitions makes it easier to change the fundamental precision (e.g., from `Single` to `Double`) or maximum array sizes across multiple dependent units by modifying only one file.
*   **Numerical Accuracy**: Enforcing `Double` precision for `REAL_AR` is critical for algorithms where accumulated floating-point errors can significantly impact results, such as iterative eigenvalue solvers or those dealing with ill-conditioned matrices.
*   **Code Clarity**: Provides clear aliases for complex data structures, improving readability.

---

## 3. Linear System Solvers

This section provides a detailed overview of the diverse methods available for solving linear systems of equations, typically represented as $AX=B$, where $A$ is the coefficient matrix, $X$ is the vector of unknowns, and $B$ is the right-hand side vector. The library offers both direct and iterative approaches, as well as specialized solvers for matrices with particular structures.

### 3.1. General Dense Linear Systems

These methods are designed for general matrices where most elements are non-zero. They solve the system in a finite number of operations, typically through matrix factorization.

#### Gauss-Jordan Method (`sysmat.pas`)

This program implements the Gauss-Jordan elimination method with full pivoting to solve linear systems ($AX=B$), calculate the inverse of a square matrix ($A^{-1}$), or compute its determinant ($\det(A)$). The method is optimized for memory usage by performing operations directly on the input matrices (`AA` and `BB`), overwriting their original content.

*   **Numerical Method:** The core idea is to transform the input matrix `AA` into an identity matrix (`I`) through a sequence of elementary row and column operations. The same transformations are applied simultaneously to the right-hand side matrix `BB`. If `BB` initially holds $B$, it becomes the solution $X$. If `BB` initially holds $I$, it becomes $A^{-1}$. Full pivoting (selecting the element with the largest absolute value in the remaining submatrix as the pivot at each step) is employed to enhance numerical stability and minimize round-off errors.
*   **Key Procedure:**
    *   `MATINV(N, M: integer; VAR AA: MAT; VAR BB: MAT1; VAR DET: DOUBLE);`
        *   **Purpose**: Computes $A^{-1}$, $\det(A)$, and solves $AX=B$.
        *   **N**: Size of the square matrix `AA`.
        *   **M**: Number of columns in matrix `BB` (representing multiple right-hand sides). If `M=0`, the procedure only computes $A^{-1}$ and $\det(A)$, treating `BB` as an auxiliary workspace.
        *   **AA** (Input/Output): On input, the coefficient matrix $A$. On output, it is replaced by its inverse $A^{-1}$.
        *   **BB** (Input/Output): On input, the right-hand side matrix $B$. On output, it is replaced by the solution matrix $X = A^{-1}B$.
        *   **DET** (Output): The determinant of the original matrix $A$.
        *   **Note**: The input matrices `AA` and `BB` are *destroyed* during the process due to in-place computation. The procedure halts with an error if the determinant is practically zero (singular matrix) or if a pivot element is too small.
*   **Example Usage (from `sysmat.txt`):**
    To solve $AX=B$ where $A$ is 3x3 and $B$ is 3x4:
    ```
    Input A:
     2.0  -1.0   1.0
     1.0   5.0  -2.0
     3.0  -2.0   3.0

    Input B:
     5.0  -1.0   1.0  2.0
     1.0   2.0   3.0  4.0
     3.0  7.556 4.0  4.0

    Output (Inverse of A, Solution Matrix X, Determinant):
    INVERSE OF MATRIX A:
     0.785714   0.071429  -0.214286
    -0.642857   0.214286   0.357143
    -1.214286   0.071429   0.785714

    SOLUTION MATRIX X:
     3.357143  -2.262000   0.142857   1.000000
    -1.928571   3.770000   1.428571   1.000000
    -3.642857   7.294000   2.142857   1.000000

    DETERMINANT= 14.000000
    ```
*   **Associated Files:** `sysmat.pas` (main program and procedures), `sysmat.txt` (explanation file).
*   **Real-world Application Ideas:**
    *   **Electrical Network Analysis:** Solving systems derived from Kirchhoff's current and voltage laws in complex circuits.
    *   **Chemical Reaction Balancing:** Determining the stoichiometric coefficients of chemical reactions by setting up a linear system where atoms must be conserved.
    *   **Robot Kinematics:** Solving inverse kinematics problems for simple robotic manipulators to determine joint angles for desired end-effector positions.

#### LU Decomposition (`lu.pas`)

This unit provides robust routines for LU decomposition, a factorization of a matrix $A$ into a lower triangular matrix $L$ and an upper triangular matrix $U$ ($A=LU$). This method is particularly efficient for solving multiple linear systems ($AX=B$) that share the same coefficient matrix $A$ but have different right-hand side vectors $B$. It also yields information to calculate the determinant.

*   **Numerical Method:** The implementation uses **Crout's algorithm** with partial pivoting (row interchanges). Partial pivoting involves swapping rows to bring the element with the largest absolute value in the current column to the pivot position. This strategy significantly enhances numerical stability and minimizes error propagation. The decomposition is performed "in-place," meaning the factors $L$ and $U$ overwrite the original matrix `A` in memory. After decomposition ($PA=LU$, where $P$ is a permutation matrix), solving $AX=B$ is achieved in two straightforward steps:
    1.  **Forward Substitution:** Solve $Ly = PB$ for the intermediate vector $y$.
    2.  **Back Substitution:** Solve $Ux = y$ for the solution vector $x$.
*   **Key Procedures:**
    *   `LUDCMP(A: pVECT; N: INTEGER; VAR INDX: pIVECT; VAR D: INTEGER; VAR CODE: INTEGER);`
        *   **Purpose**: Performs the LU decomposition of matrix `A` with partial pivoting.
        *   **A** (Input/Output): The $N \times N$ matrix. On input, the original matrix. On output, it is replaced by its LU decomposition factors (elements of $L$ and $U$).
        *   **N**: Size of the square matrix.
        *   **INDX** (Output): An integer vector (1-indexed) that records the row permutations performed during pivoting. This information is critical for `LUBKSB` and determinant calculation.
        *   **D** (Output): An integer (1 or -1) representing the sign of the permutation. This is used when calculating the determinant, where $\det(A) = D \cdot \prod(\text{diagonal elements of } U)$.
        *   **CODE** (Output): Return code: `0` for successful decomposition, `1` if the matrix is detected as singular (a pivot element is too close to zero).
    *   `LUBKSB(A: pVECT; N: INTEGER; INDX: pIVECT; VAR B: pVECT);`
        *   **Purpose**: Solves the linear system $AX=B$ using the LU decomposition factors previously computed by `LUDCMP`.
        *   **A** (Input): The LU-decomposed matrix (output from `LUDCMP`). This parameter is *not modified* by `LUBKSB`.
        *   **N**: Size of the matrix.
        *   **INDX** (Input): The row permutation vector obtained from `LUDCMP`.
        *   **B** (Input/Output): On input, the right-hand side vector $B$. On output, it is replaced by the solution vector $X$.
        *   **Note**: `A` and `INDX` remain unchanged, allowing `LUBKSB` to be called multiple times with different right-hand sides $B$ for the same decomposed matrix $A$.
*   **Matrix Storage:** Matrices passed as `pVECT` are assumed to be stored in a 1-indexed linear fashion, where the element $A[i,j]$ is accessed as `A^[i*(N+1)+j]`.
*   **Associated Files:** `lu.pas` (main unit), `lu.txt` (explanation file), `test_lu.pas` (test program for solving $AX=B$), `inv_lu.pas` (program for matrix inversion using LU).
*   **Real-world Application Ideas:**
    *   **Iterative Solvers Preconditioning:** In large-scale numerical simulations, incomplete LU (ILU) factorization is a popular preconditioning technique to accelerate the convergence of iterative solvers for sparse matrices.
    *   **Sensitivity Analysis:** When analyzing how a system's output changes with respect to variations in input parameters, LU decomposition can efficiently solve for derivatives or influence coefficients.
    *   **Real-time Optimization:** In embedded systems or real-time control, if the system matrix remains constant, pre-computing the LU decomposition allows for very fast solutions when the right-hand side changes rapidly.

#### Direct Factorization (`tdlittl.pas`)

This program provides a direct factorization method for solving linear systems $AX=B$. This method explicitly computes the $L$ and $U$ factors of $A$ (where $L$ is lower triangular with unit diagonal elements and $U$ is upper triangular). It then solves two simpler triangular systems: $LZ=B$ and $UX=Z$. It incorporates partial pivoting to ensure numerical stability by selecting the largest element in each column of the transformation matrix.

*   **Numerical Method:** The `DLITTL` procedure performs an LU decomposition (likely a variant of Doolittle's or Crout's factorization). It iteratively builds the $L$ and $U$ factors directly. At each step, it seeks the maximum element in the current column below the diagonal to use as a pivot, performing row swaps as necessary. This process overwrites the input matrix `A` and stores the $L$ and $U$ factors in a working matrix `W`. The solution is then obtained by a forward substitution step ($LZ=B$) and a backward substitution step ($UX=Z$).
*   **Key Procedure:**
    *   `DLITTL(A: pMAT; B: pVEC; N: Integer; Var X: pVEC; W: pMAT; Z: pVEC);`
        *   **Purpose**: Solves the linear system $AX=B$ using direct factorization with partial pivoting.
        *   **A** (Input/Output): The $N \times N$ input matrix. It is modified during the factorization process.
        *   **B** (Input/Output): The right-hand side vector $B$. It is also modified during forward substitution.
        *   **N**: The order (dimension) of the linear system.
        *   **X** (Output): The computed solution vector.
        *   **W** (Working Zone/Output): An $N \times N$ matrix used to store the combined $L$ and $U$ factors. Its elements are directly computed by the algorithm.
        *   **Z** (Working Zone/Output): An auxiliary $N$-length vector used to store the intermediate solution of the forward substitution step ($Z = L^{-1}B$).
        *   **Error Handling**: If a zero pivot is encountered during the factorization, it indicates that the matrix is singular, and the procedure prints a message `** DLITTL ** NO UNIQUE SOLUTION` (and does not produce a valid solution).
*   **Associated Files:** `tdlittl.pas` (main program).
*   **Real-world Application Ideas:**
    *   **Boundary Value Problems:** Solving linear systems that arise from the finite difference discretization of ordinary differential equations (ODEs) or partial differential equations (PDEs) in one dimension.
    *   **Control Systems:** Design of state-space controllers often involves solving linear systems derived from system dynamics.
    *   **Chemical Engineering:** Mass and energy balance calculations in process networks can lead to linear systems.

#### Triangularization (Gauss Elimination) (`tlinear.pas`)

This program solves a linear system $AX=B$ by transforming the coefficient matrix $A$ into an upper triangular form using basic Gauss elimination, followed by back-substitution to find the solution. This is a foundational direct method for linear systems.

*   **Numerical Method:** The program directly implements Gauss elimination. In the first phase, it performs a series of elementary row operations to eliminate elements below the main diagonal in each column, effectively reducing the matrix $A$ to an upper triangular form. The same row operations are applied to the right-hand side vector $B$. In the second phase, the transformed upper triangular system is solved for $X$ using back-substitution, starting from the last equation.
*   **Limitations:** This basic implementation of Gauss elimination does *not* explicitly include pivoting (row swaps) to handle zero or very small diagonal elements (`A[K,K]`). If a pivot element is zero, division by zero will occur. If a pivot is very small, numerical instability can lead to inaccurate results. For robustness, a more advanced implementation would incorporate partial or full pivoting.
*   **Associated Files:** `tlinear.pas` (main program), `tlinear.txt` (explanation file).
*   **Real-world Application Ideas:**
    *   **Teaching and Learning:** This program serves as an excellent educational tool to understand the fundamental mechanics of Gauss elimination and back-substitution.
    *   **Small, Well-Conditioned Systems:** For very small systems (e.g., 2x2, 3x3) or matrices known to be well-conditioned and not requiring pivoting, this simple implementation can be used.

#### Solving with Reduced Storage (`dple.pas`)

This program solves a linear system $AX=B$ using a partial pivoting algorithm specifically designed to minimize storage requirements. This is particularly advantageous for very large systems where the full matrix cannot be stored in memory, provided its rows can be generated or accessed on demand.

*   **Numerical Method:** The `dple` procedure implements the Henderson-Wassyng partial pivot algorithm. Unlike traditional methods that store the entire matrix $A$ explicitly, this approach works by calling a user-defined procedure (`rowk`) to retrieve individual rows of $A$ as needed. This "on-the-fly" row access significantly reduces the memory footprint, especially if the matrix has a specific structure (e.g., it's generated by a formula or is extremely sparse with a known pattern) that allows for implicit representation. The algorithm performs transformations to reduce the system and solves it.
*   **Key Procedure:**
    *   `dple(n: INTEGER; b, c: pVEC; VAR ierr: INTEGER);`
        *   **Purpose**: Solves the linear system $AX=B$ with reduced memory storage.
        *   **n**: The dimension (size) of the linear system.
        *   **b**: The right-hand side vector $B$.
        *   **c** (Input/Output): A vector used as internal workspace and which returns the solution vector $X$ upon completion.
        *   **ierr** (Output): Error code: `0` for success, or `k` if the system is detected as singular at step `k` (a zero pivot is encountered).
*   **Required User-Provided Procedure:**
    *   `rowk(n, k: integer; r: pVEC);`
        *   **Purpose**: This procedure defines how the system matrix $A$ is generated. The user must implement this routine to populate the vector `r` with the elements of the `k`-th row of matrix $A$. The `dple` procedure calls `rowk` whenever it needs a specific row of $A$.
*   **Real-world Application Ideas:**
    *   **Very Large-Scale Scientific Computing:** In simulations (e.g., computational physics, materials science) where the system matrix is too enormous to fit in RAM, but its elements can be computed on demand (e.g., from a kernel function or a known sparse pattern).
    *   **Graph Problems:** When the adjacency matrix of a very large graph is needed, but only row-by-row computations are feasible.
    *   **Resource-Constrained Environments:** Solving linear systems on embedded systems or older hardware with limited memory resources.

#### Singular Value Decomposition (SVD) (`tsvbksb.pas`)

This program solves linear systems $AX=B$ using the Singular Value Decomposition (SVD) method. SVD is a powerful and numerically stable factorization technique that decomposes any $m \times n$ matrix $A$ into $U \Sigma V^T$, where $U$ and $V$ are orthogonal matrices, and $\Sigma$ is a diagonal matrix containing singular values. SVD is particularly recommended for ill-conditioned or near-singular matrices, and for overdetermined or underdetermined systems where it naturally provides a least-squares solution.

*   **Numerical Method:** The SVD method involves two main steps:
    1.  **SVD Computation (`svdcmp`):** The input matrix $A$ is factored into $U$, $\Sigma$, and $V^T$. The singular values (diagonal elements of $\Sigma$, stored in vector `w`) provide insights into the matrix's rank and condition number. Small singular values indicate near-singularity.
    2.  **Solution via SVD Components (`svbksb`):** Once the SVD is computed, the solution $X$ for $AX=B$ can be found using the formula $X = V \Sigma^+ U^T B$, where $\Sigma^+$ is the pseudo-inverse of $\Sigma$ (reciprocals of non-zero singular values). A common regularization step is to set singular values below a certain threshold to zero, which effectively removes the contribution of "noisy" or redundant dimensions, leading to a more stable solution for ill-conditioned problems.
*   **Key Procedures:**
    *   `svdcmp(VAR a: pMAT; m, n: Integer; VAR w: pVEC; VAR v: pMAT);`
        *   **Purpose**: Computes the Singular Value Decomposition of matrix `a`.
        *   **a** (Input/Output): The $m \times n$ input matrix $A$. On output, it is replaced by the orthogonal matrix $U$.
        *   **m**: Number of rows in matrix `a`.
        *   **n**: Number of columns in matrix `a`.
        *   **w** (Output): A vector (length $n$) containing the singular values (diagonal elements of $\Sigma$).
        *   **v** (Output): An $n \times n$ matrix containing the orthogonal matrix $V$ (not its transpose $V^T$).
    *   `svbksb(VAR u: pMAT; VAR w: pVEC; VAR v: pMAT; m, n: Integer; b: pVEC; VAR x: pVEC);`
        *   **Purpose**: Solves the linear system $AX=B$ using the SVD components obtained from `svdcmp`.
        *   **u**, **w**, **v** (Input): The $U$, $\Sigma$, and $V$ matrices/vector obtained from `svdcmp`.
        *   **m**, **n**: Dimensions of the original matrix $A$.
        *   **b** (Input): The right-hand side vector $B$.
        *   **x** (Output): The computed solution vector $X$.
*   **Real-world Application Ideas:**
    *   **Data Fitting and Regression:** Solving overdetermined linear systems (more equations than unknowns) to find the best-fit solution in a least-squares sense, especially when the design matrix is ill-conditioned.
    *   **Image Denoising and Compression:** SVD is fundamental to Principal Component Analysis (PCA), which can be used for dimensionality reduction in image data, effectively removing noise and enabling compression.
    *   **Recommender Systems:** Algorithms like Latent Factor Models for collaborative filtering often use SVD to discover underlying patterns in user-item preference matrices.
    *   **Solving Ill-Posed Problems:** In inverse problems (e.g., medical imaging reconstruction, geophysical inversion), SVD provides a stable way to find approximate solutions even when the problem is inherently sensitive to noise.

### 3.2. Iterative Linear System Solvers

Iterative methods start with an initial guess and generate a sequence of increasingly accurate approximations that converge to the true solution. These are often preferred for very large and sparse linear systems where direct methods become computationally infeasible due to high memory or time requirements.

#### Iterative Gauss-Seidel Method (`fseidel.pas`)

This unit implements the iterative Gauss-Seidel method for solving linear systems $AX=B$. It also supports an optional relaxation coefficient ($\omega$), leading to the Successive Over-Relaxation (SOR) method, which can significantly accelerate convergence.

*   **Numerical Method:** Gauss-Seidel is an iterative technique that updates each component of the solution vector using the most recently computed values of the other components in the current iteration. For example, when calculating $x_i^{(k+1)}$, it uses $x_j^{(k+1)}$ for $j < i$ (already computed in the current iteration) and $x_j^{(k)}$ for $j > i$ (from the previous iteration). The SOR method applies a weighted average of the current Gauss-Seidel iterate and the previous iterate, using the relaxation parameter $\omega \in (0, 2)$.
    *   Convergence of Gauss-Seidel/SOR is guaranteed if the matrix $A$ satisfies certain criteria, such as being strictly diagonally dominant, or symmetric and positive definite. The unit can perform checks for row sum criterion, column sum criterion, or Schmidt-v.Mises criterion.
*   **Key Procedure:**
    *   `seidel(crit: integer; n: integer; mat: pMAT; b: pVEC; omega: double; var x: pVEC; var residu: pVEC; var iter: integer; var rc: integer);`
        *   **Purpose**: Solves the linear system $AX=B$ iteratively using the Gauss-Seidel method with relaxation.
        *   **crit**: Integer indicating which convergence criterion check to perform: `1` for row sum, `2` for column sum, `3` for Schmidt-v.Mises. Other values will skip these checks.
        *   **n**: The dimension (size) of the linear system.
        *   **mat** (Input/Output): The coefficient matrix $A$. It is transformed in-place so that its diagonal elements become 1 (by dividing each row by its diagonal element).
        *   **b** (Input/Output): The right-hand side vector $B$. Its elements are adjusted proportionally to the `mat` transformation.
        *   **omega**: The relaxation coefficient. Must be in the range `(0.0, 2.0)`. `omega = 1.0` corresponds to the standard Gauss-Seidel method without relaxation.
        *   **x** (Input/Output): On input, the starting vector for iteration (initial guess). On output, it contains the computed solution vector.
        *   **residu** (Output): The final residual vector, $B - A \cdot X$. It should be close to the zero vector for a converged solution.
        *   **iter** (Output): The number of iterations performed until convergence or maximum iterations reached.
        *   **rc** (Output): Return code: `0` for successful solution, `1` for invalid input parameters (`n < 1` or `omega` out of range), `3` if a diagonal element of `mat` vanishes (making it impossible to transform), `4` if `ITERMAX` iterations are exceeded without convergence, `11` if the row sum criterion is violated, `12` if the column sum criterion is violated, `13` if the Schmidt-v.Mises criterion is violated (indicating potential non-convergence).
*   **Associated Files:** `fseidel.pas` (main unit), `fseidel.txt` (detailed explanation of the algorithm, including convergence criteria and the SOR method), `tseidel.pas` (test program).
*   **Real-world Application Ideas:**
    *   **Finite Difference Method for PDEs:** Widely used to solve discretized elliptic partial differential equations (e.g., Laplace's equation, Poisson equation) in fields like heat transfer, fluid dynamics, and electromagnetics. The resulting matrices are often sparse and diagonally dominant.
    *   **Image Processing:** Can be applied to image reconstruction or denoising problems that can be formulated as sparse linear systems.
    *   **Power Flow Analysis:** Solving power flow equations in electrical grids.

#### Conjugate Gradient Method (`cgtst1.pas`)

This program demonstrates the Conjugate Gradient (CG) method for solving linear systems $AX=Y$ where $A$ is a **symmetric positive definite matrix**. The CG method is an iterative algorithm particularly well-suited for very large sparse systems, where direct methods become computationally infeasible.

*   **Numerical Method:** The Conjugate Gradient method is an iterative optimization algorithm that minimizes a quadratic function whose minimum corresponds to the solution of the linear system. It generates a sequence of "conjugate" search directions, ensuring that the algorithm theoretically converges to the exact solution in at most $N$ steps (where $N$ is the matrix size) in exact arithmetic. In practice, due to floating-point errors, it is stopped when the residual (error) is sufficiently small. It does not require storing the full matrix explicitly if matrix-vector products can be computed on demand.
*   **Key Procedure (within `cgtst1.pas`):**
    *   `cg_method(n: INTEGER; a: MAT; y: VEC; VAR x: VEC; VAR code: INTEGER);`
        *   **Purpose**: Solves the linear system $AX=Y$ using the Conjugate Gradient method.
        *   **n**: The size of the linear system.
        *   **a**: The coefficient matrix $A$. For symmetric matrices, typically only the upper triangle (including the diagonal) is needed. The matrix is assumed to be symmetric and positive definite.
        *   **y**: The right-hand side vector $Y$.
        *   **x** (Output): The computed solution vector. The method starts with an initial guess (often a zero vector) and refines it.
        *   **code** (Output): Return code: `0` for success, `1` for invalid input parameters (e.g., `n < 2`).
*   **Real-world Application Ideas:**
    *   **Computational Physics and Engineering:** Solving large-scale systems arising from the finite element or finite difference discretization of partial differential equations (PDEs), such as those found in heat transfer, fluid dynamics, solid mechanics, and electromagnetics.
    *   **Machine Learning:** Training optimization problems like Support Vector Machines (SVMs), or large-scale least squares problems.
    *   **Image Processing:** Iterative image reconstruction algorithms often lead to large, sparse, symmetric positive definite linear systems.
    *   **Optimal Control:** Solving the linear systems that arise in iterative methods for optimal control problems.

#### Conjugate Gradient for Sparse Systems (`tsparse.pas`)

This program provides a specific implementation of the Conjugate Gradient (CG) method tailored for solving **sparse symmetric linear systems** ($AX=B$). It highlights the efficiency gains achieved by explicitly exploiting the sparsity pattern of the matrix, avoiding the storage and processing of numerous zero elements.

*   **Numerical Method:** This implementation of the CG method is optimized for sparse matrices. Instead of requiring the full matrix $A$ in memory, it relies on user-defined subroutines (`ASUB` and `ATSUB`) to compute matrix-vector products ($AX$ and $A^TX$) on demand. These subroutines are designed to directly access and multiply only the non-zero elements of $A$, which dramatically reduces both memory consumption and computational time for very large sparse systems. The core CG algorithm then proceeds as usual, leveraging these efficient product calculations.
*   **Key Procedure:**
    *   `SPARSE(B: VECT; N: INTEGER; VAR X: VECT; VAR RSQ: REAL);`
        *   **Purpose**: Solves the sparse linear system $AX=B$ using the Conjugate Gradient method.
        *   **B**: The right-hand side vector $B$.
        *   **N**: The size (dimension) of the linear system.
        *   **X** (Input/Output): On input, an initial guess for the solution (e.g., a vector of all zeros). On output, it contains the computed solution vector $X$.
        *   **RSQ** (Output): The sum of squares of the components of the residual vector $(A \cdot X - B)$, i.e., $||AX-B||_2^2$. If this value is not small, it may indicate that the matrix is numerically singular, and the solution represents a least-squares approximation.
*   **Required User-Provided Subroutines (implementing matrix-vector products for the specific sparse matrix):**
    *   `ASUB(XIN: VECT; VAR XOUT: VECT);`
        *   **Purpose**: Calculates the product $A \cdot XIN$ and stores the result in `XOUT`. This routine must be custom-written by the user to reflect the exact non-zero pattern and values of their sparse matrix $A$.
    *   `ATSUB(XIN: VECT; VAR XOUT: VECT);`
        *   **Purpose**: Calculates the product $A^T \cdot XIN$ (transpose of $A$ times $XIN$) and stores the result in `XOUT`. For symmetric matrices, $A^T=A$, so `ATSUB` would be identical to `ASUB`.
*   **Real-world Application Ideas:**
    *   **Large-scale Network Analysis:** Solving for flows or potentials in massive networks (e.g., traffic systems, social networks, power grids) where the connectivity matrix is extremely sparse.
    *   **Computational Fluid Dynamics (CFD):** High-fidelity simulations of fluid flow often involve discretizing PDEs on millions or billions of grid points, leading to enormous sparse linear systems.
    *   **Computer Graphics (Global Illumination):** Algorithms like Radiosity, which calculate light distribution in complex scenes, often require solving very large sparse systems.
    *   **Optimization in Machine Learning:** Some machine learning algorithms, particularly those based on large datasets and sparse feature representations, lead to sparse linear systems in their optimization steps.

### 3.3. Symmetric Linear Systems

These methods are specialized for matrices where $A = A^T$. They leverage this symmetry to reduce computational cost and memory usage, and often offer enhanced numerical stability.

#### Gauss Method (`syslin.pas`)

This program solves a symmetric linear system ($AX=B$) using a modified Gauss elimination method that takes advantage of the matrix's symmetry to potentially optimize calculations. While the underlying algorithm is Gauss elimination, its adaptation for symmetric matrices can lead to efficiencies.

*   **Numerical Method:** The `RSLSG` procedure implements a form of Gauss elimination tailored for symmetric matrices. It uses row operations to transform the matrix into an upper triangular form. Due to symmetry, properties of the lower triangular part can be implicitly used. The method ensures a non-zero determinant for a unique solution.
*   **Key Procedure:**
    *   `RSLSG(eps: REAL; n: INTEGER; A: MAT; VAR it: INTEGER; VAR X: VEC);`
        *   **Purpose**: Solves a symmetric linear system.
        *   **eps**: A real value representing the desired numerical precision. Used to check if diagonal elements (pivots) are too close to zero, indicating potential singularity.
        *   **n**: The size (dimension) of the linear system.
        *   **A** (Input/Output): The extended matrix of the system, an $N \times (N+1)$ matrix where the first $N$ columns are the coefficient matrix $A$, and the last column is the right-hand side vector $B$. For symmetric matrices, typically only the upper triangular part (including the diagonal) and the right-hand side are explicitly read, as the lower triangular part can be inferred. The matrix `A` is modified in-place during the elimination process.
        *   **it** (Output): An error indicator: `1` if the system is solved successfully, `0` if the system matrix is singular (a pivot element is too close to zero).
        *   **X** (Output): The computed solution vector.
*   **Associated Files:** `syslin.pas` (main program), `matsym.dat` (example input data).
*   **Real-world Application Ideas:**
    *   **Structural Mechanics (Equilibrium Equations):** Solving for forces and displacements in structures where the stiffness matrix is symmetric.
    *   **Heat Transfer:** Analyzing steady-state temperature distributions in materials with symmetric thermal conductivity properties.
    *   **Least Squares Problems:** Solving the normal equations $A^TAx = A^Tb$, where $A^TA$ is always symmetric and positive semi-definite.

#### Cholesky Decomposition (`fcholy.pas`)

This unit provides robust and efficient procedures for solving linear systems and performing matrix inversion when the coefficient matrix is **symmetric and positive definite**. Cholesky decomposition factors such a matrix $A$ into the product of a lower triangular matrix $L$ and its transpose $L^T$ (i.e., $A = LL^T$). This method is numerically stable for positive definite matrices and generally faster than LU decomposition for this specific class of problems.

*   **Numerical Method:** The Cholesky decomposition algorithm directly computes the elements of the lower triangular matrix $L$. Once $L$ is found, solving $AX=B$ (or $LL^TX=B$) is reduced to two successive triangular solves: first $Ly = B$ (forward substitution), then $L^T x = y$ (back substitution). The method requires all principal minors of the matrix to be positive, ensuring that square roots of positive numbers are always taken.
*   **Key Procedures (from `fcholy.pas`):**
    *   `Choly(mode, n: INTEGER; VAR mat: pMAT; b: pVEC; VAR x: pVEC; VAR rc: INTEGER);`
        *   **Purpose**: A high-level interface that orchestrates the Cholesky decomposition and/or solution process.
        *   **mode**: Control parameter:
            *   `0`: Perform Cholesky decomposition of `mat` and then solve the linear system $AX=B$.
            *   `1`: Only perform the Cholesky decomposition of `mat` (useful if only factorization is needed).
            *   `2`: Solve the system $AX=B$ assuming `mat` has already been decomposed into its Cholesky factor (efficient for multiple $B$ vectors).
        *   **n**: The dimension (size) of the matrix `mat`.
        *   **mat** (Input/Output): On input, the symmetric positive definite matrix. On output, if `mode` is `0` or `1`, it is replaced by its Cholesky factor (typically the lower triangular part).
        *   **b** (Input): The right-hand side vector $B$.
        *   **x** (Output): The computed solution vector $X$ if `mode` is `0` or `2`.
        *   **rc** (Output): Return code: `0` for success, `1` for invalid input parameters (e.g., `n < 1` or null pointers), `2` if the matrix is found *not* to be positive definite (e.g., a negative number encountered during square root), `3` for an invalid `mode` value.
    *   `chodec(n: INTEGER; VAR mat: pMAT; VAR rc: INTEGER);`
        *   **Purpose**: Performs the Cholesky decomposition of the symmetric matrix `mat`.
        *   **mat** (Input/Output): The symmetric matrix. Replaced by its lower triangular Cholesky factor.
        *   **rc**: Return code (0=ok, 1=invalid n, 2=not positive definite).
    *   `chosol(n: INTEGER; lmat: pMAT; b: pVEC; VAR x: pVEC; VAR rc: INTEGER);`
        *   **Purpose**: Solves the linear system $LL^TX=B$ using the pre-computed Cholesky factor `lmat`.
        *   **lmat**: The lower triangular Cholesky factor obtained from `chodec`.
        *   **b**: Right-hand side vector.
        *   **x**: Solution vector.
        *   **rc**: Return code (0=ok, 1=improper matrix or invalid n).
*   **Matrix Storage:** Uses `pMAT` for matrices, where the underlying array is 0-indexed (`Array[0..SIZE,0..SIZE]`).
*   **Associated Files:** `fcholy.pas` (main unit), `choles.pas` (demonstrates matrix inversion and determinant using Cholesky), `tcholy.pas` (test program for solving $AX=B$).
*   **Real-world Application Ideas:**
    *   **Numerical Optimization:** Solving systems of linear equations within Newton's method for convex optimization problems, where Hessian matrices are often symmetric positive definite.
    *   **Financial Modeling:** In risk management and portfolio optimization, covariance matrices (which are symmetric positive definite) are frequently decomposed using Cholesky for various calculations.
    *   **Kalman Filtering:** In state estimation problems, the covariance matrix of the state estimate is updated using Cholesky factorization.
    *   **Spatial Statistics (Kriging):** Used in geostatistics for spatial interpolation, involving the inversion of covariance matrices.

#### SYMSOL Method (`tsymsol.pas`)

This program solves a symmetric linear system $AX=B$ using the `SYMSOL` subroutine. This method is specifically designed for symmetric matrices that are **not necessarily positive definite**, making it more general than Cholesky decomposition for symmetric systems.

*   **Numerical Method:** The `SYMSOL` procedure is based on the decomposition of the symmetric matrix $A$ into a product of three matrices: $A = U^T D U$, where $U$ is an upper triangular matrix with unit diagonal elements, and $D$ is a diagonal matrix. This is a variant of the LDLT decomposition. The method computes the factors by direct calculation, then uses forward and backward substitution to solve the system. It can handle indefinite symmetric matrices, provided they are well-balanced.
*   **Key Procedure:**
    *   `SYMSOL(N: Integer; A: pMAT; B: pVEC; Var X: pVEC; D: pVEC);`
        *   **Purpose**: Solves a symmetric linear system $AX=B$ using the $U^TDU$ factorization.
        *   **N**: The number of lines (dimension) of the linear system.
        *   **A** (Input/Output): The symmetric coefficient matrix. Typically, only the upper triangular half (including the diagonal) needs to be initialized. Its values are *overwritten* (lost) during the factorization process.
        *   **B** (Input): The right-hand side vector.
        *   **X** (Output): The computed solution vector.
        *   **D** (Working Zone/Output): A vector used to store the diagonal elements of the decomposition matrix $D$. It serves as a working area during computation.
*   **Associated Files:** `tsymsol.pas` (main program).
*   **Real-world Application Ideas:**
    *   **Constrained Optimization:** Solving KKT (Karush-Kuhn-Tucker) systems that arise from optimization problems with equality constraints. These systems often involve symmetric but indefinite matrices.
    *   **Structural Mechanics (Buckling Analysis):** In some advanced structural analyses (e.g., stability or buckling problems), the stiffness matrix can become indefinite.
    *   **Acoustics and Vibrations:** Modeling complex systems where the governing equations lead to symmetric matrices that might not always be positive definite.

### 3.4. Specialized Linear Systems

These methods are designed to take advantage of specific matrix structures (e.g., banded, tridiagonal, Vandermonde, Toeplitz) for significant gains in computational efficiency and reduced memory footprint compared to general solvers.

#### Banded Matrices (with Pivots) (`fband.pas`)

This unit provides routines for solving linear systems $A \cdot x = b$ where the coefficient matrix $A$ is **banded**. A banded matrix has non-zero elements only along a central "band" around the main diagonal, with `ld` lower co-diagonals and `ud` upper co-diagonals. The `band` procedure uses Gauss elimination with column pivot search to ensure numerical stability.

*   **Numerical Method:** The `band` procedure performs an LU factorization tailored for banded matrices. It operates only on the elements within the band, greatly reducing computations and memory. Column pivoting is incorporated to handle potential zero or small pivots, enhancing numerical stability. While pivoting helps stability, it can cause "fill-in" (generation of new non-zero elements outside the original band width), which the implementation must accommodate. The decomposed matrix is stored in a condensed format that efficiently represents the band structure.
*   **Matrix Storage:** Banded matrices are stored in a **condensed format**, typically as a 1D `pVECT`. The element $A[i, k]$ of the original $N \times N$ matrix is mapped to an index in the `pVECT` `pmat` based on its row `i` and its position relative to the main diagonal `k-i`. For example, $A[i,k]$ might be stored at `pmat^[i * (ld+ud+1) + ld + k - i]`. The actual storage details depend on the internal implementation of the unit, but this condensed form significantly reduces memory consumption for large sparse matrices.
*   **Key Procedure:**
    *   `band(mode: INTEGER; n: INTEGER; ld: INTEGER; ud: INTEGER; VAR pmat: pVECT; VAR b: pVECT; VAR perm: IVECT; VAR signd: INTEGER; VAR code: INTEGER);`
        *   **Purpose**: The main interface to factorize and/or solve banded linear systems with pivoting.
        *   **mode**: Control parameter: `0` (factor and solve), `1` (factor only), `2` (solve only).
        *   **n**: The dimension (size) of the linear system.
        *   **ld**: The number of lower co-diagonals.
        *   **ud**: The number of upper co-diagonals.
        *   **pmat** (Input/Output): The banded input matrix (in condensed form). If `mode` is `0` or `1`, it is replaced by its LU factors. If `mode` is `2`, it must contain the pre-computed LU factors.
        *   **b** (Input/Output): The right-hand side vector $B$. On output, it contains the solution vector $X$.
        *   **perm** (Output): An integer vector storing the row permutations performed during pivoting.
        *   **signd** (Output): The sign of the permutation (for determinant calculation).
        *   **code** (Output): Return code: `0` for success, `1` for invalid parameters, `3` if the matrix is numerically singular (e.g., a pivot is too small), `4` for a wrong `mode` value.
    *   **Supporting Procedures:** `banddec` (performs the factorization) and `bandsol` (solves using the factors) are internal procedures typically called by `band`.
*   **Associated Files:** `fband.pas` (main unit), `tband.pas` (test program).
*   **Real-world Application Ideas:**
    *   **Finite Difference Methods for PDEs:** Discretization of 1D or 2D Partial Differential Equations (PDEs) often yields banded matrices (e.g., heat conduction, wave propagation problems on a grid).
    *   **Spline Interpolation:** Constructing cubic splines involves solving banded linear systems.
    *   **Structural Mechanics (Beams & Trusses):** Modeling long, slender structures where connectivity between nodes is local, leading to banded stiffness matrices.

#### Banded Matrices (without Pivots) (`fbando.pas`)

This unit provides routines for solving banded linear systems without performing explicit pivot searches. This approach can be faster than methods with pivoting, but it is less numerically stable and requires the matrix to have properties that prevent zero or very small pivots (e.g., strict diagonal dominance).

*   **Numerical Method:** This unit implements Gauss elimination for banded matrices, similar to `fband.pas`, but it skips the row permutation step (pivoting). This simplification means that the diagonal elements are used directly as pivots. The factorization is performed within the condensed storage format of the banded matrix. While this avoids the overhead of pivot searching and potential fill-in outside the original band, it will fail if a pivot becomes zero or very small.
*   **Key Procedure:**
    *   `bando(mode: INTEGER; n: INTEGER; ld: INTEGER; ud: INTEGER; VAR pmat: pVECT; VAR b: pVECT; VAR code: INTEGER);`
        *   **Purpose**: The main interface to factorize and/or solve banded linear systems without pivoting.
        *   **mode**: Control parameter: `0` (factor and solve), `1` (factor only), `2` (solve only).
        *   **n, ld, ud**: System size, number of lower and upper co-diagonals, respectively.
        *   **pmat** (Input/Output): The banded input matrix (in condensed form). If `mode` is `0` or `1`, it is replaced by its LU factors. If `mode` is `2`, it must contain the pre-computed LU factors.
        *   **b** (Input/Output): The right-hand side vector $B$. On output, it contains the solution vector $X$.
        *   **code** (Output): Return code: `0` for success, `1` for invalid parameters, `2` if the matrix is numerically singular (a diagonal pivot element is too close to zero), `3` for a wrong `mode` value.
*   **Associated Files:** `fbando.pas` (main unit), `tband.pas` (test program).
*   **Real-world Application Ideas:**
    *   **Real-time Control Systems:** In applications where speed is critical and the system matrix is known to be well-conditioned and diagonally dominant.
    *   **Certain Discretized PDEs:** For problems where the discretized matrix is naturally stable (e.g., certain finite difference schemes for diffusion equations).

#### Non-Symmetric Banded Matrices (`nsbslv.pas`)

This program provides procedures for the LU factorization and subsequent solving of **non-symmetric banded linear systems**. It is a specialized adaptation of general LU decomposition techniques for matrices with a banded structure, offering efficiency over full matrix methods.

*   **Numerical Method:** The `NSBFAC` procedure performs an LU factorization with partial pivoting for a non-symmetric banded matrix. It stores the matrix in a compact banded format, operating only on the elements within the band to save memory and computations. The `NSBSLV` procedure then uses the computed LU factors and pivoting information to efficiently solve the linear system using forward and backward substitution.
*   **Matrix Storage:** The banded matrix is stored in a specific condensed format within a 2D array (`pMAT`). The mapping between the full matrix $A[i,j]$ and the banded storage $B[k,j]$ is defined such that $k = i-j+m$, where $m = ml+mu+1$. This is a common way to represent banded matrices compactly.
*   **Key Procedures:**
    *   `NSBFAC(Var B: pMAT; N, ML, MU: Integer; Var IPVT: pIVEC; Var IND:Integer);`
        *   **Purpose**: Performs LU factorization of a banded matrix `B` with partial pivoting.
        *   **B** (Input/Output): The banded matrix in its specific condensed format. On output, it is replaced by its LU elements.
        *   **N**: The size (dimension) of the matrix.
        *   **ML**: The number of lower diagonals.
        *   **MU**: The number of upper diagonals.
        *   **IPVT** (Output): An integer vector (size `N`) storing the pivoting indices used during factorization.
        *   **IND** (Output): An error flag: `0` for a non-singular matrix, `k` if the matrix may be singular at step `k` (a zero pivot is encountered).
    *   `NSBSLV(A: pMAT; N, ML, MU: Integer; IPVT: pIVEC; B: pVEC; Var X: pVEC);`
        *   **Purpose**: Solves the banded linear system using the LU factors obtained from `NSBFAC`.
        *   **A** (Input): The banded matrix containing the LU factors (output from `NSBFAC`).
        *   **N, ML, MU**: Dimensions of the system.
        *   **IPVT** (Input): The integer vector of pivoting indices from `NSBFAC`.
        *   **B** (Input): The right-hand side vector $B$.
        *   **X** (Output): The computed solution vector $X$.
*   **Associated Files:** `nsbslv.pas` (main program).
*   **Real-world Application Ideas:**
    *   **Heat Transfer and Fluid Dynamics:** Modeling convection-diffusion problems or non-symmetric discretizations of PDEs that result in non-symmetric banded matrices.
    *   **Numerical Methods for ODEs:** Some implicit methods for Ordinary Differential Equations lead to banded linear systems at each time step.

#### Tridiagonal Systems (`tridiag.pas`)

This program specifically solves **tridiagonal linear systems** ($MU=R$). A tridiagonal matrix is a very sparse type of banded matrix where non-zero elements exist only on the main diagonal and its two adjacent co-diagonals (one upper and one lower). This highly specific structure allows for extremely efficient solution algorithms.

*   **Numerical Method:** The `TRIDAG` procedure implements the Thomas algorithm (also known as the Tridiagonal Matrix Algorithm, or TDMA). This is a simplified and highly optimized form of Gaussian elimination that takes full advantage of the tridiagonal structure. It operates in $O(N)$ time (linear complexity with respect to matrix size $N$), which is significantly faster than the $O(N^3)$ required for general dense systems. The algorithm proceeds in two phases: a forward elimination phase to modify the right-hand side and the diagonal elements, and a backward substitution phase to compute the solution.
*   **Key Procedure:**
    *   `TRIDAG(A: pVECT; B: pVECT; C: pVECT; R: pVECT; VAR U: pVECT; N: INTEGER; VAR CODE: INTEGER);`
        *   **Purpose**: Solves a tridiagonal linear system.
        *   **A**: A vector representing the elements of the **lower subdiagonal**. The `i`-th element `A[i]` corresponds to $M_{i, i-1}$.
        *   **B**: A vector representing the elements of the **main diagonal**. The `i`-th element `B[i]` corresponds to $M_{i, i}$.
        *   **C**: A vector representing the elements of the **upper subdiagonal**. The `i`-th element `C[i]` corresponds to $M_{i, i+1}$.
        *   **R**: The right-hand side vector $R$.
        *   **U** (Output): The computed solution vector $U$.
        *   **N**: The size (dimension) of the linear system.
        *   **CODE** (Output): Return code: `0` for success, `1` if `B[1]` (the first pivot) is zero, `2` if a division by zero occurs during the factorization process (indicating a singular matrix).
*   **Associated Files:** `tridiag.pas` (main program).
*   **Real-world Application Ideas:**
    *   **Numerical Solution of PDEs:** The finite difference discretization of one-dimensional (1D) diffusion, heat conduction, or wave equations often results directly in tridiagonal systems.
    *   **Spline Interpolation:** The coefficients of cubic splines, used for smooth curve fitting, are typically found by solving tridiagonal systems.
    *   **Quantum Mechanics:** Solving the one-dimensional time-independent Schrödinger equation using finite difference methods often leads to tridiagonal eigenvalue problems, which are closely related to solving tridiagonal linear systems.

#### Vandermonde Systems (`tvander.pas`)

This program solves **Vandermonde linear systems**. A Vandermonde matrix has a specific structure where each column consists of powers of a given set of values (nodes). These matrices commonly arise in polynomial interpolation problems.

*   **Numerical Method:** The `Vander` procedure solves the Vandermonde system without explicitly constructing the matrix. It employs an efficient algorithm that leverages the unique structure of Vandermonde matrices. The method effectively calculates the coefficients of the interpolating polynomial. It typically involves a recursive calculation of the coefficients of a "master polynomial" whose roots are the given nodes, followed by synthetic division to find the solution vector.
*   **Key Procedure:**
    *   `Vander(X: VEC; VAR W: VEC; Q: VEC; N: Integer);`
        *   **Purpose**: Solves a Vandermonde linear system.
        *   **X**: An input vector of `N` distinct points (nodes), $x_1, x_2, \dots, x_N$. These are the base values for the powers in the matrix.
        *   **Q**: The right-hand side vector, $q_1, q_2, \dots, q_N$, representing the function values at the nodes (i.e., $f(x_i) = q_i$).
        *   **W** (Output): The solution vector, $w_1, w_2, \dots, w_N$, which represents the coefficients of the interpolating polynomial.
        *   **N**: The size (dimension) of the linear system.
*   **Real-world Application Ideas:**
    *   **Polynomial Interpolation:** The primary and most direct application. Given a set of data points $(x_i, y_i)$, the `Vander` procedure can find the unique polynomial $P(x) = w_1 + w_2 x + \dots + w_N x^{N-1}$ such that $P(x_i) = y_i$.
    *   **Numerical Differentiation and Integration:** Developing formulas for numerical differentiation or integration can involve solving Vandermonde systems.
    *   **Coding Theory:** Constructing certain error-correcting codes, such as Reed-Solomon codes, involves properties of Vandermonde matrices.

#### Toeplitz Systems (`toeplitz.pas`)

This program solves **Toeplitz linear systems**. A Toeplitz matrix is a special type of matrix where each descending diagonal from left to right has constant values. These matrices frequently appear in fields like signal processing and time series analysis.

*   **Numerical Method:** The `TOEPLZ` procedure solves Toeplitz systems using the **Levinson-Durbin algorithm**. This is a highly efficient recursive algorithm specifically designed for Toeplitz matrices. It exploits the matrix's repetitive structure to achieve a computational complexity of $O(N^2)$ (for an $N \times N$ matrix), which is significantly faster than the $O(N^3)$ required for general direct linear solvers. The algorithm builds the solution iteratively.
*   **Key Procedure:**
    *   `TOEPLZ(R: VEC1; VAR X: VEC; Y: VEC; N: Integer);`
        *   **Purpose**: Solves a Toeplitz linear system.
        *   **R**: An input vector (of size $2N-1$) that contains all the unique elements defining the Toeplitz matrix. `R[N]` represents the main diagonal, `R[N-1]` and `R[N+1]` represent the first lower and upper off-diagonals, and so on.
        *   **Y**: The right-hand side vector.
        *   **X** (Output): The computed solution vector.
        *   **N**: The size (dimension) of the linear system.
        *   **Error Handling**: The procedure prints a message `Levinson method fails: singular principal minor.` if it encounters a singular principal minor during the computation, indicating that the Toeplitz matrix is singular.
*   **Real-world Application Ideas:**
    *   **Digital Signal Processing (DSP):** Fundamental for applications like linear prediction, Wiener filtering, and autoregressive (AR) model estimation (e.g., Yule-Walker equations), where covariance matrices often exhibit a Toeplitz structure.
    *   **Time Series Analysis:** Modeling and forecasting time-dependent data, such as financial time series or climate data.
    *   **Image Restoration:** Deblurring and deconvolution problems in image processing can often be formulated as Toeplitz systems.

---

## 4. Matrix Operations

This section describes routines for fundamental matrix operations that go beyond directly solving linear systems, including matrix inversion and determinant calculation. These operations provide key insights into a matrix's properties and are building blocks for more complex algorithms.

### 4.1. Matrix Inversion

Matrix inversion involves finding a matrix $A^{-1}$ such that $AA^{-1} = A^{-1}A = I$ (the identity matrix). While often avoided in favor of solving linear systems ($AX=B$ is usually more stable and efficient than $X=A^{-1}B$), explicit matrix inversion is necessary for specific applications.

#### Inversion by LU Decomposition (`inv_lu.pas`)

This program calculates the inverse of a real square matrix using its LU decomposition. It leverages the `LUDCMP` and `LUBKSB` routines provided in the `lu.pas` unit.

*   **Numerical Method:** To find the inverse of an $N \times N$ matrix $A$, one can solve $N$ distinct linear systems: $AX_j = E_j$, where $E_j$ is the $j$-th column of the identity matrix $I$. Each solution vector $X_j$ then forms the $j$-th column of the inverse matrix $A^{-1}$. This approach is efficient with LU decomposition because $A$ needs to be decomposed only once using `LUDCMP`. Subsequently, the `LUBKSB` procedure can be called $N$ times, each time with a different column of $I$ as the right-hand side.
*   **Key Procedures:**
    *   Relies on `LUDCMP(A, N, INDX, D, rc)` from `lu.pas` for factorization.
    *   Relies on `LUBKSB(A, N, INDX, B)` from `lu.pas` for solving each column of the inverse.
    *   `MatMult(A, B: pVect; VAR C: pVect; n: Integer);` (defined within `inv_lu.pas`): Multiplies two square real matrices (stored in vectors) for verification purposes.
*   **Real-world Application Ideas:**
    *   **Covariance Matrix Inversion:** In multivariate statistics and Kalman filtering, inverting covariance matrices is a common operation.
    *   **Control Theory (Linear Quadratic Regulator):** Some formulations of the LQR problem involve explicit matrix inversion.
    *   **Numerical Optimization:** Calculating the inverse of the Hessian matrix in Newton-Raphson methods (though direct solution of the linear system is often preferred for stability).

#### Inversion by Householder's Method (`househol.pas`)

This program calculates the inverse of a real square matrix using a sequence of Householder transformations. Householder reflections are orthogonal transformations that can be used to zero out elements below the diagonal, typically transforming a matrix into an upper triangular (QR decomposition) or bidiagonal form.

*   **Numerical Method:** The algorithm applies Householder reflections to the augmented matrix $[A | I]$. Each Householder transformation is chosen to zero out elements in a column below the diagonal. By applying these transformations iteratively, $A$ is reduced to an upper triangular matrix $R$, and $I$ is simultaneously transformed into a matrix $Q$ (where $Q$ is the accumulated Householder transformations applied to $I$). The problem then becomes solving $RX = Q$ for $X$, which is a straightforward back-substitution problem since $R$ is triangular.
*   **Key Procedures (internal to `househol.pas`):**
    *   `P1000`: Performs the Householder transformations on the augmented matrix $[A | I]$ to reduce $A$ to upper triangular form.
    *   `P2000`: Solves a triangular system using back-substitution. This is called for each column of the transformed identity matrix $Q$ to find the columns of $A^{-1}$.
    *   `P3000`: Performs matrix multiplication (e.g., $A \cdot A^{-1}$) for verification.
    *   `P4000`: Calculates the determinant from the diagonal elements of the triangularized matrix.
*   **Real-world Application Ideas:**
    *   **QR Decomposition:** The Householder method is the standard way to compute the QR decomposition, which is highly stable for solving linear least squares problems.
    *   **Eigenvalue Problems:** Householder transformations are used as a preprocessing step to reduce general symmetric matrices to tridiagonal form for more efficient eigenvalue computation.
    *   **Numerical Stability in Linear Algebra:** Provides an alternative, numerically stable method for matrix inversion, particularly useful for benchmarking other inversion techniques.

#### Inversion of a Symmetric Positive Definite Matrix by Cholesky Method (`choles.pas`)

This program calculates the inverse of a symmetric positive definite matrix using Cholesky decomposition. This method is highly efficient and numerically stable for this specific class of matrices.

*   **Numerical Method:** The core idea is to first decompose the symmetric positive definite matrix $A$ into $LL^T$ using Cholesky decomposition. Then, the inverse $A^{-1}$ can be obtained by solving $N$ linear systems ($LL^T X_j = E_j$, where $E_j$ are columns of the identity matrix) or by a more direct formula involving operations on the Cholesky factor $L$. The program leverages auxiliary procedures to perform these steps.
*   **Key Procedures:**
    *   `choldc1(n:integer; Var a:pMAT; Var p:pVEC)`: The main internal method for Cholesky decomposition.
    *   `choldc(n:integer; A:pMAT; Var aa:pMAT)`: A wrapper for `choldc1` for the decomposition step.
    *   `choldcsl(n:integer; A:pMAT; Var aa:pMAT)`: Computes the inverse of the lower decomposed matrix $L$.
    *   `cholsl(n:integer; A:pMAT; Var aa:pMAT)`: Computes the full inverse of $A$ using its Cholesky decomposition.
    *   `choldet(n:integer; a:pMAT): Double`: Computes the determinant of $A$ from its Cholesky factor ($(\prod L_{ii})^2$).
    *   `Check_Matrix(n:integer;A:pMAT):Boolean`: Checks if a matrix is symmetric positive definite.
    *   `MATMULT(n:Integer;A,B:pMAT; VAR C: pMAT)`: Multiplies two matrices (for verification).
*   **Associated Files:** `choles.pas` (main program), `fcholy.pas` (unit containing Cholesky routines for solving linear systems, related to this program but this one focuses on inversion).
*   **Real-world Application Ideas:**
    *   **Gaussian Process Regression:** Inverting kernel matrices for Gaussian processes, which are symmetric positive definite.
    *   **Least Squares Optimization:** Solving the normal equations in least squares problems by explicitly inverting the Gram matrix ($A^T A$).
    *   **Statistical Inference:** Inverting Fisher information matrices or covariance matrices in statistical models.

### 4.2. Determinant Calculation

The determinant of a square matrix is a scalar value that encapsulates several fundamental properties of the matrix, including its invertibility and the scaling factor of the linear transformation it represents. The library offers various methods for its computation.

#### Determinant by Gauss Method (`deter.pas`)

This program calculates the determinant of a real square matrix using Gauss elimination with full pivoting.

*   **Numerical Method:** The core principle is that the determinant of a triangular matrix (upper or lower) is simply the product of its diagonal elements. Gauss elimination transforms a matrix into an upper triangular form through a series of elementary row and column operations.
    *   **Full Pivoting:** At each step, the largest absolute value element in the remaining submatrix is chosen as the pivot. This element is moved to the diagonal position by row and column swaps. Each swap changes the sign of the determinant, which must be tracked.
    *   The determinant is then the product of the diagonal elements of the final triangular matrix, multiplied by a sign factor that accounts for the total number of row/column exchanges.
*   **Key Procedures:**
    *   `TSRGT(eps: real; n: integer; A: pMat; VAR it: integer; VAR C: pMat; VAR Kp, Lp: pVecI);`
        *   **Purpose**: Performs upper triangularization of matrix `A` using Gauss elimination with full pivoting.
        *   **eps**: Precision used to check for numerically zero pivots.
        *   **n**: Size of the matrix.
        *   **A** (Input): The input matrix.
        *   **it** (Output): Flag: `1` if successful, `0` if the matrix is singular.
        *   **C** (Output): The matrix containing the upper triangular form (elements of $L$ and $U$ factors, and diagonal from $U$).
        *   **Kp, Lp** (Output): Integer vectors storing information about column and row permutations, respectively.
    *   `DMGT(eps: real; n: integer; VAR A: pMat): real;`
        *   **Purpose**: Calculates the determinant of matrix `A` using `TSRGT`.
        *   Returns: The determinant value. Returns `0.0` if `A` is singular.
*   **Real-world Application Ideas:**
    *   **Geometric Scaling:** In transformations (e.g., in computer graphics), the determinant indicates how much a volume (in 3D) or area (in 2D) is scaled by the transformation.
    *   **Linear Independence Check:** A non-zero determinant implies that the columns (or rows) of the matrix are linearly independent.
    *   **System Stability:** In control theory, the determinant of certain matrices related to system dynamics can be used to assess stability.

#### Determinant by LU Decomposition Method (`deter1.pas`)

This program calculates the determinant of a real square matrix by first performing its LU decomposition. This method is computationally efficient, especially if the LU decomposition is already available (e.g., from solving a linear system).

*   **Numerical Method:** For a matrix $A$ that can be decomposed as $PA=LU$ (where $P$ is a permutation matrix, $L$ is lower triangular with unit diagonal, and $U$ is upper triangular), the determinant property states that $\det(A) = \det(P) \det(L) \det(U)$. Since $\det(L)=1$ (as it has ones on its diagonal), the formula simplifies to $\det(A) = \det(P) \det(U)$. The determinant of $U$ is simply the product of its diagonal elements. The value $\det(P)$ is $\pm 1$, determined by the parity of row swaps performed during `LUDCMP` (stored in the `D` parameter).
*   **Key Procedures:**
    *   Relies on `LUDCMP(A, n, INDX, D, icode)` from `lu.pas`. `LUDCMP` performs the decomposition and directly provides `D` (the sign of the permutation) and the diagonal elements of $U$ (stored in `A^[i*(n+1)+i]`).
    *   The program then calculates the determinant as `D * product(diagonal elements of U)`.
*   **Real-world Application Ideas:**
    *   **Complement to LU Solvers:** If a system $AX=B$ is solved using LU decomposition, calculating its determinant comes at a minimal additional cost.
    *   **Condition Number Estimation:** Though more accurately done with SVD, the determinant's value can provide a rough idea of how "well-conditioned" a matrix is (matrices with determinants close to zero are ill-conditioned).

#### Determinant by a Recursive Method based on Kramer's rule (`deter2.pas`)

This program calculates the determinant of a real square matrix using a recursive function that implements cofactor expansion (often associated with Cramer's rule or Laplace expansion). This method is conceptually straightforward, directly reflecting the definition of a determinant, but it is computationally very expensive for large matrices.

*   **Numerical Method:** The `Deter` function implements the Laplace expansion. For an $N \times N$ matrix, it expands the determinant along the first row: $\det(A) = \sum_{j=1}^{N} (-1)^{1+j} A_{1j} \cdot \det(M_{1j})$, where $M_{1j}$ is the submatrix obtained by removing the first row and $j$-th column of $A$. The function recursively calls itself for the determinants of these $(N-1) \times (N-1)$ submatrices. The base case for the recursion is a $2 \times 2$ matrix, for which the determinant is calculated directly as $A_{11}A_{22} - A_{12}A_{21}$.
*   **Key Procedures:**
    *   `Deter(n: integer; A: pMat): real;`
        *   **Purpose**: The recursive function to calculate the determinant.
        *   Returns: The determinant value.
    *   `pow(n: integer): integer;`
        *   **Purpose**: Returns $-(-1)^n$ (or a similar sign factor, typically $(-1)^{row+col}$ for cofactor expansion) used for the alternating signs in the sum.
    *   `Submatrix(n, col: integer; A: pMat; VAR B: pMat);`
        *   **Purpose**: Extracts a submatrix `B` from `A` by removing the first row and the specified `col`-th column. Memory for `B` is dynamically allocated and deallocated at each recursive step to manage memory.
*   **Real-world Application Ideas:**
    *   **Educational Tool:** This method is excellent for illustrating the mathematical definition of a determinant and the concept of recursion in programming.
    *   **Small Matrix Problems:** For very small matrices (e.g., $N \le 4$), its computational cost is acceptable, and it provides a direct conceptual link to the definition. For larger matrices, its $O(N!)$ complexity makes it impractical.

#### `FindDet` Function (`tfinddet.pas`)

This program calculates the determinant of a real square matrix using a function `FindDet` that converts the matrix to upper triangular form through a series of row operations.

*   **Numerical Method:** The `FindDet` function applies elementary row operations (swapping rows, adding a multiple of one row to another) to transform the input matrix into an upper triangular form.
    *   **Property Used:** These elementary row operations do not change the determinant (or change its sign in the case of row swaps). The sign factor is tracked (variable `l`).
    *   Once the matrix is in upper triangular form, its determinant is simply the product of its diagonal elements.
    *   The function includes a basic pivot handling mechanism: if a diagonal element is zero, it attempts to swap rows to find a non-zero pivot. If no such pivot can be found, the determinant is zero.
*   **Key Function:**
    *   `FindDet(matrix: MAT; n: integer): Double;`
        *   **Purpose**: Calculates the determinant of `matrix` by triangularization.
        *   Returns: The determinant value.
*   **Real-world Application Ideas:**
    *   **Direct Invertibility Check:** Provides a straightforward way to check if a matrix is invertible (determinant is non-zero).
    *   **Educational Demonstration:** Clearly illustrates the relationship between row operations, triangular matrices, and the determinant.

---

## 5. Eigenvalue Problems

This section details various numerical methods for computing eigenvalues and eigenvectors of matrices. Eigenvalues ($\lambda$) and their corresponding eigenvectors ($x$) describe intrinsic properties of a linear transformation represented by a matrix $A$, such that $Ax = \lambda x$. They are fundamental in many scientific, engineering, and data analysis disciplines.

### 5.1. Characteristic Polynomial Calculation

The characteristic polynomial of a square matrix $A$ is defined as $P(\lambda) = \det(A - \lambda I)$, where $I$ is the identity matrix and $\lambda$ is a scalar variable. The roots of this polynomial are the eigenvalues of $A$. Calculating its coefficients can be an intermediate step to finding eigenvalues, particularly for smaller matrices or structured matrices.

#### Real Tridiagonal Matrix (`carpol.pas`)

This program calculates the coefficients of the characteristic polynomial $P(\lambda) = \det(A - \lambda I)$ for a real tridiagonal matrix $A(i,j)$.

*   **Numerical Method:** The `PCTR` procedure employs a specialized recursive relation that efficiently computes the coefficients of the characteristic polynomial specifically for tridiagonal matrices. This approach leverages the sparse structure (only three diagonals are non-zero) to significantly reduce computational complexity compared to general matrix methods. The recursion builds the polynomial by considering sub-matrices.
*   **Key Procedure:**
    *   `PCTR(n: integer; B, D: pV; VAR P: pV);`
        *   **Purpose**: Calculates the coefficients of the characteristic polynomial for a real tridiagonal matrix.
        *   **n**: The size (dimension) of the matrix.
        *   **B**: A vector storing the products of the co-diagonal terms. In the main program, this is computed as $B[k] := C[k] \cdot E[k]$, where $C[k]$ are the sub-diagonal elements and $E[k]$ are the super-diagonal elements.
        *   **D**: A vector storing the main diagonal terms of the tridiagonal matrix.
        *   **P** (Output): A vector storing the coefficients of the characteristic polynomial. `P[1]` corresponds to the coefficient of $\lambda^n$, `P[2]` to $\lambda^{n-1}$, and so on, down to `P[n+1]` for the constant term ($\lambda^0$).
*   **Real-world Application Ideas:**
    *   **Vibrational Analysis:** Used to find the natural frequencies of systems that can be modeled as chains of masses and springs, which often lead to tridiagonal matrices in their governing equations.
    *   **Quantum Mechanics (1D Systems):** Solving for energy levels in one-dimensional quantum mechanical systems with discrete potentials, where the Hamiltonian matrix can be tridiagonal.

#### Complex Square Matrix (`carpol1.pas`)

This program calculates the complex coefficients of the characteristic polynomial $P(\lambda) = \det(A - \lambda I)$ for a general square **complex** matrix $A(i,j)$.

*   **Numerical Method:** The `PCCS` procedure (Characteristic Polynomial of a Complex Square matrix) implements an algorithm that calculates the coefficients using sums involving traces of powers of the matrix, generalized for complex numbers. This method relies heavily on robust complex arithmetic operations (addition, subtraction, multiplication).
*   **Key Procedure:**
    *   `PCCS(n: integer; A: pMC; VAR P: pVC);`
        *   **Purpose**: Calculates the complex coefficients of the characteristic polynomial for a complex square matrix.
        *   **n**: The size (dimension) of the complex matrix.
        *   **A**: The input complex square matrix. It is typically a pointer to a 2D array of `Complex` records (`pMC`).
        *   **P** (Output): A vector of `Complex` numbers representing the coefficients of the characteristic polynomial.
*   **Associated Files:** `carpol1.pas` (main program), `ucomplex1.pas` (provides basic complex number arithmetic).
*   **Real-world Application Ideas:**
    *   **Quantum Computing:** Characterizing quantum gates or complex Hamiltonians, which are represented by complex matrices.
    *   **Electrical Engineering (AC Circuits):** Analyzing complex-valued system matrices in frequency domain analysis.
    *   **Control Systems:** For systems with complex-valued parameters or inputs, where stability analysis requires complex eigenvalues.

#### Real Square Matrix (`carpol2.pas`)

This program calculates the coefficients of the characteristic polynomial $P(\lambda) = \det(A - \lambda I)$ for a general **real** square matrix $A(i,j)$.

*   **Numerical Method:** The `PCMS` procedure (Characteristic Polynomial of a Real Square matrix) computes the coefficients using an algorithm based on **Newton's sums** (also known as Faddeev-Leverrier algorithm). This method iteratively calculates powers of the matrix and their traces to determine the polynomial coefficients. It's a direct approach applicable to any real square matrix.
*   **Key Procedure:**
    *   `PCMS(n: integer; A: pM; VAR P: pV);`
        *   **Purpose**: Calculates the coefficients of the characteristic polynomial for a real square matrix.
        *   **n**: The size (dimension) of the real matrix.
        *   **A**: The input real square matrix (typically a pointer to a 2D array of `Double`s, `pM`).
        *   **P** (Output): A vector of `Double`s representing the coefficients of the characteristic polynomial.
*   **Associated Files:** `carpol2.pas` (main program).
*   **Real-world Application Ideas:**
    *   **System Stability Analysis:** In linear control theory and dynamical systems, the roots of the characteristic polynomial (eigenvalues) determine the stability of the system.
    *   **Structural Mechanics:** Analyzing the dynamic behavior of structures where the system matrix might be non-symmetric.
    *   **Economic Models:** Analyzing the stability of equilibrium points in linear economic models.

#### Real Symmetric Square Matrix (Lanczos's Method) (`carpol3.pas`)

This program calculates the coefficients of the characteristic polynomial $P(\lambda) = \det(A - \lambda I)$ for a **real symmetric square matrix** $A(i,j)$ using **Lanczos's method**. This two-stage approach is particularly efficient for large symmetric matrices.

*   **Numerical Method:** Lanczos's method is an iterative algorithm that, for a symmetric matrix, constructs an orthonormal basis with respect to which the matrix is tridiagonal. This process is called **tridiagonalization**. Once the matrix is reduced to an equivalent tridiagonal form, the coefficients of its characteristic polynomial can be very efficiently calculated using the specialized `PCTR` procedure (reused from `carpol.pas`).
*   **Key Procedures:**
    *   `TDSL(eps: double; n: integer; A: pM; VAR it: integer; VAR D, E: pV; VAR Tau: pM);`
        *   **Purpose**: Tridiagonalizes a symmetric real matrix `A` using Lanczos's method.
        *   **eps**: Precision used for checks.
        *   **n**: Size of the matrix.
        *   **A**: Input symmetric matrix.
        *   **it** (Output): Flag indicating success (`1`) or failure (`0`, e.g., if a pivot becomes too small).
        *   **D** (Output): Vector storing the main diagonal elements of the resulting tridiagonal matrix.
        *   **E** (Output): Vector storing the super-diagonal elements of the resulting tridiagonal matrix.
        *   **Tau** (Working Zone): An auxiliary matrix used by the method.
    *   `PCTR(n: integer; B, D: pV; VAR P: pV);`
        *   **Purpose**: (Reused from `carpol.pas`) Calculates the coefficients of the characteristic polynomial for a tridiagonal matrix (using the `D` and `E` vectors from `TDSL`).
*   **Real-world Application Ideas:**
    *   **Large-scale Eigenvalue Problems:** For very large symmetric matrices arising in fields like quantum chemistry, materials science, or graph analysis, Lanczos's method is highly efficient for computing the characteristic polynomial or finding a subset of eigenvalues.
    *   **Model Reduction:** Can be used to create smaller, tridiagonal representations of large systems while preserving key dynamic properties.

### 5.2. Power Methods for Eigenvalues

The Power Method is an iterative algorithm primarily used for finding the dominant eigenvalue (the eigenvalue with the largest absolute value) of a matrix and its corresponding eigenvector. The Inverse Power Method is an extension used to find the eigenvalue closest to a given value, or specifically the smallest eigenvalue.

#### Greatest Eigenvalue by Power Method (`tpwm.pas`)

This program calculates the greatest eigenvalue (in absolute value, i.e., the dominant eigenvalue) of a real square matrix and its associated eigenvector using the **Power Method**.

*   **Numerical Method:** The Power Method is an iterative algorithm. It starts with an arbitrary non-zero initial vector (e.g., all ones or random values). In each iteration, this vector is multiplied by the matrix $A$. As iterations progress, the resulting vector converges towards the eigenvector corresponding to the dominant eigenvalue. The ratio of successive components (or norms) of the iterated vector converges to the dominant eigenvalue.
    *   **Limitations:** This method works best when there is a single dominant eigenvalue (its absolute value is strictly greater than all other eigenvalues). It may struggle or converge slowly if there are multiple eigenvalues with the same largest absolute value or if the dominant eigenvalue is complex.
*   **Key Procedure:**
    *   `PWM(eps, dta: Double; m, n: Integer; VAR A: MAT; VAR it: Integer; VAR gamma: double; VAR X1: VEC);`
        *   **Purpose**: Calculates the greatest eigenvalue (`gamma`) and its corresponding eigenvector (`X1`) using the Power Method.
        *   **eps**: A small positive number (epsilon) used for numerical comparisons, typically to check for near-zero values or to avoid division by zero.
        *   **dta**: The desired precision for convergence (e.g., the difference between successive approximations of the eigenvector or eigenvalue).
        *   **m**: The maximum number of iterations allowed.
        *   **n**: The size (dimension) of the square matrix `A`.
        *   **A** (Input/Output): The input real square matrix. Its content might be modified during internal computations, but its structure remains.
        *   **it** (Output): An error indicator: `1` for successful convergence, `-1` if convergence is not achieved within `m` iterations, `0` if the method cannot be applied (e.g., dominant eigenvalue is effectively zero).
        *   **gamma** (Output): The computed greatest eigenvalue (dominant eigenvalue).
        *   **X1** (Output): The computed eigenvector associated with `gamma`.
*   **Real-world Application Ideas:**
    *   **PageRank Algorithm (Google Search):** A simplified version of PageRank fundamentally relies on finding the dominant eigenvector of a matrix representing hyperlinks between web pages.
    *   **Spectral Clustering:** Identifying clusters in data by analyzing the dominant eigenvectors of similarity matrices.
    *   **Markov Chains:** Determining the long-term behavior or stationary distribution of a Markov chain.
    *   **Vibrational Analysis:** Finding the fundamental natural frequency (lowest frequency, largest period) and corresponding mode shape of a vibrating system (though often inverse power method is used for lowest).

#### Smallest Eigenvalue by Gauss and Power Methods (`tpwimgt.pas`)

This program calculates the smallest eigenvalue (in absolute value) of a real square matrix and its associated eigenvector. It achieves this by combining matrix inversion using Gauss-Jordan elimination with the Power Method.

*   **Numerical Method:** This method leverages a key property: if $\lambda$ is an eigenvalue of $A$, then $1/\lambda$ is an eigenvalue of $A^{-1}$. Therefore, the smallest eigenvalue (in absolute value) of $A$ corresponds to the largest eigenvalue (in absolute value, i.e., the dominant eigenvalue) of its inverse $A^{-1}$. The process proceeds in two main stages:
    1.  **Matrix Inversion (`MATINV`):** The input matrix $A$ is first inverted using the Gauss-Jordan elimination method. This step explicitly computes $A^{-1}$.
    2.  **Power Method on Inverse (`PWM`):** The standard Power Method (`PWM`) is then applied to the computed inverse matrix ($A^{-1}$) to find its dominant eigenvalue. This dominant eigenvalue of $A^{-1}$ is the reciprocal of the smallest eigenvalue of the original matrix $A$. The eigenvector found will be the same for both $A$ and $A^{-1}$.
*   **Key Procedures:**
    *   `MATINV(N, M: integer; VAR AA: MAT; VAR BB: MAT; VAR DET: DOUBLE);`
        *   **Purpose**: (Reused from `sysmat.pas`) Inverts the matrix `AA` and calculates its determinant. For this application, `M` is set to `0` (no right-hand side `BB` is actively used for solution).
    *   `PWM(eps, dta: Double; m, n: Integer; VAR A: MAT; VAR it: Integer; VAR gamma: double; VAR X1: VEC);`
        *   **Purpose**: (Reused from `tpwm.pas`) Calculates the greatest eigenvalue and eigenvector of the *inverted* matrix $A^{-1}$.
    *   `PWMIMGT(eps, dta: Double; m, n: Integer; VAR A: MAT; VAR it: Integer; VAR lambda: double; VAR X: VEC);`
        *   **Purpose**: Orchestrates the entire process of inversion followed by the Power Method to find the smallest eigenvalue.
        *   **lambda** (Output): The computed smallest eigenvalue (in absolute value).
*   **Limitations:** This method requires the matrix $A$ to be **invertible** (non-singular, i.e., its determinant is non-zero). If $A$ is singular or ill-conditioned, the inversion step may fail or yield inaccurate results. It also shares the limitations of the Power Method, struggling if there are multiple eigenvalues with the same smallest absolute value.
*   **Real-world Application Ideas:**
    *   **Stability Analysis:** Identifying the least stable mode in a system, which often corresponds to the eigenvalue closest to zero.
    *   **Condition Number Estimation:** The smallest eigenvalue (in absolute value) is crucial for computing the condition number of a matrix, which quantifies its sensitivity to input perturbations.
    *   **Regularization Parameter Selection:** In some ill-posed problems, the smallest eigenvalues can inform the choice of regularization parameters.

### 5.3. Jacobi's Method

Jacobi's method is a classic iterative algorithm for computing the eigenvalues and eigenvectors of a **real symmetric matrix** (and its extension, Hermitian matrices). It performs a sequence of orthogonal similarity transformations (Jacobi rotations) to iteratively make the matrix diagonal.

#### Real Symmetric Matrices (`ujacobi.pas`)

This unit provides the implementation of Jacobi's method for finding eigenvalues and eigenvectors of a real symmetric matrix. The cyclic Jacobi method iteratively applies a series of plane rotations to annihilate off-diagonal elements.

*   **Numerical Method:** The `Jacobi` procedure iteratively applies **Jacobi rotations** (2x2 orthogonal similarity transformations) to the symmetric matrix. Each rotation is chosen to set a specific off-diagonal element to zero. While one rotation may make previously zeroed elements non-zero again, the sum of squares of off-diagonal elements decreases monotonically, ensuring convergence to a diagonal matrix. The diagonal elements of the final matrix are the eigenvalues. The eigenvectors are formed by accumulating the product of all applied rotation matrices. The implementation normalizes the eigenvectors such that their largest element (in absolute value) is one.
*   **Key Procedure (from `ujacobi.pas`):**
    *   `Jacobi(Dimen: integer; MaxIter: integer; Tolerance: Float; VAR Iter: integer; VAR Error: byte);`
        *   **Purpose**: Computes eigenvalues and eigenvectors of a symmetric matrix.
        *   **Dimen**: The dimension (size) of the square symmetric matrix.
        *   **MaxIter**: The maximum number of iterations allowed for convergence.
        *   **Tolerance**: A floating-point value specifying the convergence criterion. The iteration stops when the sum of squares of all off-diagonal elements (or some other measure) falls below this tolerance.
        *   **Iter** (Output): The actual number of iterations performed.
        *   **Error** (Output): An error code: `0` for no error; `1` if `Dimen < 1`; `2` if `Tolerance` is too small; `3` if `MaxIter < 1`; `4` if the input matrix `MAT` is not symmetric; `5` if convergence is not achieved within `MaxIter`.
        *   **Input/Output via Global Pointers:** The procedure operates on globally declared pointers:
            *   `Mat`: Input symmetric matrix (must be allocated and filled by the calling program).
            *   `Eigenvalues`: Output vector storing the computed eigenvalues (unsorted).
            *   `Eigenvectors`: Output matrix storing the normalized eigenvectors (typically in rows, corresponding to eigenvalues).
        *   **Note**: The calling program is responsible for allocating (`New()`) and deallocating (`Dispose()`) these global pointer variables.
*   **Associated Files:** `ujacobi.pas` (main unit), `tujacobi.pas` (test program demonstrating usage).
*   **Real-world Application Ideas:**
    *   **Principal Component Analysis (PCA):** A cornerstone of dimensionality reduction in data science. PCA involves finding the eigenvalues and eigenvectors of the data's covariance matrix (which is symmetric).
    *   **Structural Vibrations:** Analyzing the natural frequencies and mode shapes of vibrating mechanical structures.
    *   **Quantum Chemistry:** Solving for molecular orbitals and energy levels in simple quantum systems.
    *   **Tensor Analysis:** Decomposing symmetric tensors in physics and engineering.

#### Hermitian Matrices (`tephj.pas`)

This program calculates the eigenvalues and eigenvectors of a square **Hermitian matrix** using Jacobi's Method. A Hermitian matrix is a complex square matrix that is equal to its own conjugate transpose ($A = A^H$). Its eigenvalues are always real.

*   **Numerical Method:** The `EPHJ` procedure adapts Jacobi's method for complex Hermitian matrices. Similar to the real symmetric case, it iteratively applies a sequence of **unitary similarity transformations** (complex rotations). These rotations are chosen to annihilate specific off-diagonal complex elements. The process drives the matrix towards a diagonal form, where the diagonal elements are the real eigenvalues. The product of all applied unitary transformations forms the matrix of eigenvectors.
    *   The implementation relies on basic complex arithmetic operations (addition, subtraction, multiplication, absolute value, conjugation), which are either defined internally or by a supporting unit.
*   **Key Procedure:**
    *   `EPHJ(dta: Double; M, N: Integer; A: pCMat; Var it: Integer; R: pVec; VX: pCMat);`
        *   **Purpose**: Computes the eigenvalues (`R`) and eigenvectors (`VX`) of a Hermitian matrix `A`.
        *   **dta**: The desired precision (tolerance) for convergence.
        *   **M**: The maximum number of iterations allowed.
        *   **N**: The size (dimension) of the Hermitian matrix.
        *   **A**: The input Hermitian matrix (a pointer to a 2D array of `Complex` records).
        *   **it** (Output): Flag indicating convergence: `1` for successful convergence, `-1` if convergence is not achieved within `M` iterations.
        *   **R** (Output): A vector storing the computed eigenvalues. For Hermitian matrices, these values will be real.
        *   **VX** (Output): A matrix storing the computed complex eigenvectors in its columns.
    *   **Supporting Procedure:**
        *   `NORMAL(N: Integer; R: pVec; Z: pCMat);` (internal to `tephj.pas`): Sorts the eigenvalues by absolute value in ascending order and normalizes the corresponding eigenvectors.
*   **Associated Files:** `tephj.pas` (main program), implicitly uses complex arithmetic (e.g., `ucomplex1.pas`).
*   **Real-world Application Ideas:**
    *   **Quantum Mechanics:** Hermitian operators represent observable quantities (e.g., energy, momentum, position) in quantum mechanics. Their eigenvalues correspond to the possible outcomes of measurements, and eigenvectors represent the quantum states.
    *   **Quantum Computing:** Characterizing quantum states and operations, which are often described by Hermitian or unitary matrices.
    *   **Signal Processing:** Analyzing complex-valued signals in domains like radar, sonar, and communications, where covariance matrices can be Hermitian.

### 5.4. QL Method for Symmetric Matrices

The QL algorithm is a highly efficient and numerically stable iterative method for finding all eigenvalues and eigenvectors of **symmetric matrices**. For tridiagonal matrices, it is particularly fast. For dense symmetric matrices, it is typically used as a second step after a preliminary reduction to tridiagonal form (e.g., by Householder transformations).

#### Real Symmetric Tridiagonal Matrices (`elprotd.pas`, `ttql2.pas`)

These programs compute the eigenvalues and eigenvectors of a **real symmetric tridiagonal matrix** using the QL algorithm with implicit shifts.

*   **Numerical Method:** The QL algorithm (a variant of the QR algorithm) is an iterative method that repeatedly transforms a symmetric tridiagonal matrix into a more diagonal form. It achieves this by applying a sequence of orthogonal QL decompositions ($A_k = Q_k L_k$, then $A_{k+1} = L_k Q_k$). Implicit shifts are used to accelerate convergence, especially for clustered eigenvalues. The process converges to a diagonal matrix whose elements are the eigenvalues, while the accumulated product of the orthogonal transformations forms the matrix of eigenvectors.
*   **Key Procedures:**
    *   `TQLI(Var D, E: pV; n: integer; VAR Z: pM);` (from `elprotd.pas`):
        *   **Purpose**: Implements the QL algorithm with implicit shifts for symmetric tridiagonal matrices.
        *   **D** (Input/Output): A vector containing the elements of the main diagonal of the tridiagonal matrix. On output, it stores the computed eigenvalues.
        *   **E** (Input/Output): A vector containing the elements of the subdiagonal (from $E[2]$ to $E[n]$). On output, its elements become very close to zero, indicating the diagonalization.
        *   **n**: The size (dimension) of the matrix.
        *   **Z** (Input/Output): A matrix. On input, it should be initialized as an identity matrix (if only `TQLI` is used for a raw tridiagonal matrix) or with accumulated transformations (if used after a reduction like Householder). On output, it stores the eigenvectors in its columns.
    *   `TQL2(NM, N: Integer; Var D, E: pVec; Var Z: pMat; Var IER: Integer);` (from `ttql2.pas`): A wrapper or slightly different implementation of the QL algorithm, functionally similar to `TQLI`.
        *   `NM`: First dimension of matrix `Z` (often `N` for square matrices).
        *   `IER`: Error code (0=success, L=no convergence for L-th eigenvalue).
*   **Associated Files:** `elprotd.pas` (main program), `ttql2.pas` (test program).
*   **Real-world Application Ideas:**
    *   **Vibrational Analysis:** Highly efficient for finding the natural frequencies and mode shapes of coupled mechanical oscillators that can be modeled as tridiagonal systems.
    *   **Quantum Mechanics:** Solving the Schrödinger equation for systems that can be modeled with tridiagonal Hamiltonian matrices.
    *   **Numerical Analysis of Sturm-Liouville Problems:** Discretization of such problems often leads to symmetric tridiagonal eigenvalue problems.

#### Real Symmetric Matrices (Householder Reduction + QL) (`ttred2.pas`)

This program finds the eigenvalues and eigenvectors of a general **real symmetric matrix** by employing a powerful two-stage approach: first reducing the matrix to tridiagonal form using Householder transformations, and then applying the QL method to the resulting tridiagonal matrix. This is a standard and robust method for dense symmetric matrices.

*   **Numerical Method:**
    1.  **Householder Reduction (`TRED2`):** The input symmetric matrix is transformed into an equivalent symmetric tridiagonal matrix. This is a finite (non-iterative) process that uses a sequence of orthogonal Householder reflections. This step preserves the eigenvalues of the original matrix and implicitly computes the transformation matrix needed to convert the eigenvectors of the tridiagonal form back to the original coordinate system.
    2.  **QL Algorithm (`TQL2`):** The QL algorithm (reused from `elprotd.pas` or `ttql2.pas`) is then applied to the resulting symmetric tridiagonal matrix to efficiently compute its eigenvalues and eigenvectors. The eigenvectors of the original full symmetric matrix are obtained by multiplying the eigenvectors of the tridiagonal matrix by the accumulated Householder transformation matrix.
*   **Key Procedures:**
    *   `TRED2(Var Z: pMAT; N: Integer; Var D, E: pVEC);`
        *   **Purpose**: Reduces a symmetric matrix `Z` to symmetric tridiagonal form using Householder transformations.
        *   **Z** (Input/Output): On input, the original symmetric matrix. On output, it is replaced by the accumulated Householder transformation matrix (which will be used by `TQL2`).
        *   **N**: The size (dimension) of the matrix.
        *   **D** (Output): A vector storing the main diagonal elements of the resulting tridiagonal matrix.
        *   **E** (Output): A vector storing the sub-diagonal elements of the resulting tridiagonal matrix.
    *   `TQL2(N: Integer; Var D, E: pVec; Var Z: pMat; Var IER: Integer);`
        *   **Purpose**: (Reused from `ttql2.pas`) Computes the eigenvalues and eigenvectors of the tridiagonal matrix. It takes the `D` and `E` vectors from `TRED2` and the accumulated transformation matrix `Z` to produce the final eigenvectors.
*   **Associated Files:** `ttred2.pas` (main program).
*   **Real-world Application Ideas:**
    *   **Principal Component Analysis (PCA):** A primary method for computing eigenvalues and eigenvectors of covariance matrices in high-dimensional data analysis.
    *   **Quantum Chemistry:** Calculating electronic structures of atoms and molecules where Hamiltonian matrices are symmetric.
    *   **Finite Element Analysis:** Solving eigenvalue problems in structural dynamics, heat transfer, or quantum mechanics where the discretized system matrices are symmetric.
    *   **Graph Theory:** Analyzing the eigenvalues of adjacency matrices of graphs to understand network properties.

### 5.5. Rutishauser's Method and Inverse Iteration (`elpro.pas`)

This program calculates the eigenvalues and eigenvectors of a **general real square matrix** by combining two distinct iterative methods: Rutishauser's method for eigenvalues and the Inverse Iteration method for computing their associated eigenvectors. This approach is historically significant, providing a method for non-symmetric matrices.

*   **Numerical Method:**
    1.  **Rutishauser's Method (`VAMR`):** This is an iterative algorithm for finding eigenvalues, fundamentally based on repeated similarity transformations using LU (or QR) factorization. Specifically, it involves repeatedly decomposing the matrix $A_k = L_k U_k$ and then forming a new matrix $A_{k+1} = U_k L_k$. This process (similar to the QR algorithm) drives the matrix towards an upper triangular form whose diagonal elements are the eigenvalues.
    2.  **Inverse Iteration Method (`IIM`):** Once an approximate eigenvalue $\lambda$ is known (e.g., from `VAMR`), the Inverse Iteration method is used to find its corresponding eigenvector. It works by solving the linear system $(A - \lambda I)x_{k+1} = x_k$ for a sequence of vectors $x_k$. This system is typically solved efficiently using LU decomposition, as $(A - \lambda I)$ needs to be factored only once per approximate eigenvalue. This method is highly effective for finding a specific eigenvector when a good approximation of the eigenvalue is available.
*   **Key Procedures:**
    *   `DECCRM(eps: Double; n: integer; A: MAT; VAR it: integer; VAR U, V: MAT);`
        *   **Purpose**: Performs Crout's decomposition (a variant of LU decomposition) of matrix `A`. This is a subroutine used by `VAMR`.
    *   `VAMR(eps, dta: Double; m, n: integer; A: MAT; VAR it: integer; VAR R: VEC);`
        *   **Purpose**: Calculates the eigenvalues (`R`) of a real square matrix using Rutishauser's method.
        *   **it**: Flag for convergence.
    *   `IIM(eps, dta: Double; m, n: integer; A: MAT; VAR it: integer; VAR gamma: Double; VAR X1: VEC);`
        *   **Purpose**: Calculates a real eigenvalue (`gamma`) and its associated eigenvector (`X1`) using the inverse iteration method.
    *   `EPMRI(eps, d1, d2: Double; m, n: integer; A: MAT; VAR it: integer; VAR R: VEC; VAR VX: MAT);`
        *   **Purpose**: Orchestrates the overall process: calls `VAMR` to find all eigenvalues, then iterates through them, calling `IIM` for each to find its corresponding eigenvector.
*   **Real-world Application Ideas:**
    *   **Iterative Refinement of Eigenvectors:** Useful in scenarios where an initial approximation of eigenvalues is available (e.g., from a coarser simulation), and highly accurate eigenvectors are needed.
    *   **Sensitivity Analysis of Eigenvalues:** Exploring how eigenvalues change with small perturbations in the matrix, often benefiting from the targeted nature of inverse iteration.

### 5.6. QR Algorithm for Non-Symmetric Matrices (`feigen0.pas`, `linpack.pas`, `test_hqr.pas`, `thqr.pas`)

This set of programs and units provides a comprehensive implementation of the **QR algorithm** for finding all eigenvalues and, optionally, eigenvectors of **general non-symmetric real matrices**. The QR algorithm is one of the most widely used, numerically stable, and robust algorithms for this problem.

*   **Numerical Method:** The QR algorithm is an iterative method that repeatedly applies QR decomposition to a matrix ($A_k = Q_k R_k$, then $A_{k+1} = R_k Q_k$). For non-symmetric matrices, this process generally converges to an upper triangular matrix with 1x1 blocks on the diagonal (for real eigenvalues) or 2x2 blocks on the diagonal (for complex conjugate eigenvalue pairs). The eigenvalues are then read directly from these diagonal blocks. Implicit shifts are used to accelerate convergence, especially for clustered eigenvalues. If eigenvectors are required, the orthogonal transformations ($Q_k$) are accumulated.
    The overall process typically involves three main stages:
    1.  **Balancing (`balance` in `feigen0.pas`, `Balanc` in `linpack.pas`):** The matrix is scaled (using diagonal similarity transformations) to make the norms of its rows and columns roughly equal. This improves the conditioning of the eigenvalue problem and leads to better accuracy.
    2.  **Reduction to Upper Hessenberg Form (`elmhes` in `feigen0.pas`, `ElmHes` in `linpack.pas`):** The balanced matrix is transformed into an upper Hessenberg form (upper triangular with one additional non-zero subdiagonal). This is done using Householder or elementary transformations. This reduction is a finite process that significantly reduces the computational cost of the subsequent iterative QR steps without changing the eigenvalues.
    3.  **QR Iterations (`hqr2` in `feigen0.pas`, `HQR_MR` in `linpack.pas`):** The core iterative part. A sequence of implicit QR steps is applied to the Hessenberg matrix. This iteratively transforms the matrix towards a (block) triangular form from which eigenvalues can be extracted. If eigenvectors are also desired, the accumulating orthogonal transformations (`eivec`) are updated at each step.
*   **Key Procedures/Functions (from `feigen0.pas` - full implementation with eigenvectors):**
    *   `eigen(vec, ev_norm: Boolean; n: integer; mat: Square_Matrix; VAR eivec: Square_Matrix; VAR valre, valim: Real_Vector; VAR cnt: Integer_Vector; VAR rc: Integer)`: The top-level procedure to compute all eigenvalues and, optionally, eigenvectors.
        *   **vec**: `True` to compute eigenvectors, `False` for eigenvalues only.
        *   **ev_norm**: `True` to normalize eigenvectors after computation.
        *   **mat**: Input matrix (modified in-place).
        *   **eivec**: Output matrix storing eigenvectors (columns for real, pairs of columns for complex).
        *   **valre, valim**: Output vectors for real and imaginary parts of eigenvalues.
        *   **cnt**: Iteration count for each eigenvalue.
        *   **rc**: Return code (0=OK, various error codes for sub-procedures).
    *   `balance(n, mat, scal, low, high)`: Balances the matrix.
    *   `elmhes(n, low, high, mat, perm)`: Reduces matrix to upper Hessenberg form.
    *   `elmtrans(n, low, high, mat, perm, h)`: Initializes the eigenvector matrix for Hessenberg reduction.
    *   `hqr2(vec, n, low, high, h, wr, wi, eivec, cnt, rc)`: Implements the core Francis QR algorithm (iterative part).
    *   `hqrvec(n, low, high, h, wr, wi, eivec, rc)`: Computes eigenvectors from the Hessenberg form after eigenvalues are found.
    *   `balback(n, low, high, scal, eivec)`: Reverses the balancing transformations on the eigenvectors.
    *   `norm_1(n, v, wi, rc)`: Normalizes the computed eigenvectors.
    *   `Comdiv()`, `Comabs()`: Internal helper functions for complex arithmetic, crucial for handling complex eigenvalues and their eigenvectors.
*   **`linpack.pas` (LINPACK Routines):** Provides a subset of core LINPACK library routines, including `Balanc`, `ElmHes`, and `HQR_MR` (a version of the QR algorithm focused on eigenvalues without computing eigenvectors). These can be used independently or as building blocks.
*   **Associated Files:** `feigen0.pas` (main unit for full QR with eigenvectors), `linpack.pas` (LINPACK QR routines), `test_hqr.pas` and `thqr.pas` (test programs demonstrating their usage).
*   **Real-world Application Ideas:**
    *   **Stability of Dynamical Systems:** Analyzing the eigenvalues of Jacobian matrices derived from non-linear systems to determine the stability of equilibrium points.
    *   **Control Systems Design:** Analyzing eigenvalues of state-space matrices to understand system modes, stability, and transient response.
    *   **Fluid Dynamics:** Linear stability analysis of fluid flows.
    *   **Structural Mechanics (Non-Conservative Systems):** Analyzing damped or non-symmetric structures where eigenvalues can be complex.

---
