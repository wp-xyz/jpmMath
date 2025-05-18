{------------------------------------------------------------------------*
*                                                                        *
*     Turbo Pascal Numerical Methods Toolbox                             *
*     Copyright (c) 1986, 87 by Borland International, Inc.              *
*                                                                        *
*  This unit provides the Jacobi PROCEDURE for finding eigenvalues and   *
*  eigenvectors of a real symmetric matrix.  Version with pointers for   *
*  use with Turbo Pascal for Windows by J-P Moreau (Paris).              *
*                        (www.jpmoreau.fr)                               *
*                                                                        *
*------------------------------------------------------------------------}
UNIT Ujacobi;

INTERFACE

TYPE
  Float = Double;

CONST
  TNNearlyZero = 1.2E-16;
  TNArraySize  = 50;          { Maximum size of matrix }

TYPE
  PTNV = ^TNvector;
  PTNM = ^TNmatrix;
  PTNI = ^TNIntVector;
  TNvector    = ARRAY[1..TNArraySize] OF Float;
  TNmatrix    = ARRAY[1..TNArraySize] OF TNvector;
  TNIntVector = ARRAY[1..TNArraySize] OF integer;

VAR
  Mat : PTNM;            { given square symmetric matrix }
  Eigenvalues  : PTNV;   { eigenvalues solution ARRAY    }
  Eigenvectors : PTNM;   { eigenvectors solution ARRAY   

  Nota 1 - These pointers must be initialized and disposed of after use
	   by the calling PROGRAM.

  Nota 2 - Also the optional sorting of the eigenvalues must be provided
	   by the calling PROGRAM. The normalization of the eigenvectors
	   to unity is here provided.                    }                      
                            

PROCEDURE Jacobi(Dimen        : integer;
                 MaxIter      : integer;
                 Tolerance    : Float;
	     VAR Iter         : integer;
             VAR Error        : byte);

{----------------------------------------------------------------------------}
{-                                                                          -}
{-     Input: Dimen, Mat, MaxIter, tolerance                                -}
{-    Output: Eigenvalues, Eigenvectors, Iter, Error                        -}
{-                                                                          -}
{-   Purpose: The eigensystem of a symmetric matrix can be                  -}
{-            computed much more simply than the                            -}
{-            eigensystem of an arbitrary matrix.  The                      -}
{-            cyclic Jacobi method is an iterative                          -}
{-            technique for approximating the complete                      -}
{-            eigensystem of a symmetric matrix to a given                  -}
{-            tolerance. The method consists of multiplying                 -}
{-            the matrix, Mat, by a series of rotation                      -}
{-            matrices, rot[i].  The rotation matrices are                  -}
{-            chosen so that the elements of the upper                      -}
{-            triangular part of Mat are systematically                     -}
{-            annihilated.  That is, rot[1] is chosen so                    -}
{-            that Mat[1, 1] is identically zero; rot[2] is                 -}
{-            chosen so that Mat[1, 2] is identically zero;                 -}
{-            etc.  Since each operation will probably                      -}
{-            change the value of elements annihilated IN                   -}
{-            previous operations, the method is iterative.                 -}
{-            Eventually, the matrix will be diagonal. The                  -}
{-            eigenvalues will be the elements along the                    -}
{-            diagonal of the matrix and the eigenvectors                   -}
{-            will be the rows of the matrix created by the                 -}
{-            product of all the rotation matrices rot[i].                  -}
{-                                                                          -}
{-            The user inputs the matrix, tolerance and maximum             -}
{-            number of iterations. The PROCEDURE returns the               -}
{-            eigenvalues and eigenvectors (or error code) of the           -}
{-            matrix.                                                       -}
{-                                                                          -}
{-   User-Defined Types: TNvector = ARRAY[1..TNArraySize] of real;          -}
{-                       TNmatrix = ARRAY[1..TNArraySize] of TNvector;      -}
{-                                                                          -}
{-   Global Variables:  Dimen        : integer   Dimension of square matrix -}
{-                      Mat          : TNmatrix  Square matrix              -}
{-                      MaxIter      : integer   Maximum number of          -}
{-                                               Iterations                 -}
{-                      tolerance    : real      tolerance IN answer        -}
{-                      Eigenvalues  : TNvector  Eigenvalues of Mat         -}
{-                      Eigenvectors : TNmatrix  Eigenvectors of Mat        -}
{-                      Iter         : integer   Number of iterations       -}
{-                      Error        : byte      Error code                 -}
{-                                                                          -}
{-             Errors:  0: No error                                         -}
{-                      1: Dimen < 1                                        -}
{-                      2: tolerance < TNNearlyZero                         -}
{-                      3: MaxIter < 1                                      -}
{-                      4: Mat NOT symmetric                                -}
{-                      5: Iter > MaxIter                                   -}
{-                                                                          -}
{----------------------------------------------------------------------------}

IMPLEMENTATION

PROCEDURE Jacobi{Dimen       : integer;
                 MaxIter      : integer;
                 tolerance    : Float;
             VAR Iter         : integer;
             VAR Error        : byte};

VAR
  Row, Column, Diag : integer;
  SinTheta, CosTheta : Float;
  SumSquareDiag : Float;
  Done : boolean;

PROCEDURE TestData(Dimen     : integer;
		   MaxIter   : integer;
                   tolerance : Float;
               VAR Error     : byte);

{---------------------------------------------------}
{- Input: Dimen, Mat, MaxIter, tolerance           -}
{- Output: Error                                   -}
{-                                                 -}
{- This PROCEDURE tests the input data for errors. -}
{---------------------------------------------------}

VAR
  Row, Column : integer;

BEGIN
  Error := 0;
  IF Dimen < 1 THEN
    Error := 1;
  IF tolerance <= TNNearlyZero THEN
    Error := 2;
  IF MaxIter < 1 THEN
    Error := 3;
  IF Error = 0 THEN
    FOR Row := 1 to Dimen - 1 DO
      FOR Column := Row + 1 to Dimen DO
        IF ABS(Mat^[Row, Column] - Mat^[Column, Row]) > TNNearlyZero THEN
          Error := 4   { Matrix NOT symmetric }
END; { PROCEDURE TestData }

PROCEDURE Initialize(Dimen : integer; VAR Iter : integer);

{--------------------------------------------}
{- Input: Dimen                             -}
{- Output: Iter, Eigenvectors               -}
{-                                          -}
{- This PROCEDURE initializes Iter to zero  -}
{- and Eigenvectors to the identity matrix. -}
{--------------------------------------------}

VAR
  i,j : integer;

BEGIN
  Iter := 0;
  FOR i:=1 to Dimen DO
    FOR j:= 1 to Dimen DO
      Eigenvectors^[i,j]:=0.0;
  FOR i := 1 to Dimen DO
    Eigenvectors^[i, i] := 1.0;
END; { PROCEDURE Initialize }

PROCEDURE CalculateRotation(RowRow   : Float;
                            RowCol   : Float;
                            ColCol   : Float;
                        VAR SinTheta : Float;
                        VAR CosTheta : Float);


{-----------------------------------------------------------}
{- Input: RowRow, RowCol, ColCol                           -}
{- Output: SinTheta, CosTheta                              -}
{-                                                         -}
{- This PROCEDURE calculates the sine and cosine of the    -}
{- angle Theta through which to rotate the matrix Mat.     -}
{- Given the tangent of 2-Theta, the tangent of Theta can  -}
{- be calculated WITH the quadratic formula.  The cosine   -}
{- and sine are easily calculable from the tangent. The    -}
{- rotation must be such that the Row, Column element is   -}
{- zero. RowRow is the Row,Row element; RowCol is the      -}
{- Row,Column element; ColCol is the Column,Column element -}
{- of Mat.                                                 -}
{-----------------------------------------------------------}

VAR
  TangentTwoTheta, TangentTheta, Dummy : Float;

BEGIN
  IF ABS(RowRow - ColCol) > TNNearlyZero THEN
    BEGIN
      TangentTwoTheta := (RowRow - ColCol) / (2 * RowCol);
      Dummy := Sqrt(Sqr(TangentTwoTheta) + 1);
      IF TangentTwoTheta < 0 THEN  { Choose the root nearer to zero }
        TangentTheta := -TangentTwoTheta - Dummy
      ELSE
        TangentTheta := -TangentTwoTheta + Dummy;
      CosTheta := 1 / Sqrt(1 + Sqr(TangentTheta));
      SinTheta := CosTheta * TangentTheta;
    END
  ELSE
    BEGIN
      CosTheta := Sqrt(1/2);
      IF RowCol < 0 THEN
        SinTheta := -Sqrt(1/2)
      ELSE
        SinTheta := Sqrt(1/2);
    END;
END; { PROCEDURE CalculateRotation }

PROCEDURE RotateMatrix(Dimen    : integer;
                       SinTheta : Float;
                       CosTheta : Float;
                       Row      : integer;
                       Col      : integer);

{--------------------------------------------------------------}
{- Input: Dimen, SinTheta, CosTheta, Row, Col                 -}
{- Output: Mat                                                -}
{-                                                            -}
{- This PROCEDURE rotates the matrix Mat through an angle     -}
{- Theta.  The rotation matrix is the identity matrix execept -}
{- FOR the Row,Row; Row,Col; Col,Col; and Col,Row elements.   -}
{- The rotation will make the Row,Col element of Mat          -}
{- to be zero.                                                -}
{--------------------------------------------------------------}

VAR
  CosSqr, SinSqr, SinCos : Float;
  MatRowRow, MatColCol, MatRowCol, MatRowIndex, MatColIndex : Float;

  Index : integer;

BEGIN
  CosSqr := Sqr(CosTheta);
  SinSqr := Sqr(SinTheta);
  SinCos := SinTheta * CosTheta;
  MatRowRow := Mat^[Row, Row] * CosSqr + 2 * Mat^[Row, Col] * SinCos +
	       Mat^[Col, Col] * SinSqr;
  MatColCol := Mat^[Row, Row] * SinSqr - 2 * Mat^[Row, Col] * SinCos +
	       Mat^[Col, Col] * CosSqr;
  MatRowCol := (Mat^[Col, Col] - Mat^[Row, Row]) * SinCos +
	       Mat^[Row, Col] * (CosSqr - SinSqr);

  FOR Index := 1 to Dimen DO
    IF NOT(Index IN [Row, Col]) THEN
    BEGIN
      MatRowIndex := Mat^[Row, Index] * CosTheta +
		     Mat^[Col, Index] * SinTheta;
      MatColIndex := -Mat^[Row, Index] * SinTheta +
		      Mat^[Col, Index] * CosTheta;
      Mat^[Row, Index] := MatRowIndex;
      Mat^[Index, Row] := MatRowIndex;
      Mat^[Col, Index] := MatColIndex;
      Mat^[Index, Col] := MatColIndex;
    END;
  Mat^[Row, Row] := MatRowRow;
  Mat^[Col, Col] := MatColCol;
  Mat^[Row, Col] := MatRowCol;
  Mat^[Col, Row] := MatRowCol;
END; { PROCEDURE RotateMat^rix }

PROCEDURE RotateEigenvectors(Dimen        : integer;
                             SinTheta     : Float;
                             CosTheta     : Float;
                             Row          : integer;
                             Col          : integer);
{--------------------------------------------------------------}
{- Input: Dimen, SinTheta, CosTheta, Row, Col                 -}
{- Output: Eigenvectors                                       -}
{-                                                            -}
{- This PROCEDURE rotates the Eigenvectors Matrix through an  -}
{- angle Theta.  The rotation Matrix is the identity Matrix   -}
{- except FOR the Row,Row; Row,Col; Col,Col; and Col,Row      -}
{- elements.  The Eigenvectors Matrix will be the product of  -}
{- all the rotation Matrices which operate on Mat.            -}
{--------------------------------------------------------------}
VAR
  EigenvectorsRowIndex, EigenvectorsColIndex : Float;

  Index : integer;

BEGIN
  { Transform eigenvector Mat^rix }
  FOR Index := 1 to  Dimen DO
  BEGIN
    EigenvectorsRowIndex := CosTheta * Eigenvectors^[Row, Index] +
                            SinTheta * Eigenvectors^[Col, Index];
    EigenvectorsColIndex := -SinTheta * Eigenvectors^[Row, Index] +
                             CosTheta * Eigenvectors^[Col, Index];
    Eigenvectors^[Row, Index] := EigenvectorsRowIndex;
    Eigenvectors^[Col, Index] := EigenvectorsColIndex;
  END;
END; { PROCEDURE RotateEigenvectors }

PROCEDURE NormalizeEigenvectors(Dimen : integer);
{---------------------------------------------------}
{- Input: Dimen, Eigenvectors                      -}
{- Output: Eigenvectors                            -}
{-                                                 -}
{- This PROCEDURE normalizes the eigenvectors so   -}
{- that the largest element IN each vector is one. -}
{---------------------------------------------------}
VAR
  Row : integer;
  Largest : Float;

PROCEDURE FindLargest(Dimen       : integer;
                      Eigenvector : TNvector;
		  VAR Largest     : Float);
{---------------------------------------}
{- Input: Dimen, Eigenvectors          -}
{- Output: Largest                     -}
{-                                     -}
{- This PROCEDURE returns the value of -}
{- the largest element of the vector.  -}
{---------------------------------------}
VAR
  Term : integer;

BEGIN
  Largest := Eigenvector[1];
  FOR Term := 2 to Dimen DO
    IF ABS(Eigenvector[Term]) > ABS(Largest) THEN
      Largest := Eigenvector[Term];
END; { PROCEDURE FindLargest }

PROCEDURE DivVecConst(Dimen       : integer;
                  VAR ChangingRow : TNvector;
                      Divisor     : Float);

{--------------------------------------------------------}
{- Input: Dimen, ChangingRow                            -}
{- Output: Divisor                                      -}
{-                                                      -}
{- elementary row operation - dividing by a constant    -}
{--------------------------------------------------------}

VAR
  Term : integer;

BEGIN
  FOR Term := 1 to Dimen DO
    ChangingRow[Term] := ChangingRow[Term] / Divisor;
END; { PROCEDURE DivVecConst }

BEGIN { PROCEDURE NormalizeEigenvectors }
  FOR Row := 1 to Dimen DO
  BEGIN
    FindLargest(Dimen, Eigenvectors^[Row], Largest);
    DivVecConst(Dimen, Eigenvectors^[Row], Largest);
  END;
END; { PROCEDURE NormalizeEigenvectors }


BEGIN { PROCEDURE Jacobi }
  TestData(Dimen, MaxIter, tolerance, Error);
  IF Error = 0 THEN
  BEGIN
    Initialize(Dimen, Iter);
    REPEAT
      Iter := Succ(Iter);
      SumSquareDiag := 0;
      FOR Diag := 1 to Dimen DO
	SumSquareDiag := SumSquareDiag + Sqr(Mat^[Diag, Diag]);
      Done := true;
      FOR Row := 1 to Dimen - 1  DO
        FOR Column := Row + 1 to Dimen DO
	  IF ABS(Mat^[Row, Column]) > tolerance * SumSquareDiag THEN
          BEGIN
            Done := false;
	    CalculateRotation(Mat^[Row, Row], Mat^[Row, Column],
			      Mat^[Column, Column], SinTheta, CosTheta);
	    RotateMatrix(Dimen, SinTheta, CosTheta, Row, Column);
	    RotateEigenvectors(Dimen, SinTheta, CosTheta, Row, Column);
          END;
    UNTIL Done OR (Iter > MaxIter);
    FOR Diag := 1 to Dimen DO
      Eigenvalues^[Diag] := Mat^[Diag, Diag];
    NormalizeEigenvectors(Dimen);
    IF Iter > MaxIter THEN
      Error := 5
  END;
END; { PROCEDURE Jacobi }

END. { Ujacobi }

{end of file ujacobi.pas}

