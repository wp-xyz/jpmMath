{ =====================  Unit fcholy.pas ============================= }
{   Procedures to solve a linear symmetric system by Cholesky method   }
{   See demonstration program Tcholy.pas.                              }
{ -------------------------------------------------------------------- }
{ Reference: "Numerical Algorithms with C, By Gisela Engeln-Muellges   }
{             and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].      }
{                                                                      }
{                             Pascal version by J-P Moreau, Paris.     }
{                                       (www.jpmoreau.fr)              }
{ ==================================================================== } 
UNIT Fcholy;

INTERFACE

CONST
        SIZE    = 25;     {to be increased if necessary by user}
        EPSQUAD = 1e-12;
        ZERO    = 0.0;

TYPE
        pMAT   = ^MATRIX;
        MATRIX = Array[0..SIZE,0..SIZE] of DOUBLE;
        pVEC   = ^VECTOR;
        VECTOR = Array[0..SIZE] of DOUBLE;


Procedure Choly(mode,n:INTEGER;VAR mat:pMAT;b:pVEC;VAR x:pVEC; VAR rc:INTEGER);
Procedure chodec(n:INTEGER;VAR mat:pMAT;VAR rc:INTEGER);
Procedure chosol(n:INTEGER;lmat:pMAT;b:pVEC;VAR x:pVEC;VAR rc:INTEGER);


IMPLEMENTATION
       

Procedure Choly(mode,n:INTEGER;VAR mat:pMAT;b:pVEC;VAR x:pVEC; VAR rc:INTEGER);
 {====================================================================*
 *                                                                    *
 *  The procedure cholesky solves linear systems :  mat * x = b       *
 *  for positive definite symmetric n x n matrices mat using the      *
 *  Cholesky method.                                                  *
 *                                                                    *
 *  mat must be symmetric and positive definite, or the method will   *
 *  fail. i.e. for all nonzero vectors y we must have y' * mat * y > 0*
 *                                                                    *
 *  cholesky uses only the lower triangle of mat.                     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Application:                                                     *
 *   ============                                                     *
 *                                                                    *
 *      Solve linear systems with positive definite system matrices   *
 *      efficeiently.                                                 *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Control parameter:                                               *
 *   ==================                                               *
 *      mode     integer                                              *
 *       = 0     Factor mat and solve linear system                   *
 *       = 1     Compute factorization only                           *
 *       = 2     Solve system only; for this call the factorization   *
 *               must be stored in mat from chodec. Used to solve for *
 *               many right hand sides.                               *
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer (n > 0)                                      *
 *               Dimension of mat, size of b and x                    *
 *      mat      REAL  mat[n,n];                                      *
 *                 mode = 0, 1: Matrix of the linear system           *
 *                 mode = 2   : Cholesky factor                       *
 *      b        REAL   b[n];          (for mode = 0, 2)              *
 *               Right hand side of system of equations               *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      mat      REAL   mat[n,n];       (for mode = 0, 1)             *
 *               Cholesky decomposition of input matrix mat           *
 *      x        REAL   x[n];           (for mode = 0, 2)             *
 *               solution vector                                      *
 *                                                                    *
 *   Return value rc:                                                 *
 *   ===============                                                  *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or false input parameter                       *
 *      = 2      Matrix is not positive definite                      *
 *      = 3      wrong value for mode                                 *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Procedures in use:                                               *
 *   =================                                                *
 *      chodec()  : Factorization                                     *
 *      chosol()  : Solving system                                    *
 *                                                                    *
 *====================================================================}
Begin
  rc:=3;

  if (mat = NIL) or (n < 1) then
  begin
    rc:=1;
    exit
  end;

  case mode of
  
    0: { Factor matrix and solve system ........................}
       begin chodec (n, mat, rc);
       if rc = 0 then
         chosol (n, mat, b, x, rc)
       end;

    1: { factor only ...........................................}
         chodec (n, mat, rc);

    2: { solve only ............................................}
         chosol (n, mat, b, x, rc)
  end

End;


Procedure chodec(n:INTEGER;VAR mat:pMAT;VAR rc:INTEGER);
 {====================================================================*
 *                                                                    *
 *  chodec decomposes the symmetric positive definite matrix mat.     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer (n > 0)                                      *
 *               Dimension of mat                                     *
 *      mat      REAL mat[n,n];                                       *
 *               Matrix of left hand coefficients                     *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      mat      REAL mat[n,n];                                       *
 *               Cholesky decomposition in the lower triangle         *
 *                                                                    *
 *   Return value rc:                                                 *
 *   ===============                                                  *
 *      = 0      all ok                                               *
 *      = 1      n < 1                                                *
 *      = 2      Matrix not  positive definite                        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use:                                                *
 *   ================                                                 *
 *                                                                    *
 *   From the Pascal library: sqrt()                                  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants in use:    EPSQUAD                                     *
 *   ================                                                 *
 *                                                                    *
 *====================================================================}
Var
  j, k,i:INTEGER;
  sum:DOUBLE;
Begin
  if n < 1 then
  begin
    rc:=1;                                       {n < 1  error  !  }
    exit
  end;
  if mat = NIL then
  begin
     rc:=1;
     exit
  end;

  if mat^[0][0] < EPSQUAD  then                   { mat not positive }
  begin                                           { definite         }
    rc:=2;
    exit
  end;
  mat^[0][0] := SQRT(mat^[0][0]);
  for j := 1 to n-1 do mat^[j][0] := mat^[j][0] / mat^[0][0];

  for i := 1 to n-1 do
  begin
    sum := mat^[i][i];
    for j := 0 to i-1 do  sum := sum - SQR(mat^[i][j]);

    if sum < EPSQUAD then                         { not positive definite }
    begin
      rc:=2;
      exit
    end;
    mat^[i][i] := SQRT(sum);
    for j := i + 1 to n-1 do
    begin
      sum := mat^[j][i];
      for k := 0 to i-1 do
        sum := sum - mat^[i][k] * mat^[j][k];
      mat^[j][i] := sum / mat^[i][i]
    end
  end;

  rc := 0
End;


Procedure chosol(n:INTEGER;lmat:pMAT;b:pVEC;VAR x:pVEC;VAR rc:INTEGER);
 {====================================================================*
 *                                                                    *
 *  chosol finds the solution  x  of the linear system B' *  B * x = b*
 *  for a lower triangular nonsingular matrix B as supplied in chodec.*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer  (n > 0)                                     *
 *               Dimension of lmat, size of b and x                   *
 *      lmat     REAL  lmat[n,n]                                      *
 *               lower triangular matrix as supplied by  chodec       *
 *      b        REAL   b[n]                                          *
 *               Right hand side                                      *
 *                                                                    *
 *   Output parameter:                                                *
 *   =================                                                *
 *      x        REAL   x[n]                                          *
 *               solution vector                                      *
 *                                                                    *
 *   Return value rc:                                                 *
 *   ===============                                                  *
 *      = 0      all ok                                               *
 *      = 1      improper lwer triangular matrix or  n < 1            *
 *                                                                    *
 *====================================================================}
Var
  j, k:INTEGER;
  sum:DOUBLE;
Begin
  if n < 1 then   {n<1 error !}
  begin
    rc:=1;
    exit
  end;
  if (lmat = NIL) or (b = NIL) or (x = NIL) then
  begin
    rc:=1;
    exit
  end;

  if lmat^[0][0] = ZERO then             { improper factor matrix }
  begin
    rc:=1;
    exit
  end;

  x^[0] := b^[0] / lmat^[0][0];          { update right hand side }
  for k := 1 to n-1 do
  begin
    sum := ZERO;
    for j := 0 to k-1 do
      sum := sum + lmat^[k][j] * x^[j];
    if lmat^[k][k] = ZERO then
    begin
      rc:=1;
      exit
    end;
    x^[k] := (b^[k] - sum) / lmat^[k][k]
  end;

  x^[n-1] := x^[n-1] / lmat^[n-1][n-1];   { back substitution }
  for k:=n-2 downto 0 do
  begin
    sum := ZERO;
    for j := k + 1 to n-1 do
      sum := sum + lmat^[j][k] * x^[j];
    x^[k] := (x^[k] - sum) / lmat^[k][k]
  end;

  rc := 0
End;

END.

{ --------------------- END of unit fcholy.pas ------------------------ }