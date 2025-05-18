{***********************************************************************
* This is a test program for the procedure cg_method to solve a linear *
* system                                                               *
*                       A * X  =  Y                                    *
* for a symmetric positive definite matrix A using the conjugate       *
* gradient method.                                                     *
* -------------------------------------------------------------------- *
* Ref.: "Numerical algorithms with C,  By Gisela Engeln-Muellges and   *
*        Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].               *
* -------------------------------------------------------------------- *
* SAMPLE RUN:                                                          *
*                                                                      *
* Solve the following linear system:                                   *
*                                                                      *
* 4                     size N of matrix                               *
* 5.0 -1.0 -1.0 -1.0    the upper triangle (with diagonal)             *
*      5.0 -1.0 -1.0    of the positive definite matrix A              *
*           5.0 -1.0                                                   *
*                5.0                                                   *
* 1.0                   right hand side of the linear system           *
* 1.0                                                                  *
* 1.0                                                                  *
* 1.0                                                                  *
*                                                                      *
*                                                                      *
* -----------------------------------------------------                *
*  Symmetric Linear System (Conjugate Gradient Method)                 *
* -----------------------------------------------------                *
*  Test data (Matrix and right hand side):                             *
*    5.00000  -1.00000  -1.00000  -1.00000  1.00000                    *
*   -1.00000   5.00000  -1.00000  -1.00000  1.00000                    *
*   -1.00000  -1.00000   5.00000  -1.00000  1.00000                    *
*   -1.00000  -1.00000  -1.00000   5.00000  1.00000                    *
*                                                                      *
*  Solution for the linear system:                                     *
*    0.50000   0.50000   0.50000   0.50000                             *
* -----------------------------------------------------                *
*                                                                      *
*        TPW version without dynamic allocations by J-P Moreau, Paris  *
*                            (www.jpmoreau.fr)                         *
***********************************************************************}
PROGRAM CGTST1;
USES WinCrt;

CONST   SIZE     = 25;             {maximum size of linear system}
        MACH_EPS = 1E-12;           {smallest machine real number}

TYPE
        MAT = Array[0..SIZE,0..SIZE] of REAL;      {square matrix}
        VEC = Array[0..SIZE] of REAL;                     {vector}

VAR
  a : MAT;        { the upper triangle of a positive definite real   }
                  { matrix (here stored in a vector).                }
  y : VEC;        { right hand side of the linear system             }
  x : VEC;        { solution vector of the system                    }
  n,              { size of matrix A                                 }
  i,              { row index                                        }
  j,              { column index                                     }
  fehler:INTEGER; { return value of cg_method                        }


{***********************************************************************
* CG_METHOD solves the linear system                                   *
*                         A * X = Y                                    *
* for a symmetric, positive definite matrix A via the conjugate        *
* gradient method.                                                     *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n  Size of the linear system                                         *
* a  [0..n-1,0..n-1] system matrix A. Only the upper triangle of A is  *
*    used.                                                             *
* y  [0..n-1] vector of the right hand side                            *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x  [0..n-1] vector giving the solution                               *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all is ok                                                       *
* = 1: n < 2 or other disallowed input parameters                      *
* = 2: memory exceeded (not used here)                                 *
*                                                                      *
***********************************************************************}
PROCEDURE cg_method(n:INTEGER; a:MAT; y:VEC; VAR x:VEC; VAR code:INTEGER);
LABEL  fin;
VAR
       d:VEC;         { [0..n-1] auxiliary vectors d and g            }
       g:VEC;         {                                               }
       AmalD:VEC;     { [0..n-1] auxiliary vector A * d               }
       alpha,         { coefficient                                   }
       beta,          { coefficient                                   }
       dividend,      { numerator and denominator of a fraction,      }
       divisor,       { respectively, used to compute alpha, beta     }
       hilf,          { auxiliary variables                           }
       hilf2,         {                                               }
       abstand,       { distance of two successive approximations     }
                      { for the solution vector x (taken in the       }
                      { euclidean norm)                               }
       xnorm : REAL;  { euklidean norm of x                           }
       k,i,j:INTEGER; { loop variables                                }

Begin
  if n < 2 then                           { invalid parameter?       }
  begin
    code := 1;
    exit
  end;

  {------------------------------------------------------------------}
  { start with x at the origin                                       }
  {------------------------------------------------------------------}

  for i := n - 1 downto 0 do
    x[i] := 0.0;
      
  {------------------------------------------------------------------}
  { initialize  d and g :                                            }
  { d = -g = -(a*x - y) = y (since x = 0)                            }
  {------------------------------------------------------------------}

  for i := n - 1 downto 0 do
  begin
    hilf := y[i];
    d[i] := hilf;
    g[i] := -hilf
  end;
      
  {------------------------------------------------------------------}
  { perform at most n steps of the CG Method                         }
  {------------------------------------------------------------------}

  for k := n downto 1 do
  begin

    {----------------------------------------------------------------}
    { compute new alpha:                                             }
    { alpha = -(d(transp) * g) / (d(transp) * (a * d))               }
    {----------------------------------------------------------------}

    dividend := 0.0;
    divisor  := 0.0;

    for i := n - 1 downto 0 do
    begin
      dividend := dividend + d[i] * g[i];
      hilf := 0.0;
      for j := 0 to i-1 do
        hilf := hilf + a[j,i] * d[j];
      for j := i to n-1 do
        hilf := hilf + a[i,j] * d[j];
      AmalD[i] := hilf;
      divisor := divisor + d[i] * hilf
    end;

    if divisor = 0.0 then goto fin;

    alpha := -dividend / divisor;
    
    {----------------------------------------------------------------}
    { compute the norm of x und  alpha * d  and find a new x:        }
    { x  =  x + alpha * d, then check whether x is close enough,     }
    { in order to stop the process before n complete steps           }
    {----------------------------------------------------------------}

    xnorm   := 0.0;
    abstand := 0.0;

    for i := n - 1 downto 0 do
    begin
      hilf    :=  x[i];
      xnorm   := xnorm + sqr(hilf);
      hilf2   :=  alpha * d[i];
      abstand := abstand + sqr(hilf2);
      x[i]    :=  hilf + hilf2
    end;

    if abstand < MACH_EPS * xnorm then goto fin;

    {----------------------------------------------------------------}
    { compute new g:   g  =  g + alpha * (a * d)                     }
    {----------------------------------------------------------------}

    for i := n - 1 downto 0 do
      g[i] := g[i] + alpha * AmalD[i];

    {----------------------------------------------------------------}
    { compute new beta :                                             }
    { beta = (g(transp) * (a * d)) / (d(transp) * (a * d))           }
    {----------------------------------------------------------------}

    dividend := 0.0;

    for i := n - 1 downto 0 do
      dividend := dividend + g[i] * AmalD[i];

    beta := dividend / divisor;

    {----------------------------------------------------------------}
    { compute new d :   d  =  - g + beta * d                         }
    {----------------------------------------------------------------}

    for i := n - 1 downto 0 do
      d[i] := -g[i] + beta * d[i]

  end; {k loop}

  fin:
  code := 0
End;


{main program}
BEGIN

  { ----------------------------- Write header --------------------- }
  Writeln('-----------------------------------------------------');
  Writeln(' Symmetric Linear System (Conjugate Gradient Method) ');
  Writeln('-----------------------------------------------------');

  { ----------------------------- Set inputs ----------------------- }
  n:=4;                                      { size of linear system }

  {define a matrix}
  a[0,0]:=5.0; a[0,1]:=-1.0; a[0,2]:=-1.0; a[0,3]:=-1.0;
               a[1,1]:= 5.0; a[1,2]:=-1.0; a[1,3]:=-1.0;
                             a[2,2]:= 5.0; a[2,3]:=-1.0;
                                           a[3,3]:= 5.0;
  {define y vector}
  y[0]:=1.0; y[1]:=1.0; y[2]:=1.0; y[3]:=1.0;

  { -------------- Print input data as a safeguard ---------------- }

  writeln(' Test data (Matrix and right hand side):');
  for i := 0 to n-1 do
  begin
    for j := 0 to i-1 do
      write(a[j,i]:9:5);
    for j := i to n-1 do
      write(a[i,j]:9:5);
    writeln('  ',y[i]:9:5)
  end;
      
  { ---------------- solve linear system --------------------------- }

  cg_method(n, a, y, x, fehler);               { perform  CG method  }


  { -------------------- Print results ----------------------------- }

  if fehler > 0 then                           { Error in CG method? }
    Writeln(' Error executing CG method.')
  else
  begin
    Writeln(' Solution for the linear system:');
    for i := 0 to n-1 do write(x[i]:9:5)
  end;       

  { -------------------- write footer ------------------------------ }
  Writeln;
  Writeln('-----------------------------------------------------');
  Writeln;
  { exit program }
  Readkey; DoneWinCrt

 END.

{end of file cgtst1.pas}