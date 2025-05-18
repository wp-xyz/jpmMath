{ ------------------------ UNIT brown.pas ------------------------------

************************************************************************
*                                                                      *
* Brown's method for nonlinear systems                                 *
* -------------------------------------                                *
*                                                                      *
* Programing language: Pascal                                          *
* Compiler:            Borland TPW 1.5                                 *
* Computer:            IBM Compatible PC, 400 MHz                      *
* Literature:          Brown, K. M.:                                   *
*                      A quadratically convergent Newton-like method   *
*                      based upon Gaussian elimination                 *
*                      SIAM J. Numer. Anal. Vol. 6 , (1969), p. 560    *
* Source:              QuickBASIC program                              *
* Author:              Johannes Karfusehr (FORTRAN)                    *
* Adaption:            Juergen Dietel, Computer Center, RWTH Aachen    *
*                      for C version                                   *
*                      J-P Moreau, Paris for this Pascal version       *
* Date:                9.16.2002                                       *
*                      (www.jpmoreau.fr)                               *
* -------------------------------------------------------------------- *
* Ref.: "Numerical Algorithms with C By G. Engeln-Mueller and F. Uhlig,*
*        Springer-Verlag, 1996" [BIBLI 11].                            *
***********************************************************************}
UNIT Brown;

INTERFACE

Uses Basis;

Const
      EIGHT= 8.0;
      HALF = 0.5;
      ONE  = 1.0;
      PI   = 3.1415926535;
      SIZE = 24;
      TEN  = 10.0;
      TWO  = 2.0;
      ZERO = 0.0;

Type  REAL = Double;
      MAT  = Array[0..SIZE,0..SIZE] of REAL;
      IMAT = Array[0..SIZE,0..SIZE+1] of INTEGER;
      VEC  = Array[0..SIZE] of REAL;

Var
      ihf: IMAT;
      dquot, rslin: VEC;
      hf: MAT;

Procedure brown1  { Brown's method for nonlinear systems of equations }
         (
          n:      integer;   { number of equations ...................}
          x0:     VEC;       { Starting value for iteration ..........}
          eps:    REAL;      { error bound ...........................}
          prot:   Integer;   { Protokol switch .......................}
          maxit:  integer;   { maximal number of steps ...............}
      VAR x1:     VEC;       { solution ..............................}
      VAR itanz:  integer;   { actual steps performed ................}
      VAR rc:     Integer    { error code ............................}
         );


IMPLEMENTATION


Function func0(k:Integer; x:VEC): REAL;
{***********************************************************************
*                         Test example 0                               *
***********************************************************************}
Begin
  Case k of
    0:  func0 := sqr(x[0]) - x[1] - ONE;
    1:  func0 := sqr(x[0] - TWO) + sqr(x[1] - HALF) - ONE
  End
End;

Function func1(k:Integer; x:VEC): REAL;
{***********************************************************************
*                         Test example 1                               *
***********************************************************************}
Begin
  Case k of
    0:  func1 := x[0] - x[1] + sin(x[0]) + cos(x[1]) - PI/EIGHT;
    1:  func1 := sqr(x[0]) - 4*sqr(x[1]) - Ln(sqr(x[2]));
    2:  func1 := x[2]-1 + 2*x[0] - 4*x[1]
  End
End;

Function func(k:Integer; x:VEC): REAL;
{***********************************************************************
*                         Test example 3                               *
***********************************************************************}
Begin
  Case k of
    0: func := 10*x[0] + x[1] + x[2] + x[3] -20 + sin(x[0])*sin(x[0]) + cos(x[1])*cos(x[1]);
    1: func := x[0] + 20*x[1] + x[2] + x[3] -48 + 1/(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
    2: func := (x[0]+x[1])*(x[0]+x[1]) + 30*x[2] + x[3] -97 + Ln(x[0]) + Ln(x[1]+x[2]);
    3: func := x[0] + x[1] + x[2] + 40*x[3] -166 + x[0]*x[0]
  End
End;

{ -------------------------------------------------------------------- } 
Procedure subst (
                  n:integer;
                  k:integer;
                  ihf:IMAT;
                  hf:MAT;
                  rslin:VEC;
              VAR x1:VEC
                 );
{***********************************************************************
* Solve a linear system of equations                                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* n      size of system                                                *
* k      Index for coordinates: 0,...,n-1                              *
* ihf    [0..n-1,0..n] matrix, register of row and column interchanges *
* hf     [0..n-1,0..n-1] matrix, the system matrix                     *
* rslin  [0..n-1] right nahd side vector                               *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* x1     [0..n-1] solution vector                                      *
*                                                                      *
* Constants in use :  REAL, ZERO                                       *
* ==================                                                   *
*                                                                      *
***********************************************************************}
Var
       sum: REAL;        { aux variable for finding x1[kmax]         }
       km,               { Loop counter                              }
       kmax,             { original row index                        }
       jsub,             { original column index                     }
       j:Integer;        { Loop counter                              }
Begin

  for km := k Downto 1 do
  begin
    kmax := ihf[km - 1,n];
    sum:=ZERO;
    for j := km to n-1 do
    begin
      jsub :=  ihf[km,j];
      sum  := sum + hf[km - 1,jsub] * x1[jsub]
    end;
    x1[kmax] := sum + rslin[km - 1]
  end
End;


{ -------------------------------------------------------------------- }
Procedure iter4
                (
                 n: integer;
                 epsm: REAL;
                 xalt: VEC;
             VAR x1: VEC;
             VAR sing: boolean;
             VAR rc: Integer
                );
{***********************************************************************
* Compute one approximation via  Brown's method.                       *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* fkt   Function, which evaluates the right hand side. It has          *
*       components from  0 to n - 1.                                   *
* n     number of equations                                            *
* epsm  machine constant                                               *
* xalt  previous iterate                                               *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x1    [0..n-1] vector, the new iterate for the zero                  *
* sing  error code, matrix is singular                                 *
*                                                                      *
* Return code rc:                                                      *
* ==============                                                       *
* = 0: all ok                                                          *
* = 1: error in evaluating func                                        *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* ihf, rslin, dquot, hf, nlgls, REAL, boolean, subst, ABS, FALSE,      *
* TRUE, ZERO                                                           *
***********************************************************************}
Label 10, 20;
Var
       i,          { Loop counters                                    }
       j,
       k,          { index of current function                        }
       anzahl,     { counter for those difference quotients that are 0}
       temp,       { aux variable                                     }
       kmax,       { Index of largest difference quotient             }
       jsub:Integer; { aux variable for an index in ihf               }
       hold,       { aux variable for x1[temp]                        }
       h,          { Step size for difference quotient in direction of}
                   { x1[temp]                                         }
       faktor,     { quotient  h / hold                               }
       dermax,     { maximal difference quotient magnitude            }
       sum,        { sum in rslin[k]                                  }
       f,          { value of function k at  x1                       }
       fplus:REAL; { value of function k at                           }
                   { (x1[0],..., x1[temp]+h,...,x1[n-1])              }
Begin

  { --------------------- initialize variables ----------------------- }
  temp:=0; kmax:=0;
  for j := 0 to n-1 do
  begin
    ihf[0,j] := j;
    for i:=0 to n-1 do x1[n]:=xalt[n]
  end;

  { ----------- linearize the kth coordinate function ---------------- }
  for k := 0 to n-1 do
  begin

    anzahl := 0;
    faktor := 0.001;

    for j := 0 to 2 do
    begin
      if k > 0 then
        subst(n, k, ihf, hf, rslin, x1);
      f:=func(k, x1);

      { ---- find ith diskretization size and difference quotient ---- }
      for i := k to n-1 do
      begin
        temp := ihf[k,i];
        hold := x1[temp];
        h    := faktor * hold;
        if ABS(h) <= epsm then
          h := 0.001;
        x1[temp] := hold + h;
        if k > 0 then
          subst(n, k, ihf, hf, rslin, x1);
        fplus:=func(k, x1);
        x1[temp] := hold;

        dquot[temp] := (fplus - f) / h;

        if ABS(dquot[temp]) <= epsm then
          Inc(anzahl)
        else if ABS(f / dquot[temp]) >= 1e20 then
          Inc(anzahl)
      end;

      if anzahl < n - k then
      begin
        sing := FALSE;
        goto 10;
      end
      else
      begin
        sing  :=  TRUE;
        faktor := faktor * TEN;
        anzahl :=  0
      end
    end;


10: if sing then goto 20
    else if k < n - 1 then
    begin
      kmax := ihf[k,k];

      { --- find largest magnitude difference quotient ------------ }
      dermax := ABS(dquot[kmax]);
      for i := k + 1 to n-1 do
      begin
        jsub := ihf[k,i];
        if ABS(dquot[jsub]) < dermax then
          ihf[k + 1,i] := jsub
        else
        begin
          ihf[k + 1,i] := kmax;
          kmax         := jsub
        end
      end;
      if ABS(dquot[kmax]) <= epsm then
        sing := TRUE;

      ihf[k,n] := kmax;


      if sing then goto 20
      else
      begin
        { ---------- solve kth equation for  xmax ------------------ }
        sum := ZERO;
        for j := k + 1 to n-1 do
        begin
          jsub       := ihf[k + 1,j];
          hf[k,jsub] := -dquot[jsub] / dquot[kmax];
          sum        := sum + dquot[jsub] * x1[jsub]
        end;
        rslin[k] := (sum - f) / dquot[kmax] + x1[kmax]
      end
    end
    else
    begin
      { ----- solve (n-1)th coordinate function via discrete ------- 
        ----- Newton method for one variable ----------------------- }
      if ABS(dquot[temp]) <= epsm then
        sing := TRUE
      else
      begin
        kmax     := temp;
        rslin[k] := -f / dquot[kmax] + x1[kmax]
      end
    end;
20: end;

  if Not sing then
  begin
    x1[kmax] := rslin[n - 1];                 { compute approximation }
    if n > 1 then                            { by back substitution  }
      subst(n, n - 1, ihf, hf, rslin, x1)
  end;
  rc:=0
End;

{ ------------------------------------------------------------------- }
Procedure brown1  { Brown's method for nonlinear systems of equations }
         (
          n:      Integer;   { number of equations ...................}
          x0:     VEC;       { Starting value for iteration ..........}
          eps:    REAL;      { error bound ...........................}
          prot:   Integer;   { Protokol switch .......................}
          maxit:  integer;   { maximal number of steps ...............}
      VAR x1:     VEC;       { solution ..............................}
      VAR itanz:  integer;   { actual steps performed ................}
      VAR rc:     Integer    { error code ............................}
         );
{***********************************************************************
* find a zero of a nonlinear system of n equations in n unknown using  *
* Brown's method                                                       *
*                                                                      *
* One step of Brown's method is performed by calling iter4().          *
* The iteration is continued until the admissable number of iterations *
* is reached or one of the following break-off criteria is satisfied : *
* 1. the relative change between two successive iterates is less than  *
*    eps.                                                              *
* 2. The function value has mgnitude less than MACH_EPS at the new x   *
* 3. The desired accuracy has been reached.                            *
*                                                                      *
* If specified in prot, intermediate results are tabulated, see prot   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* fkt    Function, whichdescribes the right hand side of the system;   *
*        its components are numbered 0 to n-1.                         *
* n      size of system                                                *
* x0     [0..n-1] starting vector for the iteration                    *
* eps    desired accuracy; if specified as <= 0, we set                *
*        eps = 0.8 * MACH_EPS.                                         *
* prot   Protocol flag. If TRUE we tabulate the differenz of successive*
*        iterates, the last iterate and the functional value there     *
*        after each iteration in the standard output file. If FALSE,   *
*        there is no output.                                           *
* maxit  Maximal number of iterations                                  *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x1     [0..n-1] vector, the approximate solution for the zero        *
* itanz  number of iterations executed                                 *
*                                                                      *
* Return value rc:                                                     *
* ===============                                                      *
* error code.                                                          *
* = 0: successsful iteration                                           *
* = 1: desired accuracy not achieved after maxit iterations            *
* = 2: system matrix singular                                          *
* = 3: lack of memory                                                  *
* = 4: wrong input : fkt = NULL or n < 1 or maxit < 1                  *
* = 5: error in evaluating  fkt()                                      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* ihf, rslin, dquot, hf, nlgls, REAL, boolean, MACH_EPS, FALSE,        *
* iter4, subst, IMATRIX, MATRIX, writeln, ABS, EIGHT                   *
***********************************************************************}
Label 30,35,40,50;
Var
    i,             { Loop variable                                     }
    m,             { iterations counter                                }
    krit:Integer;  { list of causes for stop of program:               }
                   { = 0: too many steps                             
                     = 1: step size too small                        
                     = 2: Function value sufficiently small          
                     = 3: desired accuracy reached                     }
    sing:boolean;  { Flag, indicating whether the linearized system    
                     is singular                                       }
    relf,          { aux variables                                     }
    fwert,         {                                                   }
    delta0,        { maximum norm of previous vector                   }
    delta1,        { maximum norm of current vector                    }
    epsm:REAL;     { machine constant                                  }
    xalt:VEC;      { [0..n-1] vector used to store old iterate         }
Begin
  relf:=ZERO; 
  { ---------------- catch wrong input ------------------------------- }
  if (n < 1) or (maxit < 1) then
  begin
    rc := 4;
    goto 50  {exit}
  end;

  if eps < MACH_EPS then           { desired accuracy rediculously }
    eps := EIGHT * MACH_EPS;       { small ? correct }

  epsm   := 100*MACH_EPS;          { initialize variables }
  sing   := FALSE;
  krit   := 0;
  delta0 := 0.01;
  for i:=0 to n-1 do xalt[n]:=x0[n];

  if prot<>0 then writeln;

  for m := 1 to maxit do           { start iteration }
  begin

    iter4(n, epsm, xalt, x1, sing, rc);
    if rc<>0 then
    begin
      rc:=5;
      goto 50 {exit}
    end;

    { ----------- if desired document each step ---------------------- }
    if prot<>0 then
    begin
      writeln(m:3,' Iteration step');
      if Not sing then
      begin
        writeln('      Difference         Component       Approximation             Function value');
        for i := 0 to n-1 do
        begin
          fwert:=func(i,x1);
          writeln(x1[i] - xalt[i], i, x1[i], fwert)
        end
      end
      else
        writeln(' Jacobi matrix singular!')
    end;

    if Not sing then
    begin

      { ----------------- test break-off criteria -------------------- }
      for i := 0 to n-1 do       { test relative change in new iterate }
      begin
        relf := (x1[i] - xalt[i]) / (xalt[i] + eps);
        if ABS(relf) >= eps then goto 30
      end;
30:   if ABS(relf) < eps then
      begin
        krit := 1;                         { step size too small }
        goto 40
      end;

      for i := 0 to n-1 do                 { check function value }
      begin
        fwert:=func(i, x1);
        if ABS(fwert) > epsm then
          goto 35
      end;

35:   if ABS(fwert) <= epsm then
      begin
        krit := 2;                       { function value small enough }
        goto 40
      end;

      delta1 := ABS(x1[0] - xalt[0]);            { compare with desired }
      for i := 1 to n-1 do                      { accuracy             }
        if delta1 < ABS(x1[i] - xalt[i]) then
          delta1 := ABS(x1[i] - xalt[i]);
      if delta1 <= 0.001 then
        if delta0 <= delta1 then
        begin
          krit := 3;                        { desired accuracy reached }
          goto 40
        end;
      delta0 := delta1;

      if m < maxit then
        for i:=0 to n-1 do xalt[i]:=x1[i]
    end
    else
      goto 40;

  end; {of m loop}


40: itanz := m;


  if sing then
    rc:= 2
  else if krit = 0 then
    rc:= 1
  else
    rc:=0;
50:End;

END.
{ -------------------------- END brown.pas -----------------