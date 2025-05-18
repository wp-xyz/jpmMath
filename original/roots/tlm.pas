{********************************************************************
*  THIS PROGRAM TESTS CODES FOR THE LEAST-SQUARES SOLUTION OF       *
*  M NONLINEAR EQUATIONS IN N VARIABLES. IT CONSISTS OF A DRIVER    *
*  AND AN INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA,       *
*  CALLS THE NONLINEAR LEAST-SQUARES SOLVER, AND FINALLY PRINTS     *
*  OUT INFORMATION ON THE PERFORMANCE OF THE SOLVER. THIS IS        *
*  ONLY A SAMPLE DRIVER, MANY OTHER DRIVERS ARE POSSIBLE. THE       *
*  INTERFACE SUBROUTINE FCN IS NECESSARY TO TAKE INTO ACCOUNT THE   *
*  FORMS OF CALLING SEQUENCES USED BY THE FUNCTION AND JACOBIAN     *
*  SUBROUTINES IN THE VARIOUS NONLINEAR LEAST-SQUARES SOLVERS.      *
*                                                                   *
*  PROCEDURES OR FUNCTIONS CALLED (SEE UNIT LM.PAS):                *
*                                                                   *
*    INITPT, SSQFCN, LMDIF1, ENORM, FCN                             *
*                                                                   *
*  From F77 program By:                                             *
*  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.        *
*  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE            *
* ----------------------------------------------------------------- *
* SAMPLE RUNS:                                                      *
* Example #1 (size 2)                                               *
* Solve following system (with initial conditions x1=1, x2=1):      *
*   X1^2 + X1 + X2^2 - 2           = 0                              *
*   X1^2 + X2 - X2^2 - 1 + log(X1) = 0                              *
*                                                                   *
* Input #example(1 to 3), size N, Number of tries: 1 2 1            *
*    PROBLEM    1    DIMENSIONS    2    2                           *
*    INITIAL L2 NORM OF THE RESIDUALS:  1.00000000000000E+0000      *
*    FINAL L2 NORM OF THE RESIDUALS..:  2.59810676458107E-0017      *
*    NUMBER OF FUNCTION EVALUATION...: 25                           *
*    NUMBER OF JACOBIAN EVALUATIONS..: 6                            *
*    EXIT PARAMETER..................: 2                            *
*    FINAL APPROXIMATE SOLUTION......:                              *
*     0.915554     0.496191                                         *
*                                                                   *
*  SUMMARY OF 1 CALL(S) TO LMDIF1:                                  *
*  NPROB   N    M   NFEV  NJEV  INFO   FINAL L2 NORM                *
*    1     2    2    25     6    2  2.59810676458107E-0017          *
*                                                                   *
* Example #2 (size 4)                                               *
* Solve following system (with initial conditions x1..x4 = 1):      *
*  10.0*x + x2 + x3  + x4 - 20.0 + Sqr(sin(x1)) + Sqr(cos(x2)) = 0  *
*  x1 + 20.0*x2 + x3 + x4 - 48.0 + one/pow^6                   = 0  *
*  Sqr(x1 + x2) + 30.0*x3 + x4 - 97.0 + log(x1) + log(x2+x3)   = 0  *
*  x1     + x2  + x3 + 40.0*x4 - 166.0 + Sqr(x1)               = 0  *
*                                                                   *
* Input #example(1 to 3), size N, Number of tries: 2 4 1            *
*    PROBLEM    2    DIMENSIONS    4    4                           *
*    INITIAL L2 NORM OF THE RESIDUALS:  1.38760694011757E+0002      *
*    FINAL L2 NORM OF THE RESIDUALS..:  2.74527174267237E-0015      *
*    NUMBER OF FUNCTION EVALUATION...: 31                           *
*    NUMBER OF JACOBIAN EVALUATIONS..: 5                            *
*    EXIT PARAMETER..................: 2                            *
*    FINAL APPROXIMATE SOLUTION......:                              *
*     1.040648     1.972398     2.745049     3.978974               *
*                                                                   *
*  SUMMARY OF 1 CALL(S) TO LMDIF1:                                  *
*  NPROB   N    M   NFEV  NJEV  INFO   FINAL L2 NORM                *
*    2     4    4    31     5    2  2.74527174267237E-0015          *
*                                                                   *
* Example #3 (size 6) - Stiff system                                *
* Solve following system (with initial conditions x1..x6 = 1):      *
*     X1 + X2 + X4 - .001 = 0                                       * 
*     X5 + X6 -55         = 0                                       *
*     X1 + X2 + X3 + 2X5 + X6 - 110.001 = 0                         *
*     X1 - 0.1X2        = 0                                         *
*     X1 - 10000 X3 X4  = 0                                         *
*     X5 - 5.5e15 X3 X6 = 0                                         *
*                                                                   *
* Input #example(1 to 3), size N, Number of tries: 3 6 1            *
*    PROBLEM    3    DIMENSIONS    6    6                           *
*    INITIAL L2 NORM OF THE RESIDUALS:  5.50000000000000E+0015      *
*    FINAL L2 NORM OF THE RESIDUALS..:  2.20108725607505E-0015      *
*    NUMBER OF FUNCTION EVALUATION...: 161                          *
*    NUMBER OF JACOBIAN EVALUATIONS..: 20                           *
*    EXIT PARAMETER..................: 2                            *
*    FINAL APPROXIMATE SOLUTION......:                              *
*     0.000083     0.000826     0.000091     0.000091               *
*    55.000000     0.000000                                         * 
*                                                                   *
*  SUMMARY OF 1 CALL(S) TO LMDIF1:                                  *
*  NPROB   N    M   NFEV  NJEV  INFO   FINAL L2 NORM                *
*    3     6    6   161    20    2  2.20108725607505E-0015          *
*                                                                   *
*                             Pascal Release By J-P Moreau, Paris.  *
*                                      (www.jpmoreau.fr)            *
********************************************************************}
PROGRAM test_lmdif;
Uses WinCrt1, LM;
Var
    i, ic, info, k, m, n, ntries: Integer;
    iwa, ma, na, nf, nj, np, nx: pIVec;
    factor, fnorm1, fnorm2, tol: Double;
    fnm, fvec, x: pVec;


Procedure initpt (n:Integer; x:pVec; nprob:Integer; factor:double);
{  **************************************************************
!  Procedure INITPT

!  THIS Procedure SPECIFIES THE STANDARD STARTING POINTS FOR THE
!  FUNCTIONS DEFINED BY Procedure SSQFCN. THE Procedure RETURNS
!  IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING POINT. FOR
!  THE 11TH FUNCTION THE STANDARD STARTING POINT IS ZERO, SO IN
!  THIS CASE, IF FACTOR IS NOT UNITY, THEN THE Procedure RETURNS
!  THE VECTOR  X(J) = FACTOR, J=1,...,N.

!  THE Procedure STATEMENT IS

!    PROCEDURE INITPT(N,X,NPROB,FACTOR)

!  WHERE

!    N IS A POSITIVE INTEGER INPUT VARIABLE.

!    X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE STANDARD
!      STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY FACTOR.

!    NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
!      NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.

!    FACTOR IS AN INPUT VARIABLE WHICH SPECIFIES THE MULTIPLE OF
!      THE STANDARD STARTING POINT.  IF FACTOR IS UNITY, NO
!      MULTIPLICATION IS PERFORMED.

!  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!  ************************************************************** }
Var
    j: Integer;
    h: double;
Const
    half = 0.5; three = 3.0; seven = 7.0; twenty = 20.0; twntf = 25.0;
    c1 = 1.2; c2 = 0.25; c3 = 0.39; c4 = 0.415; c5 = 0.02; c6 = 4000.0;
    c7 = 250.0; c8 = 0.3; c9 = 0.4; c10 = 1.5; c11 = 0.01; c12 = 1.3;
    c13 = 0.65; c14 = 0.7; c15 = 0.6; c16 = 4.5; c17 = 5.5;
Begin
{ SELECTION OF INITIAL POINT }
  For j:=1 to n do x^[j] := one;

{ COMPUTE MULTIPLE OF INITIAL POINT }

  IF factor = one then Exit;
  For j:=1 to n do x^[j] := factor*x^[j];

End; {initpt}


{main program}
Begin
  New(iwa); New(ma); New(na); New(nf); New(nj); New(np); New(nx);
  New(fnm); New(fvec); New(x);
  tol := 1E-8;
  ic := 0;
  Writeln;
  Write(' Input #example(1 to 3), size N, Number of tries: ');
  Readln(nprob, n, ntries);
  m:=n;
  factor := one;
  For k:=1 to ntries do
  begin
    ic := ic+1;
    initpt (n, x, nprob, factor);
    ssqfcn (m, n, x, fvec, nprob);
    fnorm1 := enorm(m, fvec);
    Writeln('    PROBLEM',nprob:5,'    DIMENSIONS', n:5, m:5);
    nfev := 0;
    njev := 0;
    lmdif1 (m, n, x, fvec, tol, info, iwa);
    ssqfcn (m, n, x, fvec, nprob);
    fnorm2 := enorm(m,fvec);
    np^[ic] := nprob;
    na^[ic] := n;
    ma^[ic] := m;
    nf^[ic] := nfev;
    njev := njev div n;
    nj^[ic] := njev;
    nx^[ic] := info;
    fnm^[ic] := fnorm2;
    Writeln('    INITIAL L2 NORM OF THE RESIDUALS: ',fnorm1);
    Writeln('    FINAL L2 NORM OF THE RESIDUALS..: ',fnorm2);
    Writeln('    NUMBER OF FUNCTION EVALUATION...: ',nfev);
    Writeln('    NUMBER OF JACOBIAN EVALUATIONS..: ',njev);
    Writeln('    EXIT PARAMETER..................: ',info);
    Writeln('    FINAL APPROXIMATE SOLUTION......: ');
    For i:=1 to n do write(' ',x^[i]:12:6); writeln;
    factor := ten*factor
  end;
  Writeln;
  Writeln('  SUMMARY OF ',ic,' CALL(S) TO LMDIF1:');
  Writeln('  NPROB   N    M   NFEV  NJEV  INFO   FINAL L2 NORM');
  For i:=1 to ic do
    Writeln('    ',np^[i],'     ',na^[i],'    ',ma^[i],'    ',nf^[i]:2,'    ',nj^[i]:2,'    ',nx^[i],' ',fnm^[i]);
  ReadKey;

  Dispose(iwa); Dispose(ma); Dispose(na); Dispose(nf); Dispose(nj); Dispose(np); Dispose(nx);
  Dispose(fnm); Dispose(fvec); Dispose(x);

  DoneWinCrt

END.

{end of file tlm.pas}