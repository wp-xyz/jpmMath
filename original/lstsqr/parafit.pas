{******************************************************
*    Program to demonstrate the Parafit subroutine    *
* --------------------------------------------------- *
*  Reference: BASIC Scientific Subroutines, Vol. II   *
*  By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                     *
*               Pascal version by J-P Moreau, Paris   *
*                        (www.jpmoreau.fr)            *
* --------------------------------------------------- *
* SAMPLE RUN:                                         *
*                                                     *
* Parametric least squares fit                        *
*                                                     *
* The input data are:                                 *
*                                                     *
* X( 1) =  1.000000   Y( 1) =  0.033702               *
* X( 2) =  2.000000   Y( 2) =  0.249029               *
* X( 3) =  3.000000   Y( 3) =  0.944733               *
* X( 4) =  4.000000   Y( 4) =  1.840089               *
* X( 5) =  5.000000   Y( 5) =  1.840089               *
* X( 6) =  6.000000   Y( 6) =  0.944733               *
* X( 7) =  7.000000   Y( 7) =  0.249029               *
* X( 8) =  8.000000   Y( 8) =  0.033702               *
* X( 9) =  9.000000   Y( 9) =  0.002342               *
* X(10) = 10.000000   Y(10) =  0.000084               *
*                                                     *
* The coefficients are:                               *
*  2.000000                                           *
*  4.500000                                           *
*  3.000000                                           *
*                                                     *
* The standard deviation of the fit is 0.000000       *
*                                                     *
* The number of iterations was 51                     *
******************************************************}
PROGRAM PARAFIT;
Uses WinCrt; 

CONST   SIZE = 25;

VAR
        A  : Array[1..SIZE] of DOUBLE;
        E1 : Array[1..SIZE] of DOUBLE;
        X  : Array[1..SIZE] of DOUBLE;
        Y  : Array[1..SIZE] of DOUBLE;

        a0,d,e,ee1,l1,l2,m0,m1,xx,yy : DOUBLE;
        i,l,m,n : integer;


  PROCEDURE S500;  {Function subroutine}
  begin
    yy := A[1] * EXP(-(xx - A[2]) * (xx - A[2]) / A[3])
  end;

  {Residual generation subroutine}
  PROCEDURE S200;
  var j : integer;
  Begin
    l2 := 0;
    FOR j := 1 TO n DO
    begin
      xx := X[j];
      {Obtain function}
      S500;
      l2 := l2 + (Y[j] - yy) * (Y[j] - yy)
    end;
    d := SQRT(l2 / (n - l))
  End;


{**************************************************************
*      Parametric least squares curve fit subroutine          *
* ----------------------------------------------------------- *
* This program least squares fits a function to a set of data *
* values by successively reducing the variance.  Convergence  *
* depends on the initial values and is not assured.           *
* n pairs of data values, X(i), Y(i), are given. There are l  *
* parameters, A(j), to be optimized across.                   *
* Required are initial values for the A(l) and e. Another     *
* important parameter which affects stability is e1, which is *
* initially converted to e1(l)for the first intervals.        *
* The parameters are multiplied by (1 - e1(i)) on each pass.  *
**************************************************************}
PROCEDURE Param_LS;
Label 50,100,fin;
Var i : integer;
Begin
  FOR i := 1 TO l DO e1[i] := ee1;
  m := 0;
  {Set up test residual}
  l1 := 1e6;
  {Make sweep through all parameters}
50:
  FOR i := 1 TO l DO
  begin
    a0 := A[i];
    {Get value of residual}
    A[i] := a0;
100: S200;
    {Store result in m0}
    m0 := l2;
    {Repeat for m1}
    A[i] := a0 * (1.0 - e1[i]);
    S200;
    m1 := l2;
    {Change interval size if called for
     If variance was increased, halve E1(i) }
    IF m1 > m0 THEN e1[i] := -e1[i] / 2.0;
    {If variance was reduced, increase step size by increasing E1(i)}
    IF m1 < m0 THEN e1[i] := 1.2 * e1[i];
    {If variance was increased, try to reduce it}
    IF m1 > m0 THEN A[i] := a0;
    IF m1 > m0 THEN GOTO 100
  end; {i loop}
  {End of a complete pass
   Test for convergence }
  m := m + 1;
  IF l2 < 0 THEN goto fin;
  IF ABS((l1 - l2) / l2) < e THEN goto fin;
  {If this point is reached, another pass is called for}
  l1 := l2;
  GOTO 50;

fin: End;  {Param_LS}


{main program}
BEGIN
  n := 10; l := 3;
  Writeln(' Parametric least squares fit');
  Writeln;
  Writeln(' The input data are:');
  Writeln;
  FOR i := 1 TO n DO
  begin
    X[i] := i;
    Y[i] := 2.0 * EXP(-(X[i] - 4.5) * (X[i] - 4.5) / 3.0);
    Write(' X(',i:2,') = ',X[i]:9:6);
    Writeln('   Y(',i:2,') = ',Y[i]:9:6);
  end;
  Writeln;
  e := 0.1; ee1 := 0.5; A[1] := 10.0; A[2] := 10.0; A[3] := 10.0;

  Param_LS;

  Writeln(' The coefficients are:');
  Writeln('  ',A[1]:9:6);
  Writeln('  ',A[2]:9:6);
  Writeln('  ',A[3]:9:6);
  Writeln;
  Write(' The standard deviation of the fit is ',d:11:8);
  Writeln;
  Writeln(' The number of iterations was ',m);
  Writeln;
  readkey; DoneWinCrt
END.

{End of file parafit.pas}