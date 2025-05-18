{**************************************************
* Chebyshev Approximation of a user defined real  *
* function FUNC(X) in double precision.           *
* ----------------------------------------------- *
* SAMPLE RUN:                                     *
* (Approximate derivative of sin(x) from x=0 to   *
*  x=PI).                                         *
*                                                 *
*      X          Chebyshev Eval.    COS(X)       *
*                 of derivative                   *
*  --------------------------------------------   *
*  0.00000000      0.99999707      1.00000000     *
*  0.15707963      0.98768900      0.98768834     *
*  0.31415927      0.95105644      0.95105652     *
*  0.47123890      0.89100611      0.89100652     *
*  0.62831853      0.80901694      0.80901699     *
*  0.78539816      0.70710708      0.70710678     *
*  0.94247780      0.58778552      0.58778525     *
*  1.09955743      0.45399047      0.45399050     *
*  1.25663706      0.30901672      0.30901699     *
*  1.41371669      0.15643421      0.15643447     *
*  1.57079633      0.00000000      0.00000000     *
*  1.72787596     -0.15643421     -0.15643447     *
*  1.88495559     -0.30901672     -0.30901699     *
*  2.04203522     -0.45399047     -0.45399050     *
*  2.19911486     -0.58778552     -0.58778525     *
*  2.35619449     -0.70710708     -0.70710678     *
*  2.51327412     -0.80901694     -0.80901699     *
*  2.67035376     -0.89100611     -0.89100652     *
*  2.82743339     -0.95105644     -0.95105652     *
*  2.98451302     -0.98768900     -0.98768834     *
*  3.14159265     -0.99999707     -1.00000000     *
*                                                 *
*               TPW Release By J-P Moreau, Paris. *
*                       (www.jpmoreau.fr)         *
**************************************************}
PROGRAM TESTCHDER;

Uses WinCrt, Chebyshe;

Var
    X0,X1: Double;
    COEFF, CDER: Tab;
    DX, X: Double;
    I, N,  NPTS: Integer;

BEGIN

  N:=10;
  X0:=ZERO; X1:=PI;

  {calculate Chebyshev coefficients of FUNC(X) }
  CHEBFT(X0,X1,COEFF,N);

  {calculate Chebyshev coefficients of derivative of FUNC(X) }
  CHDER(X0,X1,COEFF,CDER,N);

  NPTS:=21;
  DX:=(X1-X0)/(NPTS-1);
  X:=X0-DX;

  writeln('      X        Chebyshev Eval.     COS(X)    ');
  writeln('               of Derivative                 ');
  writeln(' --------------------------------------------');
  for I:=1 to NPTS do
  begin
    X:=X+DX;
    writeln(X:11:8, CHEBEV(X0,X1,CDER,N,X):16:8, COS(X):16:8)
  end;

  ReadKey;
  DoneWinCrt

END.

{end of file Tchebder.pas}