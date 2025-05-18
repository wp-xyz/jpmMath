{**************************************************
* Chebyshev Approximation of a user defined real  *
* function FUNC(X) in double precision.           *
* ----------------------------------------------- *
* SAMPLE RUN:                                     *
* (Approximate integral of sin(x) from x=0 to     *
*  x=PI).                                         *
*                                                 *
*      X          Chebyshev Eval.   -COS(X)+1     *
*                   of integral                   *
*  -------------------------------------------    *
*  0.00000000       0.00000000      0.00000000    *
*  0.15707963       0.01231166      0.01231166    *
*  0.31415927       0.04894348      0.04894348    *
*  0.47123890       0.10899348      0.10899348    *
*  0.62831853       0.19098301      0.19098301    *
*  0.78539816       0.29289321      0.29289322    *
*  0.94247780       0.41221474      0.41221475    *
*  1.09955743       0.54600950      0.54600950    *
*  1.25663706       0.69098301      0.69098301    *
*  1.41371669       0.84356554      0.84356553    *
*  1.57079633       1.00000000      1.00000000    *
*  1.72787596       1.15643446      1.15643447    *
*  1.88495559       1.30901699      1.30901699    *
*  2.04203522       1.45399050      1.45399050    *
*  2.19911486       1.58778526      1.58778525    *
*  2.35619449       1.70710679      1.70710678    *
*  2.51327412       1.80901699      1.80901699    *
*  2.67035376       1.89100652      1.89100652    *
*  2.82743339       1.95105651      1.95105652    *
*  2.98451302       1.98768834      1.98768834    *
*  3.14159265       2.00000000      2.00000000    *
*                                                 *
*               TPW Release By J-P Moreau, Paris. *
*                       (www.jpmoreau.fr)         *
**************************************************}
PROGRAM TESTCHEBY;

Uses WinCrt, Chebyshe;

Var
    X0,X1: Double;
    COEFF, CINT: Tab;
    DX, X: Double;
    I, N,  NPTS: Integer;

BEGIN

  N:=10;
  X0:=ZERO; X1:=PI;

  {calculate Chebyshev coefficients of FUNC(X) }
  CHEBFT(X0,X1,COEFF,N);

  {calculate Chebyshev coefficients of integral of FUNC(X) }
  CHINT(X0,X1,COEFF,CINT,N);

  NPTS:=2*N+1;
  DX:=(X1-X0)/(NPTS-1);
  X:=X0-DX;

  writeln('      X        Chebyshev Eval.   -COS(X)+1   ');
  writeln('                 of integral                 ');
  writeln(' --------------------------------------------');
  for I:=1 to NPTS do
  begin
    X:=X+DX;
    writeln(X:11:8, CHEBEV(X0,X1,CINT,N,X):16:8, (ONE-COS(X)):16:8)
  end;

  ReadKey;
  DoneWinCrt

END.

{end of file Tchebint.pas}