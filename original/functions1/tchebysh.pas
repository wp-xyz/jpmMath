{**************************************************
* Chebyshev Approximation of a user defined real  *
* function FUNC(X) in double precision.           *
* ----------------------------------------------- *
* SAMPLE RUN:                                     *
* (Approximate sin(x) from x=0 to x=PI).          *
*                                                 *
* Chebyshev coefficients (N=10):                  *
*  9.44002431536469E-0001                         *
*  4.98583921544615E-0017                         *
* -4.99403258270407E-0001                         *
* -5.16910326391409E-0017                         *
*  2.79920796175457E-0002                         *
*  4.35374934888710E-0019                         *
* -5.96695195800671E-0004                         *
*  2.25974837800291E-0017                         *
*  6.70417552415733E-0006                         *
*  3.01353146809399E-0018                         *
*      X          Chebyshev Eval.     SIN(X)      *
*  -------------------------------------------    *
*  0.00000000       0.00000005      0.00000000    *
*  0.34906585       0.34202018      0.34202014    *
*  0.69813170       0.64278757      0.64278761    *
*  1.04719755       0.86602545      0.86602540    *
*  1.39626340       0.98480773      0.98480775    *
*  1.74532925       0.98480773      0.98480775    *
*  2.09439510       0.86602545      0.86602540    *
*  2.44346095       0.64278757      0.64278761    *
*  2.79252680       0.34202018      0.34202014    *
*  3.14159265       0.00000005      0.00000000    *
*                                                 *
*               TPW Release By J-P Moreau, Paris. *
*                       (www.jpmoreau.fr)         *
**************************************************}
PROGRAM TESTCHEBY;

Uses WinCrt, Chebyshe;

Var
    X0,X1: Double; COEFF:Tab;
    DX, X: Double;
    I, N : Integer;

BEGIN

  N:=10;
  X0:=ZERO; X1:=PI;

  CHEBFT(X0,X1,COEFF,N);

  writeln(' Chebyshev coefficients (N=',N,'):');
  for I:=1 to N do
    writeln(' ',COEFF[I]);

  DX:=(X1-X0)/(N-1);
  X:=X0-DX;

  writeln('      X        Chebyshev Eval.     SIN(X)    ');
  writeln(' --------------------------------------------');
  for I:=1 to N do
  begin
    X:=X+DX;
    writeln(X:11:8, CHEBEV(X0,X1,COEFF,N,X):16:8, FUNC(X):16:8)
  end;

  ReadKey;
  DoneWinCrt

END.

{end of file Tchebysh.pas}