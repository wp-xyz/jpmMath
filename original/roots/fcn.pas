{****************************************************
*      Examples of non-linear system for UNNES      *
*                                                   *
* Release 1.1 (03/29/2007):  added example #2.      *
****************************************************}
UNIT FCN;

Interface

Uses Wincrt, Utils;

      Procedure fcn1(Var overfl:boolean; n:integer; fvec:pVec; xc:pVec);
      Procedure jacob(Var overfl:boolean; n:integer; jac:pMat; xc:pVec);

Implementation


Procedure fcn1(Var overfl:boolean; n:integer; fvec:pVec; xc:pVec);
Begin
  overfl := false;
  nfetot := nfetot + 1;

{ Nonlinear systems to solve
  Example #1:
    x1^2 + x1 + x2^2 - 2          = 0
    x1^2 + x2 - x2^2 - 1 + Ln(x1) = 0  

  Example #2:
     10*x1^2 +       + x3    + x4 - 20 + sin(x1)^2 + cos(x1)^2 = 0
     x1      + 20*x2 + x3    + x4 - 48 + 1/x1^6                = 0
     (x1+x2]^2       + 30*x3 + x4 - 97 + Ln(x1) + Ln(x2+x3)    = 0
     x1      + x2    + x3 + 40*x4 -166 + x1*x1                 = 0

  Example #3  (stiff system):
  Hiebert's 2nd Chemical Engineering Problem

  source: Hiebert; Sandia Technical Report #SAND80-0181
          Sandia National Laboratories, Albuquerque, NM (1980)

    X1 + X2 + X4 - .001 = 0
    X5 + X6 -55 = 0
    X1 + X2 + X3 + 2X5 + X6 - 110.001 = 0
    X1 - 0.1X2 = 0
    X1 - 10000 X3 X4 = 0
    X5 - 5.5e15 X3 X6 = 0

    solution: (8.264e-5, 8.264e-4, 9.091e-5, 9.091e-5, 55, 1.1e-10) }

  if n=2 then
  begin
    fvec^[1] := xc^[1]*xc^[1] + xc^[1] + xc^[2]*xc^[2] -2.0;
    fvec^[2] := xc^[1]*xc^[1] + xc^[2] - xc^[2]*xc^[2] -1.0 + Ln(xc^[1])
  end;
  if n=4 then
  begin
    fvec^[1] := 10.0*xc^[1] + xc^[2] + xc^[3] + xc^[4] - 20.0 + Sqr(sin(xc^[1])) + Sqr(cos(xc^[2]));
    fvec^[2] := xc^[1] + 20.0*xc^[2] + xc^[3] + xc^[4] - 48.0 + one/Power(xc^[1],6);
    fvec^[3] := Sqr(xc^[1] + xc^[2]) + 30.0*xc^[3] + xc^[4] - 97.0 + Ln(xc^[1]) + Ln(xc^[2]+xc^[3]);
    fvec^[4] := xc^[1] + xc^[2] + xc^[3] + 40.0*xc^[4] - 166.0 + Sqr(xc^[1])
  end;
  if n=6 then
  begin
    fvec^[1] := xc^[1] + xc^[2] + xc^[4] - 0.001;
    fvec^[2] := xc^[5] + xc^[6] - 55.0;
    fvec^[3] := xc^[1] + xc^[2] + xc^[3] + 2.0 * xc^[5] + xc^[6] - 110.001;
    fvec^[4] := xc^[1] - 0.1 * xc^[2];
    fvec^[5] := xc^[1] - 1E04 * xc^[3] * xc^[4];
    fvec^[6] := xc^[5] - 5.5E15 * xc^[3] * xc^[6]
  end;

End;


Procedure jacob(Var overfl:boolean; n:integer; jac:pMat; xc:pVec);
Var i,j: integer;
Begin
  overfl := false;
  For i:= 1 to n do
    For j:=1 to n do
     jac^[i,j] := 0.0

End;

END.