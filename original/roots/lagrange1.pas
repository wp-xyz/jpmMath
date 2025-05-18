{***************************************************
*    Root of a non-linear Equation By Lagrange     *
* ------------------------------------------------ *
* REF.: "Méthode de calcul numérique - Programmes  *
*        en Basic By Claude Nowakowski, PS1 1981". *
*                                                  *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
*  Find a root of equation: X - E SIN X            *
*                                                  *
*    Input Starting, Ending X: 1 10                *
*    Input max. error: 1E-5                        *
*                                                  *
*    Root = 2.19912286279653E+000                  *
*    Iterations: 7                                 *
*    Y = -5.40915773761579E-007                    *
*                                                  *
*                                                  *
*             Pascal Version By J-P Moreau, Paris. *
*                      (www.jpmoreau.fr)           *
***************************************************}
Program Lagrange;

Label 100, 200;

Const MAXITER = 50;

Var
    E,XB,X1,X2,YB,Y1,Y2: Double;
    K: Integer;

{User-defined non-linear function}
FUNCTION FUNC (X:Double): Double;
Begin
  FUNC := X - 2.718281828 * SIN(X)
End;


BEGIN

  Writeln;
  Write(' Input Starting, Ending X: '); Read( X1, X2);
  Write(' Input max.error: '); Read(E);
  Writeln;

  Y1 := FUNC(X1);
  Y2 := FUNC(X2);

  IF Y1 * Y2 > 0 THEN
  begin
    Writeln(' Bad choice for X1, X2!');
    goto 200
  end;

  IF ABS(Y1 * Y2) < 1E-10 THEN
  begin
    Writeln(' X1 or X2 is a root!');
    goto 200
  end;

  {main iteration loop}
  FOR K := 1 TO MAXITER DO
  begin
    XB := X2 - (X2 - X1) * Y2 / (Y2 - Y1);
    IF ABS(X1 - XB) < E THEN GOTO 100;  {normal exit}
    YB := FUNC(XB);
    IF Y1 * YB < 0 THEN
    begin
      X2 := XB; Y2 := YB
    end
    ELSE
    begin
      X1 := XB; Y1 := YB
    end
  end;

  Writeln(' No root found!');
  goto 200;

100: Writeln(' Root = ', XB);
  Writeln(' Iterations: ', K);
  Writeln(' Y = ', FUNC(XB));
  Writeln;

200: Readln
END.

{End of file Lagrange1.pas}