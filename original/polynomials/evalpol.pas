{****************************************************
*   This program evaluates a polynomial P(x) at x   *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* EVALUATION OF A POLYNOMIAL P(X) FOR X GIVEN:      *
*                                                   *
* P(X) = x3 - 5x2 + 7x + 4                          *
*                                                   *
*   3     2                                         *
*  X - 5 X + 7 X + 4                                *
*                                                   *
*  X = 1                                            *
*  P(X) = + 7                                       *
*  X = 10                                           *
*  P(X) = + 574                                     *
*  X = -5/2                                         *
*  P(X) = - 483 / 8                                 *
*  X = 1/3                                          *
*  P(X) = + 157 / 27                                *
*  X = 0                                            *
*  P(X) = + 4                                       *
*                                                   *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  EnterPolynom(), DisplayPolynom(), ReadNumber(),  *
*  EvaluatePolynom() and WriteNumber().             *
*                                                   *
*                       TPW version by J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
****************************************************}
Program EvalPol;
Uses WinCrt, Polynoms;        {see unit polynoms.pas}

Var  P  : POLYNOM;                  {polynomial type}
     a,b: NUMBER;   {real,integer or fractional type}

BEGIN

  Writeln;
  Writeln(' EVALUATION OF A POLYNOMIAL P(X) FOR X GIVEN: ');
  writeln;
  If Not EnterPolynom(' P(X) = ', P) then exit;
  DisplayPolynom(P);
  Repeat
    writeln;
    ReadNumber(' X = ', a);
    write(' P(X) = ');
    if EvaluatePolynom(P,a,b) then
      WriteNumber(b)
    else
      writeln(' Evaluation failed.')
  Until a.value=0;

  ReadKey; DoneWinCrt

END.

{end of file evalpol.pas}