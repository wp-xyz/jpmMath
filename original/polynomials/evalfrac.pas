{****************************************************
*    Evaluate a polynomial fraction for x given     *
*                                                   *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUNS:                                      *
*                                                   *
* EVALUATE A POLYNOMIAL FRACTION:                   *
*                                                   *
* Enter polynomial fraction P(x)/Q(x):              *
*                                                   *
* P(x) = x3 -5x +2                                  *
* Q(x) = x2 -1                                      *
*                                                   *
* x value: -2                                       *
* F(x) = + 4/3                                      *
* x value: 1/2                                      *
* F(x) = + 1/2                                      *
* x value: 2                                        *
* F(x) = - 0/3                                      *
* x value: 1                                        *
* Evaluation failed.                                *
* x value: 0                                        *
* F(x) = - 2                                        *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  AddNumber(), EnterPolynom(), DivNumber(),        *
*  DisplayPolynom(), MultNumber() and SetNumber().  *
*                                                   *
*                       TPW version by J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
****************************************************}
PROGRAM EVALFRAC;
Uses WinCrt, Polynoms,  Polfract;

VAR  F:FRACTION;
     x,y:NUMBER;


Function EvalPolFract(VAR F:FRACTION;x:NUMBER;VAR y:NUMBER):Boolean;
Var u,v: NUMBER;
Begin
  EvalPolFract:=FALSE;
  if Not EvaluatePolynom(F.numer,x,u) then exit;
  if Not EvaluatePolynom(F.denom,x,v) then exit;
  if Not DivNumber(u,v,y) then exit;
  EvalPolFract:=TRUE
End;


{main program}
BEGIN

  Writeln;
  Writeln(' EVALUATE A POLYNOMIAL FRACTION:');
  Writeln;
  EnterPolFract(' Enter polynomial fraction P(x)/Q(x):',F);
  writeln;
  Repeat
    ReadNumber(' x value: ',x);
    if EvalPolFract(F,x,y) then
    begin
      Write(' F(x) = '); WriteNumber(y); 
      Writeln
    end
    else
      Writeln(' Evaluation failed.');
  Until x.value=0;
  Readkey;
  DoneWinCrt

END.

{end of file evalfrac.pas}