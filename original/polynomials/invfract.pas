{------------------------------------------------------
  This utility Function allows to divide two polynomial
  fractions, by using multiply by inverted.
-------------------------------------------------------}
PROGRAM INVFRACT;
Uses WinCrt,Polynoms,Polfract;

VAR  F: FRACTION;


Function InvPolFract(VAR F:FRACTION): Boolean;
VAR  P: POLYNOM;
Begin
  InvPolFract:=FALSE;
  if (F.denom.degree=0) and (F.denom.coeff[0].value=0) then exit;
  P:=F.denom; F.denom:=F.numer; F.numer:=P;
  InvPolFract:=TRUE
End;


{main program to test inversion}
BEGIN

  Writeln;
  Writeln(' INVERSION OF A POLYNOMIAL FRACTION:');
  Writeln;
  EnterPolFract(' Enter polynomial fraction P(x)/Q(x):',F);
  writeln;

  InvPolFract(F);

  {display 1/F}
  DisplayPolFract(F);
  Writeln;
  Readkey;
  DoneWinCrt

END.

{end of file invfract.pas}