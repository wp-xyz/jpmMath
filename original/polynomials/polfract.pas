{****************************************************
*   Elementary operations on Polynomial Fractions   *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988". *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
*                                                   *
*                TPW version By J-P Moreau, Paris.  *
*                       (www.jpmoreau.fr)           *
****************************************************}
UNIT POLFRACT;

INTERFACE
Uses WinCrt,WinProcs,Polynoms;

Type

{Example of polynomial fraction: x2 +4x -5 / x2 -1}
FRACTION = Record
             numer,denom: POLYNOM
           End;

Function  EnterPolFract(tx:STRING; VAR F:FRACTION): Boolean;
Procedure DisplayPolFract(VAR F:FRACTION);


IMPLEMENTATION

Function EnterPolFract;
Var i:integer;
Begin
  EnterPolFract:=TRUE;
  for i:=1  to 3 do
  begin
    writeln(tx); fillchar(F,sizeof(F),0);
    if EnterPolynom(' P(x) = ',F.numer) then
      if EnterPolynom(' Q(x) = ',F.denom) then exit;
    MessageBeep(0)
  end;
  writeln;
  EnterPolFract:=FALSE
End;

Procedure DisplayPolFract;
Var i,x1,x2:integer;
Begin
  x1:=WhereX; DisplayPolynom(F.numer); x2:=WhereX; writeln;
  if (F.denom.degree>0) or (F.denom.coeff[0].value<>1)  then
  begin
    for i:=1 to (x2-x1) do write('-'); writeln;
    DisplayPolynom(F.denom)
  end
End;


END.

{end of file polfract.pas}