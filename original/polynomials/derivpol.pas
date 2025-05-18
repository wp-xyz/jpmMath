{****************************************************
*       Nth Derivative of a polynomial P(x)         *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* Nth DERIVATIVE OF A POLYNOMIAL:                   *
*                                                   *
* P(X) = 4x5 -3/5x3 + x2 -7/2x +11                  *
* Order of derivative: 1                            *
*                                                   *
*     4        2                                    *
* 20 X  - 9/5 X  + 2 X - 7/2                        *
*                                                   *
*                                                   *
* Nth DERIVATIVE OF A POLYNOMIAL:                   *
*                                                   *
* P(X) = 4x5 -3/5x3 + x2 -7/2x +11                  *
* Order of derivative: 3                            *
*                                                   *
*      2                                            *
* 240 X  - 18/5                                     *
*                                                   *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  EnterPolynom(), DisplayPolynom(), MultNumber(),  *
*  SetNumber().                                     *
*                                                   *
*                       TPW version by J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
****************************************************}
Program DerivPol;
Uses WinCrt, Polynoms;        {see unit polynoms.pas}

Var  P,R: POLYNOM;                  {polynomial type}
     n  : INTEGER;              {order of derivative}


Function DerivPolynom(P:POLYNOM; n:INTEGER; VAR R:POLYNOM): Boolean;
Var  i,j: INTEGER;
     u  : NUMBER;
     chv: STRING;
Begin
  DerivPolynom:=FALSE; R.degree:=P.degree-n;
  if R.degree<0 then
  begin
    R.degree:=0; if Not SetNumber(R.coeff[0], '0') then exit
  end
  else
    for i:=0 to R.degree do
    begin
      R.coeff[i]:=P.coeff[i+n];
      for j:=1 to n do
      begin
        Str(i+j,chv); if Not SetNumber(u,chv) then exit;
        if Not MultNumber(R.coeff[i], u, R.coeff[i]) then exit
      end
    end;       
  DerivPolynom:=TRUE
End;      

BEGIN

  Writeln;
  Writeln(' Nth DERIVATIVE OF A POLYNOMIAL: ');
  writeln;
  If Not EnterPolynom(' P(X) = ', P) then exit;
  write(' Order of derivative: '); readln(n);

  if DerivPolynom(P,n,R) then
    Displaypolynom(R)
  else
    writeln(' Error in calculating derivative.');
  writeln;
  ReadKey; DoneWinCrt

END.

{end of file derivpol.pas}