{****************************************************
*     Euclidian division of two polynomials         *
*             R(x) = P(x) / Q(x)                    *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* EUCLIDIAN DIVISION OF TWO POLYNOMIALS:            *
*                                                   *
* P(x) = 2x5 +2x3 -x2 +2x -1                        *
* Q(x) = x3 -1                                      *
*                                                   *
*                                                   *
* Quotient:                                         *
*                                                   *
*    2                                              *
* 2 X + 2                                           *
*                                                   *
* Remainder:                                        *
*                                                   *
*   2                                               *
*  X + 2 X + 1                                      *
*                                                   *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  AddNumber(), EnterPolynom(), DivNumber(),        *
*  DisplayPolynom(), MultNumber() and SetNumber().  *
*                                                   *
*                       TPW version by J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
****************************************************}
PROGRAM DIVPOL;
Uses WinCrt, Polynoms;

Var P,Q,H,R: POLYNOM;


  Function DivPolynom(P,Q:POLYNOM; VAR H,R:POLYNOM): Boolean;
  Var  i,j: INTEGER;
       u  : NUMBER;
  Begin
    DivPolynom:=FALSE;
    {The Q polynomial must be <> zero}
    if (Q.degree=0) and (Q.coeff[0].value=0) then exit;
    R:=P; H.degree:=P.degree - Q.degree;
    if H.degree<0 then
    begin
      H.degree:=0; if Not SetNumber(H.coeff[0],'0') then exit;
    end
    else
    begin
      for i:=H.degree downto 0 do
      begin
        if Not DivNumber(R.coeff[R.degree],Q.coeff[Q.degree],
                         H.coeff[i]) then exit;
        for j:=i to R.degree do
        begin
          if Not MultNumber(H.coeff[i],Q.coeff[j-i], u) then exit;
          u.p:=-u.p; u.value:=-u.value;
          if Not AddNumber(R.coeff[j],u, R.coeff[j]) then exit
        end;
        if R.degree > 0 then R.degree:=R.degree-1
      end;
      While (abs(R.coeff[R.degree].value) < SMALL) and
               (R.degree>0) do  R.degree:=R.degree-1
    end;
    DivPolynom:=TRUE
  End;


{main program}
BEGIN
  writeln;
  writeln(' EUCLIDIAN DIVISION OF TWO POLYNOMIALS:');
  writeln;
  if Not EnterPolynom(' P(X) = ', P) then exit;
  if Not EnterPolynom(' Q(X) = ', Q) then exit;
  writeln; writeln;
  if DivPolynom(P,Q,H,R) then
  begin
    writeln(' Quotient:');
    DisplayPolynom(H);
    writeln;
    writeln(' Remainder:');
    DisplayPolynom(R)
  end
  else
    Writeln(' Error in Euclidian division.');
  writeln;
  ReadKey; DoneWinCrt
END.

{end of file divpol.pas}