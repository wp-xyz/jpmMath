{****************************************************
*  This program multiplies a polynomial P(x) by a   *
*  polynomial Q(x)                                  *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* MULTIPLY TWO POLYNOMIALS:                         *
*                                                   *
* P(X) = x3 - 6x + 7                                *
* Q(x) = 5x5 -3x4 +x2 -3                            *
*                                                   *
*    8     7      6      5      4     3      2      *
* 5 X - 3 X - 30 X + 54 X - 21 X - 9 X  + 7 X       *
*                                                   *
* + 18 X - 21                                       *
*                                                   *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  AddNumber(), EnterPolynom(), DisplayPolynom()    *
*  and MultNumber().                                *
*                                                   *
*                       TPW version By J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
****************************************************}
PROGRAM MULTPOL;
Uses WinCrt, Polynoms;

Var P,Q,R: POLYNOM;


  { P(X) * Q(X) = R(X) }
  Function MultPolynom(P,Q:POLYNOM;VAR R:POLYNOM): Boolean;
  Var i,j, n: INTEGER;
      u     : NUMBER;
  Begin
    MultPolynom:=TRUE;
    fillchar(R,sizeof(R),0);  {set R polynomial to zero}
    {verify that P and Q are not void}
    if (P.degree=0) and (P.coeff[0].value=0) then exit;
    if (Q.degree=0) and (Q.coeff[0].value=0) then exit;
    MultPolynom:=FALSE;
    R.degree:=P.degree+Q.degree;
    if R.degree>MAXPOL then exit;  {R degree is too big}
    for n:=0 to R.degree do
    begin
      if Not SetNumber(R.coeff[n],'0') then exit;
      for i:=0 to P.degree do
      begin
        j:=n-i;
        if (j>=0) and (j<=Q.degree) then
        begin
          if Not MultNumber(P.coeff[i],Q.coeff[j],u) then exit;
          if Not AddNumber(R.coeff[n],u,R.coeff[n]) then exit
        end
      end
    end;
    MultPolynom:=TRUE
  End;


{main program}
BEGIN
  writeln;
  writeln(' MULTIPLY TWO POLYNOMIALS:');
  writeln;
  if Not EnterPolynom(' P(X) = ', P) then exit;
  if Not EnterPolynom(' Q(X) = ', Q) then exit;
  writeln;
  if MultPolynom(P,Q,R) then
    DisplayPolynom(R)
  else
    Writeln(' Error in multiplication.');
  writeln;
  ReadKey; DoneWinCrt
END.

{end of file multpol.pas}