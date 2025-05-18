{****************************************************
*  This program calculates R(x) = a*P(x) + b*Q(x),  *
*  P(x), Q(x) and R(x) being polynomials.           *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* LINEAR COMBINATION OF TWO POLYNOMIALS:            *
*                                                   *
* P(X) = x3 +5/4x2 -8                               *
* Q(x) = 2x2 -1/7                                   *
*                                                   *
* a = 1/2                                           *
* b = -1/3                                          *
*                                                   *
*      3        2                                   *
* 1/2 X - 1/24 X - 83/21                            *
*                                                   *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  AddNumber(), EnterPolynom(), DisplayPolynom(),   *
*  MultNumber() and SetNumber().                    *
*                                                   *
*                       TPW version by J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
****************************************************}
PROGRAM MULTPOL;
Uses WinCrt, Polynoms;

Var P,Q,R: POLYNOM;
    a, b : NUMBER;


  { a P(X) + b Q(X) = R(X) }
  Function CombiPolynom(P,Q:POLYNOM; a,b:NUMBER; VAR R:POLYNOM): Boolean;
  Var i  : INTEGER;
      u,v: NUMBER;
  Begin
    CombiPolynom:=FALSE;
    if Q.degree > P.degree then R.degree:=Q.degree
                           else R.degree:=P.degree;
    if R.degree > MAXPOL then exit;   {degree of R too big}
    for i:=0 to R.degree do
    begin
      if Not SetNumber(u,'0') then exit;
      if Not SetNumber(v,'0') then exit;
      if i<=P.degree then
        if Not MultNumber(a,P.coeff[i], u) then exit;
      if i<=Q.degree then
        if Not MultNumber(b,Q.coeff[i], v) then exit;
      if Not AddNumber(u,v,R.coeff[i]) then exit
    end;
    While (R.degree>0) and (abs(R.coeff[R.degree].value)<SMALL) do
      R.degree:=Pred(R.degree);
    CombiPolynom:=TRUE
  End;


{main program}
BEGIN
  writeln;
  writeln(' LINEAR COMBINATION OF TWO POLYNOMIALS:');
  writeln;
  if Not EnterPolynom(' P(X) = ', P) then exit;
  if Not EnterPolynom(' Q(X) = ', Q) then exit;
  writeln;
  ReadNumber(' a = ', a);
  ReadNumber(' b = ', b);
  if CombiPolynom(P,Q,a,b,R) then
    DisplayPolynom(R)
  else
    Writeln(' Error in linear combination.');
  writeln;
  ReadKey; DoneWinCrt
END.

{end of file combipol.pas}