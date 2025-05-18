{****************************************************
*  This program calculates R(x) = P(Q(x)),          *
*  P(x), Q(x) and R(x) being polynomials.           *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* SUBSTITUTION OF TWO POLYNOMIALS:                  *
*                                                   *
* P(x) = -5x2 +3                                    *
* Q(x) = 2x3 -x +5                                  *
*                                                   *
*                                                   *
*       6      4       3     2                      *
* - 20 X + 20 X - 100 X - 5 X + 50 X - 122          *
*                                                   *
*                                                   *
* SUBSTITUTION OF TWO POLYNOMIALS:                  *
*                                                   *
* P(x) = x4                                         *
* Q(x) = 2x2 -3                                     *
*                                                   *
*                                                   *
*     8      6       4       2                      *
* 16 X - 96 X + 216 X - 216 X + 81                  *
*                                                   *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  AddNumber(), EnterPolynom(), DisplayPolynom(),   *
*  MultNumber() and SetNumber().                    *
*                                                   *
*                       TPW Version By J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
****************************************************}
PROGRAM SUBSTPOL;
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


  { R(X) = P(Q(X)) }
  Function SubstPolynom(P,Q:POLYNOM;VAR R:POLYNOM): Boolean;
  Var i: INTEGER;
      a: NUMBER;
      V: POLYNOM;
  Begin
    SubstPolynom:=FALSE;
    fillchar(R,sizeof(R),0);   {set R polynomial to zero}
    R.degree:=P.degree * Q.degree;
    if R.degree>MAXPOL then exit;   {R degree is too big}
    R.coeff[0]:=P.coeff[0]; V.degree:=0;
    if Not SetNumber(V.coeff[0],'1') then exit;
    if Not SetNumber(a,'1') then exit;
    for i:=1 to P.degree do
    begin
      if Not MultPolynom(Q,V,V) then exit;
      if Not CombiPolynom(R,V,a,P.coeff[i],R) then exit
    end;
    SubstPolynom:=TRUE
  End;


{main program}
BEGIN
  writeln;
  writeln(' SUBSTITUTION OF TWO POLYNOMIALS:');
  writeln;
  if Not EnterPolynom(' P(X) = ', P) then exit;
  if Not EnterPolynom(' Q(X) = ', Q) then exit;
  writeln;
  if SubstPolynom(P,Q,R) then
    DisplayPolynom(R)
  else
    Writeln(' Error in polynomial substitution.');
  writeln;
  ReadKey; DoneWinCrt
END.

{end of file substpol.pas}