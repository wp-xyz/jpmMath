{****************************************************
* Division of two polynomials by increasing powers  *
*                               k                   *
*           P(x) = Q(x).H(x) + x.R(x)               *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* DIVISION OF TWO POLYNOMIALS BY INCREASING POWERS: *
*                                                   *
* P(x) = 2x2 +x                                     *
* Q(x) = x3 -x2 +1                                  *
* Order: 4                                          *
*                                                   *
*                                                   *
* Quotient:                                         *
*                                                   *
*  3    2                                           *
* X + 2X + X                                        *
*                                                   *
* Remainder:                                        *
*                                                   *
*    2                                              *
* - X - X +1                                        *
*                                                   *
*                                                   *
* DIVISION OF TWO POLYNOMIALS BY INCREASING POWERS: *
*                                                   *
* P(x) = 3/5x3 +2x -4                               *
* Q(x) = x2 -5/4                                    *
* Order: 9                                          *
*                                                   *
*                                                   *
* Quotient:                                         *
*                                                   *
*            8           7            6           5 *
* 4096/3125 X - 704/625 X + 1024/625 X - 176/125 X  *
*            4         3         2                  *
* + 256/125 X - 44/25 X + 64/25 X - 8/5 X + 16/5    *     
*                                                   *
* Remainder:                                        *
*                                                   *
*                                                   *
* - 4096/3125 X + 704/625                           *
*                                                   *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  AddNumber(), EnterPolynom(), DivNumber(),        *
*  DisplayPolynom(), MultNumber() and SetNumber().  *
*                                                   *
*                       TPW version by J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
*****************************************************
Explanations:
Given two polynomials P(x) and Q(x) and an integer k,
it exists a couple of polynomials H(x) and R(x), such
as:        Px) = Q(x).H(x) + x^k.R(x)
The degree of H(x) is < k.
Note: the constant coefficient of Q(x) must be <> 0.
----------------------------------------------------}
PROGRAM DIVPOL1;
Uses WinCrt, Polynoms;

Var P,Q,H,R: POLYNOM;
    k: INTEGER;

  {Division of two polynomials by increasing powers}
  Function DivPolynom1(P,Q:POLYNOM; k:INTEGER;VAR H,R:POLYNOM): Boolean;
  Var  i,j,mx: INTEGER;
       u     : NUMBER;
  Begin
    DivPolynom1:=FALSE;
    {The Q polynomial constant term and k must be <> zero}
    if (k<=0) or (Q.coeff[0].value=0) then exit;
    mx:=Q.degree+k-1;
    if mx<P.degree then mx:=P.degree;
    if mx>Q.degree then
      for i:=Q.degree+1 to mx do
        if Not SetNumber(Q.coeff[i],'0') then exit;
    Fillchar(R,sizeof(R), 0);
    for i:=0 to P.degree do R.coeff[i]:=P.coeff[i];
    if k>Q.degree then
      for i:=Q.degree+1 to k do
        if Not SetNumber(Q.coeff[i],'0') then exit;
    for i:=0 to Pred(k) do
    begin
      if Not DivNumber(R.coeff[i],Q.coeff[0],H.coeff[i]) then exit;
      for j:=1 to mx do
      begin
        if Not MultNumber(H.coeff[i],Q.coeff[j-i], u) then exit;
        u.p:=-u.p; u.value:=-u.value;
        if Not AddNumber(R.coeff[j], u, R.coeff[j]) then exit
      end
    end;
    H.degree:=k-1;
    While (H.degree>0) and (H.coeff[H.degree].value=0) do Dec(H.degree);
    R.degree:=mx-k;
    for i:=0 to R.degree do R.coeff[i]:=R.coeff[i+k];
    While (R.degree>0) and (R.coeff[R.degree].value=0) do Dec(R.degree);
    DivPolynom1:=TRUE
  End;


{main program}
BEGIN
  writeln;
  writeln(' DIVISION OF TWO POLYNOMIALS BY INCREASING POWERS:');
  writeln;
  if Not EnterPolynom(' P(X) = ', P) then exit;
  if Not EnterPolynom(' Q(X) = ', Q) then exit;
  Write(' Order: '); readln(k);
  writeln; writeln;
  if DivPolynom1(P,Q,k,H,R) then
  begin
    writeln(' Quotient:');
    DisplayPolynom(H);
    writeln;
    writeln(' Remainder:');
    DisplayPolynom(R)
  end
  else
    Writeln(' Error in Division by increasing powers.');
  writeln;
  ReadKey; DoneWinCrt
END.

{end of file divpol1.pas}