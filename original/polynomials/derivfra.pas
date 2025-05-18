{****************************************************
* Nth derivation of a polynomial fraction P(x)/Q(x) *
*                                                   *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUNS:                                      *
*                                                   *
* Nth DERIVATION OF A POLYNOMIAL FRACTION:          *
*                                                   *
* Enter polynomial fraction F1(x):                  *
* P(x) = 3x2 +2x +1                                 *
* Q(x) = 4x2 +5x +56                                *
*                                                   *
* Order: 1                                          *
*                                                   *
*    2                                              *
* 7 X + 328 X + 107                                 *
* -----------------                                 *
*     4      3       2                              * 
* 16 X + 40 X + 473 X + 560 X + 3136                *
*                                                   *
*                                                   *
* Nth DERIVATION OF A POLYNOMIAL FRACTION:          *
*                                                   *
* Enter polynomial fraction F1(x):                  *
* P(x) = x -1                                       *
* Q(x) = x +1                                       *
*                                                   *
* Order: 5                                          *
*                                                   *
*                                                   *
*  240                                              *
* -----------------                                 *
*   6     5      4      3      2                    * 
*  X + 6 X + 15 X + 20 X + 15 X + 6 X + 1           *
*                                                   *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  AddNumber(), EnterPolynom(), DivNumber(),        *
*  DisplayPolynom(), MultNumber() and SetNumber().  *
*                                                   *
*  Note: stack size = 45000                         *
*                                                   *
*                       TPW version by J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
*****************************************************
Explanations:
------------
This program calculates the nth derivative of the polynomial
fraction F1=P0/Q0. The first derivative is:

  F1' = (P0'Q0-P0Q0')/Q0^2 = P1/Q1

The second derivative is:

  F1" = (F1')' = (P1'Q1-P1Q1')/Q1^2 = P2/Q2

If F1(n) = Pn/Qn, then we can obtain F1(n+1) by the formula:

  F1(n+1) = (Pn'Qn-PnQn')/Qn^2 = Pn+1/Qn+1

Note: for high values of derivation order n, errors may occur:
for instance, let us seek the nth derivative of x-2/x-3.
Up to n=6, the results are correct. For n>=7, the results are
not correct (the simplification is not finished due to numerical
errors caused by successive calls of routine SimpPolFract().
--------------------------------------------------------------}
PROGRAM DERIVFRACT;
Uses WinCrt,Polynoms,Polfract;

VAR  F1,F: FRACTION;
     n: INTEGER;


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


  {Euclidian division of two polynomials}
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

  {GCD of two polynomials}
  Function GCDPolynom(P,Q:POLYNOM; VAR R:POLYNOM): Boolean;
  Var  i: INTEGER;
       V: POLYNOM;
       z: NUMBER;
       pg,pp:LONGINT;
       bb: Boolean;
       rr: REAL;
  Begin
    GCDPolynom:=FALSE;
    if (P.degree=0) and (P.coeff[P.degree].value=0) then exit;
    if (Q.degree=0) and (Q.coeff[Q.degree].value=0) then exit;
    if Q.degree>P.degree then begin R:=P; P:=Q; Q:=R end;
    R:=Q;
    While R.degree>0 do
    begin
      if Not DivPolynom(P,Q,V,R) then exit;
      While (R.degree>0) and (abs(R.coeff[R.degree].value)<SMALL) do Dec(R.degree);
      if R.degree=0 then
      begin
        if abs(R.coeff[0].value)>SMALL then Q.degree:=0
      end
      else
      begin
        P:=Q; Q:=R
      end;
      for i:=0 to Q.degree-1 do
        if Not DivNumber(Q.coeff[i],Q.coeff[Q.degree],Q.coeff[i]) then exit;
      if Not SetNumber(Q.coeff[Q.degree], '1') then exit
    end;
    R:=Q;
    if Q.degree=0 then
      bb:=SetNumber(R.coeff[0], '1')
    else
    begin
      pg:=0;
      for i:=0 to R.degree do
        if Not R.coeff[i].is_real then
          if pg=0 then pg:=R.coeff[i].p
                  else begin rr:=GCD(pg,R.coeff[i].p); pg:=Round(rr) end;
      pp:=0;
      for i:=0 to R.degree do
        if Not R.coeff[i].is_real then
          if pp=0 then pp:=R.coeff[i].q
                  else
                  begin
                    rr:=GCD(pp,R.coeff[i].q);
                    pp:=pp*R.coeff[i].q div Round(rr)
                  end;
      if pg<>0 then
      begin
        z.is_real:=FALSE; z.p:=pp;
        z.q:=pg; z.value:=pp/pg;
        for i:=0 to R.degree do
          if Not MultNumber(R.coeff[i], z, R.coeff[i]) then exit
      end
    end;
    GCDPolynom:=TRUE
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

{Simplify a polynomial fraction F(x) = P(x) / Q(x) }
Function SimpPolFract(VAR F:FRACTION;int:BOOLEAN): BOOLEAN;
Var P,R: POLYNOM;
    z:NUMBER;
    pg,pp1,pp2:REAL;
    i:INTEGER;
Begin
  SimpPolFract:=FALSE;
  if Not GCDPolynom(F.numer,F.denom,P) then exit;
  if Not DivPolynom(F.numer,P,F.numer,R) then exit;
  if Not DivPolynom(F.denom,P,F.denom,R) then exit;
  if int then   {integer coefficients asked}
  begin
    pg:=0;
    for i:=0 to F.numer.degree do
      if Not F.numer.coeff[i].is_real then
        if pg=0 then  pg:=F.numer.coeff[i].p
                else pg:=GCD(pg,F.numer.coeff[i].p);
    for i:=0 to F.denom.degree do
      if Not F.denom.coeff[i].is_real then
        if pg=0 then  pg:=F.denom.coeff[i].p
                else pg:=GCD(pg,F.denom.coeff[i].p);
    pp1:=0;
    for i:=0 to F.numer.degree do
      if Not F.numer.coeff[i].is_real then
        if pp1=0 then pp1:=F.numer.coeff[i].q
                 else pp1:=pp1*F.numer.coeff[i].q/GCD(pp1,F.numer.coeff[i].q);
    pp2:=0;
    for i:=0 to F.denom.degree do
      if Not F.denom.coeff[i].is_real then
        if pp2=0 then pp2:=F.denom.coeff[i].q
                 else pp2:=pp2*F.denom.coeff[i].q/GCD(pp2,F.denom.coeff[i].q);
    z.p:=Round(pp1*pp2/GCD(pp1,pp2)); z.q:=Round(pg);
    if z.q<>0 then
    begin
      z.value:=z.p/z.q; z.is_real:=FALSE;
      for i:=0 to F.numer.degree do
        if Not MultNumber(F.numer.coeff[i],z,F.numer.coeff[i]) then exit;
      for i:=0 to F.denom.degree do
        if Not MultNumber(F.denom.coeff[i],z,F.denom.coeff[i]) then exit
    end
  end;
  SimpPolFract:=TRUE;
End;

Function DerivPolFract(F:FRACTION;n:INTEGER;VAR G:FRACTION): Boolean;
Var P,Q:POLYNOM;i:integer;a,b:NUMBER;
Begin
  DerivPolFract:=FALSE;
  if Not SetNumber(a,'1') then exit;   {a=1}
  if Not SetNumber(b,'-1') then exit;  {b=-1}
  for i:=1 to n do
  begin
    if Not MultPolynom(F.denom,F.denom,G.denom) then exit;
    if Not DerivPolynom(F.numer,1,P) then exit;
    if Not DerivPolynom(F.denom,1,Q) then exit;
    if Not MultPolynom(F.denom,P,P) then exit;
    if Not MultPolynom(F.numer,Q,Q) then exit;
    if Not CombiPolynom(P,Q,a,b,G.numer) then exit;
    G.denom:=F.denom;
    if Not SimpPolFract(G,FALSE) then exit;
    if Not MultPolynom(G.denom,F.denom,G.denom) then exit;
    if Not SimpPolFract(G,FALSE) then exit;
    F:=G
  end;
  DerivPolFract:=TRUE
End;


{main program}
BEGIN

  Writeln;
  Writeln(' Nth DERIVATION OF A POLYNOMIAL FRACTION:');
  Writeln;
  EnterPolFract(' Enter polynomial fraction F1(x)):',F1);
  writeln;
  write(' Order: '); readln(n);

  if DerivPolFract(F1,n,F) then
    {display result F}
    DisplayPolFract(F)
  else
    writeln(' Error in derivation.');
  Writeln;
  Readkey;
  DoneWinCrt

END.

{end of file derivfra.pas}