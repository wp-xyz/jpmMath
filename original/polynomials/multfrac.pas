{****************************************************
*        Multiply two polynomial fractions          *
*                                                   *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUNS:                                      *
*                                                   *
* MULTIPLY TWO POLYNOMIAL FRACTIONS:                *
*                                                   *
* Enter polynomial fraction F1(x):                  *
*                                                   *
* P(x) = x -1                                       *
* Q(x) = x2 +5x -7                                  *
*                                                   *
* Enter polynomial fraction F2(x):                  *
*                                                   *
* P(x) = x +3                                       *
* Q(x) = x2 -1                                      *
*                                                   *
*                                                   *
*  X + 3                                            *
* -----------                                       *
*   3    2                                          * 
*  X + 6X - 2X -7                                   *
*                                                   *
* ------------------------------------------------- *
* Functions used (of unit Polynoms):                *
*                                                   *
*  AddNumber(), EnterPolynom(), DivNumber(),        *
*  DisplayPolynom(), MultNumber() and SetNumber().  *
*                                                   *
*  Note: stack size = 45000                         *
*                                                   *
*                       TPW version By J-P Moreau.  *
*                           (www.jpmoreau.fr)       *
****************************************************}
PROGRAM MULTFRACT;
Uses WinCrt,Polynoms,Polfract;

VAR  F,F1,F2: FRACTION;


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

Function MultPolFract(F1,F2:FRACTION;VAR F:FRACTION): Boolean;
Begin
  MultPolFract:=FALSE;
  if Not SimpPolFract(F1,FALSE) then exit;
  if Not SimpPolFract(F2,FALSE) then exit;
  F.numer:=F1.numer;
  F1.numer:=F2.numer;
  F2.numer:=F.numer;
  if Not SimpPolFract(F1,FALSE) then exit;
  if Not SimpPolFract(F2,FALSE) then exit;
  if Not MultPolynom(F1.numer,F2.numer,F.numer) then exit;
  if Not MultPolynom(F1.denom,F2.denom,F.denom) then exit;
  MultPolFract:=TRUE
End;


{main program}
BEGIN

  Writeln;
  Writeln(' MULTIPLICATION OF TWO POLYNOMIAL FRACTIONS:');
  Writeln;
  EnterPolFract(' Enter polynomial fraction F1(x)):',F1);
  writeln;
  EnterPolFract(' Enter polynomial fraction F2(x)):',F2);
  writeln;

  if MultPolFract(F1,F2,F) then
    {display result F}
    DisplayPolFract(F)
  else
    writeln(' Error in multiplication.');
  Writeln;
  Readkey;
  DoneWinCrt

END.

{end of file multfrac.pas}