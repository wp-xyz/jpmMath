{****************************************************
*          GCD and SCM of two polynomials           *
*                                                   *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* GCD AND SCM OF TWO POLYNOMIALS:                   *
*                                                   *
* P(x) = x5 +5x4 +7x3 +5x2 +x -1                    *
* Q(x) = x4 +4x3 -7x +2                             *
*                                                   *
*                                                   *
* Greatest common divisor:                          *
*                                                   *
*  2                                                *
* X + 3X - 1                                        *
*                                                   *
*                                                   *
* Smallest common multiple:                         *
*                                                   *
*  7    6     5    4   3      2                     *
* X + 6X + 10X + 2X -8X - 10 X -3X +2               *
*                                                   *
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
 The Greatest common divisor (GCD) and the Smallest Common
 multiple (SCM) of two polynomials are defined but for a
 coefficient. That is to say, if a polynomial R(x) is the
 GCD of P(x) and Q(x), any polynomial µ*R(x) is also a GCD.
 To find the GCD, after permutating P(x) and Q(x), if the
 degree of P is lower than the degree of Q, we first perform
 the Euclidian division  P/Q: P=Q.H1+R1 (see program Divpol.pas).
 if R1<>0, we then divide Q by R1: Q=R1.H2+R2.
 if R2<>0, we continue until Rn=0. Rn-1 is then the GCD of
 P(x) and Q(x).
 To find the SCM of P(x) and Q(x), we use the relation:
    P(x).Q(x)=GCD(P,Q).SCM(P,Q).  
----------------------------------------------------------------}
PROGRAM GCDPOL;
Uses WinCrt, Polynoms;

Var  P,Q,R: POLYNOM;


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

  {SCM of two polynomials}
  Function SCMPolynom(P,Q:POLYNOM; VAR R:POLYNOM): Boolean;
  Var  U,V: POLYNOM;
  Begin
    SCMPolynom:=FALSE;
    if Not GCDPolynom(P,Q,U) then exit;
    if Not MultPolynom(P,Q,V) then exit;
    if Not DivPolynom(V,U,R,V) then exit;
    SCMPolynom:=TRUE
  End;


{main program}
BEGIN
  writeln;
  writeln(' GCD AND SCM OF TWO POLYNOMIALS:');
  writeln;
  if Not EnterPolynom(' P(X) = ', P) then exit;
  if Not EnterPolynom(' Q(X) = ', Q) then exit;
  writeln; writeln;
  if GCDPolynom(P,Q,R) then
  begin
    writeln(' Greatest common divisor:');
    DisplayPolynom(R)
  end
  else
    Writeln(' Error in looking for GCD.');
  writeln;
  writeln;
  if SCMPolynom(P,Q,R) then
  begin
    writeln(' Smallest common multiple:');
    DisplayPolynom(R)
  end
  else
    Writeln(' Error in looking for SCM.');
  writeln;
  ReadKey; DoneWinCrt
END.

{end of file gcdpol.pas}