{*****************************************************************
*      Purpose: This program computes the beta function          * 
*               B(p,q) for p > 0 and q > 0 using                 *
*               subroutine BETA                                  *
*      Input :  p  --- Parameter  ( p > 0 )                      *
*               q  --- Parameter  ( q > 0 )                      *
*      Output:  BT --- B(p,q)                                    *
*      Examples:                                                 *
*                p       q           B(p,q)                      *
*              ---------------------------------                 *
*               1.5     2.0     .2666666667D+00                  *
*               2.5     2.0     .1142857143D+00                  *
*               1.5     3.0     .1523809524D+00                  *
* -------------------------------------------------------------- *
* REFERENCE: "Fortran Routines for Computation of Special        *
*             Functions jin.ece.uiuc.edu/routines/routines.html" *
*                                                                *
*                              TPW Release By J-P Moreau, Paris. *
*                                      (www.jpmoreau.fr)         *
*****************************************************************}
Program mbeta;

var
  BT, P, Q: double;

function GAMMA(X: double): double; Forward;

function BETA(P, Q: double): double;

{       ==========================================
!       Purpose: Compute the beta function B(p,q)
!       Input :  p  --- Parameter  ( p > 0 )
!                q  --- Parameter  ( q > 0 )
!       Output:  B(p,q)
!       Routine called: GAMMA for computing â(x)
!       ========================================== }
var
  GP, GQ, GPQ, PPQ: double;
begin
  GP := GAMMA(P);
  GQ := GAMMA(Q);
  PPQ := P + Q;
  GPQ := GAMMA(PPQ);
  BETA := (GP * GQ / GPQ)
end;


function GAMMA(X: double): double;

{       ==================================================
!       Purpose: Compute gamma function gamma(x)
!       Input :  x  --- Argument of gamma(x)
!                       ( x is not equal to 0,-1,-2)
!       Output:  gamma(x)
!       ================================================== }
var
  GA, GR, R, Z: double;
  G: Array[1..26] of double;
  K, M, M1: integer;
begin
  if X = int(X) then
    if (X > 0.0) then
    begin
      GA := 1.0;
      M1 := Round(X-1);
      for K:=2 to M1 do
        GA := GA * K;
    end else
      GA := 1E+100
  else
  begin
    if abs(X) > 1.0 then
    begin
      Z := abs(X);
      M := Round(Z);
      R := 1.0;
      for K := 1 to M do
        R := R * (Z - K);
      Z := Z - M;
    end
    else
      Z := X;

    G[1] := 1.0;
    G[2] := 0.5772156649015329;
    G[3] := -0.6558780715202538;
    G[4] := -0.420026350340952e-1;
    G[5] := 0.1665386113822915;
    G[6] := -0.421977345555443e-1;
    G[7] := -0.96219715278770e-2;
    G[8] := 0.72189432466630e-2;
    G[9] := -0.11651675918591e-2;
    G[10] := -0.2152416741149e-3;
    G[11] := 0.1280502823882e-3;
    G[12] := -0.201348547807e-4;
    G[13] := -0.12504934821e-5;
    G[14] := 0.11330272320e-5;
    G[15] := -0.2056338417e-6;
    G[16] := 0.61160950e-8;
    G[17] := 0.50020075e-8;
    G[18] := -0.11812746e-8;
    G[19] := 0.1043427e-9;
    G[20] := 0.77823e-11;
    G[21] := -0.36968e-11;
    G[22] := 0.51e-12;
    G[23] := -0.206e-13;
    G[24] := -0.54e-14;
    G[25] := 0.14e-14;
    G[26] := 0.1e-15;

    GR := G[26];
    for K := 25 downto 1 do
      GR := GR * Z + G[K];
    GA := 1.0 / (GR * Z);
    if abs(X) > 1.0 then
    begin
      GA := GA * R;
      if X < 0.0 then
        GA := -PI / (X * GA * sin(PI*X));
    end
  end;
  Result := GA;
end;

begin
  WriteLn;
  WriteLn('    p       q                B(p,q)        ');
  WriteLn('  -----------------------------------------');
  P := 1.5;
  Q := 2.0;
  BT := BETA(P, Q);
  WriteLn('  ', P:5:1, '   ', Q:5:1, '     ', BT);
  P := 2.5;
  Q := 2.0;
  BT := BETA(P, Q);
  WriteLn('  ', P:5:1, '   ', Q:5:1, '     ', BT);
  P := 1.5;
  Q := 3.0;
  BT := BETA(P, Q);
  WriteLn('  ', P:5:1, '   ', Q:5:1, '     ', BT);
  WriteLn;
  Write('Press ENTER to Close...');
  ReadLn;
end.

{ end of file mbeta.pas }
