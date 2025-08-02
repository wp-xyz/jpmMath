unit jpmAiry;

{$mode ObjFPC}{$H+}

interface

uses
  jpmTypes;

procedure AiryA(X: Float; out AI, BI, AD, BD: Float);
procedure AiryB(X: Float; out AI, BI, AD, BD: Float);
procedure AiryZero(NT, KF: integer; out XA, XB, XC, XD: TFloatArray);

implementation

{ Calculate y power x }
function Power(y, x: Float): Float;
begin
  if y <= 0 then
    Power := 0.0
  else
    Power := Exp(x*Ln(y))
end;

{===============================================================================
!  Purpose: Compute Bessel functions Jv(x) and Yv(x),
!           and modified Bessel functions Iv(x) and Kv(x),
!           and their derivatives with v=1/3,2/3
!  Input :  x --- Argument of Jv(x),Yv(x),Iv(x) and Kv(x) ( x <> 0 )
!  Output:  VJ1 --- J1/3(x)
!           VJ2 --- J2/3(x)
!           VY1 --- Y1/3(x)
!           VY2 --- Y2/3(x)
!           VI1 --- I1/3(x)
!           VI2 --- I2/3(x)
!           VK1 --- K1/3(x)
!           VK2 --- K2/3(x)
! =============================================================================}
procedure AJYIK(X: Float; out VJ1, VJ2, VY1, VY2, VI1, VI2, VK1, VK2: Float);
var
  A0,B0,CK,GN1,GN2,GP1,GP2,QX,PX,R,RP,RP2,RQ,SK,UJ1,UJ2,UU0,VV,VJL,VL,VV0,X2,XK: Float;
  C0,GN,PV1,PV2,SUM,VIL,VSL: Float;
  K, K0, L: Integer;
begin
  if X = 0.0 then
  begin
    VJ1 := 0.0;
    VJ2 := 0.0;
    VY1 := -1.0E+100;
    VY2 := 1.0E+100;
    VI1 := 0.0;
    VI2 := 0.0;
    VK1 := -1.0E+100;
    VK2 := -1.0E+100;
    exit;
  end;

  RP2 := 0.63661977236758;
  GP1 := 0.892979511569249;
  GP2 := 0.902745292950934;
  GN1 := 1.3541179394264;
  GN2 := 2.678938534707747;
  VV0 := 0.444444444444444;
  UU0 := 1.1547005383793;
  X2 := X*X;
  K0 := 12;
  if (X >= 35.0) then K0 := 10;
  if (X >= 50.0) then K0 := 8;
  if (X <= 12.0) then
  for L := 1 to 2 do
  begin
    VL := L/3.0;
    VJL := 1.0;
    R := 1.0;
    for K := 1 to 40 do
    begin
      R := -0.25 * R * X2 / (K*(K + VL));
      VJL := VJL + R;
      if (abs(R) < 1.0E-15) then break;
    end;
    A0 := Power(0.5*X, VL);
    if (L = 1) then VJ1 := A0 / GP1 * VJL;
    if (L = 2) then VJ2 := A0 / GP2 * VJL
  end
  else
    for L := 1 to 2 do
    begin
      VV := VV0 * L * L;
      PX := 1.0;
      RP := 1.0;
      for K := 1 to K0 do
      begin
        RP := -0.78125E-2 * RP * (VV - Power(4.0*K-3.0, 2.0)) * (VV - Power(4.0*K-1.0, 2.0)) / (K * (2.0*K - 1.0) * X2);
        PX := PX + RP;
      end;
      QX := 1.0;
      RQ := 1.0;
      for K := 1 to K0 do
      begin
        RQ := -0.78125E-2 * RQ * (VV - Power(4.0*K-1.0, 2.0)) * (VV - Power(4.0*K+1.0, 2.0)) / (K * (2.0*K + 1.0) * X2);
        QX := QX + RQ;
      end;

      QX := 0.125 * (VV-1.0) * QX / X;
      XK := X - (0.5*L/3.0 + 0.25) * PI;
      A0 := sqrt(RP2/ X);
      CK := cos(XK);
      SK := sin(XK);

      if (L = 1) then
      begin
        VJ1 := A0 * (PX*CK - QX*SK);
        VY1 := A0 * (PX*SK + QX*CK)
      end
      else if (L = 2) then
      begin
        VJ2 := A0 * (PX*CK - QX*SK);
        VY2 := A0 * (PX*SK + QX*CK)
      end
    end;

  if (X <= 12.0) then
  begin
    for L := 1 to 2 do
    begin
      VL := L / 3.0;
      VJL := 1.0;
      R := 1.0;
      for K := 1 to 40 do
      begin
        R := -0.25 * R * X2 / (K*(K - VL));
        VJL := VJL + R;
        if (abs(R) < 1.0E-15) then break;
      end;
      B0 := Power(2.0/X, VL);
      if (L = 1) then UJ1 := B0 * VJL / GN1;
      if (L = 2) then UJ2 := B0 * VJL / GN2
    end;
    PV1 := PI / 3.0;
    PV2 := PI / 1.5;
    VY1 := UU0 * (VJ1*cos(PV1) - UJ1);
    VY2 := UU0 * (VJ2*cos(PV2) - UJ2)
  end;

  if (X <= 18.0) then
    for L := 1 to 2 do
    begin
      VL := L / 3.0;
      VIL := 1.0;
      R := 1.0;
      for K := 1 to 40 do
      begin
        R := 0.25 * R * X2 / (K*(K + VL));
        VIL := VIL + R;
        if (abs(R) < 1.0E-15) then break;
      end;
      A0 := Power(0.5*X, VL);
      if (L = 1) then VI1 := A0 / GP1 * VIL;
      if (L = 2) then VI2 := A0 / GP2 * VIL
    end
  else
  begin
    C0 := exp(X) / sqrt(2.0 * PI * X);
    for L := 1 to 2 do
    begin
      VV := VV0 * L * L;
      VSL := 1.0;
      R := 1.0;
      for K := 1 to K0 do
      begin
        R := -0.125 * R * (VV - power(2.0*K-1.0, 2.0)) / (K*X);
        VSL := VSL + R;
      end;
      if (L = 1) then VI1 := C0 * VSL;
      if (L = 2) then VI2 := C0 * VSL
    end
  end;

  if (X <= 9.0) then
    for L := 1 to 2 do
    begin
      VL := L / 3.0;
      if (L = 1) then GN := GN1;
      if (L = 2) then GN := GN2;
      A0 := Power(2.0/X, VL) / GN;
      SUM := 1.0;
      R := 1.0;
      for K := 1 to 60 do
      begin
        R := 0.25 * R * X2 / (K*(K - VL));
        SUM := SUM + R;
        if (abs(R) < 1.0E-15) then break;
      end;
      if (L = 1) then VK1 := 0.5 * UU0 * PI * (SUM*A0 - (VI1));
      if (L = 2) then VK2 := 0.5 * UU0 * PI * (SUM*A0 - (VI2));
    end
  else
  begin
    C0 := exp(-X) * sqrt(0.5 * PI / X);
    for L := 1 to 2 do
    begin
      VV := VV0 * L * L;
      SUM := 1.0;
      R := 1.0;
      for K := 1 to K0 do
      begin
        R := 0.125 * R * (VV - power(2.0*K-1.0, 2.0)) / (K*X);
        SUM := SUM + R;
      end;
      if (L = 1) then VK1 := C0 * SUM;
      if (L = 2) then VK2 := C0 * SUM;
    end
  end;
end;

{===============================================================================
!  Purpose: Compute Airy functions and their derivatives
!  Input:   x  --- Argument of Airy function
!  Output:  AI --- Ai(x)
!           BI --- Bi(x)
!           AD --- Ai'(x)
!           BD --- Bi'(x)
!  Routine called:
!           AJYIK for computing Jv(x), Yv(x), Iv(x) and Kv(x) with v=1/3 and 2/3
! =============================================================================}
procedure AiryA(X: Float; out AI, BI, AD, BD: Float);
var
  C1,C2,PIR,SR3,VI1,VI2,VJ1,VJ2,VK1,VK2,VY1,VY2,XA,XQ,Z: Float;
begin
  XA := abs(X);
  PIR := 0.318309886183891;
  C1  := 0.355028053887817;
  C2  := 0.258819403792807;
  SR3 := 1.732050807568877;

  Z := Power(XA, 1.5) / 1.5;
  XQ := sqrt(XA);

  AJYIK(Z, VJ1, VJ2, VY1, VY2, VI1, VI2, VK1, VK2);

  if (X = 0.0) then
  begin
    AI := C1;
    BI := SR3*C1;
    AD := -C2;
    BD := SR3*C2
  end
  else if (X > 0.0) then
  begin
    AI := PIR * XQ / SR3 * VK1;
    BI := XQ * (PIR*VK1 + 2.0/SR3*VI1);
    AD := -XA / SR3*PIR*VK2;
    BD := XA * (PIR*VK2 + 2.0/SR3*VI2)
  end
  else
  begin
    AI :=  0.5*XQ * (VJ1 - VY1/SR3);
    BI := -0.5*XQ * (VJ1/SR3 + VY1);
    AD :=  0.5*XA * (VJ2 + VY2/SR3);
    BD :=  0.5*XA * (VJ2/SR3 - VY2)
  end;
end;

{===============================================================================
!  Purpose: Compute Airy functions and their derivatives
!  Input:   x  --- Argument of Airy function
!  Output:  AI --- Ai(x)
!           BI --- Bi(x)
!           AD --- Ai'(x)
!           BD --- Bi'(x)
! =============================================================================}
procedure AiryB(X: Float; out AI, BI, AD, BD: Float);

type
  VEC1 = Array[1..40] of Float;
var
  CK, DK: VEC1;
  C1,C2,DF,DG,EPS,FX,GX,R,RP,SAD,SAI,SR3,XA,XAR,XE,XF,XM,XQ,XR1: Float;
  SBD,SBI,SDA,SDB,SSA,SSB,XCS,XP1,XR2,XSS: Float;
  K,KM: integer;
begin
  EPS := 1E-15;
  C1 := 0.355028053887817;
  C2 := 0.258819403792807;
  SR3 := 1.732050807568877;
  XA := abs(X);
  XQ := sqrt(XA);
  if (X > 0.0) then XM := 5.0;
  if (X <= 0.0) then XM := 8.0;
  if (X = 0.0) then
  begin
    AI := C1;
    BI := SR3 * C1;
    AD := -C2;
    BD := SR3 * C2;
    Exit;
  end;

  if (XA <= XM) then
  begin
    FX := 1.0;
    R := 1.0;
    for K := 1 to 40 do
    begin
      R := R*X / (3.0*K) * X / (3.0*K - 1.0) * X;
      FX := FX + R;
      if (abs(R/FX) < EPS) then break;
    end;

    GX := X;
    R := X;
    for K := 1 to 40 do
    begin
      R := R * X / (3.0*K) * X / (3.0*K + 1.0) * X;
      GX := GX + R;
      if (abs(R/GX) < EPS) then break;
    end;

    AI := C1*FX - C2*GX;
    BI := SR3 * (C1*FX + C2*GX);
    DF := 0.5*X*X;
    R := DF;
    for K := 1 to 40 do
    begin
      R := R * X / (3.0*K) * X / (3.0*K + 2.0) * X;
      DF := DF + R;
      if (abs(R/DF) < EPS) then break;
    end;

    DG := 1.0;
    R := 1.0;
    for K := 1 to 40 do
    begin
      R := R * X / (3.0*K) * X / (3.0*K - 2.0) * X;
      DG := DG + R;
      if (abs(R/DG) < EPS) then break;
    end;

    AD := C1*DF - C2*DG;
    BD := SR3 * (C1*DF + C2*DG);
  end else
  begin
    XE := XA * XQ / 1.5;
    XR1 := 1.0 / XE;
    XAR := 1.0 / XQ;
    XF := sqrt(XAR);
    RP := 0.5641895835477563;
    R := 1.0;
    for K := 1 to 40 do
    begin
      R := R * (6.0*K - 1.0) / 216.0 * (6.0*K - 3.0) / K * (6.0*K - 5.0) / (2.0*K - 1.0);
      CK[K] := R;
      DK[K] := -(6.0*K + 1.0) / (6.0*K - 1.0) * CK[K];
    end;
    KM := Round(24.5 - XA);
    if (XA < 6.0) then KM := 14;
    if (XA > 15.0) then KM := 10;
    if (X > 0.0) then
    begin
      SAI := 1.0;
      SAD := 1.0;
      R := 1.0;
      for K := 1 to KM do
      begin
        R := -R * XR1;
        SAI := SAI + CK[K]*R;
        SAD := SAD + DK[K]*R;
      end;

      SBI := 1.0;
      SBD := 1.0;
      R := 1.0;
      for K := 1 to KM do
      begin
        R := R * XR1;
        SBI := SBI + CK[K]*R;
        SBD := SBD + DK[K]*R
      end;
      XP1 := exp(-XE);
      AI := 0.5 * RP * XF * XP1 * SAI;
      BI := RP * XF / XP1 * SBI;
      AD := -0.5 * RP / XF * XP1 * SAD;
      BD := RP / XF / XP1 * SBD;
    end else
    begin
      XCS := cos(XE + PI/4.0);
      XSS := sin(XE + PI/4.0);

      SSA := 1.0;
      SDA := 1.0;
      R := 1.0;
      XR2 := 1.0 / (XE*XE);
      for K := 1 to KM do
      begin
        R := -R * XR2;
        SSA := SSA + CK[2*K]*R;
        SDA := SDA + DK[2*K]*R
      end;

      SSB := CK[1] * XR1;
      SDB := DK[1] * XR1;
      R := XR1;
      for K := 1 to KM do
      begin
        R := -R * XR2;
        SSB := SSB + CK[2*K + 1] * R;
        SDB := SDB + DK[2*K + 1] * R
      end;

      AI := RP * XF * (XSS*SSA - XCS*SSB);
      BI := RP * XF * (XCS*SSA + XSS*SSB);
      AD := -RP / XF * (XCS*SDA + XSS*SDB);
      BD := RP / XF * (XSS*SDA - XCS*SDB)
    end
  end;
end;

{===============================================================================
!  Purpose: Compute the first NT zeros of Airy functions Ai(x) and Ai'(x),
!           a and a', and the associated values of Ai(a') and Ai'(a);
!           and the first NT zeros of Airy functions Bi(x) and Bi'(x), b and b',
!           and the associated values of Bi(b') and Bi'(b)
!  Input :  NT    --- Total number of zeros
!           KF    --- Function code
!                       KF=1 for Ai(x) and Ai'(x)
!                       KF=2 for Bi(x) and Bi'(x)
!  Output:  XA(m) --- a, the m-th zero of Ai(x) or b, the m-th zero of Bi(x)
!           XB(m) --- a', the m-th zero of Ai'(x) or b', the m-th zero of Bi'(x)
!           XC(m) --- Ai(a') or Bi(b')
!           XD(m) --- Ai'(a) or Bi'(b)
!           ( m --- Serial number of zeros, 0 .. NT-1 )
!  Routine called: AiryB for computing Airy functions and their derivatives
! =============================================================================}
procedure AiryZero(NT, KF: integer; out XA, XB, XC, XD: TFloatArray);
var
  I: integer;
  RT, RT0, U, U1, X: Float;
  AI, BI, AD, BD: Float;
begin
  SetLength(XA, NT);
  SetLength(XB, NT);
  SetLength(XC, NT);
  SetLength(XD, NT);
  for I := 1 to NT do
  begin
    if (KF = 1) then
    begin
      U := 3.0 * PI * (4.0*I - 1) / 8.0;
      U1 := 1 / (U*U);
      RT0 := -Power(U*U, 1.0/3.0) * ((((-15.5902*U1 + 0.929844)*U1 - 0.138889)*U1 + 0.10416667)*U1 + 1.0);
    end else
    if (KF = 2) then
      if (I = 1) then
        RT0 := -1.17371
      else
      begin
        U := 3.0 * PI * (4.0*I - 3.0) / 8.0;
        U1 := 1.0 / (U*U);
        RT0 := -Power(U*U, 1.0/3.0) * ((((-15.5902*U1 + 0.929844)*U1 - 0.138889)*U1 + 0.10416667)*U1 + 1.0);
      end;

    repeat
      X := RT0;

      AiryB(X, AI, BI, AD, BD);

      if (KF = 1) then RT := RT0 - AI/AD;
      if (KF = 2) then RT := RT0 - BI/BD;
      if (abs((RT-RT0)/RT) > 1E-9) then
        RT0 := RT
      else
      begin
        XA[I-1] := RT;
        if (KF=1) then XD[I-1] := AD;
        if (KF=2) then XD[I-1] := BD;
        break;
      end
    until false;
  end;

  for I := 1 to NT do
  begin
    if KF=1 then
      if (I=1) then
        RT0 := -1.01879
      else
      begin
        U := 3.0 * PI * (4.0*I - 3.0) / 8.0;
        U1 := 1 / (U*U);
        RT0 := -Power(U*U, 1.0/3.0)*((((15.0168*U1 - 0.873954)*U1 + 0.121528)*U1 - 0.145833)*U1 + 1.0)
      end
    else
    if (KF=2) then
      if (I=1) then
        RT0 := -2.29444
      else
      begin
        U := 3.0 * PI * (4.0*I - 1.0) / 8.0;
        U1 := 1.0 / (U*U);
        RT0 := -Power(U*U, 1.0/3.0)*((((15.0168*U1 - 0.873954)*U1 + 0.121528)*U1 - 0.145833)*U1 + 1.0)
      end;

    repeat
      X := RT0;

      AiryB(X, AI, BI, AD, BD);

      if (KF = 1) then RT := RT0 - AD / (AI*X);
      if (KF = 2) then RT := RT0 - BD / (BI*X);
      if (abs((RT-RT0)/RT) > 1E-9) then
        RT0 := RT
      else
      begin
        XB[I-1] := RT;
        if (KF = 1) then XC[I-1] := AI;
        if (KF = 2) then XC[I-1] := BI;
        break;
      end;
    until false;
  end;
end;

end.

