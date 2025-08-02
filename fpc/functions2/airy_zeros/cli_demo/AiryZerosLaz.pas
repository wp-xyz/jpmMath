{*******************************************************************
!*      Purpose: This program computes the first NT zeros of Airy  * 
!*               functions Ai(x) and Ai'(x), and the associated    *
!*               values of Ai(a') and Ai'(a), and the first NT     *
!*               zeros of Airy functions Bi(x) and Bi'(x), and     *
!*               the associated values of Bi(b') and Bi'(b) using  * 
!*               subroutine AIRYZO                                 *
!*      Input :  NT    --- Total number of zeros                   *
!*               KF    --- Function code                           *
!*                         KF=1 for Ai(x) and Ai'(x)               *
!*                         KF=2 for Bi(x) and Bi'(x)               *
!*      Output:  XA(m) --- a, the m-th zero of Ai(x) or            *
!*                         b, the m-th zero of Bi(x)               *
!*               XB(m) --- a', the m-th zero of Ai'(x) or          *
!*                         b', the m-th zero of Bi'(x)             *
!*               XC(m) --- Ai(a') or Bi(b')                        *
!*               XD(m) --- Ai'(a) or Bi'(b)                        *
!*                         ( m --- Serial number of zeros )        * 
!*      Example: NT=5                                              *
!*                                                                 *
!*      m         a            Ai'(a)         a'          Ai(a')   *
!*     ----------------------------------------------------------- *
!*      1    -2.33810741     .70121082   -1.01879297    .53565666  *
!*      2    -4.08794944    -.80311137   -3.24819758   -.41901548  *
!*      3    -5.52055983     .86520403   -4.82009921    .38040647  *
!*      4    -6.78670809    -.91085074   -6.16330736   -.35790794  *
!*      5    -7.94413359     .94733571   -7.37217726    .34230124  *
!*                                                                 *
!*      m         b            Bi'(b)         b'          Bi(b')   *
!*     ----------------------------------------------------------- *
!*      1    -1.17371322     .60195789   -2.29443968   -.45494438  *
!*      2    -3.27109330    -.76031014   -4.07315509    .39652284  *
!*      3    -4.83073784     .83699101   -5.51239573   -.36796916  *
!*      4    -6.16985213    -.88947990   -6.78129445    .34949912  *
!*      5    -7.37676208     .92998364   -7.94017869   -.33602624  *
!*                                                                 *
!* --------------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special Func-   *
!*             tions jin.ece.uiuc.edu/routines/routines.html".     *
!*                                                                 *
!*                               TPW Release By J-P Moreau, Paris. *
!*                                       (www.jpmoreau.fr)         *
!******************************************************************}
Program AiryZerosLaz;

Procedure AIRYB(X: double; out AI,BI,AD,BD: double); Forward;

{ Calculate y power x }
function Power(y, x: Double): Double;
begin
  if y <= 0 then
    Power := 0.0
  else
    Power := Exp(x*Ln(y))
end;

type
  TDoubleArray = array of double;

procedure AIRYZO(NT, KF: integer; out XA, XB, XC, XD: TDoubleArray);

{       ========================================================
!       Purpose: Compute the first NT zeros of Airy functions
!                Ai(x) and Ai'(x), a and a', and the associated
!                values of Ai(a') and Ai'(a); and the first NT
!                zeros of Airy functions Bi(x) and Bi'(x), b and
!                b', and the associated values of Bi(b') and
!                Bi'(b)
!       Input :  NT    --- Total number of zeros
!                KF    --- Function code
!                          KF=1 for Ai(x) and Ai'(x)
!                          KF=2 for Bi(x) and Bi'(x)
!       Output:  XA(m) --- a, the m-th zero of Ai(x) or
!                          b, the m-th zero of Bi(x) 
!                XB(m) --- a', the m-th zero of Ai'(x) or
!                          b', the m-th zero of Bi'(x)
!                XC(m) --- Ai(a') or Bi(b')
!                XD(m) --- Ai'(a) or Bi'(b)
!                          ( m --- Serial number of zeros )
!       Routine called: AIRYB for computing Airy functions and
!                       their derivatives
!       ======================================================= }
//label 10, 20;
var
  I: integer;
  RT,RT0,U,U1,X: double;
  AI,BI,AD,BD: double;
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
    end
    else
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

      AIRYB(X, AI, BI, AD, BD);

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
           
      AIRYB(X,AI,BI,AD,BD);
           
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

Procedure AIRYB(X:double; out AI, BI, AD, BD: double);

{       =======================================================
!       Purpose: Compute Airy functions and their derivatives
!       Input:   x  --- Argument of Airy function
!       Output:  AI --- Ai(x)
!                BI --- Bi(x)
!                AD --- Ai'(x)
!                BD --- Bi'(x)
!       ======================================================= }
type
  VEC1 = Array[1..40] of double;
var
  CK, DK: VEC1;
  C1,C2,DF,DG,EPS,FX,GX,R,RP,SAD,SAI,SR3,XA,XAR,XE,XF,XM,XQ,XR1: double;
  SBD,SBI,SDA,SDB,SSA,SSB,XCS,XP1,XR2,XSS: double;
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

var
  K, KF, NT: integer;
  XA, XB, XC, XD: TDoubleArray;

begin
  WriteLn;
  WriteLn(' KF=1 for Ai(x) and Ai''(x); KF=2 for Bi(x) and Bi''(x)');
  WriteLn(' NT is the number of the zeros.');
  Write(' Please enter NT: ');
  ReadLn(NT);

  Writeln;
  for KF := 1 to 2 do
  begin
    if (KF=1) then
      WriteLn('  m        a         Ai''(a)         a''         Ai(a'')')
    else if (KF=2) then
      WriteLn('  m        b         Bi''(b)         b''         Bi(b'')');
    WriteLn('----------------------------------------------------------');
        
    AIRYZO(NT, KF, XA, XB, XC, XD);
        
    for K := 0 to NT-1 do
    begin
      WriteLn(K:3, ' ', XA[K]:12:8, ' ', XD[K]:12:8, ' ' ,XB[K]:12:8, ' ', XC[K]:12:8)
    end;
    WriteLn;
  end;
  Write('Press ENTER to close...');
  ReadLn;
end.

{end of file mairyzo.pas}
