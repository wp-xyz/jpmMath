{*************************************************************************
*    Procedures and Functions used By programs TZBESJ, TZBESK, TZBESY    *
*    (Evalute Bessel Functions with complex argument, 1st to 3rd kind)   *
* ---------------------------------------------------------------------- *
* Reference:  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES, 1983.       *
*                                                                        *
*                      Pascal Release By J-P Moreau, Paris (07/01/2005). *
*                                     (www.jpmoreau.fr)                  *
*************************************************************************}
UNIT CBESS3;

INTERFACE

Uses Complex, CBess0, CBess00, CBess1;

  Procedure ZBUNK(ZR, ZI, FNU:double; KODE, MR, N:Integer; Var YR, YI:VEC;
                  Var NZ: Integer; TOL, ELIM, ALIM:double);

  Procedure ZACON(ZR, ZI, FNU:double; KODE, MR, N:integer; Var YR, YI:VEC; Var NZ:Integer;
                  RL, FNUL, TOL, ELIM, ALIM:Double);


IMPLEMENTATION

  Procedure ZUNK1(ZR, ZI, FNU: double; KODE, MR, N: integer; Var YR, YI:VEC;
                  Var NZ:integer; TOL, ELIM, ALIM:double); Forward;

  Procedure ZUNK2(ZR, ZI, FNU: double; KODE, MR, N: integer; Var YR, YI:VEC;
                  Var NZ:integer; TOL, ELIM, ALIM:double); Forward;

Procedure ZBUNK(ZR, ZI, FNU:double; KODE, MR, N:Integer; Var YR, YI:VEC;
                Var NZ: Integer; TOL, ELIM, ALIM:double);
{***BEGIN PROLOGUE  ZBUNK
!***REFER TO  ZBESK,ZBESH
!
!     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
!     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
!
!***ROUTINES CALLED  ZUNK1,ZUNK2
!***END PROLOGUE  ZBUNK
!     COMPLEX Y,Z }
Label 10, Return;
Var
      AX, AY: Double;
Begin
      NZ := 0;
      AX := ABS(ZR)*1.7321;
      AY := ABS(ZI);
      IF AY > AX THEN GOTO 10;
{-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3.LE.ARG(Z).LE.PI/3
!----------------------------------------------------------------------}
      ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM);
      GOTO RETURN;
{-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M:=+I OR -I
!     AND HPI:=PI/2
!----------------------------------------------------------------------}
10:   ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM);
Return: End; {ZBUNK}


Procedure ZUNK1(ZR, ZI, FNU: double; KODE, MR, N: integer; Var YR, YI:VEC;
                Var NZ:integer; TOL, ELIM, ALIM:double);
{***BEGIN PROLOGUE  ZUNK1
!***REFER TO  ZBESK
!
!     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSION.
!     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ:=-1 MEANS AN OVERFLOW WILL OCCUR
!
!***ROUTINES CALLED  ZKSCL,ZS1S2,ZUCHK,ZUNIK,D1MACH,ZABS
!***END PROLOGUE  ZUNK1
!     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO,
!     C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR}
Label 10,20,30,40,50,60,70,75,80,90,95,100,120,160,170,172,175,180,
      200,210,220,230,250,255,260,270,275,280,290,300,Return;
Type  pMAT16 = ^MAT16;
      MAT16 = Array[1..16,1..3] of Double;
Var
      ANG, APHI, ASC, ASCLE, CKI, CKR, CONEI, CONER, CRSC, CSCL, CSGNI,
      CSPNI, CSPNR, CSR, C1I, C1R, C2I, C2M, C2R, FMR, FN, FNF, PHIDI,
      PHIDR, RAST, RAZR, RS1, RZI, RZR, SGN, STI, STR, SUMDI, SUMDR,
      S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZET1DI, ZET1DR, ZET2DI, ZET2DR,
      ZRI, ZRR: Double;
      I, IB, IFLAG, IFN, II, IL, INU, IUF, K, KDFLG, KFLAG, KK, NW, INITD,
      IC, IPARD, J, M: Integer;
      BRY, SUMR, SUMI, ZETA1R, ZETA1I, ZETA2R, ZETA2I, CYR, CYI: VEC;
      INIT:array[1..3] of Integer;
      CWRKR, CWRKI: pMAT16;
      CSSR, CSRR, PHIR, PHII: VEC;
      TMPR,TMPI:VEC16;
Begin
      New(CWRKR); New(CWRKI);
      ZEROR:=0.0; ZEROI:=0.0; CONER:=1.0; CONEI:=0.0;
      KDFLG := 1;
      NZ := 0;
{-----------------------------------------------------------------------
!     EXP(-ALIM):=EXP(-ELIM)/TOL:=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!----------------------------------------------------------------------}
      CSCL := 1.0/TOL;
      CRSC := TOL;
      CSSR[1] := CSCL;
      CSSR[2] := CONER;
      CSSR[3] := CRSC;
      CSRR[1] := CRSC;
      CSRR[2] := CONER;
      CSRR[3] := CSCL;
      BRY[1] := 1000*D1MACH(1)/TOL;
      BRY[2] := 1.0/BRY[1];
      BRY[3] := D1MACH(2);
      ZRR := ZR;
      ZRI := ZI;
      IF ZR >= 0.0 THEN GOTO 10;
      ZRR := -ZR;
      ZRI := -ZI;
10:   J := 2;
      For I:=1 to N do
      begin
{-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J := 3 - J
!----------------------------------------------------------------------}
        J := 3 - J;
        FN := FNU + 1.0*(I-1);
        INIT[J] := 0;
        For II:=1 to 16 do
        begin
          TMPR[II]:=CWRKR^[II,J];
          TMPI[II]:=CWRKI^[II,J]
        end;
        ZUNIK(ZRR, ZRI, FN, 2, 0, TOL, INIT[J], PHIR[J], PHII[J],
              ZETA1R[J], ZETA1I[J], ZETA2R[J], ZETA2I[J], SUMR[J], SUMI[J],
              TMPR, TMPI);
        For II:=1 to 16 do
        begin
          CWRKR^[II,J]:=TMPR[II];
          CWRKI^[II,J]:=TMPI[II]
        end;
        IF KODE = 1 THEN GOTO 20;
        STR := ZRR + ZETA2R[J];
        STI := ZRI + ZETA2I[J] ;
        RAST := FN/ZABS(STR,STI);
        STR := STR*RAST*RAST;
        STI := -STI*RAST*RAST;
        S1R := ZETA1R[J] - STR;
        S1I := ZETA1I[J] - STI;
        GOTO 30;
20:     S1R := ZETA1R[J] - ZETA2R[J];
        S1I := ZETA1I[J] - ZETA2I[J];
30:     RS1 := S1R;
{-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!----------------------------------------------------------------------}
        IF ABS(RS1) > ELIM THEN GOTO 60;
        IF KDFLG = 1 THEN KFLAG := 2;
        IF ABS(RS1) < ALIM THEN GOTO 40;
{-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!----------------------------------------------------------------------}
        APHI := ZABS(PHIR[J],PHII[J]);
        RS1 := RS1 + Ln(APHI);
        IF ABS(RS1) > ELIM THEN GOTO 60;
        IF KDFLG = 1 THEN KFLAG := 1;
        IF RS1 < 0.0 THEN GOTO 40;
        IF KDFLG = 1 THEN KFLAG := 3;
{-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!----------------------------------------------------------------------}
40:     S2R := PHIR[J]*SUMR[J] - PHII[J]*SUMI[J];
        S2I := PHIR[J]*SUMI[J] + PHII[J]*SUMR[J];
        STR := EXP(S1R)*CSSR[KFLAG];
        S1R := STR*COS(S1I);
        S1I := STR*SIN(S1I);
        STR := S2R*S1R - S2I*S1I;
        S2I := S1R*S2I + S2R*S1I;
        S2R := STR;
        IF KFLAG <> 1 THEN GOTO 50;
        ZUCHK(S2R, S2I, NW, BRY[1], TOL);
        IF NW <> 0 THEN GOTO 60;
50:     CYR[KDFLG] := S2R;
        CYI[KDFLG] := S2I;
        YR[I] := S2R*CSRR[KFLAG];
        YI[I] := S2I*CSRR[KFLAG];
        IF KDFLG = 2 THEN GOTO 75;
        KDFLG := 2;
        GOTO 70;
60:     IF RS1 > 0.0 THEN GOTO 300;
{-----------------------------------------------------------------------
!     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!----------------------------------------------------------------------}
        IF ZR < 0.0 THEN GOTO 300;
        KDFLG := 1;
        YR[I]:=ZEROR;
        YI[I]:=ZEROI;
        NZ:=NZ+1;
        IF I = 1 THEN GOTO 70;
        IF (YR[I-1] = ZEROR) AND (YI[I-1] = ZEROI) THEN GOTO 70;
        YR[I-1]:=ZEROR;
        YI[I-1]:=ZEROI;
        NZ:=NZ+1;
70:   end;
      I := N;
75:   RAZR := 1.0/ZABS(ZRR,ZRI);
      STR := ZRR*RAZR;
      STI := -ZRI*RAZR;
      RZR := (STR+STR)*RAZR;
      RZI := (STI+STI)*RAZR;
      CKR := FN*RZR;
      CKI := FN*RZI;
      IB := I + 1;
      IF N < IB THEN GOTO 160;
{-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
!     ON UNDERFLOW.
!----------------------------------------------------------------------}
      FN := FNU + 1.0*(N-1);
      IPARD := 1;
      IF MR <> 0 THEN IPARD := 0;
      INITD := 0;
      For II:=1 to 16 do
      begin
        TMPR[II]:=CWRKR^[II,3];
        TMPI[II]:=CWRKI^[II,3]
      end;
      ZUNIK(ZRR, ZRI, FN, 2, IPARD, TOL, INITD, PHIDR, PHIDI, ZET1DR,
            ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI, TMPR, TMPI);
      For II:=1 to 16 do
      begin
        CWRKR^[II,3]:=TMPR[II];
        CWRKI^[II,3]:=TMPI[II]
      end;
      IF KODE = 1 THEN GOTO 80;
      STR := ZRR + ZET2DR;
      STI := ZRI + ZET2DI;
      RAST := FN/ZABS(STR,STI);
      STR := STR*RAST*RAST;
      STI := -STI*RAST*RAST;
      S1R := ZET1DR - STR;
      S1I := ZET1DI - STI;
      GOTO 90;
80:   S1R := ZET1DR - ZET2DR;
      S1I := ZET1DI - ZET2DI;
90:   RS1 := S1R;
      IF ABS(RS1) > ELIM THEN GOTO 95;
      IF ABS(RS1) < ALIM THEN GOTO 100;
{-----------------------------------------------------------------------
!     REFINE ESTIMATE AND TEST
!----------------------------------------------------------------------}
      APHI := ZABS(PHIDR,PHIDI);
      RS1 := RS1+Ln(APHI);
      IF ABS(RS1) < ELIM THEN GOTO 100;
95:   IF ABS(RS1) > 0.0 THEN GOTO 300;
{-----------------------------------------------------------------------
!     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!----------------------------------------------------------------------}
      IF ZR < 0.0 THEN GOTO 300;
      NZ := N;
      For I:=1 to N do
      begin
        YR[I] := ZEROR;
        YI[I] := ZEROI
      end;
      GOTO RETURN;
{-----------------------------------------------------------------------
!     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
!----------------------------------------------------------------------}
100:  S1R := CYR[1];
      S1I := CYI[1];
      S2R := CYR[2];
      S2I := CYI[2];
      C1R := CSRR[KFLAG];
      ASCLE := BRY[KFLAG];
      For I:=IB to N do
      begin
        C2R := S2R;
        C2I := S2I;
        S2R := CKR*C2R - CKI*C2I + S1R;
        S2I := CKR*C2I + CKI*C2R + S1I;
        S1R := C2R;
        S1I := C2I;
        CKR := CKR + RZR;
        CKI := CKI + RZI;
        C2R := S2R*C1R;
        C2I := S2I*C1R;
        YR[I] := C2R;
        YI[I] := C2I;
        IF KFLAG >= 3 THEN GOTO 120;
        STR := ABS(C2R);
        STI := ABS(C2I);
        C2M := DMAX(STR,STI);
        IF C2M <= ASCLE THEN GOTO 120;
        KFLAG := KFLAG + 1;
        ASCLE := BRY[KFLAG];
        S1R := S1R*C1R;
        S1I := S1I*C1R;
        S2R := C2R;
        S2I := C2I;
        S1R := S1R*CSSR[KFLAG];
        S1I := S1I*CSSR[KFLAG];
        S2R := S2R*CSSR[KFLAG];
        S2I := S2I*CSSR[KFLAG];
        C1R := CSRR[KFLAG];
120:  end;
160:  IF MR = 0 THEN GOTO RETURN;
{-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
!----------------------------------------------------------------------}
      NZ := 0;
      FMR := 1.0*MR;
      SGN := -SIGN(PI,FMR);
{-----------------------------------------------------------------------
!     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
!----------------------------------------------------------------------}
      CSGNI := SGN;
      INU := Round(FNU);
      FNF := FNU - 1.0*INU;
      IFN := INU + N - 1;
      ANG := FNF*SGN;
      CSPNR := COS(ANG);
      CSPNI := SIN(ANG);
      IF (IFN Mod 2) = 0 THEN GOTO 170;
      CSPNR := -CSPNR;
      CSPNI := -CSPNI;
170:  ASC := BRY[1];
      IUF := 0;
      KK := N;
      KDFLG := 1;
      IB := IB - 1;
      IC := IB - 1;
      For K:=1 to N do
      begin
        FN := FNU + 1.0*(KK-1);
{-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!----------------------------------------------------------------------}
        M:=3;
        IF N > 2 THEN GOTO 175;
172:    INITD := INIT[J];
        PHIDR := PHIR[J];
        PHIDI := PHII[J];
        ZET1DR := ZETA1R[J];
        ZET1DI := ZETA1I[J];
        ZET2DR := ZETA2R[J];
        ZET2DI := ZETA2I[J];
        SUMDR := SUMR[J];
        SUMDI := SUMI[J];
        M := J;
        J := 3 - J;
        GOTO 180;
175:    IF (KK = N) AND (IB < N) THEN GOTO 180;
        IF (KK = IB) OR (KK = IC) THEN GOTO 172;
        INITD := 0;
180:    For II:=1 to 16 do
        begin
          TMPR[II]:=CWRKR^[II,M];
          TMPI[II]:=CWRKI^[II,M]
        end;
        ZUNIK(ZRR, ZRI, FN, 1, 0, TOL, INITD, PHIDR, PHIDI,
              ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI,
              TMPR,TMPI);
        For II:=1 to 16 do
        begin
          CWRKR^[II,M]:=TMPR[II];
          CWRKI^[II,M]:=TMPI[II]
        end;
        IF KODE = 1 THEN GOTO 200;
        STR := ZRR + ZET2DR;
        STI := ZRI + ZET2DI;
        RAST := FN/ZABS(STR,STI);
        STR := STR*RAST*RAST;
        STI := -STI*RAST*RAST;
        S1R := -ZET1DR + STR;
        S1I := -ZET1DI + STI;
        GOTO 210;
200:    S1R := -ZET1DR + ZET2DR;
        S1I := -ZET1DI + ZET2DI;
{-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!----------------------------------------------------------------------}
210:    RS1 := S1R;
        IF ABS(RS1) > ELIM THEN GOTO 260;
        IF KDFLG = 1 THEN IFLAG := 2;
        IF ABS(RS1) < ALIM THEN GOTO 220;
{-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!----------------------------------------------------------------------}
        APHI := ZABS(PHIDR,PHIDI);
        RS1 := RS1 + Ln(APHI);
        IF ABS(RS1) > ELIM THEN GOTO 260;
        IF KDFLG = 1 THEN IFLAG := 1;
        IF RS1 < 0.0 THEN GOTO 220;
        IF KDFLG = 1 THEN IFLAG := 3;
220:    STR := PHIDR*SUMDR - PHIDI*SUMDI;
        STI := PHIDR*SUMDI + PHIDI*SUMDR;
        S2R := -CSGNI*STI;
        S2I := CSGNI*STR;
        STR := EXP(S1R)*CSSR[IFLAG];
        S1R := STR*COS(S1I);
        S1I := STR*SIN(S1I);
        STR := S2R*S1R - S2I*S1I;
        S2I := S2R*S1I + S2I*S1R;
        S2R := STR;
        IF IFLAG <> 1 THEN GOTO 230;
        ZUCHK(S2R, S2I, NW, BRY[1], TOL);
        IF NW = 0 THEN GOTO 230;
        S2R := ZEROR;
        S2I := ZEROI;
230:    CYR[KDFLG] := S2R;
        CYI[KDFLG] := S2I;
        C2R := S2R;
        C2I := S2I;
        S2R := S2R*CSRR[IFLAG];
        S2I := S2I*CSRR[IFLAG];
{-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I:=1,N
!----------------------------------------------------------------------}
        S1R := YR[KK];
        S1I := YI[KK];
        IF KODE = 1 THEN GOTO 250;
        ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF);
        NZ := NZ + NW;
250:    YR[KK] := S1R*CSPNR - S1I*CSPNI + S2R;
        YI[KK] := CSPNR*S1I + CSPNI*S1R + S2I;
        KK := KK - 1;
        CSPNR := -CSPNR;
        CSPNI := -CSPNI;
        IF (C2R <> 0.0) OR (C2I <> 0.0) THEN GOTO 255;
        KDFLG := 1;
        GOTO 270;
255:    IF KDFLG = 2 THEN GOTO 275;
        KDFLG := 2;
        GOTO 270;
260:    IF RS1 > 0.0 THEN GOTO 300;
        S2R := ZEROR;
        S2I := ZEROI;
        GOTO 230;
270:  end; {K loop}
      K := N;
275:  IL := N - K;
      IF IL = 0 THEN GOTO RETURN;
{-----------------------------------------------------------------------
!     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
!     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
!----------------------------------------------------------------------}
      S1R := CYR[1];
      S1I := CYI[1];
      S2R := CYR[2];
      S2I := CYI[2];
      CSR := CSRR[IFLAG];
      ASCLE := BRY[IFLAG];
      FN := 1.0*(INU+IL);
      For I:=1 to IL do
      begin
        C2R := S2R;
        C2I := S2I;
        S2R := S1R + (FN+FNF)*(RZR*C2R-RZI*C2I);
        S2I := S1I + (FN+FNF)*(RZR*C2I+RZI*C2R);
        S1R := C2R;
        S1I := C2I;
        FN := FN - 1.0;
        C2R := S2R*CSR;
        C2I := S2I*CSR;
        CKR := C2R;
        CKI := C2I;
        C1R := YR[KK];
        C1I := YI[KK];
        IF KODE = 1 THEN GOTO 280;
        ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF);
        NZ := NZ + NW;
280:    YR[KK] := C1R*CSPNR - C1I*CSPNI + C2R;
        YI[KK] := C1R*CSPNI + C1I*CSPNR + C2I;
        KK := KK - 1;
        CSPNR := -CSPNR;
        CSPNI := -CSPNI;
        IF IFLAG >= 3 THEN GOTO 290;
        C2R := ABS(CKR);
        C2I := ABS(CKI);
        C2M := DMAX(C2R,C2I);
        IF C2M <= ASCLE THEN GOTO 290;
        IFLAG := IFLAG + 1;
        ASCLE := BRY[IFLAG];
        S1R := S1R*CSR;
        S1I := S1I*CSR;
        S2R := CKR;
        S2I := CKI;
        S1R := S1R*CSSR[IFLAG];
        S1I := S1I*CSSR[IFLAG];
        S2R := S2R*CSSR[IFLAG];
        S2I := S2I*CSSR[IFLAG];
        CSR := CSRR[IFLAG];
290:  end;
      GOTO RETURN;
300:  NZ := -1;
Return:Dispose(CWRKR); Dispose(CWRKI);
End; {ZUNK1}


Procedure ZUNK2(ZR, ZI, FNU:double; KODE, MR, N:Integer; Var YR, YI:VEC;
                Var NZ:Integer; TOL, ELIM, ALIM:Double);
{***BEGIN PROLOGUE  ZUNK2
!***REFER TO  ZBESK
!
!     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
!     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND:=(3-MR)/2, MR:=+1 OR
!     -1. HERE ZN:=ZR*I OR -ZR*I WHERE ZR:=Z IF Z IS IN THE RIGHT
!     HALF PLANE OR ZR:=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
!     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ:=-1 MEANS AN OVERFLOW WILL OCCUR
!
!***ROUTINES CALLED  ZAIRY,ZKSCL,ZS1S2,ZUCHK,ZUNHJ,D1MACH,ZABS
!***END PROLOGUE  ZUNK2
!     COMPLEX AI,ARG,ARGD,ASUM,ASUMD,BSUM,BSUMD,CFN,CI,CIP,CK,CONE,CRSC,
!     CR1,CR2,CS,CSCL,CSGN,CSPN,CSR,CSS,CY,CZERO,C1,C2,DAI,PHI,PHID,RZ,
!     S1,S2,Y,Z,ZB,ZETA1,ZETA1D,ZETA2,ZETA2D,ZN,ZR }
Label 10,20,30,40,50,60,70,80,85,90,100,105,120,130,172,175,180,190,210,
      220,230,240,250,255,270,280,290,295,300,310,320,Return;
Var
      AARG, AIC, AII, AIR, ANG, APHI, ARGDI, ARGDR, ASC, ASCLE, ASUMDI,
      ASUMDR, BSUMDI, BSUMDR,CAR, CKI, CKR, CONEI, CONER, CRSC, CR1I, CR1R,
      CR2I, CR2R, CSCL, CSGNI, CSI, CSPNI, CSPNR, CSR, C1I, C1R, C2I, C2M,
      C2R, DAII, DAIR, FMR, FN, FNF, HPI, PHIDI, PHIDR, PTI, PTR, RAST,
      RAZR, RS1, RZI, RZR, SAR, SGN, STI, STR, S1I, S1R, S2I, S2R, YY,
      ZBI, ZBR, ZEROI, ZEROR, ZET1DI, ZET1DR, ZET2DI, ZET2DR, ZNI, ZNR,
      ZRI, ZRR: Double;
      I, IB, IFLAG, IFN, IL, IN0, INU, IUF, K, KDFLG, KFLAG, KK, NAI,
      NDAI, NW, IDUM, J, IPARD, IC: Integer;
      BRY, ASUMR, ASUMI, BSUMR, BSUMI, PHIR, PHII, ARGR, ARGI, ZETA1R,
      ZETA1I, ZETA2R, ZETA2I, CYR, CYI, CIPR, CIPI, CSSR, CSRR: VEC;
Begin

      ZEROR:=0.0; ZEROI:=0.0; CONER:=1.0; CONEI:=0.0;
      CR1R:= 1.0;  CR1I:= 1.73205080756887729;
      CR2R:=-0.5; CR2I:=-8.66025403784438647E-01;
      HPI:=1.57079632679489662; 
      AIC:=1.26551212348464539;
 
      CIPR[1]:= 1.0; CIPI[1]:=0.0; CIPR[2]:=0.0; CIPI[2]:=-1.0;
      CIPR[3]:=-1.0; CIPI[3]:=0.0; CIPR[4]:=0.0; CIPI[4]:= 1.0;

      KDFLG := 1;
      NZ := 0;
{-----------------------------------------------------------------------
!     EXP(-ALIM):=EXP(-ELIM)/TOL:=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!----------------------------------------------------------------------}
      CSCL := 1.0/TOL;
      CRSC := TOL;
      CSSR[1] := CSCL;
      CSSR[2] := CONER;
      CSSR[3] := CRSC;
      CSRR[1] := CRSC;
      CSRR[2] := CONER;
      CSRR[3] := CSCL;
      BRY[1] := 1000.0*D1MACH(1)/TOL;
      BRY[2] := 1.0/BRY[1];
      BRY[3] := D1MACH(2);
      ZRR := ZR;
      ZRI := ZI;
      IF ZR >= 0.0 THEN GOTO 10;
      ZRR := -ZR;
      ZRI := -ZI;
10:   YY := ZRI;
      ZNR := ZRI;
      ZNI := -ZRR;
      ZBR := ZRR;
      ZBI := ZRI;
      INU := Round(FNU);
      FNF := FNU - 1.0*INU;
      ANG := -HPI*FNF;
      CAR := COS(ANG);
      SAR := SIN(ANG);
      C2R := HPI*SAR;
      C2I := -HPI*CAR;
      KK := (INU Mod 4) + 1;
      STR := C2R*CIPR[KK] - C2I*CIPI[KK];
      STI := C2R*CIPI[KK] + C2I*CIPR[KK];
      CSR := CR1R*STR - CR1I*STI;
      CSI := CR1R*STI + CR1I*STR;
      IF YY > 0.0 THEN GOTO 20;
      ZNR := -ZNR;
      ZBI := -ZBI;
{-----------------------------------------------------------------------
!     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
!     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
!     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
!----------------------------------------------------------------------}
20:   J := 2;
      For I:=1 to N do
      begin
{-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J := 3 - J
!----------------------------------------------------------------------}
        J := 3 - J;
        FN := FNU + 1.0*(I-1);
        ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR[J], PHII[J], ARGR[J], ARGI[J],
              ZETA1R[J], ZETA1I[J], ZETA2R[J], ZETA2I[J], ASUMR[J],
              ASUMI[J], BSUMR[J], BSUMI[J]);
        IF KODE = 1 THEN GOTO 30;
        STR := ZBR + ZETA2R[J];
        STI := ZBI + ZETA2I[J];
        RAST := FN/ZABS(STR,STI);
        STR := STR*RAST*RAST;
        STI := -STI*RAST*RAST;
        S1R := ZETA1R[J] - STR;
        S1I := ZETA1I[J] - STI;
        GOTO 40;
30:     S1R := ZETA1R[J] - ZETA2R[J];
        S1I := ZETA1I[J] - ZETA2I[J];
{-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!----------------------------------------------------------------------}
40:     RS1 := S1R;
        IF ABS(RS1) > ELIM THEN GOTO 70;
        IF KDFLG = 1 THEN KFLAG := 2;
        IF ABS(RS1) < ALIM THEN GOTO 50;
{-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!----------------------------------------------------------------------}
        APHI := ZABS(PHIR[J],PHII[J]);
        AARG := ZABS(ARGR[J],ARGI[J]);
        RS1 := RS1 + Ln(APHI) - 0.25*Ln(AARG) - AIC;
        IF ABS(RS1) > ELIM THEN GOTO 70;
        IF KDFLG = 1 THEN KFLAG := 1;
        IF RS1 < 0.0 THEN GOTO 50;
        IF KDFLG = 1 THEN KFLAG := 3;
{-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!----------------------------------------------------------------------}
50:     C2R := ARGR[J]*CR2R - ARGI[J]*CR2I;
        C2I := ARGR[J]*CR2I + ARGI[J]*CR2R;
        ZAIRY(C2R, C2I, 0, 2, AIR, AII, NAI, IDUM);
        ZAIRY(C2R, C2I, 1, 2, DAIR, DAII, NDAI, IDUM);
        STR := DAIR*BSUMR[J] - DAII*BSUMI[J];
        STI := DAIR*BSUMI[J] + DAII*BSUMR[J];
        PTR := STR*CR2R - STI*CR2I;
        PTI := STR*CR2I + STI*CR2R;
        STR := PTR + (AIR*ASUMR[J]-AII*ASUMI[J]);
        STI := PTI + (AIR*ASUMI[J]+AII*ASUMR[J]);
        PTR := STR*PHIR[J] - STI*PHII[J];
        PTI := STR*PHII[J] + STI*PHIR[J];
        S2R := PTR*CSR - PTI*CSI;
        S2I := PTR*CSI + PTI*CSR;
        STR := EXP(S1R)*CSSR[KFLAG];
        S1R := STR*COS(S1I);
        S1I := STR*SIN(S1I);
        STR := S2R*S1R - S2I*S1I;
        S2I := S1R*S2I + S2R*S1I;
        S2R := STR;
        IF KFLAG <> 1 THEN GOTO 60;
        ZUCHK(S2R, S2I, NW, BRY[1], TOL);
        IF NW <> 0 THEN GOTO 70;
60:     IF YY <= 0.0 THEN S2I := -S2I;
        CYR[KDFLG] := S2R;
        CYI[KDFLG] := S2I;
        YR[I] := S2R*CSRR[KFLAG];
        YI[I] := S2I*CSRR[KFLAG];
        STR := CSI;
        CSI := -CSR;
        CSR := STR;
        IF KDFLG = 2 THEN GOTO 85;
        KDFLG := 2;
        GOTO 80;
70:     IF RS1 > 0.0 THEN GOTO 320;
{-----------------------------------------------------------------------
!     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!----------------------------------------------------------------------}
        IF ZR < 0.0 THEN GOTO 320;
        KDFLG := 1;
        YR[I]:=ZEROR;
        YI[I]:=ZEROI;
        NZ:=NZ+1;
        STR := CSI;
        CSI :=-CSR;
        CSR := STR;
        IF I = 1 THEN GOTO 80;
        IF (YR[I-1] = ZEROR) AND (YI[I-1] = ZEROI) THEN GOTO 80;
        YR[I-1]:=ZEROR;
        YI[I-1]:=ZEROI;
        NZ:=NZ+1;
80:   end; {I loop}
      I := N;
85:   RAZR := 1.0/ZABS(ZRR,ZRI);
      STR := ZRR*RAZR;
      STI := -ZRI*RAZR;
      RZR := (STR+STR)*RAZR;
      RZI := (STI+STI)*RAZR;
      CKR := FN*RZR;
      CKI := FN*RZI;
      IB := I + 1;
      IF N < IB THEN GOTO 180;
{-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
!     ON UNDERFLOW.
!----------------------------------------------------------------------}
      FN := FNU + 1.0*(N-1);
      IPARD := 1;
      IF MR <> 0 THEN IPARD := 0;
      ZUNHJ(ZNR, ZNI, FN, IPARD, TOL, PHIDR, PHIDI, ARGDR, ARGDI,
            ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR, ASUMDI, BSUMDR, BSUMDI);
      IF KODE = 1 THEN GOTO 90;
      STR := ZBR + ZET2DR;
      STI := ZBI + ZET2DI;
      RAST := FN/ZABS(STR,STI);
      STR := STR*RAST*RAST;
      STI := -STI*RAST*RAST;
      S1R := ZET1DR - STR;
      S1I := ZET1DI - STI;
      GOTO 100;
90:   S1R := ZET1DR - ZET2DR;
      S1I := ZET1DI - ZET2DI;
100:  RS1 := S1R;
      IF ABS(RS1) > ELIM THEN GOTO 105;
      IF ABS(RS1) < ALIM THEN GOTO 120;
{-----------------------------------------------------------------------
!     REFINE ESTIMATE AND TEST
!----------------------------------------------------------------------}
      APHI := ZABS(PHIDR,PHIDI);
      RS1 := RS1+Ln(APHI);
      IF ABS(RS1) < ELIM THEN GOTO 120;
105:  IF RS1 > 0.0 THEN GOTO 320;
{-----------------------------------------------------------------------
!     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!----------------------------------------------------------------------}
      IF ZR < 0.0 THEN GOTO 320;
      NZ := N;
      For I:=1 to N do
      begin
        YR[I] := ZEROR;
        YI[I] := ZEROI
      end;
      GOTO RETURN;
120:  S1R := CYR[1];
      S1I := CYI[1];
      S2R := CYR[2];
      S2I := CYI[2];
      C1R := CSRR[KFLAG];
      ASCLE := BRY[KFLAG];
      For I:=IB to N do
      begin
        C2R := S2R;
        C2I := S2I;
        S2R := CKR*C2R - CKI*C2I + S1R;
        S2I := CKR*C2I + CKI*C2R + S1I;
        S1R := C2R;
        S1I := C2I;
        CKR := CKR + RZR;
        CKI := CKI + RZI;
        C2R := S2R*C1R;
        C2I := S2I*C1R;
        YR[I] := C2R;
        YI[I] := C2I;
        IF KFLAG >= 3 THEN GOTO 130;
        STR := ABS(C2R);
        STI := ABS(C2I);
        C2M := DMAX(STR,STI);
        IF C2M <= ASCLE THEN GOTO 130;
        KFLAG := KFLAG + 1;
        ASCLE := BRY[KFLAG];
        S1R := S1R*C1R;
        S1I := S1I*C1R;
        S2R := C2R;
        S2I := C2I;
        S1R := S1R*CSSR[KFLAG];
        S1I := S1I*CSSR[KFLAG];
        S2R := S2R*CSSR[KFLAG];
        S2I := S2I*CSSR[KFLAG];
        C1R := CSRR[KFLAG];
130:  end; {I loop}
180:  IF MR = 0 THEN GOTO RETURN;
{-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
!----------------------------------------------------------------------}
      NZ := 0;
      FMR := 1.0*MR;
      SGN := -SIGN(PI,FMR);
{-----------------------------------------------------------------------
!     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
!----------------------------------------------------------------------}
      CSGNI := SGN;
      IF YY <= 0.0 THEN CSGNI := -CSGNI;
      IFN := INU + N - 1;
      ANG := FNF*SGN;
      CSPNR := COS(ANG);
      CSPNI := SIN(ANG);
      IF (IFN Mod 2) = 0 THEN GOTO 190;
      CSPNR := -CSPNR;
      CSPNI := -CSPNI;
{-----------------------------------------------------------------------
!     CS:=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
!     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
!     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
!     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
!----------------------------------------------------------------------}
190:  CSR := SAR*CSGNI;
      CSI := CAR*CSGNI;
      IN0 := (IFN Mod 4) + 1;
      C2R := CIPR[IN0];
      C2I := CIPI[IN0];
      STR := CSR*C2R + CSI*C2I;
      CSI := -CSR*C2I + CSI*C2R;
      CSR := STR;
      ASC := BRY[1];
      IUF := 0;
      KK := N;
      KDFLG := 1;
      IB := IB - 1;
      IC := IB - 1;
      For K:=1 to N do
      begin
        FN := FNU + 1.0*(KK-1);
{-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!----------------------------------------------------------------------}
        IF N > 2 THEN GOTO 175;
172:    PHIDR := PHIR[J];
        PHIDI := PHII[J];
        ARGDR := ARGR[J];
        ARGDI := ARGI[J];
        ZET1DR := ZETA1R[J];
        ZET1DI := ZETA1I[J];
        ZET2DR := ZETA2R[J];
        ZET2DI := ZETA2I[J];
        ASUMDR := ASUMR[J];
        ASUMDI := ASUMI[J];
        BSUMDR := BSUMR[J];
        BSUMDI := BSUMI[J];
        J := 3 - J;
        GOTO 210;
175:    IF (KK = N) AND (IB < N) THEN GOTO 210;
        IF (KK = IB) OR (KK = IC) THEN GOTO 172;
        ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIDR, PHIDI, ARGDR,
              ARGDI, ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR,
              ASUMDI, BSUMDR, BSUMDI);
210:    IF KODE = 1 THEN GOTO 220;
        STR := ZBR + ZET2DR;
        STI := ZBI + ZET2DI;
        RAST := FN/ZABS(STR,STI);
        STR := STR*RAST*RAST;
        STI := -STI*RAST*RAST;
        S1R := -ZET1DR + STR;
        S1I := -ZET1DI + STI;
        GOTO 230;
220:    S1R := -ZET1DR + ZET2DR;
        S1I := -ZET1DI + ZET2DI;
{-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!----------------------------------------------------------------------}
230:    RS1 := S1R;
        IF ABS(RS1) > ELIM THEN GOTO 280;
        IF KDFLG = 1 THEN IFLAG := 2;
        IF ABS(RS1) < ALIM THEN GOTO 240;
{-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!----------------------------------------------------------------------}
        APHI := ZABS(PHIDR,PHIDI);
        AARG := ZABS(ARGDR,ARGDI);
        RS1 := RS1 + Ln(APHI) - 0.25*Ln(AARG) - AIC;
        IF ABS(RS1) > ELIM THEN GOTO 280;
        IF KDFLG = 1 THEN IFLAG := 1;
        IF RS1 < 0.0 THEN GOTO 240;
        IF KDFLG = 1 THEN IFLAG := 3;
240:    ZAIRY(ARGDR, ARGDI, 0, 2, AIR, AII, NAI, IDUM);
        ZAIRY(ARGDR, ARGDI, 1, 2, DAIR, DAII, NDAI, IDUM);
        STR := DAIR*BSUMDR - DAII*BSUMDI;
        STI := DAIR*BSUMDI + DAII*BSUMDR;
        STR := STR + (AIR*ASUMDR-AII*ASUMDI);
        STI := STI + (AIR*ASUMDI+AII*ASUMDR);
        PTR := STR*PHIDR - STI*PHIDI;
        PTI := STR*PHIDI + STI*PHIDR;
        S2R := PTR*CSR - PTI*CSI;
        S2I := PTR*CSI + PTI*CSR;
        STR := EXP(S1R)*CSSR[IFLAG];
        S1R := STR*COS(S1I);
        S1I := STR*SIN(S1I);
        STR := S2R*S1R - S2I*S1I;
        S2I := S2R*S1I + S2I*S1R;
        S2R := STR;
        IF IFLAG <> 1 THEN GOTO 250;
        ZUCHK(S2R, S2I, NW, BRY[1], TOL);
        IF NW = 0 THEN GOTO 250;
        S2R := ZEROR;
        S2I := ZEROI;
250:    IF YY <= 0.0 THEN S2I := -S2I;
        CYR[KDFLG] := S2R;
        CYI[KDFLG] := S2I;
        C2R := S2R;
        C2I := S2I;
        S2R := S2R*CSRR[IFLAG];
        S2I := S2I*CSRR[IFLAG];
{-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I:=1,N
!----------------------------------------------------------------------}
        S1R := YR[KK];
        S1I := YI[KK];
        IF KODE = 1 THEN GOTO 270;
        ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF);
        NZ := NZ + NW;
270:    YR[KK] := S1R*CSPNR - S1I*CSPNI + S2R;
        YI[KK] := S1R*CSPNI + S1I*CSPNR + S2I;
        KK := KK - 1;
        CSPNR := -CSPNR;
        CSPNI := -CSPNI;
        STR := CSI;
        CSI := -CSR;
        CSR := STR;
        IF (C2R <> 0.0) OR (C2I <> 0.0) THEN GOTO 255;
        KDFLG := 1;
        GOTO 290;
255:    IF KDFLG = 2 THEN GOTO 295;
        KDFLG := 2;
        GOTO 290;
280:    IF RS1 > 0.0 THEN GOTO 320;
        S2R := ZEROR;
        S2I := ZEROI;
        GOTO 250;
290:  end; {K loop} 
      K := N;
295:  IL := N - K;
      IF IL = 0 THEN GOTO RETURN;
{-----------------------------------------------------------------------
!     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
!     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
!----------------------------------------------------------------------}
      S1R := CYR[1];
      S1I := CYI[1];
      S2R := CYR[2];
      S2I := CYI[2];
      CSR := CSRR[IFLAG];
      ASCLE := BRY[IFLAG];
      FN := 1.0*(INU+IL);
      For I:=1 to IL do
      begin
        C2R := S2R;
        C2I := S2I;
        S2R := S1R + (FN+FNF)*(RZR*C2R-RZI*C2I);
        S2I := S1I + (FN+FNF)*(RZR*C2I+RZI*C2R);
        S1R := C2R;
        S1I := C2I;
        FN := FN - 1.0;
        C2R := S2R*CSR;
        C2I := S2I*CSR;
        CKR := C2R;
        CKI := C2I;
        C1R := YR[KK];
        C1I := YI[KK];
        IF KODE = 1 THEN GOTO 300;
        ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF);
        NZ := NZ + NW;
300:    YR[KK] := C1R*CSPNR - C1I*CSPNI + C2R;
        YI[KK] := C1R*CSPNI + C1I*CSPNR + C2I;
        KK := KK - 1;
        CSPNR := -CSPNR;
        CSPNI := -CSPNI;
        IF IFLAG >= 3 THEN GOTO 310;
        C2R := ABS(CKR);
        C2I := ABS(CKI);
        C2M := DMAX(C2R,C2I);
        IF C2M <= ASCLE THEN GOTO 310;
        IFLAG := IFLAG + 1;
        ASCLE := BRY[IFLAG];
        S1R := S1R*CSR;
        S1I := S1I*CSR;
        S2R := CKR;
        S2I := CKI;
        S1R := S1R*CSSR[IFLAG];
        S1I := S1I*CSSR[IFLAG];
        S2R := S2R*CSSR[IFLAG];
        S2I := S2I*CSSR[IFLAG];
        CSR := CSRR[IFLAG];
310:  end; {I loop}
      GOTO RETURN;
320:  NZ := -1;
Return: End; {ZUNK2}


Procedure ZACON(ZR, ZI, FNU:double; KODE, MR, N:integer; Var YR, YI:VEC; Var NZ:Integer;
                RL, FNUL, TOL, ELIM, ALIM:Double);
{***BEGIN PROLOGUE  ZACON
!***REFER TO  ZBESK,ZBESH
!
!     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP)):=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP:=PI*MR*CMPLX(0.0,1.0)
!
!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE
!
!***ROUTINES CALLED  ZBINU,ZBKNU,ZS1S2,D1MACH,ZABS,ZMLT
!***END PROLOGUE  ZACON
!     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST,
!     S1,S2,Y,Z,ZN }
Label 10,20,30,40,50,60,70,80,90,Return;
Var
      ARG, ASCLE, AS2, AZN, BSCLE, CKI, CKR, CONEI, CONER, CPN, CSCL, CSCR,
      CSGNI, CSGNR, CSPNI, CSPNR, CSR, C1I, C1M, C1R, C2I, C2R, FMR, FN,
      PTI, PTR, RAZN, RZI, RZR, SC1I, SC1R, SC2I, SC2R, SGN, SPN, STI, STR,
      S1I, S1R, S2I, S2R, YY, ZEROI, ZEROR, ZNI, ZNR: Double;
      I, INU, IUF, KFLAG, NN, NW: Integer;
      CYR, CYI, CSSR, CSRR, BRY: VEC;
Begin

      ZEROR:=0.0; ZEROI:=0.0; CONER:=1.0; CONEI:=0.0;
      NZ := 0;
      ZNR := -ZR;
      ZNI := -ZI;
      NN := N;
      ZBINU(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, FNUL, TOL, ELIM, ALIM);
      IF NW < 0 THEN GOTO 90;
{-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!----------------------------------------------------------------------}
      NN := IMIN(2,N);
      ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM);
      IF NW <> 0 THEN GOTO 90;
      S1R := CYR[1];
      S1I := CYI[1];
      FMR := 1.0*MR;
      SGN := -SIGN(PI,FMR);
      CSGNR := ZEROR;
      CSGNI := SGN;
      IF KODE = 1 THEN GOTO 10;
      YY := -ZNI;
      CPN := COS(YY);
      SPN := SIN(YY);
      ZMLT(CSGNR, CSGNI, CPN, SPN, CSGNR, CSGNI);
{-----------------------------------------------------------------------
!     CALCULATE CSPN:=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!----------------------------------------------------------------------}
10:   INU := Round(FNU);
      ARG := (FNU-1.0*INU)*SGN;
      CPN := COS(ARG);
      SPN := SIN(ARG);
      CSPNR := CPN;
      CSPNI := SPN;
      IF (INU Mod 2) = 0 THEN GOTO 20;
      CSPNR := -CSPNR;
      CSPNI := -CSPNI;
20:   IUF := 0;
      C1R := S1R;
      C1I := S1I;
      C2R := YR[1];
      C2I := YI[1];
      ASCLE := 1000.0*D1MACH(1)/TOL;
      IF KODE = 1 THEN GOTO 30;
      ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF);
      NZ := NZ + NW;
      SC1R := C1R;
      SC1I := C1I;
30:   ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI);
      ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI);
      YR[1] := STR + PTR;
      YI[1] := STI + PTI;
      IF N = 1 THEN GOTO RETURN;
      CSPNR := -CSPNR;
      CSPNI := -CSPNI;
      S2R := CYR[2];
      S2I := CYI[2];
      C1R := S2R;
      C1I := S2I;
      C2R := YR[2];
      C2I := YI[2];
      IF KODE = 1 THEN GOTO 40;
      ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF);
      NZ := NZ + NW;
      SC2R := C1R;
      SC2I := C1I;
40:   ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI);
      ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI);
      YR[2] := STR + PTR;
      YI[2] := STI + PTI;
      IF N = 2 THEN GOTO RETURN;
      CSPNR := -CSPNR;
      CSPNI := -CSPNI;
      AZN := ZABS(ZNR,ZNI);
      RAZN := 1.0/AZN;
      STR := ZNR*RAZN;
      STI := -ZNI*RAZN;
      RZR := (STR+STR)*RAZN;
      RZI := (STI+STI)*RAZN;
      FN := FNU + 1.0;
      CKR := FN*RZR;
      CKI := FN*RZI;
{-----------------------------------------------------------------------
!     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
!----------------------------------------------------------------------}
      CSCL := 1.0/TOL;
      CSCR := TOL;
      CSSR[1] := CSCL;
      CSSR[2] := CONER;
      CSSR[3] := CSCR;
      CSRR[1] := CSCR;
      CSRR[2] := CONER;
      CSRR[3] := CSCL;
      BRY[1] := ASCLE;
      BRY[2] := 1.0/ASCLE;
      BRY[3] := D1MACH(2);
      AS2 := ZABS(S2R,S2I);
      KFLAG := 2;
      IF AS2 > BRY[1] THEN GOTO 50;
      KFLAG := 1;
      GOTO 60;
50:   IF AS2 < BRY[2] THEN GOTO 60;
      KFLAG := 3;
60:   BSCLE := BRY[KFLAG];
      S1R := S1R*CSSR[KFLAG];
      S1I := S1I*CSSR[KFLAG];
      S2R := S2R*CSSR[KFLAG];
      S2I := S2I*CSSR[KFLAG];
      CSR := CSRR[KFLAG];
      For I:=3 to N do
      begin
        STR := S2R;
        STI := S2I;
        S2R := CKR*STR - CKI*STI + S1R;
        S2I := CKR*STI + CKI*STR + S1I;
        S1R := STR;
        S1I := STI;
        C1R := S2R*CSR;
        C1I := S2I*CSR;
        STR := C1R;
        STI := C1I;
        C2R := YR[I];
        C2I := YI[I];
        IF KODE = 1 THEN GOTO 70;
        IF IUF < 0 THEN GOTO 70;
        ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF);
        NZ := NZ + NW;
        SC1R := SC2R;
        SC1I := SC2I;
        SC2R := C1R;
        SC2I := C1I;
        IF IUF <> 3 THEN GOTO 70;
        IUF := -4;
        S1R := SC1R*CSSR[KFLAG];
        S1I := SC1I*CSSR[KFLAG];
        S2R := SC2R*CSSR[KFLAG];
        S2I := SC2I*CSSR[KFLAG];
        STR := SC2R;
        STI := SC2I;
70:     PTR := CSPNR*C1R - CSPNI*C1I;
        PTI := CSPNR*C1I + CSPNI*C1R;
        YR[I] := PTR + CSGNR*C2R - CSGNI*C2I;
        YI[I] := PTI + CSGNR*C2I + CSGNI*C2R;
        CKR := CKR + RZR;
        CKI := CKI + RZI;
        CSPNR := -CSPNR;
        CSPNI := -CSPNI;
        IF KFLAG >= 3 THEN GOTO 80;
        PTR := ABS(C1R);
        PTI := ABS(C1I);
        C1M := DMAX(PTR,PTI);
        IF C1M <= BSCLE THEN GOTO 80;
        KFLAG := KFLAG + 1;
        BSCLE := BRY[KFLAG];
        S1R := S1R*CSR;
        S1I := S1I*CSR;
        S2R := STR;
        S2I := STI;
        S1R := S1R*CSSR[KFLAG];
        S1I := S1I*CSSR[KFLAG];
        S2R := S2R*CSSR[KFLAG];
        S2I := S2I*CSSR[KFLAG];
        CSR := CSRR[KFLAG];
80:   end; {I loop}
      GOTO RETURN;
90:   NZ := -1;
      IF NW = -2 THEN NZ:=-2;
Return: End; {ZACON}

END.

{end of file Cbess3.pas}