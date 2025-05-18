{*************************************************************************
*    Procedures and Functions used By programs TZBESJ, TZBESK, TZBESY    *
*    (Evalute Bessel Functions with complex argument, 1st to 3rd kind)   *
* ---------------------------------------------------------------------- *
* Reference:  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES, 1983.       *
*                                                                        *
*                      Pascal Release By J-P Moreau, Paris (07/01/2005). *
*                                     (www.jpmoreau.fr)                  *
*************************************************************************}
UNIT CBESS2;    {ZBUNI}

INTERFACE

Uses Complex, CBess0, CBess00;

  Procedure ZBUNI(ZR, ZI, FNU:double; KODE, N:Integer; Var YR, YI:VEC;
                  Var NZ, NUI, NLAST:Integer; Var FNUL, TOL, ELIM, ALIM:double);


IMPLEMENTATION

  Procedure ZUNI1(ZR, ZI, FNU:double; KODE, N:integer; var YR, YI:VEC; var NZ, NLAST:integer;
                  var FNUL, TOL, ELIM, ALIM:double); Forward;
  Procedure ZUNI2(ZR, ZI, FNU:double; KODE, N:integer; var YR, YI:VEC; var NZ, NLAST:integer;
                  var FNUL, TOL, ELIM, ALIM:double); Forward;


Procedure ZBUNI(ZR, ZI, FNU:double; KODE, N:Integer; Var YR, YI:VEC;
                Var NZ, NUI, NLAST:Integer; Var FNUL, TOL, ELIM, ALIM:double);
{***BEGIN PROLOGUE  ZBUNI
!***REFER TO  ZBESI,ZBESK
!
!     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z) >
!     FNUL AND FNU+N-1 < FNUL. THE ORDER IS INCREASED FROM
!     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
!     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2.
!
!***ROUTINES CALLED  ZUNI1,ZUNI2,ZABS,D1MACH
!***END PROLOGUE  ZBUNI
!     COMPLEX CSCL,CSCR,CY,RZ,ST,S1,S2,Y,Z }
Label 10,20,21,25,30,40,50,60,70,80,90, Return;
Var
      AX, AY, CSCLR, CSCRR, DFNU, FNUI, GNU, RAZ, RZI, RZR, STI, STR,
      S1I, S1R, S2I, S2R, ASCLE, C1R, C1I, C1M:Double;
      I, IFLAG, IFORM, K, NL, NW:integer;
      CYR, CYI, BRY: VEC;
Begin
      NZ := 0;
      AX := ABS(ZR)*1.7321;
      AY := ABS(ZI);
      IFORM := 1;
      IF AY > AX THEN IFORM := 2;
      IF NUI = 0 THEN GOTO 60;
      FNUI := 1.0*NUI;
      DFNU := FNU + 1.0*(N-1);
      GNU := DFNU + FNUI;
      IF IFORM = 2 THEN GOTO 10;
{-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3.LE.ARG(Z).LE.PI/3
!----------------------------------------------------------------------}
      ZUNI1(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL, ELIM, ALIM);
      GOTO 20;
{-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M:=+I OR -I
!     AND HPI:=PI/2
!----------------------------------------------------------------------}
10:   ZUNI2(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL, ELIM, ALIM);
20:   IF NW < 0 THEN GOTO 50;
      IF NW <> 0 THEN GOTO 90;
      STR := ZABS(CYR[1],CYI[1]);
{----------------------------------------------------------------------
!     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
!---------------------------------------------------------------------}
      BRY[1]:=1000*D1MACH(1)/TOL;
      BRY[2] := 1.0/BRY[1];
      BRY[3] := BRY[2];
      IFLAG := 2;
      ASCLE := BRY[2];
      CSCLR := 1.0;
      IF STR > BRY[1] THEN GOTO 21;
      IFLAG := 1;
      ASCLE := BRY[1];
      CSCLR := 1.0/TOL;
      GOTO 25;
21:   IF STR < BRY[2] THEN GOTO 25;
      IFLAG := 3;
      ASCLE:=BRY[3];
      CSCLR := TOL;
25:   CSCRR := 1.0/CSCLR;
      S1R := CYR[2]*CSCLR;
      S1I := CYI[2]*CSCLR;
      S2R := CYR[1]*CSCLR;
      S2I := CYI[1]*CSCLR;
      RAZ := 1.0/ZABS(ZR,ZI);
      STR := ZR*RAZ;
      STI := -ZI*RAZ;
      RZR := (STR+STR)*RAZ;
      RZI := (STI+STI)*RAZ;
      For I:=1 to NUI do
      begin
        STR := S2R;
        STI := S2I;
        S2R := (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R;
        S2I := (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I;
        S1R := STR;
        S1I := STI;
        FNUI := FNUI - 1.0;
        IF IFLAG >= 3 THEN GOTO 30;
        STR := S2R*CSCRR;
        STI := S2I*CSCRR;
        C1R := ABS(STR);
        C1I := ABS(STI);
        C1M := DMAX(C1R,C1I);
        IF C1M <= ASCLE THEN GOTO 30;
        IFLAG := IFLAG+1;
        ASCLE := BRY[IFLAG];
        S1R := S1R*CSCRR;
        S1I := S1I*CSCRR;
        S2R := STR;
        S2I := STI;
        CSCLR := CSCLR*TOL;
        CSCRR := 1.0/CSCLR;
        S1R := S1R*CSCLR;
        S1I := S1I*CSCLR;
        S2R := S2R*CSCLR;
        S2I := S2I*CSCLR;
30:   end;
      YR[N] := S2R*CSCRR;
      YI[N] := S2I*CSCRR;
      IF N = 1 THEN GOTO RETURN;
      NL := N - 1;
      FNUI := 1.0*NL;
      K := NL;
      For I:=1 to NL do
      begin
        STR := S2R;
        STI := S2I;
        S2R := (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R;
        S2I := (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I;
        S1R := STR;
        S1I := STI;
        STR := S2R*CSCRR;
        STI := S2I*CSCRR;
        YR[K] := STR;
        YI[K] := STI;
        FNUI := FNUI - 1.0;
        K := K - 1;
        IF IFLAG >= 3 THEN GOTO 40;
        C1R := ABS(STR);
        C1I := ABS(STI);
        C1M := DMAX(C1R,C1I);
        IF  C1M <= ASCLE THEN GOTO 40;
        IFLAG := IFLAG+1;
        ASCLE := BRY[IFLAG];
        S1R := S1R*CSCRR;
        S1I := S1I*CSCRR;
        S2R := STR;
        S2I := STI;
        CSCLR := CSCLR*TOL;
        CSCRR := 1.0/CSCLR;
        S1R := S1R*CSCLR;
        S1I := S1I*CSCLR;
        S2R := S2R*CSCLR;
        S2I := S2I*CSCLR;
40:   end;
      GOTO RETURN;
50:   NZ := -1;
      IF NW = -2 THEN NZ:=-2;
      GOTO RETURN;
60:   IF IFORM = 2 THEN GOTO 70;
{-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3.LE.ARG(Z).LE.PI/3
!----------------------------------------------------------------------}
      ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, ELIM, ALIM);
      GOTO 80;
{-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M:=+I OR -I
!     AND HPI:=PI/2
!----------------------------------------------------------------------}
70:   ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL, ELIM, ALIM);
80:   IF NW < 0 THEN GOTO 50;
      NZ := NW;
      GOTO RETURN;
90:   NLAST := N;
Return: End; {ZBUNI}


Procedure ZUNI1(ZR, ZI, FNU:double; KODE, N:integer; var YR, YI:VEC; var NZ, NLAST:integer;
                var FNUL, TOL, ELIM, ALIM:double);
{***BEGIN PROLOGUE  ZUNI1
!***REFER TO  ZBESI,ZBESK
!
!     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
!     EXPANSION FOR I(FNU,Z) IN -PI/3 <= ARG Z <= PI/3.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST <> 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I):=CZERO FOR I=NLAST+1,N.
!
!***ROUTINES CALLED  ZUCHK,ZUNIK,ZUOIK,D1MACH,ZABS
!***END PROLOGUE  ZUNI1
!     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1,
!    *S2,Y,Z,ZETA1,ZETA2 }
Label 10,20,30,40,50,60,70,90,100,110,120,130,Return;
Var
      APHI, ASCLE, CONEI, CONER, CRSC, CSCL, C1R, C2I, C2M, C2R, FN,
      PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI, SUMR, S1I, S1R,
      S2I, S2R, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R: Double;
      I, IFLAG, INIT, K, M, ND, NN, NUF, NW: Integer;
      CWRKR, CWRKI: VEC16;
      BRY, CSSR, CSRR, CYR, CYI: VEC;
Begin

      ZEROR:=0.0; ZEROI:=0.0; CONER:=1.0; CONEI:=0.0;

      NZ := 0;
      ND := N;
      NLAST := 0;
{-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
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
{-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!----------------------------------------------------------------------}
      FN := DMAX(FNU,1.0);
      INIT := 0;
      ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R,
            ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI);
      IF KODE = 1 THEN GOTO 10;
      STR := ZR + ZETA2R;
      STI := ZI + ZETA2I;
      RAST := FN/ZABS(STR,STI);
      STR := STR*RAST*RAST;
      STI := -STI*RAST*RAST;
      S1R := -ZETA1R + STR;
      S1I := -ZETA1I + STI;
      GOTO 20;
10:   S1R := -ZETA1R + ZETA2R;
      S1I := -ZETA1I + ZETA2I;
20:   RS1 := S1R;
      IF ABS(RS1) > ELIM THEN GOTO 130;
30:   NN := IMIN(2,ND);
      For I:=1 to NN do
      begin
        FN := FNU + 1.0*(ND-I);
        INIT := 0;
        ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R,
              ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI);
        IF KODE = 1 THEN GOTO 40;
        STR := ZR + ZETA2R;
        STI := ZI + ZETA2I;
        RAST := FN/ZABS(STR,STI);
        STR := STR*RAST*RAST;
        STI := -STI*RAST*RAST;
        S1R := -ZETA1R + STR;
        S1I := -ZETA1I + STI + ZI;
        GOTO 50;
40:     S1R := -ZETA1R + ZETA2R;
        S1I := -ZETA1I + ZETA2I;
{-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!----------------------------------------------------------------------}
50:     RS1 := S1R;
        IF ABS(RS1) > ELIM THEN GOTO 110;
        IF I = 1 THEN IFLAG := 2;
        IF ABS(RS1) < ALIM THEN GOTO 60;
{-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!----------------------------------------------------------------------}
        APHI := ZABS(PHIR,PHII);
        RS1 := RS1 + Ln(APHI);
        IF ABS(RS1) > ELIM THEN GOTO 110;
        IF I = 1 THEN IFLAG := 1;
        IF RS1 < 0.0 THEN GOTO 60;
        IF I = 1 THEN IFLAG := 3;
{-----------------------------------------------------------------------
!     SCALE S1 IF CABS(S1) < ASCLE
!----------------------------------------------------------------------}
60:     S2R := PHIR*SUMR - PHII*SUMI;
        S2I := PHIR*SUMI + PHII*SUMR;
        STR := EXP(S1R)*CSSR[IFLAG];
        S1R := STR*COS(S1I);
        S1I := STR*SIN(S1I);
        STR := S2R*S1R - S2I*S1I;
        S2I := S2R*S1I + S2I*S1R;
        S2R := STR;
        IF IFLAG <> 1 THEN GOTO 70;
        ZUCHK(S2R, S2I, NW, BRY[1], TOL);
        IF NW <> 0 THEN GOTO 110;
70:     CYR[I] := S2R;
        CYI[I] := S2I;
        M := ND - I + 1;
        YR[M] := S2R*CSRR[IFLAG];
        YI[M] := S2I*CSRR[IFLAG]
      end; {i loop}
      IF ND <= 2 THEN GOTO 100;
      RAST := 1.0/ZABS(ZR,ZI);
      STR := ZR*RAST;
      STI := -ZI*RAST;
      RZR := (STR+STR)*RAST;
      RZI := (STI+STI)*RAST;
      BRY[2] := 1.0/BRY[1];
      BRY[3] := D1MACH(2);
      S1R := CYR[1];
      S1I := CYI[1];
      S2R := CYR[2];
      S2I := CYI[2];
      C1R := CSRR[IFLAG];
      ASCLE := BRY[IFLAG];
      K := ND - 2;
      FN := 1.0*K;
      For I:=3 to ND do
      begin
        C2R := S2R;
        C2I := S2I;
        S2R := S1R + (FNU+FN)*(RZR*C2R-RZI*C2I);
        S2I := S1I + (FNU+FN)*(RZR*C2I+RZI*C2R);
        S1R := C2R;
        S1I := C2I;
        C2R := S2R*C1R;
        C2I := S2I*C1R;
        YR[K] := C2R;
        YI[K] := C2I;
        K := K - 1;
        FN := FN - 1.0;
        IF IFLAG >= 3 THEN GOTO 90;
        STR := ABS(C2R);
        STI := ABS(C2I);
        C2M := DMAX(STR,STI);
        IF C2M <= ASCLE THEN GOTO 90;
        IFLAG := IFLAG + 1;
        ASCLE := BRY[IFLAG];
        S1R := S1R*C1R;
        S1I := S1I*C1R;
        S2R := C2R;
        S2I := C2I;
        S1R := S1R*CSSR[IFLAG];
        S1I := S1I*CSSR[IFLAG];
        S2R := S2R*CSSR[IFLAG];
        S2I := S2I*CSSR[IFLAG];
        C1R := CSRR[IFLAG];
90:   end; {i loop}
100:  GOTO RETURN;
{-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!----------------------------------------------------------------------}
110:  IF RS1 > 0.0 THEN GOTO 120;
      YR[ND] := ZEROR;
      YI[ND] := ZEROI;
      NZ := NZ + 1;
      ND := ND - 1;
      IF ND = 0 THEN GOTO 100;
      ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM);
      IF NUF < 0 THEN GOTO 120;
      ND := ND - NUF;
      NZ := NZ + NUF;
      IF ND = 0 THEN GOTO 100;
      FN := FNU + 1.0*(ND-1);
      IF FN >= FNUL THEN GOTO 30;
      NLAST := ND;
      GOTO RETURN;
120:  NZ := -1;
      GOTO RETURN;
130:  IF RS1 > 0.0 THEN GOTO 120;
      NZ := N;
      For I:=1 to N do
      begin
        YR[I] := ZEROR;
        YI[I] := ZEROI
      end;
Return: End; {ZUNI1}


Procedure ZUNI2(ZR, ZI, FNU:double; KODE, N:integer; var YR, YI:VEC; var NZ, NLAST:integer;
                var FNUL, TOL, ELIM, ALIM:double);
{***BEGIN PROLOGUE  ZUNI2
!***REFER TO  ZBESI,ZBESK
!
!     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
!     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
!     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST <> 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N.
!
!***ROUTINES CALLED  ZAIRY,ZUCHK,ZUNHJ,ZUOIK,D1MACH,ZABS
!***END PROLOGUE  ZUNI2
!     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS,
!     CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN
!Note: variable IN ==> IN0 (IN is a reserved word in Pascal). }
Label 10,20,30,40,50,60,70,80, 100,110,120,130,140,150,Return;
Var
      AARG, AIC, AII, AIR, ANG, APHI, ARGI, ARGR, ASCLE, ASUMI, ASUMR,
      BSUMI, BSUMR, CIDI, CONEI, CONER, CRSC, CSCL, C1R, C2I, C2M, C2R,
      DAII, DAIR, FN, HPI, PHII, PHIR, RAST, RAZ, RS1, RZI, RZR, STI,
      STR, S1I, S1R, S2I, S2R, ZBI, ZBR, ZEROI, ZEROR, ZETA1I, ZETA1R,
      ZETA2I, ZETA2R, ZNI, ZNR: Double;
      I, IFLAG, IN0, INU, J, K, NAI, ND, NDAI, NN, NUF, NW, IDUM: Integer;
      BRY, CIPR, CIPI, CSSR, CSRR, CYR, CYI: VEC;
Begin
      ZEROR:=0.0; ZEROI:=0.0; CONER:=1.0; CONEI:=0.0;
      CIPR[1]:= 1.0; CIPI[1]:=0.0; CIPR[2]:=0.0; CIPI[2]:= 1.0;
      CIPR[3]:=-1.0; CIPI[3]:=0.0; CIPR[4]:=0.0; CIPI[4]:=-1.0;

      HPI:=1.57079632679489662; AIC:=1.265512123484645396;

      NZ := 0;
      ND := N;
      NLAST := 0;
{-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
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
{----------------------------------------------------------------------
!     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
!----------------------------------------------------------------------}
      ZNR := ZI;
      ZNI := -ZR;
      ZBR := ZR;
      ZBI := ZI;
      CIDI := -CONER;
      INU := Round(FNU);
      ANG := HPI*(FNU-1.0*INU);
      C2R := COS(ANG);
      C2I := SIN(ANG);
      IN0 := INU + N - 1;
      IN0 := (IN0 Mod 4) + 1;
      STR := C2R*CIPR[IN0] - C2I*CIPI[IN0];
      C2I := C2R*CIPI[IN0] + C2I*CIPR[IN0];
      C2R := STR;
      IF ZI > 0.0 THEN GOTO 10;
      ZNR := -ZNR;
      ZBI := -ZBI;
      CIDI := -CIDI;
      C2I := -C2I;
{-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!----------------------------------------------------------------------}
10:   FN := DMAX(FNU,1.0);
      ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
            ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI);
      IF KODE = 1 THEN GOTO 20;
      STR := ZBR + ZETA2R;
      STI := ZBI + ZETA2I;
      RAST := FN/ZABS(STR,STI);
      STR := STR*RAST*RAST;
      STI := -STI*RAST*RAST;
      S1R := -ZETA1R + STR;
      S1I := -ZETA1I + STI;
      GOTO 30;
20:   S1R := -ZETA1R + ZETA2R;
      S1I := -ZETA1I + ZETA2I;
30:   RS1 := S1R;
      IF ABS(RS1) > ELIM THEN GOTO 150;
40:   NN := IMIN(2,ND);
      For I:=1 to NN do
      begin
        FN := FNU + 1.0*(ND-I);
        ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI,
              ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI);
        IF KODE = 1 THEN GOTO 50;
        STR := ZBR + ZETA2R;
        STI := ZBI + ZETA2I;
        RAST := FN/ZABS(STR,STI);
        STR := STR*RAST*RAST;
        STI := -STI*RAST*RAST;
        S1R := -ZETA1R + STR;
        S1I := -ZETA1I + STI + ABS(ZI);
        GOTO 60;
50:     S1R := -ZETA1R + ZETA2R;
        S1I := -ZETA1I + ZETA2I;
{-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!----------------------------------------------------------------------}
60:     RS1 := S1R;
        IF ABS(RS1) > ELIM THEN GOTO 120;
        IF I = 1 THEN IFLAG := 2;
        IF ABS(RS1) < ALIM THEN GOTO 70;
{-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!----------------------------------------------------------------------}
        APHI := ZABS(PHIR,PHII);
        AARG := ZABS(ARGR,ARGI);
        RS1 := RS1 + Ln(APHI) - 0.25*Ln(AARG) - AIC;
        IF ABS(RS1) > ELIM THEN GOTO 120;
        IF I = 1 THEN IFLAG := 1;
        IF RS1 < 0.0 THEN GOTO 70;
        IF I = 1 THEN IFLAG := 3;
{-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!----------------------------------------------------------------------}
70:     ZAIRY(ARGR, ARGI, 0, 2, AIR, AII, NAI, IDUM);
        ZAIRY(ARGR, ARGI, 1, 2, DAIR, DAII, NDAI, IDUM);
        STR := DAIR*BSUMR - DAII*BSUMI;
        STI := DAIR*BSUMI + DAII*BSUMR;
        STR := STR + (AIR*ASUMR-AII*ASUMI);
        STI := STI + (AIR*ASUMI+AII*ASUMR);
        S2R := PHIR*STR - PHII*STI;
        S2I := PHIR*STI + PHII*STR;
        STR := EXP(S1R)*CSSR[IFLAG];
        S1R := STR*COS(S1I);
        S1I := STR*SIN(S1I);
        STR := S2R*S1R - S2I*S1I;
        S2I := S2R*S1I + S2I*S1R;
        S2R := STR;
        IF IFLAG <> 1 THEN GOTO 80;
        ZUCHK(S2R, S2I, NW, BRY[1], TOL);
        IF NW <> 0 THEN GOTO 120;
80:     IF ZI <= 0.0 THEN S2I := -S2I;
        STR := S2R*C2R - S2I*C2I;
        S2I := S2R*C2I + S2I*C2R;
        S2R := STR;
        CYR[I] := S2R;
        CYI[I] := S2I;
        J := ND - I + 1;
        YR[J] := S2R*CSRR[IFLAG];
        YI[J] := S2I*CSRR[IFLAG];
        STR := -C2I*CIDI;
        C2I := C2R*CIDI;
        C2R := STR
      end; {I loop}
      IF ND <= 2 THEN GOTO 110;
      RAZ := 1.0/ZABS(ZR,ZI);
      STR := ZR*RAZ;
      STI := -ZI*RAZ;
      RZR := (STR+STR)*RAZ;
      RZI := (STI+STI)*RAZ;
      BRY[2] := 1.0/BRY[1];
      BRY[3] := D1MACH(2);
      S1R := CYR[1];
      S1I := CYI[1];
      S2R := CYR[2];
      S2I := CYI[2];
      C1R := CSRR[IFLAG];
      ASCLE := BRY[IFLAG];
      K := ND - 2;
      FN := 1.0*K;
      For I:=3 to ND do
      begin
        C2R := S2R;
        C2I := S2I;
        S2R := S1R + (FNU+FN)*(RZR*C2R-RZI*C2I);
        S2I := S1I + (FNU+FN)*(RZR*C2I+RZI*C2R);
        S1R := C2R;
        S1I := C2I;
        C2R := S2R*C1R;
        C2I := S2I*C1R;
        YR[K] := C2R;
        YI[K] := C2I;
        K := K - 1;
        FN := FN - 1.0;
        IF IFLAG >= 3 THEN GOTO 100;
        STR := ABS(C2R);
        STI := ABS(C2I);
        C2M := DMAX(STR,STI);
        IF C2M <= ASCLE THEN GOTO 100;
        IFLAG := IFLAG + 1;
        ASCLE := BRY[IFLAG];
        S1R := S1R*C1R;
        S1I := S1I*C1R;
        S2R := C2R;
        S2I := C2I;
        S1R := S1R*CSSR[IFLAG];
        S1I := S1I*CSSR[IFLAG];
        S2R := S2R*CSSR[IFLAG];
        S2I := S2I*CSSR[IFLAG];
        C1R := CSRR[IFLAG];
100:  end;
110:  GOTO RETURN;
120:  IF RS1 > 0.0 THEN GOTO 140;
{-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!----------------------------------------------------------------------}
      YR[ND] := ZEROR;
      YI[ND] := ZEROI;
      NZ := NZ + 1;
      ND := ND - 1;
      IF ND = 0 THEN GOTO 110;
      ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM);
      IF NUF < 0 THEN GOTO 140;
      ND := ND - NUF;
      NZ := NZ + NUF;
      IF ND = 0 THEN GOTO 110;
      FN := FNU + 1.0*(ND-1);
      IF FN < FNUL THEN GOTO 130;
      FN := CIDI;
      J := NUF + 1;
      K := (J Mod 4) + 1;
      S1R := CIPR[K];
      S1I := CIPI[K];
      IF FN < 0.0 THEN S1I := -S1I;
      STR := C2R*S1R - C2I*S1I;
      C2I := C2R*S1I + C2I*S1R;
      C2R := STR;
      GOTO 40;
130:  NLAST := ND;
      GOTO RETURN;
140:  NZ := -1;
      GOTO RETURN;
150:  IF RS1 < 0.0 THEN GOTO 140;
      NZ := N;
      For I:=1 to N do
      begin
        YR[I] := ZEROR;
        YI[I] := ZEROI
      end;
Return:End; {ZUNI2}

END.

{end of file cbess2.pas}