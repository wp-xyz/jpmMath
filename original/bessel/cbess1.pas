{*************************************************************************
*    Procedures and Functions used By programs TZBESJ, TZBESK, TZBESY    *
*    (Evalute Bessel Functions with complex argument, 1st to 3rd kind)   *
* ---------------------------------------------------------------------- *
* Reference:  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES, 1983.       *
*                                                                        *
*                      Pascal Release By J-P Moreau, Paris (07/01/2005). *
*                                      (www.jpmoreau.fr)                 *
*************************************************************************}
UNIT CBESS1;    {ZBINU}

INTERFACE

Uses Complex, CBess0, CBess00, CBess2;

  Procedure ZBINU(ZR, ZI, FNU:Double; KODE, N: Integer; Var CYR, CYI: VEC;
                  Var NZ:Integer; Var RL, FNUL, TOL, ELIM, ALIM: Double);


IMPLEMENTATION

Procedure ZWRSK(ZRR, ZRI, FNU:double; KODE, N:Integer; Var YR, YI:VEC;
                Var NZ: Integer; Var CWR, CWI: VEC; Var TOL, ELIM, ALIM:double);
                Forward;

Procedure ZBINU(ZR, ZI, FNU:Double; KODE, N: Integer; Var CYR, CYI: VEC;
                Var NZ:Integer; Var RL, FNUL, TOL, ELIM, ALIM: Double);
{***BEGIN PROLOGUE  ZBINU
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZAIRY,ZBIRY

!   ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE

!***ROUTINES CALLED  ZABS,ZASYI,ZBUNI,ZMLRI,ZSERI,ZUOIK,ZWRSK
!***END PROLOGUE  ZBINU }
Label 10,20,30,40,50,60,70,80, 100,110,120,130, Return;
Var
      AZ, DFNU, ZEROI, ZEROR: double;
      I, INW, NLAST, NN, NUI, NW: Integer;
      CWR, CWI: VEC;
Begin

      ZEROR:=0.0; ZEROI:=0.0;

      NZ := 0;
      AZ := ZABS(ZR,ZI);
      NN := N;
      DFNU := FNU + 1.0*(N-1);
      IF AZ <= 2.0 THEN GOTO 10;
      IF AZ*AZ*0.25 > DFNU+1.0 THEN GOTO 20;
{-----------------------------------------------------------------------
!     POWER SERIES
!----------------------------------------------------------------------}
10:   ZSERI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM);
      INW := ABS(NW);
      NZ := NZ + INW;
      NN := NN - INW;
      IF NN = 0 THEN GOTO RETURN;
      IF NW >= 0 THEN GOTO 120;
      DFNU := FNU + 1.0*(NN-1);
20:   IF AZ < RL THEN GOTO 40;
      IF DFNU <= 1.0 THEN GOTO 30;
      IF AZ+AZ < DFNU*DFNU THEN GOTO 50;
{-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z
!----------------------------------------------------------------------}
30:   ZASYI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, RL, TOL, ELIM, ALIM);
      IF NW < 0 THEN GOTO 130;
      GOTO 120;
40:   IF DFNU <= 1.0 THEN GOTO 70;
{-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
!----------------------------------------------------------------------}
50:   ZUOIK(ZR, ZI, FNU, KODE, 1, NN, CYR, CYI, NW, TOL, ELIM, ALIM);
      IF NW < 0 THEN GOTO 130;
      NZ := NZ + NW;
      NN := NN - NW;
      IF NN = 0 THEN GOTO RETURN;
      DFNU := FNU+1.0*(NN-1);
      IF DFNU > FNUL THEN GOTO 110;
      IF AZ > FNUL THEN GOTO 110;
60:   IF AZ > RL THEN GOTO 80;
{-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES
!----------------------------------------------------------------------}
70:   ZMLRI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL);
      IF NW < 0 THEN GOTO 130;
      GOTO 120;
{-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
!----------------------------------------------------------------------}
80:   ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM, ALIM);
      IF NW >= 0 THEN GOTO 100;
      NZ := NN;
      For I:=1 to NN do
      begin
        CYR[I] := ZEROR;
        CYI[I] := ZEROI
      end;
      GOTO RETURN;
100:  IF NW > 0 THEN GOTO 130;
      ZWRSK(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL, ELIM, ALIM);
      IF NW < 0 THEN GOTO 130;
      GOTO 120;
{-----------------------------------------------------------------------
!     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
!----------------------------------------------------------------------}
110:  NUI := Round(FNUL-DFNU) + 1;
      NUI := IMAX(NUI,0);
      ZBUNI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL,TOL,ELIM,ALIM);
      IF NW < 0 THEN GOTO 130;
      NZ := NZ + NW;
      IF NLAST = 0 THEN GOTO 120;
      NN := NLAST;
      GOTO 60;
120:  GOTO RETURN;
130:  NZ := -1;
      IF NW = -2 THEN NZ:=-2;
Return: End; {ZBINU}


Procedure ZWRSK(ZRR, ZRI, FNU:double; KODE, N:Integer; Var YR, YI:VEC;
                Var NZ: Integer; Var CWR, CWI: VEC; Var TOL, ELIM, ALIM:double);
{***BEGIN PROLOGUE  ZWRSK
!***REFER TO  ZBESI,ZBESK
!
!     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
!     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
!
!***ROUTINES CALLED  D1MACH,ZBKNU,ZRATI,ZABS
!***END PROLOGUE  ZWRSK
!     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR }
Label 10,20,30,50, Return;
Var
      ACT, ACW, ASCLE, CINUI, CINUR, CSCLR, CTI, CTR, C1I, C1R,
      C2I, C2R, PTI, PTR, RACT, STI, STR: Double;
      I, NW: Integer;
Begin
{-----------------------------------------------------------------------
!     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
!     Y(I):=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
!     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
!----------------------------------------------------------------------}
      NZ := 0;
      ZBKNU(ZRR, ZRI, FNU, KODE, 2, CWR, CWI, NW, TOL, ELIM, ALIM);
      IF NW <> 0 THEN GOTO 50;
      ZRATI(ZRR, ZRI, FNU, N, YR, YI, TOL);
{-----------------------------------------------------------------------
!     RECUR FORWARD ON I(FNU+1,Z) := R(FNU,Z)*I(FNU,Z),
!     R(FNU+J-1,Z):=Y(J),  J:=1,...,N
!----------------------------------------------------------------------}
      CINUR := 1.0;
      CINUI := 0.0;
      IF KODE = 1 THEN GOTO 10;
      CINUR := COS(ZRI);
      CINUI := SIN(ZRI);
{-----------------------------------------------------------------------
!     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
!     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
!     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
!     THE RESULT IS ON SCALE.
!----------------------------------------------------------------------}
10:   ACW := ZABS(CWR[2],CWI[2]);
      ASCLE := 1000*D1MACH(1)/TOL;
      CSCLR := 1.0;
      IF ACW > ASCLE THEN GOTO 20;
      CSCLR := 1.0/TOL;
      GOTO 30;
20:   ASCLE := 1.0/ASCLE;
      IF ACW < ASCLE THEN GOTO 30;
      CSCLR := TOL;
30:   C1R := CWR[1]*CSCLR;
      C1I := CWI[1]*CSCLR;
      C2R := CWR[2]*CSCLR;
      C2I := CWI[2]*CSCLR;
      STR := YR[1];
      STI := YI[1];
{-----------------------------------------------------------------------
!     CINU:=CINU*(CONJG(CT)/CABS(CT))*(1.0D0/CABS(CT) PREVENTS
!     UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
!----------------------------------------------------------------------}
      PTR := STR*C1R - STI*C1I;
      PTI := STR*C1I + STI*C1R;
      PTR := PTR + C2R;
      PTI := PTI + C2I;
      CTR := ZRR*PTR - ZRI*PTI;
      CTI := ZRR*PTI + ZRI*PTR;
      ACT := ZABS(CTR,CTI);
      RACT := 1.0/ACT;
      CTR := CTR*RACT;
      CTI := -CTI*RACT;
      PTR := CINUR*RACT;
      PTI := CINUI*RACT;
      CINUR := PTR*CTR - PTI*CTI;
      CINUI := PTR*CTI + PTI*CTR;
      YR[1] := CINUR*CSCLR;
      YI[1] := CINUI*CSCLR;
      IF N = 1 THEN GOTO RETURN;
      For I:=2 to N do
      begin
        PTR := STR*CINUR - STI*CINUI;
        CINUI := STR*CINUI + STI*CINUR;
        CINUR := PTR;
        STR := YR[I];
        STI := YI[I];
        YR[I] := CINUR*CSCLR;
        YI[I] := CINUI*CSCLR
      end;
      GOTO RETURN;
50:   NZ := -1;
      IF NW = -2 THEN NZ:=-2;
Return: End; {ZWRSK}

END.

{end of file Cbess1.pas}