{*************************************************************************
*    Procedures and Functions used By programs TZBESJ, TZBESK, TZBESY    *
*    (Evalute Bessel Functions with complex argument, 1st to 3rd kind)   *
* ---------------------------------------------------------------------- *
* Reference:  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES, 1983.       *
*                                                                        *
*                      Pascal Release By J-P Moreau, Paris (07/01/2005). *
*                                     (www.jpmoreau.fr)                  *
*************************************************************************}
UNIT CBESS00;   {procedures ZAIRY, ZBKNU}

INTERFACE

Uses Complex, Cbess0;

  Procedure ZAIRY(ZR, ZI: double; ID, KODE: Integer; Var AIR, AII:Double;
                  Var NZ, IERR: Integer);

  Procedure ZBKNU(ZR, ZI, FNU:double; KODE, N:Integer; Var YR, YI: VEC;
                  Var NZ: Integer; TOL, ELIM, ALIM:Double);


IMPLEMENTATION

Procedure ZACAI(ZR, ZI, FNU:double; KODE, MR, N:Integer; Var YR, YI:VEC;
                Var NZ: Integer; Var RL, TOL, ELIM, ALIM:Double); Forward;

Procedure ZKSCL(ZRR,ZRI,FNU:Double; N:Integer; Var YR,YI:VEC; Var NZ:Integer;
                Var RZR,RZI,ASCLE,TOL,ELIM: Double); Forward;

Procedure ZAIRY(ZR, ZI: double; ID, KODE: Integer; Var AIR, AII:Double;
                Var NZ, IERR: Integer);
{***BEGIN PROLOGUE  ZAIRY
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  830501   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
!***DESCRIPTION
!
!                      ***A DOUBLE PRECISION ROUTINE***
!         ON KODE=1, ZAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
!         ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY. ON
!         KODE=2, A SCALING OPTION CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*
!         DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
!         -PI/3 < ARG(Z) < PI/3 AND THE EXPONENTIAL GROWTH IN
!         PI/3 < ABS(ARG(Z)) < PI, WHERE ZTA:=(2/3)*Z*CSQRT(Z).
!
!         WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
!         THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
!         FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
!         DEFINTIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
!         MATHEMATICAL FUNCTIONS (REF. 1).
!
!         INPUT      ZR,ZI ARE DOUBLE PRECISION
!           ZR,ZI  - Z=CMPLX(ZR,ZI)
!           ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
!           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!                    KODE= 1  RETURNS
!                             AI=AI(Z)                 ON ID=0 OR
!                             AI=DAI(Z)/DZ             ON ID=1
!                        = 2  RETURNS
!                             AI=CEXP(ZTA)*AI(Z)       ON ID=0 OR
!                             AI=CEXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
!                             ZTA=(2/3)*Z*CSQRT(Z)
!
!         OUTPUT     AIR,AII ARE DOUBLE PRECISION
!           AIR,AII- COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND
!                    KODE
!           NZ     - UNDERFLOW INDICATOR
!                    NZ= 0   , NORMAL RETURN
!                    NZ= 1   , AI=CMPLX(0.0,0.0) DUE TO UNDERFLOW IN
!                              -PI/3 < ARG(Z) < PI/3 ON KODE=1.
!           IERR   - ERROR FLAG
!                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!                    IERR=1, INPUT ERROR   - NO COMPUTATION
!                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
!                            TOO LARGE ON KODE:=1
!                    IERR=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
!                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
!                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
!                    IERR=4, CABS(Z) TOO LARGE  - NO COMPUTATION
!                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
!                            REDUCTION
!                    IERR=5, ERROR              - NO COMPUTATION,
!                            ALGORITHM TERMINATION CONDITION NOT MET
!
!***LONG DESCRIPTION
!
!         AI AND DAI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE K BESSEL
!         FUNCTIONS BY
!
!            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
!                           C=1.0/(PI*SQRT(3.0))
!                           ZTA=(2/3)*Z^(3/2)
!
!         WITH THE POWER SERIES FOR CABS(Z) <= 1.
!
!         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
!         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
!         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
!         THE MAGNITUDE OF ZETA=(2/3)*Z^1.5 EXCEEDS U1=SQRT(0.5/UR),
!         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
!         FLAG IERR:=3 IS TRIGGERED WHERE UR:=DMAX(D1MACH(4),1.0D-18) IS
!         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
!         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2:=0.5/UR, THEN
!         ALL SIGNIFICANCE IS LOST AND IERR:=4. IN ORDER TO USE THE INT
!         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
!         LARGEST INTEGER, U3:=I1MACH(9). THUS, THE MAGNITUDE OF ZETA
!         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
!         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
!         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
!         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
!         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
!         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
!         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
!         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
!         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
!         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
!         MACHINES.
!
!         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
!         BESSEL FUNCTION CAN BE EXPRESSED BY P*10^S WHERE P=MAX(UNIT
!         ROUNDOFF,1E-18) IS THE NOMINAL PRECISION AND 10^S REPRE-
!         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
!         ELEMENTARY FUNCTIONS. HERE, S:=MAX(1,ABS(LOG10(CABS(Z))),
!         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S:=MAX(1,ABS(EXPONENT OF
!         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
!         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
!         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
!         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10^K LARGER
!         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
!         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
!         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
!         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
!         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
!         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
!         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
!         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
!         OR -PI/2+P.
!
!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
!                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
!                 COMMERCE, 1955.
!
!               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
!
!               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
!                 1018, MAY, 1985
!
!               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
!                 MATH. SOFTWARE, 1986
!
!***ROUTINES CALLED  ZACAI,ZBKNU,ZEXP,ZSQRT,I1MACH,D1MACH
!***END PROLOGUE  ZAIRY
!     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3 }
Label 40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,260,270,280,Return;
Var
      AA, AD, AK, ALIM, ATRM, AZ, AZ3, BK,
      CC, CK, COEF, CONEI, CONER, CSQI, CSQR, C1, C2, DIG,
      DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR,
      S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TTH, ZEROI,
      ZEROR, ZTAI, ZTAR, Z3I, Z3R, ALAZ, BB: Double;
      IFLAG, K, K1, K2, MR, NN: Integer;
      CYI,CYR: VEC;
Begin

      TTH := 6.66666666666666667E-01;
      C1  := 3.55028053887817240E-01;
      C2  := 2.58819403792806799E-01; 
      COEF:= 1.83776298473930683E-01;
      ZEROR:=0.0; ZEROI:=0.0; CONER:=1.0; CONEI:=0.0;
      IERR := 0;
      NZ:=0;
      IF (ID < 0) OR (ID >1) THEN IERR:=1;
      IF (KODE < 1) OR (KODE > 2) THEN IERR:=1;
      IF IERR <> 0 THEN GOTO RETURN;
      AZ := ZABS(ZR,ZI);
      TOL := DMAX(D1MACH(4),1E-18);
      FID := 1.0*ID;
      IF AZ > 1.0 THEN GOTO 70;
{-----------------------------------------------------------------------
!     POWER SERIES FOR CABS(Z).LE.1.
!----------------------------------------------------------------------}
      S1R := CONER;
      S1I := CONEI;
      S2R := CONER;
      S2I := CONEI;
      IF AZ < TOL THEN GOTO 170;
      AA := AZ*AZ;
      IF AA < TOL/AZ THEN GOTO 40;
      TRM1R := CONER;
      TRM1I := CONEI;
      TRM2R := CONER;
      TRM2I := CONEI;
      ATRM := 1.0;
      STR := ZR*ZR - ZI*ZI;
      STI := ZR*ZI + ZI*ZR;
      Z3R := STR*ZR - STI*ZI;
      Z3I := STR*ZI + STI*ZR;
      AZ3 := AZ*AA;
      AK := 2.0 + FID;
      BK := 3.0 - FID - FID;
      CK := 4.0 - FID;
      DK := 3.0 + FID + FID;
      D1 := AK*DK;
      D2 := BK*CK;
      AD := DMIN(D1,D2);
      AK := 24.0 + 9.0*FID;
      BK := 30.0 - 9.0*FID;
      For K:=1 to 25 do
      begin
        STR := (TRM1R*Z3R-TRM1I*Z3I)/D1;
        TRM1I := (TRM1R*Z3I+TRM1I*Z3R)/D1;
        TRM1R := STR;
        S1R := S1R + TRM1R;
        S1I := S1I + TRM1I;
        STR := (TRM2R*Z3R-TRM2I*Z3I)/D2;
        TRM2I := (TRM2R*Z3I+TRM2I*Z3R)/D2;
        TRM2R := STR;
        S2R := S2R + TRM2R;
        S2I := S2I + TRM2I;
        ATRM := ATRM*AZ3/AD;
        D1 := D1 + AK;
        D2 := D2 + BK;
        AD := DMIN(D1,D2);
        IF ATRM < TOL*AD THEN GOTO 40;
        AK := AK + 18.0;
        BK := BK + 18.0
      end;
40:   IF ID = 1 THEN GOTO 50;
      AIR := S1R*C1 - C2*(ZR*S2R-ZI*S2I);
      AII := S1I*C1 - C2*(ZR*S2I+ZI*S2R);
      IF KODE = 1 THEN GOTO RETURN;
      ZSQRT(ZR, ZI, STR, STI);
      ZTAR := TTH*(ZR*STR-ZI*STI);
      ZTAI := TTH*(ZR*STI+ZI*STR);
      ZEXP(ZTAR, ZTAI, STR, STI);
      PTR := AIR*STR - AII*STI;
      AII := AIR*STI + AII*STR;
      AIR := PTR;
      GOTO RETURN;
50:   AIR := -S2R*C2;
      AII := -S2I*C2;
      IF AZ <= TOL THEN GOTO 60;
      STR := ZR*S1R - ZI*S1I;
      STI := ZR*S1I + ZI*S1R;
      CC := C1/(1.0+FID);
      AIR := AIR + CC*(STR*ZR-STI*ZI);
      AII := AII + CC*(STR*ZI+STI*ZR);
60:   IF KODE = 1 THEN GOTO RETURN;
      ZSQRT(ZR, ZI, STR, STI);
      ZTAR := TTH*(ZR*STR-ZI*STI);
      ZTAI := TTH*(ZR*STI+ZI*STR);
      ZEXP(ZTAR, ZTAI, STR, STI);
      PTR := STR*AIR - STI*AII;
      AII := STR*AII + STI*AIR;
      AIR := PTR;
      GOTO RETURN;
{-----------------------------------------------------------------------
!     CASE FOR CABS(Z).GT.1.0
!----------------------------------------------------------------------}
70:   FNU := (1.0+FID)/3.0;
{-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM).LT.EXP(-ALIM):=EXP(-ELIM)/TOL    AND
!     EXP(ELIM).GT.EXP(ALIM):=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG := NUMBER OF BASE 10 DIGITS IN TOL := 10**(-DIG).
!----------------------------------------------------------------------}
      K1 := I1MACH(15);
      K2 := I1MACH(16);
      R1M5 := D1MACH(5);
      K := IMIN(ABS(K1),ABS(K2));
      ELIM := 2.303*(K*R1M5-3.0);
      K1 := I1MACH(14) - 1;
      AA := R1M5*K1;
      DIG := DMIN(AA,18.0);
      AA := AA*2.303;
      ALIM := ELIM + DMAX(-AA,-41.45);
      RL := 1.2*DIG + 3.0;
      ALAZ := Ln(AZ);
{-----------------------------------------------------------------------
!     TEST FOR PROPER RANGE
!----------------------------------------------------------------------}
      AA:=0.5/TOL;
      BB:=0.5*I1MACH(9);
      AA:=DMIN(AA,BB);
      AA:=Power(AA,TTH);
      IF AZ > AA THEN GOTO 260;
      AA:=SQRT(AA);
      IF AZ > AA THEN IERR:=3;
      ZSQRT(ZR, ZI, CSQR, CSQI);
      ZTAR := TTH*(ZR*CSQR-ZI*CSQI);
      ZTAI := TTH*(ZR*CSQI+ZI*CSQR);
{-----------------------------------------------------------------------
!     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
!----------------------------------------------------------------------}
      IFLAG := 0;
      SFAC := 1.0;
      AK := ZTAI;
      IF ZR >= 0.0 THEN GOTO 80;
      BK := ZTAR;
      CK := -ABS(BK);
      ZTAR := CK;
      ZTAI := AK;
80:   IF ZI <> 0.0 THEN GOTO 90;
      IF ZR > 0.0 THEN GOTO 90;
      ZTAR := 0.0;
      ZTAI := AK;
90:   AA := ZTAR;
      IF (AA >= 0.0) AND (ZR > 0.0) THEN GOTO 110;
      IF KODE = 2 THEN GOTO 100;
{-----------------------------------------------------------------------
!     OVERFLOW TEST
!----------------------------------------------------------------------}
      IF AA > -ALIM THEN GOTO 100;
      AA := -AA + 0.25*ALAZ;
      IFLAG := 1;
      SFAC := TOL;
      IF AA > ELIM THEN GOTO 270;
{-----------------------------------------------------------------------
!     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE:=2
!----------------------------------------------------------------------}
100:  MR := 1;
      IF ZI < 0.0 THEN MR := -1;
      ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL, ELIM, ALIM);
      IF NN < 0 THEN GOTO 280;
      NZ := NZ + NN;
      GOTO 130;
110:  IF KODE = 2 THEN GOTO 120;
{-----------------------------------------------------------------------
!     UNDERFLOW TEST
!----------------------------------------------------------------------}
      IF AA < ALIM THEN GOTO 120;
      AA := -AA - 0.25*ALAZ;
      IFLAG := 2;
      SFAC := 1.0/TOL;
      IF AA < -ELIM THEN GOTO 210;
120:  ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM, ALIM);
130:  S1R := CYR[1]*COEF;
      S1I := CYI[1]*COEF;
      IF IFLAG <> 0 THEN GOTO 150;
      IF ID = 1 THEN GOTO 140;
      AIR := CSQR*S1R - CSQI*S1I;
      AII := CSQR*S1I + CSQI*S1R;
      GOTO RETURN;
140:  AIR := -(ZR*S1R-ZI*S1I);
      AII := -(ZR*S1I+ZI*S1R);
      GOTO RETURN;
150:  S1R := S1R*SFAC;
      S1I := S1I*SFAC;
      IF ID = 1 THEN GOTO 160;
      STR := S1R*CSQR - S1I*CSQI;
      S1I := S1R*CSQI + S1I*CSQR;
      S1R := STR;
      AIR := S1R/SFAC;
      AII := S1I/SFAC;
      GOTO RETURN;
160:  STR := -(S1R*ZR-S1I*ZI);
      S1I := -(S1R*ZI+S1I*ZR);
      S1R := STR;
      AIR := S1R/SFAC;
      AII := S1I/SFAC;
      GOTO RETURN;
170:  AA := 1000*D1MACH(1);
      S1R := ZEROR;
      S1I := ZEROI;
      IF ID = 1 THEN GOTO 190;
      IF AZ <= AA THEN GOTO 180;
      S1R := C2*ZR;
      S1I := C2*ZI;
180:  AIR := C1 - S1R;
      AII := -S1I;
      GOTO RETURN;
190:  AIR := -C2;
      AII := 0.0;
      AA := SQRT(AA);
      IF AZ <= AA THEN GOTO 200;
      S1R := 0.5*(ZR*ZR-ZI*ZI);
      S1I := ZR*ZI;
200:  AIR := AIR + C1*S1R;
      AII := AII + C1*S1I;
      GOTO RETURN;
210:  NZ := 1;
      AIR := ZEROR;
      AII := ZEROI;
      GOTO RETURN;
270:  NZ := 0;
      IERR:=2;
      GOTO RETURN;
280:  IF NN = -1 THEN GOTO 270;
      NZ:=0;
      IERR:=5;
      GOTO RETURN;
260:  IERR:=4;
      NZ:=0;
Return: End; {ZAIRY}


Procedure ZACAI(ZR, ZI, FNU:double; KODE, MR, N:Integer; Var YR, YI:VEC;
                Var NZ: Integer; Var RL, TOL, ELIM, ALIM:Double);
{***BEGIN PROLOGUE  ZACAI
!***REFER TO  ZAIRY
!
!     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
!     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
!     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
!     IS CALLED FROM ZAIRY.
!
!***ROUTINES CALLED  ZASYI,ZBKNU,ZMLRI,ZSERI,ZS1S2,D1MACH,ZABS
!***END PROLOGUE  ZACAI
!     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY }
Label 10,20,30,40,50,60,70,80,Return;
Var
      ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR, CSPNI, C1R, C1I, C2R, C2I,
      DFNU, FMR, SGN, YY, ZNR, ZNI: Double;
      INU, IUF, NN, NW: Integer;
      CYR, CYI: VEC;
Begin
      NZ := 0;
      ZNR := -ZR;
      ZNI := -ZI;
      AZ := ZABS(ZR,ZI);
      NN := N;
      DFNU := FNU + 1.0*(N-1);
      IF AZ <= 2.0 THEN GOTO 10;
      IF AZ*AZ*0.25 > DFNU+1.0 THEN GOTO 20;
{-----------------------------------------------------------------------
!     POWER SERIES FOR THE I FUNCTION
!----------------------------------------------------------------------}
10:   ZSERI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM);
      GOTO 40;
20:   IF AZ < RL THEN GOTO 30;
{-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
!----------------------------------------------------------------------}
      ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM, ALIM);
      IF NW < 0 THEN GOTO 80;
      GOTO 40;
{-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
!----------------------------------------------------------------------}
30:   ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL);
      IF NW < 0 THEN GOTO 80;
{-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!----------------------------------------------------------------------}
40:   ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM);
      IF NW <> 0 THEN GOTO 80;
      FMR := 1.0*MR;
      SGN := -SIGN(PI,FMR);
      CSGNR := 0.0;
      CSGNI := SGN;
      IF KODE = 1 THEN GOTO 50;
      YY := -ZNI;
      CSGNR := -CSGNI*SIN(YY);
      CSGNI := CSGNI*COS(YY);
{-----------------------------------------------------------------------
!     CALCULATE CSPN:=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!----------------------------------------------------------------------}
50:   INU := Round(FNU);
      ARG := (FNU-1.0*INU)*SGN;
      CSPNR := COS(ARG);
      CSPNI := SIN(ARG);
      IF (INU Mod 2) = 0 THEN GOTO 60;
      CSPNR := -CSPNR;
      CSPNI := -CSPNI;
60:   C1R := CYR[1];
      C1I := CYI[1];
      C2R := YR[1];
      C2I := YI[1];
      IF KODE = 1 THEN GOTO 70;
      IUF := 0;
      ASCLE := 1000*D1MACH(1)/TOL;
      ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF);
      NZ := NZ + NW;
70:   YR[1] := CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I;
      YI[1] := CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R;
      GOTO RETURN;
80:   NZ := -1;
      IF NW = -2 THEN NZ:=-2;
Return:End; {ZACAI}


Procedure ZBKNU(ZR, ZI, FNU:double; KODE, N:Integer; Var YR, YI: VEC;
                Var NZ: Integer; TOL, ELIM, ALIM:Double);
{***BEGIN PROLOGUE  ZBKNU
!***REFER TO  ZBESI,ZBESK,ZAIRY,ZBESH
!
!     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
!
!***ROUTINES CALLED  DGAMLN,I1MACH,D1MACH,ZKSCL,ZSHCH,ZUCHK,ZABS,ZDIV,
!                    ZEXP,ZLOG,ZMLT,ZSQRT
!***END PROLOGUE  ZBKNU }
Label 10,30,40,50,60,70,80,90,100,110,120,130,140,160,170,180,200,210,215,220,
      225,230,240,250,260,261,262,263,264, 270,280,290,300,310,Return;
Var
      AA, AK, ASCLE, A1, A2, BB, BK, CAZ,
      CBI, CBR, CCHI, CCHR, CKI, CKR, COEFI, COEFR, CONEI, CONER,
      CRSCR, CSCLR, CSHI, CSHR, CSI, CSR, CTWOI, CTWOR,
      CZEROI, CZEROR, CZI, CZR, DNU, DNU2, DPI, ETEST, FC, FHS,
      FI, FK, FKS, FMUI, FMUR, FPI, FR, G1, G2, HPI, PI, PR, PTI,
      PTR, P1I, P1R, P2I, P2M, P2R, QI, QR, RAK, RCAZ, RTHPI, RZI,
      RZR, R1, S, SMUI, SMUR, SPI, STI, STR, S1I, S1R, S2I, S2R, TM,
      TTH, T1, T2, ELM, CELMR, ZDR, ZDI, AS, ALAS, HELIM: Double;
      I, IFLAG, INU, K, KFLAG, KK, KMAX, KODED, IDUM, J, IC, INUB, NW: Integer;
      CC, CSSR, CSRR, BRY, CYR, CYI: VEC;
{     COMPLEX: Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH,
      CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK}
Begin

      KMAX:=30;
      CZEROR:=0.0; CZEROI:=0.0; CONER:=1.0; CONEI:=0.0;
      CTWOR:=2.0; CTWOI:=0.0; R1:=2.0;

      DPI:=3.14159265358979324; RTHPI:=1.25331413731550025;
      SPI:=1.90985931710274403; HPI:=1.57079632679489662;
      FPI:=1.89769999331517738; TTH:=6.66666666666666666E-01;

      CC[1]:= 5.77215664901532861E-01; CC[2]:=-4.20026350340952355E-02;
      CC[3]:=-4.21977345555443367E-02; CC[4]:= 7.21894324666309954E-03;
      CC[5]:=-2.15241674114950973E-04; CC[6]:=-2.01348547807882387E-05;
      CC[7]:= 1.13302723198169588E-06; CC[8]:= 6.11609510448141582E-09;

      CAZ := ZABS(ZR,ZI);
      CSCLR := 1.0/TOL;
      CRSCR := TOL;
      CSSR[1] := CSCLR;
      CSSR[2] := 1.0;
      CSSR[3] := CRSCR;
      CSRR[1] := CRSCR;
      CSRR[2] := 1.0;
      CSRR[3] := CSCLR;
      BRY[1] := 1000*D1MACH(1)/TOL;
      BRY[2] := 1.0/BRY[1];
      BRY[3] := D1MACH(2);
      NZ := 0;
      IFLAG := 0;
      KODED := KODE;
      RCAZ := 1.0/CAZ;
      STR := ZR*RCAZ;
      STI := -ZI*RCAZ;
      RZR := (STR+STR)*RCAZ;
      RZI := (STI+STI)*RCAZ;
      INU := Round(FNU+0.5);
      DNU := FNU - 1.0*INU;
      IF ABS(DNU) =  0.5 THEN GOTO 110;
      DNU2 := 0.0;
      IF ABS(DNU) > TOL THEN DNU2 := DNU*DNU;
      IF CAZ > R1 THEN GOTO 110;
{-----------------------------------------------------------------------
!     SERIES FOR CABS(Z).LE.R1
!----------------------------------------------------------------------}
      FC := 1.0;
      ZLOG(RZR, RZI, SMUR, SMUI, IDUM);
      FMUR := SMUR*DNU;
      FMUI := SMUI*DNU;
      ZSHCH(FMUR, FMUI, CSHR, CSHI, CCHR, CCHI);
      IF DNU = 0.0 THEN GOTO 10;
      FC := DNU*DPI;
      FC := FC/SIN(FC);
      SMUR := CSHR/DNU;
      SMUI := CSHI/DNU;
10:   A2 := 1.0 + DNU;
{-----------------------------------------------------------------------
!     GAM(1-Z)*GAM(1+Z):=PI*Z/SIN(PI*Z), T1:=1/GAM(1-DNU), T2:=1/GAM(1+DNU)
!----------------------------------------------------------------------}
      T2 := EXP(-DGAMLN(A2,IDUM));
      T1 := 1.0/(T2*FC);
      IF ABS(DNU) > 0.1 THEN GOTO 40;
{-----------------------------------------------------------------------
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
!----------------------------------------------------------------------}
      AK := 1.0;
      S := CC[1];
      For K:=2 to 8 do
      begin
        AK := AK*DNU2;
        TM := CC[K]*AK;
        S := S + TM;
        IF ABS(TM) < TOL THEN GOTO 30
      end;
30:   G1 := -S;
      GOTO 50;
40:   G1 := (T1-T2)/(DNU+DNU);
50:   G2 := (T1+T2)*0.5;
      FR := FC*(CCHR*G1+SMUR*G2);
      FI := FC*(CCHI*G1+SMUI*G2);
      ZEXP(FMUR, FMUI, STR, STI);
      PR := 0.5*STR/T2;
      PI := 0.5*STI/T2;
      ZDIV(0.5, 0.0, STR, STI, PTR, PTI);
      QR := PTR/T1;
      QI := PTI/T1;
      S1R := FR;
      S1I := FI;
      S2R := PR;
      S2I := PI;
      AK := 1.0;
      A1 := 1.0;
      CKR := CONER;
      CKI := CONEI;
      BK := 1.0 - DNU2;
      IF (INU > 0) OR (N >1) THEN GOTO 80;
{-----------------------------------------------------------------------
!     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N:=1
!----------------------------------------------------------------------}
      IF CAZ < TOL THEN GOTO 70;
      ZMLT(ZR, ZI, ZR, ZI, CZR, CZI);
      CZR := 0.25*CZR;
      CZI := 0.25*CZI;
      T1 := 0.25*CAZ*CAZ;
60:   FR := (FR*AK+PR+QR)/BK;
      FI := (FI*AK+PI+QI)/BK;
      STR := 1.0/(AK-DNU);
      PR := PR*STR;
      PI := PI*STR;
      STR := 1.0/(AK+DNU);
      QR := QR*STR;
      QI := QI*STR;
      STR := CKR*CZR - CKI*CZI;
      RAK := 1.0/AK;
      CKI := (CKR*CZI+CKI*CZR)*RAK;
      CKR := STR*RAK;
      S1R := CKR*FR - CKI*FI + S1R;
      S1I := CKR*FI + CKI*FR + S1I;
      A1 := A1*T1*RAK;
      BK := BK + AK + AK + 1.0;
      AK := AK + 1.0;
      IF A1 > TOL THEN GOTO 60;
70:   YR[1] := S1R;
      YI[1] := S1I;
      IF KODED = 1 THEN GOTO RETURN;
      ZEXP(ZR, ZI, STR, STI);
      ZMLT(S1R, S1I, STR, STI, YR[1], YI[1]);
      GOTO RETURN;
{-----------------------------------------------------------------------
!     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
!----------------------------------------------------------------------}
80:   IF CAZ < TOL THEN GOTO 100;
      ZMLT(ZR, ZI, ZR, ZI, CZR, CZI);
      CZR := 0.25*CZR;
      CZI := 0.25*CZI;
      T1 := 0.25*CAZ*CAZ;
90:   FR := (FR*AK+PR+QR)/BK;
      FI := (FI*AK+PI+QI)/BK;
      STR := 1.0/(AK-DNU);
      PR := PR*STR;
      PI := PI*STR;
      STR := 1.0/(AK+DNU);
      QR := QR*STR;
      QI := QI*STR;
      STR := CKR*CZR - CKI*CZI;
      RAK := 1.0/AK;
      CKI := (CKR*CZI+CKI*CZR)*RAK;
      CKR := STR*RAK;
      S1R := CKR*FR - CKI*FI + S1R;
      S1I := CKR*FI + CKI*FR + S1I;
      STR := PR - FR*AK;
      STI := PI - FI*AK;
      S2R := CKR*STR - CKI*STI + S2R;
      S2I := CKR*STI + CKI*STR + S2I;
      A1 := A1*T1*RAK;
      BK := BK + AK + AK + 1.0;
      AK := AK + 1.0;
      IF A1 > TOL THEN GOTO 90;
100:  KFLAG := 2;
      A1 := FNU + 1.0;
      AK := A1*ABS(SMUR);
      IF AK > ALIM THEN KFLAG := 3;
      STR := CSSR[KFLAG];
      P2R := S2R*STR;
      P2I := S2I*STR;
      ZMLT(P2R, P2I, RZR, RZI, S2R, S2I);
      S1R := S1R*STR;
      S1I := S1I*STR;
      IF KODED = 1 THEN GOTO 210;
      ZEXP(ZR, ZI, FR, FI);
      ZMLT(S1R, S1I, FR, FI, S1R, S1I);
      ZMLT(S2R, S2I, FR, FI, S2R, S2I);
      GOTO 210;
{-----------------------------------------------------------------------
!     IFLAG:=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG:=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED:=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
!----------------------------------------------------------------------}
110:  ZSQRT(ZR, ZI, STR, STI);
      ZDIV(RTHPI, CZEROI, STR, STI, COEFR, COEFI);
      KFLAG := 2;
      IF KODED = 2 THEN GOTO 120;
      IF ZR > ALIM THEN GOTO 290;

      STR := EXP(-ZR)*CSSR[KFLAG];
      STI := -STR*SIN(ZI);
      STR := STR*COS(ZI);
      ZMLT(COEFR, COEFI, STR, STI, COEFR, COEFI);
120:  IF ABS(DNU) = 0.5 THEN GOTO 300;
{-----------------------------------------------------------------------
!     MILLER ALGORITHM FOR CABS(Z).GT.R1
!----------------------------------------------------------------------}
      AK := COS(DPI*DNU);
      AK := ABS(AK);
      IF AK = CZEROR THEN GOTO 300;
      FHS := ABS(0.25-DNU2);
      IF FHS = CZEROR THEN GOTO 300;
{-----------------------------------------------------------------------
!     COMPUTE R2:=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE TO
!     DETERMINE THE BACKWARD INDEX K. R2:=F(E) IS A STRAIGHT LINE ON
!     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E):=B**(1-I1MACH(14)):=
!     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
!----------------------------------------------------------------------}
      T1 := 1.0*(I1MACH(14)-1);
      T1 := T1*D1MACH(5)*3.321928094;
      T1 := DMAX(T1,12.0);
      T1 := DMIN(T1,60.0);
      T2 := TTH*T1 - 6.0;
      IF  ZR <> 0.0 THEN GOTO 130;
      T1 := HPI;
      GOTO 140;
130:  T1 := ArcTan(ZI/ZR);
      T1 := ABS(T1);
140:  IF T2 > CAZ THEN GOTO 170;
{-----------------------------------------------------------------------
!     FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2
!----------------------------------------------------------------------}
      ETEST := AK/(DPI*CAZ*TOL);
      FK := CONER;
      IF ETEST < CONER THEN GOTO 180;
      FKS := CTWOR;
      CKR := CAZ + CAZ + CTWOR;
      P1R := CZEROR;
      P2R := CONER;
      For I:=1 to KMAX do
      begin
        AK := FHS/FKS;
        CBR := CKR/(FK+CONER);
        PTR := P2R;
        P2R := CBR*P2R - P1R*AK;
        P1R := PTR;
        CKR := CKR + CTWOR;
        FKS := FKS + FK + FK + CTWOR;
        FHS := FHS + FK + FK;
        FK := FK + CONER;
        STR := ABS(P2R)*FK;
        IF ETEST < STR THEN GOTO 160
      end;
      GOTO 310;
160:  FK := FK + SPI*T1*SQRT(T2/CAZ);
      FHS := ABS(0.25-DNU2);
      GOTO 180;
{-----------------------------------------------------------------------
!     COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2
!----------------------------------------------------------------------}
170:  A2 := SQRT(CAZ);
      AK := FPI*AK/(TOL*SQRT(A2));
      AA := 3.0*T1/(1.0+CAZ);
      BB := 14.7*T1/(28.0+CAZ);
      AK := (Ln(AK)+CAZ*COS(AA)/(1.0+0.008*CAZ))/COS(BB);
      FK := 0.12125*AK*AK/CAZ + 1.5;
{-----------------------------------------------------------------------
!     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
!----------------------------------------------------------------------}
180:  K := Round(FK);
      FK := 1.0*K;
      FKS := FK*FK;
      P1R := CZEROR;
      P1I := CZEROI;
      P2R := TOL;
      P2I := CZEROI;
      CSR := P2R;
      CSI := P2I;
      For I:=1 to K do
      begin
        A1 := FKS - FK;
        AK := (FKS+FK)/(A1+FHS);
        RAK := 2.0/(FK+CONER);
        CBR := (FK+ZR)*RAK;
        CBI := ZI*RAK;
        PTR := P2R;
        PTI := P2I;
        P2R := (PTR*CBR-PTI*CBI-P1R)*AK;
        P2I := (PTI*CBR+PTR*CBI-P1I)*AK;
        P1R := PTR;
        P1I := PTI;
        CSR := CSR + P2R;
        CSI := CSI + P2I;
        FKS := A1 - FK + CONER;
        FK := FK - CONER
      end;
{-----------------------------------------------------------------------
!     COMPUTE (P2/CS):=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR BETTER
!     SCALING
!----------------------------------------------------------------------}
      TM := ZABS(CSR,CSI);
      PTR := 1.0/TM;
      S1R := P2R*PTR;
      S1I := P2I*PTR;
      CSR := CSR*PTR;
      CSI := -CSI*PTR;
      ZMLT(COEFR, COEFI, S1R, S1I, STR, STI);
      ZMLT(STR, STI, CSR, CSI, S1R, S1I);
      IF (INU > 0) OR (N > 1) THEN GOTO 200;
      ZDR := ZR;
      ZDI := ZI;
      IF IFLAG = 1 THEN GOTO 270;
      GOTO 240;
{-----------------------------------------------------------------------
!     COMPUTE P1/P2:=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR SCALING
!----------------------------------------------------------------------}
200:  TM := ZABS(P2R,P2I);
      PTR := 1.0/TM;
      P1R := P1R*PTR;
      P1I := P1I*PTR;
      P2R := P2R*PTR;
      P2I := -P2I*PTR;
      ZMLT(P1R, P1I, P2R, P2I, PTR, PTI);
      STR := DNU + 0.5 - PTR;
      STI := -PTI;
      ZDIV(STR, STI, ZR, ZI, STR, STI);
      STR := STR + 1.0;
      ZMLT(STR, STI, S1R, S1I, S2R, S2I);
{-----------------------------------------------------------------------
!     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
!     SCALING NEAR EXPONENT EXTREMES ON KFLAG:=1 OR KFLAG:=3
!----------------------------------------------------------------------}
210:  STR := DNU + 1.0;
      CKR := STR*RZR;
      CKI := STR*RZI;
      IF N = 1 THEN INU := INU - 1;
      IF INU > 0 THEN GOTO 220;
      IF N > 1 THEN GOTO 215;
      S1R := S2R;
      S1I := S2I;
215:  ZDR := ZR;
      ZDI := ZI;
      IF IFLAG = 1 THEN GOTO 270;
      GOTO 240;
220:  INUB := 1;
      IF IFLAG = 1 THEN GOTO 261;
225:  P1R := CSRR[KFLAG];
      ASCLE := BRY[KFLAG];
      For I:=INUB to INU do
      begin
        STR := S2R;
        STI := S2I;
        S2R := CKR*STR - CKI*STI + S1R;
        S2I := CKR*STI + CKI*STR + S1I;
        S1R := STR;
        S1I := STI;
        CKR := CKR + RZR;
        CKI := CKI + RZI;
        IF KFLAG >= 3 THEN GOTO 230;
        P2R := S2R*P1R;
        P2I := S2I*P1R;
        STR := ABS(P2R);
        STI := ABS(P2I);
        P2M := DMAX(STR,STI);
        IF P2M <= ASCLE THEN GOTO 230;
        KFLAG := KFLAG + 1;
        ASCLE := BRY[KFLAG];
        S1R := S1R*P1R;
        S1I := S1I*P1R;
        S2R := P2R;
        S2I := P2I;
        STR := CSSR[KFLAG];
        S1R := S1R*STR;
        S1I := S1I*STR;
        S2R := S2R*STR;
        S2I := S2I*STR;
        P1R := CSRR[KFLAG];
230:  end;
      IF N <> 1 THEN GOTO 240;
      S1R := S2R;
      S1I := S2I;
240:  STR := CSRR[KFLAG];
      YR[1] := S1R*STR;
      YI[1] := S1I*STR;
      IF N = 1 THEN GOTO RETURN;
      YR[2] := S2R*STR;
      YI[2] := S2I*STR;
      IF N = 2 THEN GOTO RETURN;
      KK := 2;
250:  KK := KK + 1;
      IF KK > N THEN GOTO RETURN;
      P1R := CSRR[KFLAG];
      ASCLE := BRY[KFLAG];
      For I:=KK to N do
      begin
        P2R := S2R;
        P2I := S2I;
        S2R := CKR*P2R - CKI*P2I + S1R;
        S2I := CKI*P2R + CKR*P2I + S1I;
        S1R := P2R;
        S1I := P2I;
        CKR := CKR + RZR;
        CKI := CKI + RZI;
        P2R := S2R*P1R;
        P2I := S2I*P1R;
        YR[I] := P2R;
        YI[I] := P2I;
        IF  KFLAG >= 3 THEN GOTO 260;
        STR := ABS(P2R);
        STI := ABS(P2I);
        P2M := DMAX(STR,STI);
        IF P2M <= ASCLE THEN GOTO 260;
        KFLAG := KFLAG + 1;
        ASCLE := BRY[KFLAG];
        S1R := S1R*P1R;
        S1I := S1I*P1R;
        S2R := P2R;
        S2I := P2I;
        STR := CSSR[KFLAG];
        S1R := S1R*STR;
        S1I := S1I*STR;
        S2R := S2R*STR;
        S2I := S2I*STR;
        P1R := CSRR[KFLAG];
260:  end;
      GOTO RETURN;
{-----------------------------------------------------------------------
!     IFLAG:=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
!----------------------------------------------------------------------}
261:  HELIM := 0.5*ELIM;
      ELM := EXP(-ELIM);
      CELMR := ELM;
      ASCLE := BRY[1];
      ZDR := ZR;
      ZDI := ZI;
      IC := -1;
      J := 2;
      For I:=1 to INU do
      begin
        STR := S2R;
        STI := S2I;
        S2R := STR*CKR-STI*CKI+S1R;
        S2I := STI*CKR+STR*CKI+S1I;
        S1R := STR;
        S1I := STI;
        CKR := CKR+RZR;
        CKI := CKI+RZI;
        AS := ZABS(S2R,S2I);
        ALAS := Ln(AS);
        P2R := -ZDR+ALAS;
        IF P2R < -ELIM THEN GOTO 263;
        ZLOG(S2R,S2I,STR,STI,IDUM);
        P2R := -ZDR+STR;
        P2I := -ZDI+STI;
        P2M := EXP(P2R)/TOL;
        P1R := P2M*COS(P2I);
        P1I := P2M*SIN(P2I);
        ZUCHK(P1R,P1I,NW,ASCLE,TOL);
        IF NW <> 0 THEN GOTO 263;
        J := 3 - J;
        CYR[J] := P1R;
        CYI[J] := P1I;
        IF IC = I-1 THEN GOTO 264;
        IC := I;
        GOTO 262;
263:    IF ALAS < HELIM THEN GOTO 262;
        ZDR := ZDR-ELIM;
        S1R := S1R*CELMR;
        S1I := S1I*CELMR;
        S2R := S2R*CELMR;
        S2I := S2I*CELMR;
262:  end;
      IF N <> 1 THEN GOTO 270;
      S1R := S2R;
      S1I := S2I;
      GOTO 270;
264:  KFLAG := 1;
      INUB := I+1;
      S2R := CYR[J];
      S2I := CYI[J];
      J := 3 - J;
      S1R := CYR[J];
      S1I := CYI[J];
      IF INUB <= INU THEN GOTO 225;
      IF N <> 1 THEN GOTO 240;
      S1R := S2R;
      S1I := S2I;
      GOTO 240;
270:  YR[1] := S1R;
      YI[1] := S1I;
      IF N = 1 THEN GOTO 280;
      YR[2] := S2R;
      YI[2] := S2I;
280:  ASCLE := BRY[1];
      ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM);
      INU := N - NZ;
      IF INU <= 0 THEN GOTO RETURN;
      KK := NZ + 1;
      S1R := YR[KK];
      S1I := YI[KK];
      YR[KK] := S1R*CSRR[1];
      YI[KK] := S1I*CSRR[1];
      IF INU = 1 THEN GOTO RETURN;
      KK := NZ + 2;
      S2R := YR[KK];
      S2I := YI[KK];
      YR[KK] := S2R*CSRR[1];
      YI[KK] := S2I*CSRR[1];
      IF INU = 2 THEN GOTO RETURN;
      T2 := FNU + 1.0*(KK-1);
      CKR := T2*RZR;
      CKI := T2*RZI;
      KFLAG := 1;
      GOTO 250;
{-----------------------------------------------------------------------
!     SCALE BY DEXP(Z), IFLAG := 1 CASES
!----------------------------------------------------------------------}
290:  KODED := 2;
      IFLAG := 1;
      KFLAG := 2;
      GOTO 120;
{-----------------------------------------------------------------------
!     FNU:=HALF ODD INTEGER CASE, DNU:=-0.5
!----------------------------------------------------------------------}
300:  S1R := COEFR;
      S1I := COEFI;
      S2R := COEFR;
      S2I := COEFI;
      GOTO 210;
310:  NZ:=-2;
Return: End; {ZBKNU}


Procedure ZKSCL(ZRR,ZRI,FNU:Double; N:Integer; Var YR,YI:VEC; Var NZ:Integer;
                Var RZR,RZI,ASCLE,TOL,ELIM: Double);
{***BEGIN PROLOGUE  ZKSCL
!***REFER TO  ZBESK
!
!     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
!     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
!     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
!
!***ROUTINES CALLED  ZUCHK,ZABS,ZLOG
!***END PROLOGUE  ZKSCL
!     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM }
Label 10,20,25,30,40,45,Return;
Var
      ACS, AS, CKI, CKR, CSI, CSR, FN, STR, S1I, S1R, S2I, S2R,
      ZEROI, ZEROR, ZDR, ZDI, CELMR, ELM, HELIM, ALAS: Double;
      I, IC, IDUM, KK, NN, NW: Integer;
      CYR, CYI: VEC;
Begin

      ZEROR:=0.0; ZEROI:=0.0;
      NZ := 0;
      IC := 0;
      NN := IMIN(2,N);
      For I:=1 to NN do
      begin
        S1R := YR[I];
        S1I := YI[I];
        CYR[I] := S1R;
        CYI[I] := S1I;
        AS := ZABS(S1R,S1I);
        ACS := -ZRR + Ln(AS);
        NZ := NZ + 1;
        YR[I] := ZEROR;
        YI[I] := ZEROI;
        IF ACS < -ELIM THEN GOTO 10;
        ZLOG(S1R, S1I, CSR, CSI, IDUM);
        CSR := CSR - ZRR;
        CSI := CSI - ZRI;
        STR := EXP(CSR)/TOL;
        CSR := STR*COS(CSI);
        CSI := STR*SIN(CSI);
        ZUCHK(CSR, CSI, NW, ASCLE, TOL);
        IF NW <> 0 THEN GOTO 10;
        YR[I] := CSR;
        YI[I] := CSI;
        IC := I;
        NZ := NZ - 1;
10:   end;
      IF N = 1 THEN GOTO RETURN;
      IF IC > 1 THEN  GOTO 20;
      YR[1] := ZEROR;
      YI[1] := ZEROI;
      NZ := 2;
20:   IF N = 2 THEN GOTO RETURN;
      IF NZ = 0 THEN GOTO RETURN;
      FN := FNU + 1.0;
      CKR := FN*RZR;
      CKI := FN*RZI;
      S1R := CYR[1];
      S1I := CYI[1];
      S2R := CYR[2];
      S2I := CYI[2];
      HELIM := 0.5*ELIM;
      ELM := EXP(-ELIM);
      CELMR := ELM;
      ZDR := ZRR;
      ZDI := ZRI;

{     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
      S2 GETS LARGER THAN EXP(ELIM/2)  }

      For I:=3 to N do
      begin
        KK := I;
        CSR := S2R;
        CSI := S2I;
        S2R := CKR*CSR - CKI*CSI + S1R;
        S2I := CKI*CSR + CKR*CSI + S1I;
        S1R := CSR;
        S1I := CSI;
        CKR := CKR + RZR;
        CKI := CKI + RZI;
        AS := ZABS(S2R,S2I);
        ALAS := Ln(AS);
        ACS := -ZDR + ALAS;
        NZ := NZ + 1;
        YR[I] := ZEROR;
        YI[I] := ZEROI;
        IF ACS < -ELIM THEN GOTO 25;
        ZLOG(S2R, S2I, CSR, CSI, IDUM);
        CSR := CSR - ZDR;
        CSI := CSI - ZDI;
        STR := EXP(CSR)/TOL;
        CSR := STR*COS(CSI);
        CSI := STR*SIN(CSI);
        ZUCHK(CSR, CSI, NW, ASCLE, TOL);
        IF NW <> 0 THEN GOTO 25;
        YR[I] := CSR;
        YI[I] := CSI;
        NZ := NZ - 1;
        IF IC = KK-1 THEN GOTO 40;
        IC := KK;
        GOTO 30;
25:     IF ALAS < HELIM THEN GOTO 30;
        ZDR := ZDR - ELIM;
        S1R := S1R*CELMR;
        S1I := S1I*CELMR;
        S2R := S2R*CELMR;
        S2I := S2I*CELMR;
30:   end;
      NZ := N;
      IF IC = N THEN NZ:=N-1;
      GOTO 45;
40:   NZ := KK - 2;
45:   For I:=1 to NZ do
      begin
        YR[I] := ZEROR;
        YI[I] := ZEROI
      end;
Return:End; {ZKSCL}

END.

{end of file Cbess00.pas}