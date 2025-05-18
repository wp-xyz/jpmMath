{****************************************************************
* EVALUATE A K-BESSEL FUNCTION OF COMPLEX ARGUMENT (THIRD KIND) *
* ------------------------------------------------------------- *
* SAMPLE RUN:                                                   *
* (Evaluate K0 to K4 for argument Z=(1.0,2.0) ).                *
*                                                               *
* zr(0) =  -0.242345                                            *
* zi(0) =  -0.176267                                            *
* zr(1) =  -0.300362                                            *
* zi(1) =  -0.151186                                            *
* zr(2) =  -0.483439                                            *
* zi(2) =   0.003548                                            *
* zr(3) =  -0.681436                                            *
* zi(3) =   0.625155                                            *
* zr(4) =   0.199208                                            *
* zi(4) =   2.389181                                            *
* NZ = 0                                                        *
* Error code: 0                                                 *
*                                                               *
* ------------------------------------------------------------- *
* Ref.: From Numath Library By Tuan Dang Trong in Fortran 77    *
*       [BIBLI 18].                                             *
*                                                               *
*                        TPW Release 1.0 By J-P Moreau, Paris   *
*                                 (www.jpmoreau.fr)             *
*****************************************************************
Note: Used files: CBess0,CBess00,CBess1,CBess3,Complex.
------------------------------------------------------- }
PROGRAM TEST_ZBESK;

Uses WinCrt, Complex, CBess0, CBess00, CBess3;

Var
    zr,zi: Double;
    cyr, cyi: VEC;
    i,ierr,n,nz: Integer;


Procedure ZBESK(ZR, ZI, FNU:double; KODE, N:integer;
                Var CYR, CYI: VEC; Var NZ, IERR:integer);
{***BEGIN PROLOGUE  ZBESK
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  830501   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION OF THE SECOND KIND,
!             BESSEL FUNCTION OF THE THIRD KIND
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!***DESCRIPTION
!
!                      ***A DOUBLE PRECISION ROUTINE***
!
!         ON KODE=1, CBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
!         BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE
!         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z.NE.CMPLX(0.0,0.0)
!         IN THE CUT PLANE -PI < ARG(Z) <= PI. ON KODE=2, CBESK
!         RETURNS THE SCALED K FUNCTIONS,
!
!         CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,
!
!         WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
!         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
!         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
!         FUNCTIONS (REF. 1).
!
!         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
!           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
!                    -PI.LT.ARG(Z).LE.PI
!           FNU    - ORDER OF INITIAL K FUNCTION, FNU.GE.0.0D0
!           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
!           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!                    KODE= 1  RETURNS
!                             CY(I)=K(FNU+I-1,Z), I=1,...,N
!                        = 2  RETURNS
!                             CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
!
!         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
!           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
!                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
!                    CY(I)=K(FNU+I-1,Z), I=1,...,N OR
!                    CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
!                    DEPENDING ON KODE
!           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
!                    NZ= 0   , NORMAL RETURN
!                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
!                              TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0),
!                              I=1,...,N WHEN X >= 0. WHEN X < 0,
!                              NZ STATES ONLY THE NUMBER OF UNDERFLOWS
!                              IN THE SEQUENCE.
!
!           IERR   - ERROR FLAG
!                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!                    IERR=1, INPUT ERROR   - NO COMPUTATION
!                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS
!                            TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
!                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
!                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
!                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
!                            ACCURACY
!                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
!                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
!                            CANCE BY ARGUMENT REDUCTION
!                    IERR=5, ERROR              - NO COMPUTATION,
!                            ALGORITHM TERMINATION CONDITION NOT MET
!
!***LONG DESCRIPTION
!
!         EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS
!         DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X.GE.0.0. FORWARD
!         RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
!         HALF PLANE BY THE RELATION
!
!         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
!         MP=MR*PI*I, MR=+1 OR -1, RE(Z) > 0, I^2=-1
!
!         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
!
!         FOR LARGE ORDERS, FNU > FNUL, THE K FUNCTION IS COMPUTED
!         BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.
!
!         FOR NEGATIVE ORDERS, THE FORMULA
!
!                       K(-FNU,Z) = K(FNU,Z)
!
!         CAN BE USED.
!
!         CBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS
!         AVAILABLE.
!
!         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
!         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
!         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
!         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
!         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
!         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
!         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
!         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
!         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
!         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
!         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
!         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
!         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
!         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
!         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
!         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
!         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
!         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
!         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
!
!         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
!         BESSEL FUNCTION CAN BE EXPRESSED BY P*10^S WHERE P=MAX(UNIT
!         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10^S REPRE-
!         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
!         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
!         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
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
!                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
!
!               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983.
!
!               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
!                 1018, MAY, 1985
!
!               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
!                 MATH. SOFTWARE, 1986
!
!***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,ZABS,I1MACH,D1MACH
!***END PROLOGUE  ZBESK
!
!     COMPLEX CY,Z }
Label 50,60,70,80,90,100, 180,200,260,Return;
Var
      AA, ALIM, ALN, ARG, AZ, DIG, ELIM, FN, FNUL, RL, R1M5, TOL, UFL, BB: Double;
      K, K1, K2, MR, NN, NUF, NW: Integer;
Begin
{***FIRST EXECUTABLE STATEMENT  ZBESK }
      IERR := 0;
      NZ:=0;
      IF (ZI = 0.0) AND (ZR = 0.0) THEN IERR:=1;
      IF FNU < 0.0 THEN IERR:=1;
      IF (KODE < 1) OR (KODE > 2) THEN IERR:=1;
      IF N < 1 THEN IERR:=1;
      IF IERR <> 0 THEN GOTO RETURN;   {bad parameter(s) }
      NN := N;
{-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM).LT.EXP(-ALIM):=EXP(-ELIM)/TOL    AND
!     EXP(ELIM).GT.EXP(ALIM):=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG := NUMBER OF BASE 10 DIGITS IN TOL := 10**(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
!----------------------------------------------------------------------}
      TOL := DMAX(D1MACH(4),1E-18);
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
      FNUL := 10.0 + 6.0*(DIG-3.0);
      RL := 1.2*DIG + 3.0;
{-----------------------------------------------------------------------
!     TEST FOR PROPER RANGE
!----------------------------------------------------------------------}
      AZ := ZABS(ZR,ZI);
      FN := FNU + 1.0*(NN-1);
      AA := 0.5/TOL;
      BB:=0.5*I1MACH(9);
      AA := DMIN(AA,BB);
      IF AZ > AA THEN GOTO 260;
      IF FN > AA THEN GOTO 260;
      AA := SQRT(AA);
      IF AZ > AA THEN IERR:=3;
      IF FN > AA THEN IERR:=3;
{-----------------------------------------------------------------------
!     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
!----------------------------------------------------------------------}
      UFL := EXP(-ELIM);
      IF AZ < UFL THEN GOTO 180;
      IF FNU > FNUL THEN GOTO 80;
      IF FN <= 1.0 THEN GOTO 60;
      IF FN > 2.0 THEN GOTO 50;
      IF AZ > TOL THEN GOTO 60;
      ARG := 0.5*AZ;
      ALN := -FN*Ln(ARG);
      IF ALN > ELIM THEN GOTO 180;
      GOTO 60;
50:   ZUOIK(ZR, ZI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM, ALIM);
      IF NUF < 0 THEN GOTO 180;
      NZ := NZ + NUF;
      NN := NN - NUF;
{-----------------------------------------------------------------------
!     HERE NN:=N OR NN:=0 SINCE NUF:=0,NN, OR -1 ON RETURN FROM CUOIK
!     IF NUF:=NN, THEN CY(I):=CZERO FOR ALL I
!----------------------------------------------------------------------}
      IF NN = 0 THEN GOTO 100;
60:   IF ZR < 0.0 THEN GOTO 70;
{-----------------------------------------------------------------------
!     RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0.
!----------------------------------------------------------------------}
      ZBKNU(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM);
      IF NW < 0 THEN GOTO 200;
      NZ:=NW;
      GOTO RETURN;
{-----------------------------------------------------------------------
!     LEFT HALF PLANE COMPUTATION
!     PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2.
!----------------------------------------------------------------------}
70:   IF NZ <> 0 THEN GOTO 180;
      MR := 1;
      IF ZI < 0.0 THEN MR := -1;
      ZACON(ZR, ZI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL, TOL, ELIM, ALIM);
      IF NW < 0 THEN GOTO 200;
      NZ:=NW;
      GOTO RETURN;
{-----------------------------------------------------------------------
!     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
!----------------------------------------------------------------------}
80:   MR := 0;
      IF ZR >= 0.0 THEN GOTO 90;
      MR := 1;
      IF ZI < 0.0 THEN MR := -1;
90:   ZBUNK(ZR, ZI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM, ALIM);
      IF NW < 0 THEN GOTO 200;
      NZ := NZ + NW;
      GOTO RETURN;
100:  IF ZR < 0.0 THEN GOTO 180;
      GOTO RETURN;
180:  NZ := 0;
      IERR:=2;
      GOTO RETURN;
200:  IF NW = -1 THEN GOTO 180;
      NZ:=0;
      IERR:=5;
      GOTO RETURN;
260:  NZ:=0;
      IERR:=4;
Return: End; {ZBESK}


{main program}
BEGIN

  n:=5;
  zr:=1.0; zi:=2.0;

  ZBESK(zr,zi,0.0,1,n,cyr,cyi,nz,ierr);

  writeln;
  for i:=1 to n do
  begin
    writeln(' zr(',i-1,') = ',cyr[i]:10:6);
    writeln(' zi(',i-1,') = ',cyi[i]:10:6)
  end;
  writeln(' NZ = ', NZ);
  writeln(' Error code: ', ierr);
  writeln;
  ReadKey;
  DoneWinCrt
  
END.

{end of file tzbesk.pas}