{****************************************************************
* EVALUATE A I-BESSEL FUNCTION OF COMPLEX ARGUMENT (MODIFIED    *
* FIRST KIND)                                                   *
* ------------------------------------------------------------- *
* SAMPLE RUN:                                                   *
* (Evaluate I0 to I4 for argument Z=(1.0,2.0) ).                *
*                                                               *
* zr(0) =   0.187854                                            *
* zi(0) =   0.646169                                            *
* zr(1) =  -0.079933                                            *
* zi(1) =   0.790623                                            *
* zr(2) =  -0.412672                                            *
* zi(2) =   0.265974                                            *
* zr(3) =  -0.175353                                            *
* zi(3) =  -0.082431                                            *
* zr(4) =  -0.004414                                            *
* zi(4) =  -0.055957                                            *
* NZ = 0                                                        *
* Error code: 0                                                 *
*                                                               *
* ------------------------------------------------------------- *
* Ref.: From Numath Library By Tuan Dang Trong in Fortran 77.   *
*                                                               *
*                        TPW Release 1.0 By J-P Moreau, Paris   *
*                                  (www.jpmoreau.fr)            *   
*****************************************************************
Note: Used files: CBess0,CBess00,CBess1,CBess2,Complex.
------------------------------------------------------- }
PROGRAM TEST_ZBESI;

Uses WinCrt, Complex, CBess1 {for ZBINU};

Var
  zr,zi: double;
  cyr, cyi: VEC;
  i,n,nz,ierr: integer;


Procedure ZBESI(ZR, ZI, FNU:double; KODE, N: integer; Var CYR, CYI: VEC;
                var NZ, IERR:integer);
{***BEGIN PROLOGUE  ZBESI
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  830501   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  I-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION OF THE FIRST KIND
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE I-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!***DESCRIPTION
!
!                    ***A DOUBLE PRECISION ROUTINE***
!         ON KODE=1, ZBESI COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
!         BESSEL FUNCTIONS CY(J)=I(FNU+J-1,Z) FOR REAL, NONNEGATIVE
!         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z IN THE CUT PLANE
!         -PI.LT.ARG(Z).LE.PI. ON KODE=2, ZBESI RETURNS THE SCALED
!         FUNCTIONS
!
!         CY(J)=EXP(-ABS(X))*I(FNU+J-1,Z)   J = 1,...,N , X=REAL(Z)
!
!         WITH THE EXPONENTIAL GROWTH REMOVED IN BOTH THE LEFT AND
!         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
!         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
!         (REF. 1).
!
!         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
!           ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI.LT.ARG(Z).LE.PI
!           FNU    - ORDER OF INITIAL I FUNCTION, FNU.GE.0.0D0
!           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!                    KODE= 1  RETURNS
!                             CY(J)=I(FNU+J-1,Z), J=1,...,N
!                        = 2  RETURNS
!                             CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X)), J=1,...,N
!           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
!
!         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
!           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
!                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
!                    CY(J)=I(FNU+J-1,Z)  OR
!                    CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X))  J=1,...,N
!                    DEPENDING ON KODE, X=REAL(Z)
!           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
!                    NZ= 0   , NORMAL RETURN
!                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET TO ZERO
!                              TO UNDERFLOW, CY(J)=CMPLX(0.0D0,0.0D0)
!                              J = N-NZ+1,...,N
!           IERR   - ERROR FLAG
!                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!                    IERR=1, INPUT ERROR   - NO COMPUTATION
!                    IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z) TOO
!                            LARGE ON KODE=1
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
!         THE COMPUTATION IS CARRIED OUT BY THE POWER SERIES FOR
!         SMALL CABS(Z), THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z),
!         THE MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN AND A
!         NEUMANN SERIES FOR IMTERMEDIATE MAGNITUDES, AND THE
!         UNIFORM ASYMPTOTIC EXPANSIONS FOR I(FNU,Z) AND J(FNU,Z)
!         FOR LARGE ORDERS. BACKWARD RECURRENCE IS USED TO GENERATE
!         SEQUENCES OR REDUCE ORDERS WHEN NECESSARY.
!
!         THE CALCULATIONS ABOVE ARE DONE IN THE RIGHT HALF PLANE AND
!         CONTINUED INTO THE LEFT HALF PLANE BY THE FORMULA
!
!         I(FNU,Z*EXP(M*PI)) = EXP(M*PI*FNU)*I(FNU,Z)  REAL(Z).GT.0.0
!                       M = +I OR -I,  I**2=-1
!
!         FOR NEGATIVE ORDERS,THE FORMULA
!
!              I(-FNU,Z) = I(FNU,Z) + (2/PI)*SIN(PI*FNU)*K(FNU,Z)
!
!         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
!         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
!         INTEGER,THE MAGNITUDE OF I(-FNU,Z)=I(FNU,Z) IS A LARGE
!         NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
!         K(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
!         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
!         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
!         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
!         LARGE MEANS FNU.GT.CABS(Z).
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
!         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
!         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
!         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
!         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
!         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
!         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
!         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
!         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
!         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
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
!***ROUTINES CALLED  ZBINU,I1MACH,D1MACH
!***END PROLOGUE  ZBESI
!     COMPLEX CONE,CSGN,CW,CY,CZERO,Z,ZN }
Label 40, 120, 130, 260, Return;
Var
      AA, ALIM, ARG, CONEI, CONER, CSGNI, CSGNR, DIG, ELIM,
      FNUL, PI, RL, R1M5, STR, TOL, ZNI, ZNR, AZ, BB, FN: Double;
      I, INU, K, K1, K2, NN: Integer;
Begin
{***FIRST EXECUTABLE STATEMENT  ZBESI }
      CONER:=1.0; CONEI:=0.0; 
      IERR := 0;
      NZ:=0;
      IF FNU < 0.0 Then IERR:=1;
      IF (KODE < 1) OR (KODE > 2) Then IERR:=1;
      IF N < 1 Then IERR:=1;
      IF IERR <> 0 Then goto Return;
{-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10^(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
!----------------------------------------------------------------------}
      TOL := DMAX(D1MACH(4),1E-18);
      K1 := I1MACH(15);
      K2 := I1MACH(16);
      R1M5 := D1MACH(5);
      K := IMIN(ABS(K1), ABS(K2));
      ELIM := 2.303*(K*R1M5-3.0);
      K1 := I1MACH(14) - 1;
      AA := R1M5*K1;
      DIG := DMIN(AA,18.0);
      AA := AA*2.303;
      ALIM := ELIM + DMAX(-AA,-41.45);
      RL := 1.2*DIG + 3.0;
      FNUL := 10.0 + 6.0*(DIG-3.0);
{-----------------------------------------------------------------------
!     TEST FOR PROPER RANGE
!----------------------------------------------------------------------}
      AZ := ZABS(ZR,ZI);
      FN := FNU+1.0*(N-1);
      AA := 0.5/TOL;
      BB:=I1MACH(9)*0.5;
      AA := DMIN(AA,BB);
      IF AZ > AA Then GOTO 260;
      IF FN > AA Then GOTO 260;
      AA := SQRT(AA);
      IF AZ > AA Then IERR:=3;
      IF FN > AA Then IERR:=3;
      ZNR := ZR;
      ZNI := ZI;
      CSGNR := CONER;
      CSGNI := CONEI;
      IF ZR >= 0.0 Then GOTO 40;
      ZNR := -ZR;
      ZNI := -ZI;
{-----------------------------------------------------------------------
!     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!----------------------------------------------------------------------}
      INU := Round(FNU);
      ARG := (FNU-1.0*INU)*PI;
      IF ZI < 0.0 Then ARG := -ARG;
      CSGNR := COS(ARG);
      CSGNI := SIN(ARG);
      IF (INU MOD 2) = 0 Then GOTO 40;
      CSGNR := -CSGNR;
      CSGNI := -CSGNI;
40:   ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL, ELIM, ALIM);
      IF NZ < 0 Then GOTO 120;
      IF ZR >= 0.0 Then goto RETURN;
{-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
!----------------------------------------------------------------------}
      NN := N - NZ;
      IF NN = 0 Then goto RETURN;
      For I:=1 to NN do
      begin
        STR := CYR[I]*CSGNR - CYI[I]*CSGNI;
        CYI[I] := CYR[I]*CSGNI + CYI[I]*CSGNR;
        CYR[I] := STR;
        CSGNR := -CSGNR;
        CSGNI := -CSGNI
      end;
      goto RETURN;
120:  IF NZ = -2 Then GOTO 130;
      NZ := 0;
      IERR:=2;
      goto RETURN;
130:  NZ:=0;
      IERR:=5;
      goto RETURN;
260:  NZ:=0;
      IERR:=4;
Return: End; {ZBESI}


{main program}
Begin
  n:=5;
  zr:=1.0; zi:=2.0;

  ZBESI(zr,zi,0,1,n,cyr,cyi,nz,ierr);

  writeln;
  for i:=1 to n do
  begin
    writeln(' zr(', i-1, ') = ', cyr[i]:10:6);
    writeln(' zi(', i-1, ') = ', cyi[i]:10:6)
  end;
  writeln(' NZ = ', NZ);
  writeln(' Error code: ', ierr);
  writeln;
  ReadKey;
  DoneWinCrt
  
END.

{end of file tzbesi.pas}