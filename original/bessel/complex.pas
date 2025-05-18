{ ---------------------------------------------------------------------
!      Utility subroutines used by any program from Numath library
!      with (not intrinsic) complex numbers z = (zr,zi).
! ---------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                               TPW Release 1.0 By J-P Moreau, Paris
!                                        (www.jpmoreau.fr)
! --------------------------------------------------------------------}
UNIT COMPLEX;

INTERFACE

  Const NMAX = 10;

  Type  VEC = Array[1..NMAX] of double;
  Type  VEC16 = Array[1..16] of double;

  FUNCTION  ZABS(ZR, ZI:double): double;
  Procedure ZSQRT(AR, AI:double; var BR, BI:double);
  Procedure ZEXP(AR, AI:double; var BR, BI:double);
  Procedure ZMLT(AR, AI, BR, BI:double; var CR, CI:double);
  Procedure ZDIV(AR, AI, BR, BI:double; var CR, CI:double);
  Procedure ZLOG(AR, AI:double; var BR, BI:double; var IERR:integer);

  FUNCTION  D1MACH(I:Integer): double;
  FUNCTION  I1MACH(I:Integer):LongInt;

  Function DMAX(a,b:Double):Double;
  Function DMIN(a,b:Double):Double;
  Function IMAX(a,b:integer):integer;
  Function IMIN(a,b:integer):integer;


IMPLEMENTATION

FUNCTION ZABS(ZR, ZI:double): double;
{***BEGIN PROLOGUE  ZABS
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY

!     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
!     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)

!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  ZABS }
Label 10, 20, 30;
Var   U, V, Q, S: double;
Begin
      U := ABS(ZR);
      V := ABS(ZI);
      S := U + V;
{---------------------------------------------------------------------
!     S*1.0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
!     TRUE FLOATING ZERO
!--------------------------------------------------------------------}
      S := S*1.0;
      IF S = 0.0 THEN GOTO 20;
      IF U > V THEN GOTO 10;
      Q := U/V;
      ZABS := V*SQRT(1.0+Q*Q);
      GOTO 30;
10:   Q := V/U;
      ZABS := U*SQRT(1.0+Q*Q);
      GOTO 30;
20:   ZABS := 0.0;
30:End; {ZABS}

Procedure ZSQRT(AR, AI:double; var BR, BI:double);
{***BEGIN PROLOGUE  ZSQRT
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
!
!     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
!
!***ROUTINES CALLED  ZABS
!***END PROLOGUE  ZSQRT}
Label 10,20,30,40,50,60,70,return;
Var ZM, DTHETA, DPI, DRT: double;
Begin
      DRT:=7.071067811865475244008443621E-0001;
      DPI:=3.141592653589793238462643383E+0000;
      ZM := ZABS(AR,AI);
      ZM := SQRT(ZM);
      IF AR = 0.0 THEN GOTO 10;
      IF AI = 0.0 THEN GOTO 20;
      DTHETA := ARCTAN(AI/AR);
      IF DTHETA <= 0.0 THEN GOTO 40;
      IF AR < 0.0 THEN DTHETA := DTHETA - DPI;
      GOTO 50;
10:   IF AI > 0.0 THEN GOTO 60;
      IF AI < 0.0 THEN GOTO 70;
      BR := 0.0;
      BI := 0.0;
      GOTO RETURN;
20:   IF AR > 0.0 THEN GOTO 30;
      BR := 0.0;
      BI := SQRT(ABS(AR));
      GOTO RETURN;
30:   BR := SQRT(AR);
      BI := 0.0;
      GOTO RETURN;
40:   IF AR < 0.0 THEN DTHETA := DTHETA + DPI;
50:   DTHETA := DTHETA*0.5;
      BR := ZM*COS(DTHETA);
      BI := ZM*SIN(DTHETA);
      GOTO RETURN;
60:   BR := ZM*DRT;
      BI := ZM*DRT;
      GOTO RETURN;
70:   BR := ZM*DRT;
      BI := -ZM*DRT;
return: End; {ZSQRT}

Procedure ZEXP(AR, AI:double; var BR, BI:double);
{***BEGIN PROLOGUE  ZEXP
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
!
!     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B:=EXP(A)
!
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  ZEXP }
Var ZM, CA, CB: double;
Begin
      ZM := EXP(AR);
      CA := ZM*COS(AI);
      CB := ZM*SIN(AI);
      BR := CA;
      BI := CB
End; {ZEXP}

Procedure ZMLT(AR, AI, BR, BI:double; var CR, CI:double);
{***BEGIN PROLOGUE  ZMLT
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
!
!     DOUBLE PRECISION COMPLEX MULTIPLY, C:=A*B.
!
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  ZMLT }
Var CA, CB: double;
Begin
      CA := AR*BR - AI*BI;
      CB := AR*BI + AI*BR;
      CR := CA;
      CI := CB
End; {ZMLT}

Procedure ZDIV(AR, AI, BR, BI:double; var CR, CI:double);
{***BEGIN PROLOGUE  ZDIV
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY
!
!     DOUBLE PRECISION COMPLEX DIVIDE C:=A/B.

!***ROUTINES CALLED  ZABS
!***END PROLOGUE  ZDIV }
Var BM, CA, CB, CC, CD: double;
Begin
      BM := 1.0/ZABS(BR,BI);
      CC := BR*BM;
      CD := BI*BM;
      CA := (AR*CC+AI*CD)*BM;
      CB := (AI*CC-AR*CD)*BM;
      CR := CA;
      CI := CB
End; {ZDIV}

Procedure ZLOG(AR, AI:double; var BR, BI:double; var IERR:integer);
{***BEGIN PROLOGUE  ZLOG
!***REFER TO  ZBESH,ZBESI,ZBESJ,ZBESK,ZBESY,ZAIRY,ZBIRY

!     DOUBLE PRECISION COMPLEX LOGARITHM B:=CLOG(A)
!     IERR:=0,NORMAL RETURN      IERR:=1, Z:=CMPLX(0.0,0.0)
!***ROUTINES CALLED  ZABS
!***END PROLOGUE  ZLOG }
Label 10,20,30,40,50,60,return;
Var ZM, DTHETA, DPI, DHPI: double;
Begin
      DPI := 3.141592653589793238462643383E+0000;
      DHPI:= 1.570796326794896619231321696E+0000;

      IERR:=0;
      IF AR = 0.0 THEN GOTO 10;
      IF AI = 0.0 THEN GOTO 20;
      DTHETA := ArcTan(AI/AR);
      IF DTHETA <= 0.0 THEN GOTO 40;
      IF AR < 0.0 THEN DTHETA := DTHETA - DPI;
      GOTO 50;
10:   IF AI = 0.0 THEN GOTO 60;
      BI := DHPI;
      BR := Ln(ABS(AI));
      IF AI < 0.0 THEN BI := -BI;
      GOTO RETURN;
20:   IF AR > 0.0 THEN GOTO 30;
      BR := Ln(ABS(AR));
      BI := DPI;
      GOTO RETURN;
30:   BR := Ln(AR);
      BI := 0.0;
      GOTO RETURN;
40:   IF AR < 0.0 THEN DTHETA := DTHETA + DPI;
50:   ZM := ZABS(AR,AI);
      BR := Ln(ZM);
      BI := DTHETA;
      GOTO RETURN;
60:   IERR:=1;
return:End; {ZLOG}

FUNCTION D1MACH(I:Integer): double;
{***BEGIN PROLOGUE  D1MACH
!***DATE WRITTEN   750101   (YYMMDD)
!***REVISION DATE  860501   (YYMMDD)
!***CATEGORY NO.  R1
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  FOX, P. A., (BELL LABS)
!           HALL, A. D., (BELL LABS)
!           SCHRYER, N. L., (BELL LABS)
!***PURPOSE  RETURN DOUBLE PRECISION MACHINE DEPENDENT CONSTANTS.
!***DESCRIPTION

!     D1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
!     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
!     SUBPROGRAM WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
!     AS FOLLOWS, FOR EXAMPLE

!          D = D1MACH(I)

!     WHERE I=1,...,5.  THE (OUTPUT) VALUE OF D ABOVE IS
!     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
!     VARIOUS VALUES OF I ARE DISCUSSED BELOW.

!  DOUBLE-PRECISION MACHINE CONSTANTS
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  D1MACH( 5) = LOG10(B)
!***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
!                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
!***ROUTINES CALLED  XERROR
!***END PROLOGUE  D1MACH }

Var   DMACH: Array[1..5] of double;

Begin
{***FIRST EXECUTABLE STATEMENT  D1MACH }
      IF (I < 1) OR (I > 5) THEN
        Writeln(' D1MACH -- I OUT OF BOUNDS');

      {For IBM PC or APOLLO}
      DMACH[1] := 2.22559E-308;
      DMACH[2] := 1.79728E308;
      DMACH[3] := 1.11048E-16;
      DMACH[4] := 2.22096E-16;
      DMACH[5] := 0.301029995663981198;

      D1MACH := DMACH[I]

End; {D1MACH}

FUNCTION I1MACH(I:Integer):LongInt;
{***BEGIN PROLOGUE  I1MACH
!***DATE WRITTEN   750101   (YYMMDD)
!***REVISION DATE  890313   (YYMMDD)
!***CATEGORY NO.  R1
!***KEYWORDS  LIBRARY=SLATEC,TYPE=INTEGER(I1MACH-I),MACHINE CONSTANTS
!***AUTHOR  FOX, P. A., (BELL LABS)
!           HALL, A. D., (BELL LABS)
!           SCHRYER, N. L., (BELL LABS)
!***PURPOSE  Return integer machine dependent constants.
!***DESCRIPTION

!     I1MACH can be used to obtain machine-dependent parameters
!     for the local machine environment.  It is a function
!     subroutine with one (input) argument, and can be called
!     as follows, for example

!          K = I1MACH(I)

!     where I=1,...,16.  The (output) value of K above is
!     determined by the (input) value of I.  The results for
!     various values of I are discussed below.

!  I/O unit numbers.
!    I1MACH( 1) = the standard input unit.
!    I1MACH( 2) = the standard output unit.
!    I1MACH( 3) = the standard punch unit.
!    I1MACH( 4) = the standard error message unit.

!  Words.
!    I1MACH( 5) = the number of bits per integer storage unit.
!    I1MACH( 6) = the number of characters per integer storage unit.

!  Integers.
!    assume integers are represented in the S-digit, base-A form

!               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )

!               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!    I1MACH( 7) = A, the base.
!    I1MACH( 8) = S, the number of base-A digits.
!    I1MACH( 9) = A**S - 1, the largest magnitude.

!  Floating-Point Numbers.
!    Assume floating-point numbers are represented in the T-digit,
!    base-B form
!               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )

!               where 0 .LE. X(I) .LT. B for I=1,...,T,
!               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!    I1MACH(10) = B, the base.

!  Single-Precision
!    I1MACH(11) = T, the number of base-B digits.
!    I1MACH(12) = EMIN, the smallest exponent E.
!    I1MACH(13) = EMAX, the largest exponent E.

!  Double-Precision
!    I1MACH(14) = T, the number of base-B digits.
!    I1MACH(15) = EMIN, the smallest exponent E.
!    I1MACH(16) = EMAX, the largest exponent E.

!  To alter this function for a particular environment,
!  the desired set of DATA statements should be activated by
!  removing the ! from column 1.  Also, the values of
!  I1MACH(1) - I1MACH(4) should be checked for consistency
!  with the local operating system.

!***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
!                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
!                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  I1MACH }
Var IMACH: Array[1..16] of LongInt;
Begin
{***FIRST EXECUTABLE STATEMENT  D1MACH }
      IF (I < 1) OR (I > 16) THEN
        Writeln(' I1MACH -- I OUT OF BOUNDS');

      {For IBM PC or APOLLO}
      IMACH[ 1] :=    5;
      IMACH[ 2] :=    6;
      IMACH[ 3] :=    0;
      IMACH[ 4] :=    0;
      IMACH[ 5] :=   32;
      IMACH[ 6] :=    4;
      IMACH[ 7] :=    2;
      IMACH[ 8] :=   31;
      IMACH[ 9] := 2147483647;
      IMACH[10] :=    2;
      IMACH[11] :=   24;
      IMACH[12] := -125;
      IMACH[13] :=  127;
      IMACH[14] :=   53;
      IMACH[15] :=-1021;
      IMACH[16] := 1023;

      I1MACH := IMACH[i]

End; {I1MACH}

Function DMAX(a,b:Double):Double;
Begin
  if a>=b then DMAX:=a
          else DMAX:=b
End;

Function DMIN(a,b:Double):Double;
Begin
  if a<=b then DMIN:=a
          else DMIN:=b
End;

Function IMAX(a,b:integer):integer;
Begin
  if a>=b then IMAX:=a
          else IMAX:=b
End;

Function IMIN(a,b:integer):integer;
Begin
  if a<=b then IMIN:=a
          else IMIN:=b
End;

END.

{ end of file Complex.pas }