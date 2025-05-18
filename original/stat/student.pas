{************************************************************
*                STUDENT  T-PROBABILITY  LAW                *
* --------------------------------------------------------- *
* SAMPLE RUN:                                               *
* (Calculate Student T-probability (lower tail and upper    *
*  tail) for T=0.257).                                      *
*                                                           *
*  X= 0.2570000                                             *
*  PROB1= 0.6000294                                         *
*  ERROR=0                                                  *
*  X= 0.2570000                                             *
*  PROB2= 0.3999705                                         *
*  ERROR=0                                                  *
*  PROB1+PROB2= 0.9999998                                   *
*                                                           *
* --------------------------------------------------------- *
* Ref.:"JOURNAL OF APPLIED STATISTICS (1968) VOL.17, P.189, *
*       & VOL.19, NO.1".                                    *
*                                                           *
*                     Pascal Release By J-P Moreau, Paris.  *
*                               (www.jpmoreau.fr)           *
************************************************************}
PROGRAM STUDENT;
Uses WinCrt;

Const
      ZERO = 0.0;
      HALF = 0.5;
      ONE = 1.0;
      TWO = 2.0;
      FOUR = 4.0;

Var
      X,XNU,PROB1,PROB2: REAL;
      NU, ERROR: INTEGER;


    FUNCTION Power(x:REAL;n:Integer):REAL;
    {Calculate x power n}
    var i,m : integer;
        result :real;
    begin
      result := ONE;
      if n=0 then
      begin
        power:=result;
        exit;
      end;
      m:=  n;
      if n<0 then m:=-n;
      for i:=1 to m do result :=x*result;
      Power :=result;
      if n<0 then Power:=ONE/result;
    end;


    FUNCTION PROBST(T:REAL; IDF:Integer; Var IFAULT:Integer): REAL;
{ ---------------------------------------------------------------------
!        ALGORITHM AS 3  APPL. STATIST. (1968) VOL.17, P.189
!
!        STUDENT T PROBABILITY (LOWER TAIL)
! --------------------------------------------------------------------}
    Label 20, 30, Return;
    Const G1:REAL = ONE/PI;
    Var
       A, B, C, F, S, FK, ZSQRT, ZATAN: REAL;
       IM2,IOE,K,KS:Integer;
    Begin
      IFAULT := 1;
      PROBST := ZERO;
      IF IDF < 1 then goto RETURN;
      IFAULT := 0;
      F := ONE*IDF;
      A := T / SQRT(F);
      B := F / (F + T * T);
      IM2 := IDF - 2;
      IOE := IDF MOD 2;
      S := ONE;
      C := ONE;
      F := ONE;
      KS := 2 + IOE;
      FK := KS;
      IF IM2 < 2 THEN GOTO 20;
      K:=KS;
      While K <= IM2 do
      begin
        C := C * B * (FK - ONE) / FK;
        S := S + C;
        IF S = F THEN GOTO 20;
        F := S;
        FK := FK + TWO;
        Inc(K,2)
      end;
20:   IF IOE = 1 THEN GOTO 30;
      PROBST := HALF + HALF * A * SQRT(B) * S;
      GOTO RETURN;
30:   IF IDF = 1 THEN S := ZERO;
      PROBST := HALF + (A * B * S + ARCTAN(A)) * G1;
Return: End;


    FUNCTION STUDNT (T, DOFF:REAL; Var IFAULT:Integer): REAL;
{ ----------------------------------------------------------------
!     ALGORITHM AS 27  APPL. STATIST. VOL.19, NO.1
!
!     Calculate the upper tail area under Student's t-distribution
!
!     Translated from Algol by Alan Miller
! ---------------------------------------------------------------}
    Label Return;
    Var
      V, X, TT: REAL;
      A1, A2, A3, A4, A5, B1, B2, C1, C2, C3, C4, C5, D1, D2: REAL;
      E1, E2, E3, E4, E5, F1, F2, G1, G2, G3, G4, G5, H1, H2: REAL;
      I1, I2, I3, I4, I5, J1, J2: REAL;
      POS: Boolean;

    Begin

      A1:=0.09979441; A2:=-0.581821; A3:=1.390993; A4:=-1.222452; A5:=2.151185;
      B1:=5.537409; B2:=11.42343;
      C1:=0.04431742; C2:=-0.2206018; C3:=-0.03317253; C4:=5.679969; C5:=-12.96519;
      D1:=5.166733; D2:=13.49862;
      E1:=0.009694901; E2:=-0.1408854; E3:=1.88993; E4:=-12.75532; E5:=25.77532;
      F1:=4.233736; F2:=14.3963;
      G1:=-9.187228E-5; G2:=0.03789901; G3:=-1.280346; G4:=9.249528; G5:=-19.08115;
      H1:=2.777816; H2:=16.46132;
      I1:=5.79602E-4; I2:=-0.02763334; I3:=0.4517029; I4:=-2.657697; I5:=5.127212;
      J1:=0.5657187; J2:=21.83269;

{     Check that number of degrees of freedom > 4 }

      IF DOFF < TWO THEN
      begin
	IFAULT := 1;
	STUDNT := - ONE;
	goto RETURN
      end;

      IF DOFF <= FOUR THEN
	IFAULT := Round(DOFF)
      ELSE
	IFAULT := 0;

{     Evaluate series }

      V := ONE / DOFF;
      POS := (T >= ZERO);
      TT := ABS(T);
      X := (ONE + TT * (((A1 + V * (A2 + V * (A3 + V * (A4 + V * A5)))) / (ONE - V * (B1 - V * B2))) +
           TT * (((C1 + V * (C2 + V * (C3 + V * (C4 + V * C5)))) / (ONE - V * (D1 - V * D2))) +
           TT * (((E1 + V * (E2 + V * (E3 + V * (E4 + V * E5)))) / (ONE - V * (F1 - V * F2))) +
           TT * (((G1 + V * (G2 + V * (G3 + V * (G4 + V * G5)))) / (ONE - V * (H1 - V * H2))) +
           TT * ((I1 + V * (I2 + V * (I3 + V * (I4 + V * I5)))) /
           (ONE - V * (J1 - V * J2))))))));
      X:=HALF*Power(X,-8);
      IF (POS) THEN
	STUDNT := X
      ELSE
	STUDNT := ONE - X;

RETURN: End;


{main program}
BEGIN

  X:=0.257;
  NU:=19;

  PROB1:=PROBST(X,NU,ERROR);
  writeln;
  writeln('  X=', X:10:7);
  writeln('  PROB1=', PROB1:10:7);
  writeln('  ERROR=', ERROR);

  X:=0.257;
  XNU:=19.0;

  PROB2:=STUDNT(X,XNU,ERROR);
  writeln('  X=', X:10:7);
  writeln('  PROB2=', PROB2:10:7);
  writeln('  ERROR=', ERROR);

  writeln('  PROB1+PROB2=', PROB1+PROB2:10:7);
  writeln;

  ReadKey;
  DoneWinCrt

END.

{end of file Student.pas}