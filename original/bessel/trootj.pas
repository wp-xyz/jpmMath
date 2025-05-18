{ --------------------------------------------------------------------------
! Program to calculate the first zeroes (root abscissas) of the first 
! kind Bessel function of integer order N using the subroutine ROOTJ.
! --------------------------------------------------------------------------
! SAMPLE RUN:
!
! (Calculate the first 10 zeroes of 1st kind Bessel function of order 2).
! 
! Zeroes of Bessel Function of order:  2
! 
! Number of calculated zeroes: 10 
!
! Table of root abcissas (4 items per line)
!  5.1356223E+0000  8.4172442E+0000  1.1619841E+0001  1.4795952E+0001
!  1.7959819E+0001  2.1116997E+0001  2.4270112E+0001  2.7420574E+0001
!  3.0569204E+0001  3.3716520E+0001
! 
! Table of error codes (4 items per line)
!   0   0   0   0
!   0   0   0   0
!   0   0
!
! --------------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                               TPW Release 1.0 By J-P Moreau, Paris
!                                        (www.jpmoreau.fr)
! -------------------------------------------------------------------------}
PROGRAM TROOTJ;

Uses Wincrt;

Type
     Tab10  = Array[1..10] of Double;
     ITab10 = Array[1..10] of Integer;

Var
     JZ: Tab10;
     IE: ITab10;
     I,J,N,NR: Integer;

    Procedure SECANT(N,NITMX:Integer;TOL:Double;Var ZEROJ:Double;Var IER:Integer); Forward;
    Function  BESSJ(N:Integer;X:Double): Double; Forward;
    Function  BESSJ0(X:Double): Double; Forward;
    Function  BESSJ1(X:Double): Double; Forward;


    Function Power(y,x:double): double;
    Begin
      IF x<0 THEN EXIT;
      Power:=Exp(x*Ln(y))
    End;

    Function Sign(a,b:Double): Double;
    Begin
      if b>=0.0 then Sign:=ABS(a)
                else Sign:=-ABS(a)
    End;

{ ---------------------------------------------------------------------}
    PROCEDURE ROOTJ(N,NK:Integer; Var JZERO:Tab10; Var IER:ITab10);
{ ----------------------------------------------------------------------
!     CALCULATE THE FIRST NK ZEROES OF BESSEL FUNCTION J(N,X)
!
!     INPUTS:
!       N    ORDER OF FUNCTION J (INTEGER >= 0)                  I*4
!       NK   NUMBER OF FIRST ZEROES  (INTEGER > 0)               I*4
!     OUTPUTS:
!       JZERO(NK)  TABLE OF FIRST ZEROES (ABCISSAS)              R*8
!       IER(NK)    TABLE OF ERROR CODES (MUST BE ZEROES)         I*4
!
!     REFERENCE :
!     ABRAMOWITZ M. & STEGUN IRENE A.
!     HANDBOOK OF MATHEMATICAL FUNCTIONS
! ---------------------------------------------------------------------}
    Var
      ZEROJ,B0,B1,B2,B3,B5,B7,T0,T1,T3,T5,T7,FN,FK: Double;
      C1,C2,C3,C4,C5,F1,F2,F3,TOL,ERRJ: Double;
      IERROR,K,NITMX: Integer;
    Begin  
      TOL:=1E-8; NITMX:=10;
      C1:=1.8557571; C2:=1.033150; C3:=0.00397; C4:=0.0908; C5:=0.043;
      FN := 1.0*N;

{     FIRST ZERO }
      IF N=0 THEN
      begin

        ZEROJ := C1+C2-C3-C4+C5;
{       WRITE(*,'(1X,A,I5,E15.6)') 'K=1,N=0,ZEROJ',K,ZEROJ }

        SECANT(N,NITMX,TOL,ZEROJ,IERROR);
        
	IER[1]:=IERROR;
        JZERO[1]:=ZEROJ
      end
      ELSE
      begin
        F1 := Power(FN,(1.0/3.0));
        F2 := F1*F1*FN;
        F3 := F1*FN*FN;
        ZEROJ := FN+C1*F1+(C2/F1)-(C3/FN)-(C4/F2)+(C5/F3);
        SECANT(N,NITMX,TOL,ZEROJ,IERROR);
        IER[1]:=IERROR;
        JZERO[1]:=ZEROJ

      END;

      T0 := 4.0*FN*FN;
      T1 := T0-1.0;
      T3 := 4.0*T1*(7.0*T0-31.0);
      T5 := 32.0*T1*((83.0*T0-982.0)*T0+3779.0);
      T7 := 64.0*T1*(((6949.0*T0-153855.0)*T0+1585743.0)*T0-6277237.0);

{    OTHER ZEROES }
      For K := 2 to NK do
      begin
        JZERO[K] := 0.0;
        FK := 1.0*K;

{    MAC MAHON'S SERIES FOR K>>N }

        B0 := (FK+0.5*FN-0.25)*PI;
        B1 := 8.0*B0;
        B2 := B1*B1;
        B3 := 3.0*B1*B2;
        B5 := 5.0*B3*B2;
        B7 := 7.0*B5*B2;

        ZEROJ := B0-(T1/B1)-(T3/B3)-(T5/B5)-(T7/B7);

        ERRJ:=ABS(BESSJ(N,ZEROJ));

{    IMPROVE SOLUTION USING PROCEDURE SECANT }

        IF ERRJ > TOL THEN SECANT(N,NITMX,TOL,ZEROJ,IERROR);
        JZERO[K]:=ZEROJ;
        IER[K]:=IERROR;
      end
    End;
{ ------------------------------------------------------------------------------ }
    PROCEDURE SECANT(N,NITMX:Integer;TOL:Double;Var ZEROJ:Double;Var IER:Integer);
    Label 5,10,15,20;
    Var
        P0,P1,Q0,Q1,DP,P: Double;
        C: array[1..2] of Double;
        IT,NEV,NTRY: Integer;
    Begin
      C[1]:=0.95; C[2]:=0.999;
      NTRY:=1;
      IER:=0;

5:    P0 := C[NTRY]*ZEROJ;

      P1 := ZEROJ;
      NEV := 2;
      Q0 := BESSJ(N,P0);
      Q1 := BESSJ(N,P1);
{     WRITE(*,'(1X,A,I5,4E15.6)') 'NTRY,P0,Q0,P1,Q1',NTRY,P0,Q0,P1,Q1 }
      For IT := 1 to NITMX do
      begin
        IF Q1=Q0 THEN GOTO 15;
        P := P1-Q1*(P1-P0)/(Q1-Q0);
        DP := P-P1;
{       WRITE(*,'(1X,A,I5,4E15.6)') 'IT,P,DP',IT,P,DP }
        IF IT=1 THEN GOTO 10;
        IF ABS(DP) < TOL THEN GOTO 20;

10:     NEV := NEV+1;
        P0 := P1;
        Q0 := Q1;
        P1 := P;
        Q1 := BESSJ(N,P1)
      end;
15:   NTRY:=NTRY+1;
      IF NTRY <= 2 THEN GOTO 5;
      IER := NTRY;
20:   ZEROJ := P
    End;

{ ---------------------------------------------------------------------- }
    FUNCTION BESSJ(N:Integer;X:Double): Double;

{     THIS FUNCTION RETURNS THE VALUE OF THE FIRST KIND BESSEL FUNCTION 
      OF ORDER N, INTEGER FOR ANY REAL X. WE USE HERE THE CLASSICAL 
      RECURRENT FORMULA, WHEN  X > N. FOR X < N, THE MILLER'S ALGORITHM
      IS USED TO AVOID OVERFLOWS.
      REFERENCE :
      C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
      MATHEMATICAL TABLES, VOL.5, 1962.  }
    Label  EndProc;
    Const
      IACC = 40; BIGNO = 1E10; BIGNI = 1E-10;
    Var
      TOX,BJM,BJ,BJP,SUM,TMP: Double;
      J,JSUM,M: Integer;
    Begin
      IF N=0 THEN
      begin
        BESSJ := BESSJ0(X);
        goto EndProc
      end;
      IF N=1 THEN
      begin
        BESSJ := BESSJ1(X);
        goto EndProc
      end;
      IF X=0.0 THEN
      begin
        BESSJ := 0.0;
        goto EndProc 
      end;
      TOX := 2.0/X;
      IF X > N THEN
      begin
        BJM := BESSJ0(X);
        BJ  := BESSJ1(X);
        For J := 1 to N-1 do
        begin
          BJP := J*TOX*BJ-BJM;
          BJM := BJ;
          BJ  := BJP;
          BESSJ := BJ
        end
      end
      ELSE
      begin
        M := 2*((N+Round(SQRT(IACC*N))) div 2);
        TMP := 0.0;
        JSUM := 0;
        SUM := 0.0;
        BJP := 0.0;
        BJ  := 1.0;
        For J := M Downto 1 do
        begin
          BJM := J*TOX*BJ-BJP;
          BJP := BJ;
          BJ  := BJM;
          IF ABS(BJ) > BIGNO THEN
          begin
            BJ  := BJ*BIGNI;
            BJP := BJP*BIGNI;
            TMP := TMP * BIGNI;
            SUM := SUM * BIGNI
          end;
          IF JSUM <> 0 THEN SUM := SUM+BJ;
          JSUM := 1-JSUM;
          IF J = N THEN TMP := BJP;
          SUM := 2.0*SUM-BJ
        end;
        BESSJ := TMP/SUM
      end;
EndProc: End;
{ ---------------------------------------------------------------------- }
    FUNCTION BESSJ0(X:Double): Double;
      
{     THIS FUNCTION RETURNS THE VALUE OF THE FIRST KIND BESSEL FUNCTION 
      OF ORDER 0 FOR ANY REAL X. WE USE HERE THE POLYNOMIAL APPROXIMATION
      BY SERIES OF CHEBYSHEV POLYNOMIALS FOR 0<X<8 AND 0<8/X<1.
      REFERENCES :
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.  }
    Label 1, EndFunc;
    Var
        AX,FR,FS,Z,FP,FQ,XX: Double;
        Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6: Double;
    Begin

      P1:=1.0; P2:=-0.1098628627E-2; P3:=0.2734510407E-4;
      P4:=-0.2073370639E-5; P5:=0.2093887211E-6;
      Q1:=-0.1562499995E-1; Q2:=0.1430488765E-3; Q3:=-0.6911147651E-5;
      Q4:=0.7621095161E-6; Q5:=-0.9349451520E-7;
      R1:=57568490574.0; R2:=-13362590354.0; R3:=651619640.7;
      R4:=-11214424.18; R5:=77392.33017; R6:=-184.9052456;
      S1:=57568490411.0; S2:=1029532985.0; S3:=9494680.718;
      S4:=59272.64853; S5:=267.8532712; S6:=1.0;

      IF X=0.0 THEN GOTO 1;
      AX := ABS(X);
      IF AX < 8 THEN
      begin
        Y := X*X;
        FR := R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
        FS := S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
        BESSJ0 := FR/FS
      end
      ELSE
      begin
        Z := 8./AX;
        Y := Z*Z;
        XX := AX-0.785398164;
        FP := P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
        FQ := Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
        BESSJ0 := SQRT(0.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      end;
      goto EndFunc;
1:    BESSJ0 := 1.0;

EndFunc: End;
{ ---------------------------------------------------------------------- }
    FUNCTION BESSJ1(X:Double): Double;
      
{     THIS FUNCTION RETURNS THE VALUE OF THE FIRST KIND BESSEL FUNCTION 
      OF ORDER 0 FOR ANY REAL X. WE USE HERE THE POLYNOMIAL APPROXIMATION
      BY SERIES OF CHEBYSHEV POLYNOMIALS FOR 0<X<8 AND 0<8/X<1.
      REFERENCES :
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.  }
    Var
        AX,FR,FS,Z,FP,FQ,XX: Double;
        Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6: Double;
    Begin

      P1:=1.0; P2:=0.183105E-2; P3:=-0.3516396496E-4;
      P4:=0.2457520174E-5; P5:=-0.240337019E-6; P6:=0.636619772;
      Q1:=0.04687499995; Q2:=-0.2002690873E-3; Q3:=0.8449199096E-5;
      Q4:=-0.88228987E-6; Q5:=0.105787412E-6;
      R1:=72362614232.0; R2:=-7895059235.0; R3:=242396853.1;
      R4:=-2972611.439; R5:=15704.48260; R6:=-30.16036606;
      S1:=144725228442.0; S2:=2300535178.0; S3:=18583304.74;
      S4:=99447.43394; S5:=376.9991397; S6:=1.0;

      AX := ABS(X);
      IF AX < 8 THEN
      begin
        Y := X*X;
        FR := R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
        FS := S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
        BESSJ1 := X*(FR/FS)
      end
      ELSE
      begin
        Z := 8./AX;
        Y := Z*Z;
        XX := AX-2.35619491;
        FP := P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
        FQ := Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
        BESSJ1 := SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*Sign(S6,X)
      end
    End;

{main program}
BEGIN

  N:=2;
  NR:=10;

  ROOTJ(N,NR,JZ,IE);

  writeln;
  writeln(' Zeroes of Bessel Function of order: ', N);
  writeln;
  writeln(' Number of calculated zeroes: ', NR);
  writeln;
  writeln(' Table of root abcissas (4 items per line)');
  j:=1;
  for i:=1 to NR do
  begin
    write(' ',JZ[i]:16); Inc(j);
    if j MOD 5 = 0 then
    begin
      writeln; j:=1
    end
  end;
  writeln; writeln;
  writeln(' Table of error codes (4 items per line)');
  j:=1;
  for i:=1 to NR do
  begin
    write('   ',IE[i]); Inc(j);
    if j MOD 5 = 0 then
    begin
      writeln; j:=1
    end
  end;
  writeln;
  ReadKey; DoneWinCrt

END.

{ end of file trootj.pas }