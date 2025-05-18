{*************************************************************
*     CALCULATE THE Kth ZERO OF THE DERIVATIVE OF BESSEL     *
*     FUNCTION OF ORDER N, J(N,X)                            *
* ---------------------------------------------------------- *
* SAMPLE RUN:                                                *
* (Calculate the 10th zero of the derivative of J(1,X) ).    *
*                                                            *
*  X=  3.06019229735146E+0001                                *
*                                                            *
*  BESSJP(X)= -1.75243069973916E-0016                        *
*                                                            *
* ---------------------------------------------------------- *
* Reference:   From Numath Library By Tuan Dang Trong in     *
*              Fortran 77 [BIBLI 18].                        *
*                                                            *
*                          TPW Release By J-P Moreau, Paris  *
*                                 (www.jpmoreau.fr)          *
*************************************************************}
PROGRAM TEST_ZEROJP;

Uses WinCrt;

VAR
     RES: Double; K,N: Integer;


    FUNCTION Power(y,x:Double):Double;
    Begin
      IF x<0 THEN EXIT;
      Power:=Exp(x*Ln(y))
    End;

    Function Sign(a,b:Double): Double;
    Begin
      if b<0 then Sign:=-ABS(a)
             else Sign:=ABS(a)
    End;


    FUNCTION BESSJ (N:Integer;X:Double):Double; Forward;
    FUNCTION BESSJ0(X:Double):Double; Forward;
    FUNCTION BESSJ1(X:Double):Double; Forward;
    FUNCTION BESSJP(N:Integer;X:Double):Double; Forward;


    FUNCTION ZEROJP(N,K:Integer): Double;
{--------------------------------------------------------------------
!     CALCULATE THE Kth ZERO OF THE DERIVATIVE OF BESSEL FUNCTION
!     OF ORDER N, J(N,X) 
!--------------------------------------------------------------------
!     CALLING MODE:
!       RES = ZEROJP(N,K)
!
!     INPUTS:
!       N    ORDER OF BESSEL FUNCTION J (INTEGER >= 0)            I*4
!       K    RANK OF ZERO (INTEGER > 0)                           I*4
!     OUTPUT:
!       ZEROJP                                                    R*8
!     REFERENCE:
!     ABRAMOWITZ M. & STEGUN IRENE A.
!     HANDBOOK OF MATHEMATICAL FUNCTIONS
!--------------------------------------------------------------------}
    Label 10,20,25, 40,50;
    Var
      B0,B1,B2,B3,B5,B7,T0,T1,T3,T5,T7,FN,FK: Double;
      C1,C2,C3,C4,F1,F2,P,DP,P0,P1,Q0,Q1,TOL: Double;
      ZJP: Double;
      IER,IT,NEV,NITMX: Integer;
      IMPROV: Boolean;
    Begin

      TOL:=1.e-7; NITMX:=15;
      C1:=0.8086165; C2:=0.072490; C3:=0.05097; C4:=0.0094;
      IMPROV:=TRUE;
      FN := 1.0*N;
      FK := 1.0*K;

      IF K > 1 THEN GOTO 10;

{     if N = 0 ET K = 1 }

      IF N = 0 THEN
      begin
        ZEROJP:= 0.0;
        goto 50;
      end

{     TCHEBYCHEV'S SERIES FOR K <= N }

      ELSE
      begin
        F1 := Power(FN,1.0/3.0);
        F2 := F1*F1*FN;
        ZJP := FN+C1*F1+(C2/F1)-(C3/FN)+(C4/F2);
        GOTO 20
      end;

{     MAC MAHON'S SERIES FOR K >> N }

10:   B0 := (FK+0.5*FN-0.75)*PI;
      B1 := 8.0*B0;
      B2 := B1*B1;
      B3 := 3.0*B1*B2;
      B5 := 5.0*B3*B2;
      B7 := 7.0*B5*B2;
      T0 := 4.0*FN*FN;
      T1 := T0 + 3.0;
      T3 := 4.0*((7.0*T0+82.0)*T0-9.0);
      T5 := 32.0*(((83.0*T0+2075.0)*T0-3039.0)*T0+3537.0);
      T7 := 64.0*((((6949.0*T0+296492.0)*T0-1248002.0)*T0
                      +7414380.0)*T0-5853627.0);

      ZJP := B0-(T1/B1)-(T3/B3)-(T5/B5)-(T7/B7);

20:   IF IMPROV THEN
      begin
{     IMPROVE SOLUTION BY SECANT METHOD WHEN K > N
      AND IMPROV = TRUE  }
        P0 := 0.9*ZJP;
        P1 := ZJP;
        IER := 0;
        NEV := 2;
        Q0 := BESSJP(N,P0);
        Q1 := BESSJP(N,P1);
        For IT := 1 to NITMX do
        begin
          P := P1-Q1*(P1-P0)/(Q1-Q0);
          DP := P-P1;
          IF IT = 1 THEN GOTO 25;
          IF ABS(DP) < TOL THEN GOTO 40;
25:       NEV := NEV+1;
          P0 := P1;
          Q0 := Q1;
          P1 := P ;
          Q1 := BESSJP(N,P1);
        end;
        IER := 1;
        WRITELN(' ** ZEROJP ** NITMX EXCEEDED.');
        ZEROJP:=ZJP;
        goto 50;
40:     ZEROJP := P
      end;
50: End;

    FUNCTION BESSJP(N:Integer; X:Double): Double;
{ ----------------------------------------------------------------------
!    NAME  :  BESSJP
!    DATE  :  06/01/1982
!    IV    :  1
!    IE    :  1
!    AUTHOR:  DANG TRONG TUAN
! ......................................................................
!
!    FIRST DERIVATIVE OF FIRST KIND BESSEL FUNCTION OF ORDER N, FOR REAL X
!
!                          MODULE BESSJP
! ......................................................................
!
!    THIS SUBROUTINE CALCULATES THE FIRST DERIVATIVE OF FIRST KIND BESSEL 
!    FUNCTION OF ORDER N, FOR REAL X.
!
! ......................................................................
!
!  I  VARIABLE DIMENSION/TYPE  DESCRIPTION  (INPUTS)
!        N       I*4           ORDER OF FUNCTION
!        X       R*8           ABSCISSA OF FUNCTION BESSJP(N,X)
!
!  O  VARIABLE,DIMENSION/TYPE  DESCRIPTION  (OUTPUT)
!
!      BESSJP    R*8           FUNCTION EVALUATION AT X
!.......................................................................
!    CALLED SUBROUTINE                                                  
!
!      BESSJ     FIRST KIND BESSEL FUNCTION                            
!
! ----------------------------------------------------------------------}
    Begin
      IF N = 0 THEN
        BESSJP:=-BESSJ(1,X)
      ELSE IF X = 0.0 THEN
        X:=1.e-30
      ELSE
        BESSJP:=BESSJ(N-1,X)-(1.0*N/X)*BESSJ(N,X);
    End;


    FUNCTION BESSJ(N:Integer; X:Double): Double;
{-----------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE 1ST KIND BESSEL FUNCTION OF ORDER
!     N, INTEGER FOR ANY REAL X. THE CLASSICAL RECURRENCE FORMULA IS USED
!     HERE, WHEN X IS > N, FOR X < N, THE MILLER'S ALGORITHM ALLOWS
!     AVOIDING OVERFLOWS.
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!----------------------------------------------------------------------}
    Label 10;
    Const IACC = 40; BIGNO = 1e10;  BIGNI = 1e-10;
    Var TOX,BJM,BJ,BJP,BSJ,SUM: Double;
        J,JSUM,M:Integer;
    Begin
      IF N = 0 THEN
      begin
        BESSJ := BESSJ0(X);
        goto 10
      end;
      IF N = 1 THEN
      begin
        BESSJ := BESSJ1(X);
        goto 10
      end;
      IF X = 0.0 THEN
      begin
        BESSJ := 0.0;
        goto 10
      end;
      TOX := 2.0/X;
      IF X > 1.0*N THEN
      begin
        BJM := BESSJ0(X);
        BJ  := BESSJ1(X);
        For J := 1 to N-1 do
        begin
          BJP := J*TOX*BJ-BJM;
          BJM := BJ;
          BJ  := BJP
        end;
        BESSJ := BJ
      end
      ELSE
      begin
        M := Round(2*((N+Int(SQRT(IACC*N)))/2));
        BSJ := 0.0;
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
            BSJ := BSJ*BIGNI;
            SUM := SUM*BIGNI
          end;
          IF JSUM <> 0 then SUM := SUM+BJ;
          JSUM := 1-JSUM;
          IF J = N then  BSJ := BJP
        end;
        SUM := 2.0*SUM-BJ;
        BESSJ := BSJ/SUM
      end;
10: End;

    FUNCTION BESSJ0(X:Double):Double;
{----------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE 1ST KIND BESSEL FUNCTION OF ORDER
!     ZERO FOR ANY REAL X. THE CHEBYSHEV'S POLYNOMIAL SERIES IS USED 
!     FOR 0<X<8 AND 0<8/X<1.
!     REFERENCE:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
!----------------------------------------------------------------------}
    Label 1,2;
    Var AX,FR,FS,Z,FP,FQ,XX:Double;
        Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6: Double;
        Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6: Double;
    Begin

      P1:=1.0; P2:=-0.1098628627e-02; P3:=0.2734510407e-04;
      P4:=-0.2073370639e-05; P5:=0.2093887211e-06;
      Q1:=-0.1562499995e-01; Q2:=0.1430488765e-03; Q3:=-0.6911147651e-05;
      Q4:=0.7621095161e-06; Q5:=-0.9349451520e-07;
      R1:=57568490574.0; R2:=-13362590354.0; R3:=651619640.7;
      R4:=-11214424.18; R5:=77392.33017; R6:=-184.9052456;
      S1:=57568490411.0; S2:=1029532985.0; S3:=9494680.718;
      S4:=59272.64853; S5:=267.8532712; S6:=1.0;

      IF X = 0.0 then GOTO 1;
      AX := ABS(X);
      IF AX < 8.0 THEN
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
      goto 2;
1:    BESSJ0 := 1.0;
2:  End;

    FUNCTION BESSJ1(X:Double): Double;
{----------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE 1ST KIND BESSEL FUNCTION OF ORDER
!     ONE FOR ANY REAL X. THE CHEBYSHEV'S POLYNOMIAL SERIES IS USED 
!     FOR 0<X<8 AND 0<8/X<1.
!     REFERENCE:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
!----------------------------------------------------------------------}
    Var AX,FR,FS,Z,FP,FQ,XX: Double;
        Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6: Double;
        Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6: Double;
    Begin

      P1:=1.0; P2:=0.183105e-02; P3:=-0.3516396496e-04; 
      P4:=0.2457520174e-05; P5:=-0.240337019e-06; P6:=0.636619772;
      Q1:=0.04687499995; Q2:=-0.2002690873e-03; Q3:=0.8449199096e-05;
      Q4:=-0.88228987e-06; Q5:=0.105787412e-06;
      R1:=72362614232.0; R2:=-7895059235.0; R3:=242396853.1;
      R4:=-2972611.439; R5:=15704.48260; R6:=-30.16036606;
      S1:=144725228442.0; S2:=2300535178.0; S3:=18583304.74;
      S4:=99447.43394; S5:=376.9991397; S6:=1.0;

      AX := ABS(X);
      IF AX < 8.0 THEN
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

  N := 1;
  K := 10;
  
  RES := ZEROJP(N,K);


  writeln;
  writeln('  X= ', RES);
  writeln;
  writeln('  BESSJP(X)= ', BESSJP(N,RES));
  writeln;

  ReadKey;
  DoneWinCrt

END.

{end of file tzerojp.pas}