{ ---------------------------------------------------------------------
! Program to calculate the first kind modified Bessel function
! of integer order N, for any REAL X, using the function BESSI(N,X).
! ---------------------------------------------------------------------
! SAMPLE RUN:
!
! (Calculate Bessel function for N=2, X=0.75).
! 
! Bessel function of order 2 for X =  0.7500:
!
!      Y =  7.36668780479450E-0002
!
! ---------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                               TPW Release 1.0 By J-P Moreau, Paris
!                                        (www.jpmoreau.fr) 
! --------------------------------------------------------------------}
PROGRAM TBESSI;

Uses WinCrt;

Var
  X, Y: Double;
  N: Integer;

  Function BESSI0(X:Double): Double; Forward;
  Function BESSI1(X:Double): Double; Forward;

{ --------------------------------------------------------------------- }
      FUNCTION BESSI(N:Integer; X:Double): Double;
{
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
}
LABEL endlabel;
CONST IACC = 40; BIGNO = 1.e10; BIGNI = 1.e-10;
VAR   TOX, BIM, BI, BIP, BSI: Double;
      J, M: Integer;
Begin
      IF N=0 THEN
      begin
        BESSI := BESSI0(X);
        goto endlabel
      end;
      IF N=1 THEN
      begin
        BESSI := BESSI1(X);
        goto endlabel
      end;
      IF X=0.0 THEN
      begin
        BESSI:=0.0;
        goto endlabel
      end;
      TOX := 2.0/X;
      BIP := 0.0;
      BI  := 1.0;
      BSI := 0.0;
      M := 2*((N+Round(SQRT(IACC*N))));
      For J := M Downto 1 do
      begin
        BIM := BIP+J*TOX*BI;
        BIP := BI;
        BI  := BIM;
        IF ABS(BI) > BIGNO THEN
        begin
          BI  := BI*BIGNI;
          BIP := BIP*BIGNI;
          BSI := BSI*BIGNI
        end;
        IF J=N THEN BSI := BIP
      end;
      BESSI := BSI*BESSI0(X)/BI;
endlabel: End;
{ ----------------------------------------------------------------------
  Auxiliary Bessel functions for N=0, N=1  }
    FUNCTION BESSI0(X:Double): Double;
    VAR
      Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX: Double;
    Begin
      P1:=1.0; P2:=3.5156229; P3:=3.0899424; P4:=1.2067429;
      P5:=0.2659732; P6:=0.360768e-1; P7:=0.45813e-2;
      Q1:=0.39894228; Q2:=0.1328592e-1; Q3:=0.225319e-2;
      Q4:=-0.157565e-2; Q5:=0.916281e-2; Q6:=-0.2057706e-1;
      Q7:=0.2635537e-1; Q8:=-0.1647633e-1; Q9:=0.392377e-2;
      IF ABS(X) < 3.75 THEN
      begin
        Y:=(X/3.75)*(X/3.75);
        BESSI0:=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      end
      ELSE
      begin
        AX:=ABS(X);
        Y:=3.75/AX;
        BX:=EXP(AX)/SQRT(AX);
        AX:=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        BESSI0:=AX*BX
      end
    End;
{ --------------------------------------------------------------------- }
    FUNCTION BESSI1(X:Double): Double;
    VAR
      Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX: Double;
    Begin
      P1:=0.5; P2:=0.87890594; P3:=0.51498869; P4:=0.15084934;
      P5:=0.2658733e-1; P6:=0.301532e-2; P7:=0.32411e-3;
      Q1:=0.39894228; Q2:=-0.3988024e-1; Q3:=-0.362018e-2;
      Q4:=0.163801e-2; Q5:=-0.1031555e-1; Q6:=0.2282967e-1;
      Q7:=-0.2895312e-1; Q8:=0.1787654e-1; Q9:=-0.420059e-2;
      IF ABS(X) < 3.75 THEN
      begin
        Y:=(X/3.75)*(X/3.75);
        BESSI1:=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      end
      ELSE
      begin
        AX:=ABS(X);
        Y:=3.75/AX;
        BX:=EXP(AX)/SQRT(AX);
        AX:=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        BESSI1:=AX*BX
      end
    End;
{ --------------------------------------------------------------------- }

BEGIN
  N:=2;
  X:=0.75;

  Y := BESSI(N,X);

  writeln;
  writeln(' Bessel Function of order ', N,' for X=', X:8:4,':');
  writeln;
  writeln('      Y = ', Y);

  readkey;
  Donewincrt

END.

{ end of file tbessi.pas}