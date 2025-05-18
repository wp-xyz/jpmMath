{ ---------------------------------------------------------------------
! Program to integrate a user-defined function f(x) from x1 to x2 by
! the QANC8 subroutine with control of absolute and relative precisions
! ---------------------------------------------------------------------
! SAMPLE RUN:
!
!  (Integrate function cos(x) - 2 sin(x) from x:=0 to x:=1,
!   with a precision <= 1e-10).
! 
!  Integral Value = -7.79244034558251E-0002
!
!  Estimated error =  1.69979160539506E-0016
!
!  Error code =  0.00000000000000E+000
!
!  Number of function evaluations = 33
!
! ---------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                               TPW Release 1.0 By J-P Moreau, Paris
!                                       (www.jpmoreau.fr)
! --------------------------------------------------------------------}
PROGRAM TQANC8;

Uses WinCrt;

VAR
    AERROR,CODE,ERROR,RERROR,X1,X2,VAL: Double;
    NBRF: Integer;


    FUNCTION FCT(X:Double): Double;
    Begin
  
      FCT:=COS(X)-2*SIN(X)

    End;

    Function IPower(x,n:integer): Integer;
    {x power n, n must be >= 0 else error}
    var i,result : integer;
    begin
      result :=1;
      if n=0 then
      begin
        IPower:=result;
        exit;
      end
      else
        for i:=1 to n do
          result := x * result;
      IPower := result
    end;

    Function MAX(a,b:Double): Double;
    Begin
      if a>=b then MAX:=a
              else MAX:=b
    End;

{ ---------------------------------------------------------------------}
    PROCEDURE QANC8 (A,B,AERR,RERR:Double; VAR RES,ERR: Double; VAR
                     NBF: Integer; VAR FLG:Double);
{
!     INTEGRATE A REAL FUNCTION FCT(X) FROM X:=A TO X:=B, WITH
!     GIVEN ABSOLUTE AND RELATIVE PRECISIONS, AERR, RERR. 

!     INPUTS:
!     FCT     EXTERNAL USER-DEFINED FUNCTION FOR ANY X VALUE
!             IN INTERVAL (A,B)
!     A,B     LIMITS OF INTERVAL
!     AERR,RERR   RESPECTIVELY ABSOLUTE ERROR AND RELATIVE ERROR
!                 REQUIRED BY USER
!
!     OUTPUTS:
!     RES     VALUE OF INTEGRAL
!     ERR     ESTIMATED ERROR
!     NBF     NUMBER OF NECESSARY FCT(X) EVALUATIONS
!     FLG     INDICATOR
!             := 0.0       CORRECT RESULT
!             := NNN.RRR   NO CONVERGENCE DU TO A SINGULARITY.
!             THE SINGULAR POINT ABCISSA IS GIVEN BY FORMULA:
!             XS := B-.RRR*(B-A)

!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL, INC.
! ----------------------------------------------------------------------}
    Label 30, 50, 60, 62, 70, 72, 75, 80, 82, EndProc;
    Var
      QR: Array[1..31] of REAL;
      F, X: Array[1..16] of REAL;
      FS, XS: Array[1..8,1..30] of REAL;
      COR,F0,PAS,PAS1,QP,SUM,W0,W1,W2,W3,W4,X0: Double;
      ERR1,QD,QL,QN,TEMP,TOL1: Double;
      I,J,L,LMIN,LMAX,LOUT,NFIN,NIM,NMAX: Integer;
    Begin
      LMIN := 1;
      LMAX := 30;
      LOUT := 6;
      NMAX := 5000;
      NFIN := NMAX-8*(LMAX-LOUT+IPower(2,LOUT+1));
      W0   :=   3956.0/14175.0;
      W1   :=  23552.0/14175.0;
      W2   :=  -3712.0/14175.0;
      W3   :=  41984.0/14175.0;
      W4   := -18160.0/14175.0;
      FLG  := 0.0;
      RES  := 0.0;
      COR  := 0.0;
      ERR  := 0.0;
      SUM  := 0.0;
      NBF  := 0;
      if A=B then goto EndProc;
      L := 0;
      NIM := 1;
      X0  := A;
      X[16] := B;
      QP  := 0.0;
      F0   := FCT(X0);
      PAS1  := (B-A)/16.0;
      X[8]  := (X0+X[16])*0.5;
      X[4] := (X0+X[8])*0.5;
      X[12] := (X[8]+X[16])*0.5;
      X[2]  := (X0+X[4])*0.5;
      X[6] := (X[4]+X[8])*0.5;
      X[10]:= (X[8]+X[12])*0.5;
      X[14] := (X[12]+X[16])*0.5;
      J:=2;
      While J<=16 do
      begin
        F[J] := FCT(X[J]);
        Inc(J,2)
      end;
      NBF := 9;
30:   X[1] := (X0+X[2])*0.5;
      F[1] := FCT(X[1]);
      J:=3;
      While J<=15 do
      begin
        X[J] := (X[J-1]+X[J+1])*0.5;
        F[J] := FCT(X[J]);
        Inc(J,2)
      end;
      NBF := NBF+8;
      PAS := (X[16]-X0)/16.0;
      QL  := (W0*(F0+F[8])+W1*(F[1]+F[7])+W2*(F[2]+F[6])+W3*(F[3]+F[5])+W4*F[4])*PAS;
      QR[L+1] := (W0*(F[8]+F[16])+W1*(F[9]+F[15])+W2*(F[10]+F[14])+W3*(F[11]+F[13])+W4*F[12])*PAS;
      QN := QL + QR[L+1];
      QD := QN - QP;
      SUM := SUM + QD;
      ERR1 := ABS(QD)/1023.0;
      TOL1 := MAX(AERR,RERR*ABS(SUM))*(PAS/PAS1);
      IF L < LMIN THEN GOTO 50;
      IF L >= LMAX THEN GOTO 62;
      IF NBF > NFIN THEN GOTO 60;
      IF ERR1 <= TOL1 THEN GOTO 70;
50:   NIM := 2*NIM;
      L := L+1;
      For I := 1 to 8 do
      begin
        FS[I,L] := F[I+8];
        XS[I,L] := X[I+8]
      end;
      QP := QL;
      For I := 1 to 8 do
      begin
        F[18-2*I] := F[9-I];
        X[18-2*I] := X[9-I]
      end;
      GOTO 30;

60:   NFIN := 2*NFIN;
      LMAX := LOUT;
      FLG := FLG + (B-X0)/(B-A);
      GOTO 70;
62:   FLG := FLG + 1.0;
70:   RES := RES + QN;
      ERR := ERR + ERR1;
      COR := COR + QD/1023.0;
72:   IF NIM=2*(NIM DIV 2) THEN GOTO 75;
      NIM := NIM DIV 2;
      L := L-1;
      GOTO 72;
75:   NIM := NIM+1;
      IF L <= 0 THEN GOTO 80;
      QP := QR[L];
      X0 := X[16];
      F0 := F[16];
      For I := 1 to 8 do
      begin
        F[2*I] := FS[I,L];
        X[2*I] := XS[I,L]
      end;
      GOTO 30;

80:   RES := RES + COR;
      IF ERR=0.0 THEN GOTO EndProc;
82:   TEMP := ABS(RES) + ERR;
      IF TEMP<>ABS(RES) THEN GOTO EndProc;
      ERR := 2.0*ERR;
      GOTO 82;
EndProc: End;


{main program}
BEGIN
  X1:=0.0;
  X2:=1.0;
  AERROR:=1E-9;
  RERROR:=1E-10;

  QANC8(X1,X2,AERROR,RERROR,VAL,ERROR,NBRF,CODE);

  writeln;
  writeln(' Integral Value = ', VAL);
  writeln;
  writeln(' Estimated error = ', ERROR);
  writeln;
  writeln(' Error Code = ', CODE);
  writeln;
  writeln(' Number of function evaluations = ', NBRF);
  ReadKey;
  DoneWinCrt
     
END.

{ end of file tqanc8.pas}