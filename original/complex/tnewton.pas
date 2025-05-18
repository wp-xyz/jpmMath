{***************************************************
*   Find all roots of a complex polynomial using   *
*   Newton's iterative formulation.                *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
* (Find roots of complex polynomial:               *
*  (2.5 + I) X3 + (7.5 - 12 I) X2 + (-3.75 + 0.4 I *
*  ) X + (4.25 - I)  )                             *
*                                                  *
*  Error code = 0                                  *
*                                                  *
*  Complex roots are:                              *
*  2.59967402504241E-0001 -3.99661390631072E-0001  *
* -1.11845742991334E-0001  6.68006959856486E-0001  *
* -1.07915614227153E+0000  4.90406822387803E+0000  *
*                                                  *
* ------------------------------------------------ *
* Ref.: From Numath Library By Tuan Dang Trong in  *
*       Fortran 77 [BIBLI 18].                     *
*                                                  *
*                TPW Release By J-P Moreau, Paris. *
*                        (www.jpmoreau.fr)         *
***************************************************}   
PROGRAM TEST_NEWTON;

Uses WinCrt;   {Borland only}

Const  NMAX = 50;

Type
     pVEC = ^VEC;
     VEC = Array[1..NMAX] of Double;

Var
     I,IERR, N  : Integer; {IERR: error code, N: order of polynomial}
     A,B,RR,RI  : pVEC;    {A: real parts of complex coefficients
                            B: imaginary parts of complex coefficients
                            (stored in decreasing powers)
                            RR: real parts of roots
                            RI: imaginary parts of roots}
     W1,W2,W3,W4: pVEC;    {working zones}



     Procedure CNEWTON (N:Integer; A,B: pVEC; Var X,Y:pVEC; Var IER:Integer;
                        G,H,T,F:pVEC);
{---------------------------------------------------------------------
!     CALCULATE ALL THE ROOTS OF A COMPLEX POLYNOMIAL USING
!     NEWTON'S ITERATIVE FORMULATION
!
!     INPUTS:
!     N        ORDER OF POLYNOMIAL
!     A,B      TABLES OF DIMENSION N+1, RESPECTIVELY REAL,
!              IMAGINARY PART OF POLYNOMIAL COMPLEX COEFFICIENTS 
!              C := A+I*B, STORED IN DECREASING POWERS.
!     OUTPUTS:
!     X,Y      TABLES OF DIMENSION N, RESPECTIVELY REAL,
!              IMAGINARY PART OF COMPLEX ROOTS FOUND,  Z := X+I*Y
!     IER      ERROR CODE := 0, CONVERGENCE OK
!                         := I, NO CONVERGENCE FOR I-TH ROOT
!     WORKING ZONES:
!     G,H,T,F  TABLES OF DIMENSION N
!----------------------------------------------------------------------
!     REFERENCE:
!     E.DURAND -SOLUTIONS NUMERIQUES DES EQUATIONS ALGEBRIQUES,
!               TOME 1, MASSON & CIE, 1960, PP.260-266
!----------------------------------------------------------------------}
    Label 5,7,10,15,30,40,45,55,60,70,75;
    Var
        I,IT1,IT2,K,K1,K2,L,M: Integer;
        EPS,EPS1,P,Q,R,TOL,X0,X1,Y0,Y1: Double;
        CABS,DX,DY,GD,GX,GY,R1,R2,R3,X2,Y2: Double;
    Begin
      IT1:=25; IT2:=100; TOL:=1e-2; EPS:=0.0;
      IER := 0;
      IF EPS <> 0.0 THEN GOTO 7;
      EPS := 1.0;
5:    EPS := 0.50*EPS;
      EPS1 := EPS+1.0;
      IF EPS1 > 1.0 THEN GOTO 5;
7:    L := N;
      X0 := 1.0;
      Y0 := 1.0;
      G^[1] := A^[1];
      T^[1] := A^[1];
      H^[1] := B^[1];
      F^[1] := B^[1];
10:   IF L = 1 THEN GOTO 55;
      M := L+1;
      I := N-L+1;
      K1 := 0;
      K2 := 0;
      X1 := X0;
      Y1 := Y0;
15:   For K := 2 to M do
      begin
        G^[K] := A^[K]+X1*G^[K-1]-Y1*H^[K-1];
        H^[K] := B^[K]+Y1*G^[K-1]+X1*H^[K-1]
      end;
      For K := 2 to L do
      begin
        T^[K] := G^[K]+X1*T^[K-1]-Y1*F^[K-1];
        F^[K] := H^[K]+Y1*T^[K-1]+X1*F^[K-1]
      end;
      GD := T^[L]*T^[L]+F^[L]*F^[L];
      GX := G^[M]*T^[L]+H^[M]*F^[L];
      GY := H^[M]*T^[L]-G^[M]*F^[L];
      DX := -GX/GD;
      DY := -GY/GD;
      X2 := X1+DX;
      Y2 := Y1+DY;
      R1 := (ABS(DX)+ABS(DY))/(ABS(X2)+ABS(Y2));
      K2 := K2+1;
      X1 := X2;
      Y1 := Y2;
      IF R1 < TOL THEN GOTO 30;
      IF K2 <= IT2 THEN GOTO 15;
      IER := I;
      GOTO 60;
30:   K1 := K1+1;
      IF K1 <= IT1 THEN GOTO 15;
      X^[I] := X2;
      Y^[I] := Y2;
      IF ABS(X2) <= EPS THEN X2 := EPS;
      R2 := ABS(Y2/X2);
      IF R2 < 0.1 THEN GOTO 40;
      IF R2*EPS > 1.0 THEN X^[I] := 0.0;
      X0 := X2;
      Y0 := Y2;
      GOTO 45;
40:   IF R2 <= EPS THEN Y^[I] := 0.0;
      X0 := X2+Y2;
      Y0 := X2-Y2;
45:   For K := 1 to L do
      begin
        A^[K] := G^[K];
        B^[K] := H^[K]
      end;
      L := L-1;
      GOTO 10;
55:   R3 := G^[1]*G^[1]+H^[1]*H^[1];
      X^[N] := -(G^[1]*G^[2]+H^[1]*H^[2])/R3;
      Y^[N] :=  (H^[1]*G^[2]-G^[1]*H^[2])/R3;
      IF ABS(X^[N]) <= EPS THEN  X^[N] := EPS;
      R2 := ABS(Y^[N]/X^[N]);
      IF R2*EPS > 1.0 THEN  X^[N] := 0.0;
      IF R2 <= EPS THEN  Y^[N] := 0.0;
60:   For L:=1 to N do
      begin
        P:=X^[L];
        Q:=Y^[L];
        R:=P*P+Q*Q;
        IF L = 1 THEN GOTO 70;
        For I:=L Downto 2 do
        begin
          CABS:=Sqr(X^[I-1])+Sqr(Y^[I-1]);
          IF R >= CABS THEN GOTO 75;
          X^[I]:=X^[I-1];
          Y^[I]:=Y^[I-1]
        end;
70:     I:=1;
75:     X^[I]:=P;
        Y^[I]:=Q
      end
    End; {CNEWTON}

{main program}
BEGIN

  N := 3;   {order of polynomial}

  {allocate memory for Tables}
  New(A); New(B); New(RR); New(RI);
  New(W1); New(W2); New(W3); New(W4);

  {define vectors A, B of coefficients}
  A^[1]:= 2.5;  B^[1]:=1.0;
  A^[2]:= 7.5;  B^[2]:=-12.0;
  A^[3]:=-3.75; B^[3]:=0.4;
  A^[4]:= 4.25; B^[4]:=-1.0;

{ Example #2 (real roots are -3, 1, 2)
  A^[1]:= 1.;   B^[1]:=0.0;
  A^[2]:= 0.;   B^[2]:=0.0;
  A^[3]:=-7.;   B^[3]:=0.0;
  A^[4]:= 6.;   B^[4]:=0.0;  }

  {call main procedure}
  CNEWTON(N,A,B,RR,RI,IERR,W1,W2,W3,W4);

  {print results}
  writeln;
  writeln(' Error code = ', IERR);
  writeln;
  writeln(' Complex roots are:');

  For I:=1 to N do
    writeln(' ', RR^[I], ' ', RI^[I]);
  writeln;

  {ending section}
  ReadKey;
  {free memory}
  Dispose(A); Dispose(B); Dispose(RR); Dispose(RI);
  Dispose(W1); Dispose(W2); Dispose(W3); Dispose(W4);

  {Borland only: close application window}
  DoneWinCrt

END.

{end of file tnewton.pas}