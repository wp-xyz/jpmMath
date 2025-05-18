{****************************************************
*   Find all roots of a complex polynomial using    *
*   Laguerre formulation in complex domain.         *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
* (Find roots of complex polynomial:                *
*  (2.5 + I) X3 + (7.5 - 12 I) X2 + (-3.75 + 0.4 I  *
*  ) X + (4.25 - I)  )                              *
*                                                   *
*  Error code =         0                           *
*                                                   *
*  Complex roots are:                               *
*                                                   *
* (-1.07915614173380E+0000, 4.90406822631543E+0000) *
* (-1.11845747804128E-0001, 6.68006961189352E-0001) *
* ( 2.59967402738909E-0001,-3.99661390407884E-0001) *
*                                                   *
* ------------------------------------------------- *
* Ref.: From Numath Library By Tuan Dang Trong in   *
*       Fortran 77 [BIBLI 18].                      *
*                                                   *
*                TPW Release By J-P Moreau, Paris.  *
*                        (www.jpmoreau.fr)          *
****************************************************}   
PROGRAM TEST_LAGUERRE;

Uses WinCrt;

Const N=3;  

Type
     Complex = Array[1..2] of Double;
     CTab = Array[0..N] of Complex;

Var
     A, RAC, W1, W2, W3: CTab;
     I, IMP, ITMAX: Integer;
     EPS,EPS2: Double;
     CZERO: Complex;


    {*** Utility procedures for complex numbers ***} 

    {ABSOLUTE VALUE OF A COMPLEX NUMBER Z=X+I*Y }
    FUNCTION CABS(Z:COMPLEX): Double;
    Label 10, 20;
    Var
        U,V,Q,S: Double;
    Begin
      U := ABS(Z[1]);
      V := ABS(Z[2]);
      S := U+V;
{--------------------------------------------------------------------
      S*1.0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
      TRUE FLOATING ZERO
--------------------------------------------------------------------}
      S := S*1.0;
      if S = 0.0 then goto 20;
      if U > V  then goto 10;
      Q := U/V;
      CABS:=V*sqrt(1.0+Q*Q);
      exit;
10:   Q := V/U;
      CABS:=U*sqrt(1.0+Q*Q);
      exit;
20:   CABS:=0.0
    End;

    {z3=z1+z2}
    Procedure CAdd(z1,z2:Complex; Var z3:Complex);
    Begin
      z3[1]:=z1[1]+z2[1];
      z3[2]:=z1[2]+z2[2]
    End;

    {let z1=z}
    Procedure CAssign(z:Complex; Var z1:Complex);
    Begin
      z1[1]:=z[1]; z1[2]:=z[2]
    End;

    {Z=Z1/Z2}
    Procedure CDIV(Z1,Z2:Complex; Var Z:Complex);
    Var D:double;
    Begin
      D:=Z2[1]*Z2[1]+Z2[2]*Z2[2];
      if D<1e-12 then exit;
      Z[1]:=(Z1[1]*Z2[1]+Z1[2]*Z2[2])/D;
      Z[2]:=(Z1[2]*Z2[1]-Z1[1]*Z2[2])/D
    End;

    {z1=exp(z) }
    Procedure CEXP(Z:Complex; Var Z1:Complex);
    Var tmp:double;
    Begin
      if exp(Z[1]) > 1e16 then tmp:=1e16 else tmp:=exp(Z[1]);
      Z1[1]:=tmp*cos(Z[2]);
      Z1[2]:=tmp*sin(Z[2])
    End;

    {Z=Z1*Z2}
    Procedure CMUL(Z1,Z2:Complex; Var Z:Complex);
    Begin
      Z[1]:=Z1[1]*Z2[1] - Z1[2]*Z2[2];
      Z[2]:=Z1[1]*Z2[2] + Z1[2]*Z2[1]
    End;

    {z=x+iy }
    Procedure CMPLX(Var z:Complex; x,y:Double);
    Begin
      z[1]:=x; z[2]:=y
    End; 

    {print a complex number}
    Procedure CPrint(z:Complex);
    Begin
      writeln(' (',z[1],',',z[2],')')
    End;

    Procedure CSQRT(z:Complex; Var z1:Complex);
{   SQUARE ROOT OF A COMPLEX NUMBER  A+I*B = SQRT(X+I*Y) }
    Var
        X,Y,A,B: Double;
    Begin
      X:=z[1]; Y:=z[2];
      IF (X = 0.0) AND (Y = 0.0) THEN
      begin
        A:=0.0;
        B:=0.0
      end
      ELSE
      begin
        A:=SQRT(ABS(X)+CABS(z)*0.5);
        IF X >= 0.0 THEN
          B:=Y/(A+A)
        ELSE
          IF Y < 0.0 THEN
            B:=-A
          ELSE
          begin
            B:=A;
            A:=Y/(B+B)
          end
      end;
      CMPLX(z1,A,B)
    End;

    {z3=z1-z2}
    Procedure CSUB(z1,z2:Complex; Var z3:Complex);
    Begin
      z3[1]:=z1[1]-z2[1];
      z3[2]:=z1[2]-z2[2]
    End;


    Procedure CLAGUE(N: Integer; A: CTab; ITMAX: Integer; EPS,EPS2: Double;
                     Var IMP: Integer; Var X: CTab; B, C, D: CTab);
{=================================================================
!     ROOTS OF A COMPLEX COEFFICIENTS POLYNOMIAL
!     BY LAGUERRE FORMULA IN COMPLEX DOMAIN
!=================================================================
!     CALLING MODE:
!       CLAGUE(N,A,ITMAX,EPS,EPS2,IMP,X,B,C,D);
!     INPUTS:
!     N : ORDER OF POLYNOMIAL
!     A : TABLE OF SIZE (0:N) OF COMPLEX COEFFICIENTS STORED IN
!         DECREASING ORDER OF X POWERS
!     ITMAX: MAXIMUN NUMBER OF ITERATIONS FOR EACH ROOT
!     EPS:   MINIMAL RELATIVE ERROR
!     EPS2:  MAXIMAL RELATIVE ERROR
!     X(0):  APPROXIMATE VALUE OF FIRST ROOT (=0. GENERALLY)
!     OUTPUTS:
!     IMP:  FLAG = 0  CONVERGENCE WITH AT LEAST EPS2 PRECISION
!                = 1  NO CONVERGENCE
!     X:  TABLE OF SIZE (0:N) OF FOUND ROOTS
!         X(I),I=1,N
!     WORKING ZONE:
!     W:  TABLE OF SIZE (0:3*N), here divided into B, C, D.
! ----------------------------------------------------------------
!     REFERENCE:
!     E.DURAND. SOLUTIONS NUMERIQUES DES EQUATIONS ALGEBRIQUES
!               TOME I, MASSON & CIE PAGES 269-270
!=================================================================}
    Label 1,2,3,4, 10 {return};
    Var
      XK,XR,F,FP,FS,H,DEN,SQ,D2,TMP1,TMP2,TMP3: Complex;
      EPS1,TEST: Double;
      I,IK,IT: Integer;
    Begin

      IMP:=0;
      CAssign(A[0],B[0]);
      CAssign(B[0],C[0]);
      CAssign(C[0],D[0]);

      IF N = 1 THEN
      begin
        {X[1]:=-A[1]/A[0]}
        CDIV(A[1],A[0],X[1]);
        X[1][1]:=-X[1][1];  {change sign of X[1]}
        X[1][2]:=-X[1][2];
        goto 10; {return}
      end;
      
      CAssign(X[0],XK);
1:    IK:=0;
      EPS1:=EPS;
      IT:=0;

2:    For I:=1 to N do
      begin
        {B[I]:=A[I]+XK*B[I-1];}
        CMUL(XK,B[I-1],TMP1);
        CADD(A[I],TMP1,B[I]);
        IF I <= N-1 THEN
        begin
          {C[I]:=B[I)+XK*C[I-1) }
          CMUL(XK,C[I-1],TMP1);
          CADD(B[I],TMP1,C[I]);
        end;
        IF I <= N-2 THEN
        begin
          {D[I]:=C[I)+XK*D[I-1) }
          CMUL(XK,D[I-1],TMP1);
          CADD(C[I],TMP1,D[I]);
        end
      end;

      CAssign(B[N],F);
      CAssign(C[N-1],FP);
      CAssign(D[N-2],FS);
      {H:=((N-1)*FP)^2-N*(N-1)*F*FS }
      CMPLX(TMP1,1.0*(N-1),0.0); CMUL(TMP1,FP,TMP2);
      CMUL(TMP2,TMP2,TMP1); {TMP1=((N-1)*FP)^2 }
      CMPLX(TMP2,1.0*N*(N-1),0.0); CMUL(TMP2,F,TMP3);
      CMUL(TMP3,FS,TMP2);   {TMP2=N*(N-1)*F*FS }
      CSUB(TMP1,TMP2,H);
      IF CABS(H) = 0.0 THEN
      begin
        For I:=N Downto 1 do
        begin
          {X[I]:=-A[1]/(N*A[0]) }
          CMPLX(TMP1,1.0*N,0.0); CMUL(TMP1,A[0],TMP2);
          CDIV(A[1],TMP2,X[I]);
          X[I][1]:=-X[I][1];
          X[I][2]:=-X[I][2]
        end;
        goto 10 {return}
      end;
      IF CABS(FP) = 0.0 THEN 
        CSQRT(H,DEN)
      ELSE
      begin
        CSQRT(H,SQ);
        CADD(FP,SQ,DEN);
        CSUB(FP,SQ,D2);
        IF CABS(DEN) < CABS(D2) THEN CAssign(D2,DEN)
      end;
      IF CABS(DEN) = 0.0 THEN
      begin
        XK[1]:=XK[1] + 0.1;
        XK[2]:=XK[2] + 0.1;
        GOTO 2
      end;
      IK:=IK+1;
      {XR:=XK-N*F/DEN }
      CMPLX(TMP1,1.0*N,0.0); CMUL(TMP1,F,TMP2);
      CDIV(TMP2,DEN,TMP1);
      CSUB(XK,TMP1,XR);
      CSUB(XR,XK,TMP1);
      TEST:=CABS(TMP1);
      CAssign(XR,XK);
      IF TEST < EPS1*CABS(XR) THEN GOTO 3;
{     WRITELN(IK,H[1],H[2],DEN[1],DEN[2],XK[1],XK[2],TEST); }
      IF IK <= ITMAX THEN GOTO 2;
      EPS1:=EPS1*10.0;
      IK:=0;
      IF EPS1 < EPS2 THEN GOTO 2;
      IMP:=1;
      goto 10;

3:    CAssign(XR,X[N]);
      N:=N-1;
      IF N = 1 THEN GOTO 4;
      IF N <= 0 THEN goto 10;
      For I:=0 to N do CAssign(B[I],A[I]);
      GOTO 1;

4:    {X[N]:=-B[1]/B[0] }
      CDIV(B[1],B[0],X[N]);
      X[N][1]:=-X[N][1];
      X[N][2]:=-X[N][2];

10: End;


{main program}
BEGIN

  CMPLX(CZERO,0.0,0.0);

{ Example #1 }
  CMPLX(A[0],2.5,1.0);
  CMPLX(A[1],7.5,-12.0);
  CMPLX(A[2],-3.75,0.4);
  CMPLX(A[3],4.25,-1.0);

{ Example #2 (real roots are -3, 1, 2)
  CMPLX(A[0],1.0,0.0);
  CMPLX(A[1],0d0,0d0);
  CMPLX(A[2],-7.0,0.0);
  CMPLX(A[3],6.0,0.0);
}

  ITMAX:=10;                  {Maximum number of iterations}
  EPS:=1e-8; EPS2:=1e-6;      {Minimum, maximum relative error}
  CAssign(CZERO,RAC[0]);      {Approximate value of 1st root}

  {call complex Laguerre procedure}
  CLAGUE(N,A,ITMAX,EPS,EPS2,IMP,RAC,W1,W2,W3);

  {print results}
  writeln;
  writeln(' Error code =', IMP);
  writeln;
  writeln(' Complex roots are:');
  writeln;
  For I:=1 to N do CPrint(RAC[I]);
  writeln;

  ReadKey;
  DoneWinCrt

END.

{end of file tclague.pas}