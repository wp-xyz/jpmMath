{------------------------------------------------------
! Conjugate Gradient Method for a Sparse Linear System
!------------------------------------------------------
! Ref.: "Numerical Recipes, Cambridge University Press
!        1986 (chapter 2.10)" [BIBLI 08].
!
!        Pascal Demo1 Version By J-P Moreau, Paris.
!                     (www.jpmoreau.fr)
!-------------------------------------------------------
! SAMPLE RUN:
!
! Solve sparse linear system (size=10):
!
!    2x1 -x10     =  0
!    2x2 -x9 -x10 =  0
!    2x3 -x8 -x9  =  0
!    2x4 -x7 -x6  =  0
!    2x5 -x6 -x7  =  0
!    -x5 +2x6     = 11
!    -x4 -x5 +2x7 =  0
!    -x3 -x4 +2x8 =  0
!    -x2 -x3 +2x9 =  0
!    -x1 -x2 +2x10 = 0
!
!
! Structure of sparse matrix A:
!
! X        X
!  X      XX
!   X    XX
!    X  XX
!     XXX
!     XX
!    XX X
!   XX   X
!  XX     X
! XX       X
!
! Any key to continue...
!
!
! Number of iterations: 11
!
! Solution vector:
!
!  1.0000  3.0000  5.0000  7.0000  9.0000
! 10.0000  8.0000  6.0000  4.0000  2.0000
!
!
! RSQ =   1.37619393585021E-0011
!
!
! Product A.X:
!
!  0.0000  0.0000  0.0000  0.0000  0.0000
! 11.0000  0.0000  0.0000  0.0000  0.0000
!
!
!------------------------------------------------------}
PROGRAM TSPARSE;
Uses WinCrt;

CONST N=10;

TYPE VECT = Array[1..N] of REAL;

VAR A: Array[1..N,1..N] of REAL;

    B, B1, X: VECT;
    RSQ: REAL;
    i,j: INTEGER;


Procedure ASUB(X:VECT; VAR V:VECT); Forward;
Procedure ATSUB(X:VECT; VAR V:VECT); Forward;


PROCEDURE SPARSE(B:VECT; N:INTEGER; VAR X:VECT; VAR RSQ:REAL);
{-------------------------------------------------------------------
!Solves the linear system A.x = b for the vector X of length N,
!given the right-hand vector B, and given two subroutines,
!ASUB(XIN,XOUT) and ATSUB(XIN,XOUT), which respectively calculate
!A.x and AT.x (AT for Transpose of A) for x given as first argument,
!returning the result in their second argument.
!These subroutines should take every advantage of the sparseness
!of the matrix A. On input, X should be set to a first guess of the
!desired solution (all zero components is fine). On output, X is
!the solution vector, and RSQ is the sum of the squares of the
!components of the residual vector A.x - b. If this is not small,
!then the matrix is numerically singular and the solution represents
!a least-squares approximation. 
!------------------------------------------------------------------}
LABEL 1,fin;
CONST EPS=1.E-6;
{maximum anticipated N, and r.m.s. accuracy desired}
VAR G,H,XI,XJ: VECT;
    EPS2,BSQ,RP,ANUM,ADEN,GAM,GG,DGG:REAL;
    IRST,ITER,J:INTEGER;
BEGIN
   EPS2:=N*EPS*EPS;
   IRST:=0;
1: IRST:=IRST+1;
   ASUB(X,XI);
   RP:=0.0;
   BSQ:=0.0;
   for J:=1 to N do
   begin
     BSQ:=BSQ+B[J]*B[J];
     XI[J]:=XI[J]-B[J];
     RP:=RP+XI[J]*XI[J]
   end;
   ATSUB(XI,G);
   for J:=1 to N do
   begin
     G[J]:=-G[J];
     H[J]:=G[J]
   end;
   for ITER:=1 to 10*N do   {main loop}
   begin
     ASUB(H,XI);
     ANUM:=0.0;
     ADEN:=0.0;
     for J:=1 to N do
     begin
       ANUM:=ANUM+G[J]*H[J];
       ADEN:=ADEN+XI[J]*XI[J]
     end;
     if ADEN=0.0 then writeln(' Very singular matrix.');
     ANUM:=ANUM/ADEN;
     for J:=1 to N do
     begin
       XI[J]:=X[J];
       X[J]:=X[J]+ANUM*H[J]
     end;
     ASUB(X,XJ);
     RSQ:=0.0;
     for J:=1 to N do
     begin
       XJ[J]:=XJ[J]-B[J];
       RSQ:=RSQ+XJ[J]*XJ[J]
     end;
     if (RSQ=RP) or (RSQ<=BSQ*EPS2) then
     begin
       writeln(' Number of iterations: ',ITER);
       goto fin  {normal return}
     end;
     if RSQ > RP then
     begin
       for J:=1 to N do X[J]:=XI[J];
       if IRST>=3  then goto fin;  {return if too many restarts}
       goto 1
     end;
     RP:=RSQ;
     ATSUB(XJ,XI);  {compute gradient for next iteration}
     GG:=0.0;
     DGG:=0.0;
     for J:=1 to N do
     begin
       GG:=GG+G[J]*G[J];
       DGG:=DGG+(XI[J]+G[J])*XI[J]
     end;
     if GG=0.0 then
     begin
       writeln(' Number of iterations:', ITER);
       goto fin  {rare but normal return}
     end;
     GAM:=DGG/GG;
     for J:=1 to N do
     begin
       G[J]:=-XI[J];
       H[J]:=G[J]+GAM*H[J]
     end
   end; {main loop}
   writeln(' Too many iterations.');
fin:END; 

{In these versions of ASUB1 and ATSUB1,
 we do not take any advantage of
 sparseness!  On the contrary,
 ASUB and ATSUB take into account the
 sparseness of the matrix A.

 N.B. For big values of N, the difference
 in computation time may be important. }

Procedure ASUB1(X:VECT; VAR V:VECT);
Var I,J:INTEGER;
Begin
  for I:=1 to N do
  begin
    V[I]:=0.0;
    for J:=1 to N do
      V[I]:=V[I]+A[I,J]*X[J];
  end
End;

Procedure ASUB(X:VECT; VAR V:VECT);
Var I,J:INTEGER;
Begin
  V[1]:=A[1,1]*X[1]+A[1,10]*X[10];
  For I:=2 to 5 do
    V[i]:=A[i,i]*X[i]+A[I,N-I+1]*X[N-I+1]+A[I,N-I+2]*X[N-I+2];
  V[6]:=A[6,6]*X[6]+A[6,5]*X[5];
  For I:=7 to 10 do
    V[i]:=A[i,i]*X[i]+A[I,N-I+1]*X[N-I+1]+A[I,N-I+2]*X[N-I+2]
End;

Procedure ATSUB1(X:VECT; VAR V:VECT);
Var I,J:INTEGER;
Begin
  for I:=1 to N do
  begin
    V[I]:=0.0;
    for J:=1 to N do
      V[I]:=V[I]+A[J,I]*X[J];
  end
End;

Procedure ATSUB(X:VECT; VAR V:VECT);
Var I,J:INTEGER;
Begin
  V[1]:=A[1,1]*X[1]+A[10,1]*X[10];
  For I:=2 to 5 do
    V[i]:=A[i,i]*X[i]+A[N-I+1,I]*X[N-I+1]+A[N-I+2,I]*X[N-I+2];
  V[6]:=A[6,6]*X[6]+A[5,6]*X[5];
  For I:=7 to 10 do
    V[i]:=A[i,i]*X[i]+A[N-I+1,I]*X[N-I+1]+A[N-I+2,I]*X[N-I+2];
End;

{show structure of sparse matrix A}
Procedure ShowMat;
Var i,j:INTEGER;
Begin
  writeln;
  for i:=1 to N do
  begin
    write(' ');
    for j:=1 to N do
      if A[i,j]<>0 then write('X') else write(' ');
    writeln
  end;
  writeln
End;

{main program}
BEGIN
{matrix A}
  For i:=1 to N do
    for j:=1 to N do
      A[i,j]:=0.0;
  A[1,1]:=2.0; A[1,10]:=-1.0;
  A[2,2]:=2.0; A[2,9]:=-1.0; A[2,10]:=-1.0;
  A[3,3]:=2.0; A[3,8]:=-1.0; A[3,9]:=-1.0;
  A[4,4]:=2.0; A[4,7]:=-1.0; A[4,8]:=-1.0;
  A[5,5]:=2.0; A[5,6]:=-1.0; A[5,7]:=-1.0;
  A[6,5]:=-1.0; A[6,6]:=2.0;
  A[7,4]:=-1.0; A[7,5]:=-1.0; A[7,7]:=2.0;
  A[8,3]:=-1.0; A[8,4]:=-1.0; A[8,8]:=2.0;
  A[9,2]:=-1.0; A[9,3]:=-1.0; A[9,9]:=2.0;
  A[10,1]:=-1.0; A[10,2]:=-1.0; A[10,10]:=2.0;
{vector B}
  For i:=1 to N do B[i]:=0.0; B[6]:=11.0;
{initial guess for X}
  for i:=1 to N do X[i]:=1.0;

  writeln;
  writeln(' Structure of sparse matrix A:');
  ShowMat;
  Write(' Any key to continue... ');
  Readkey;

  Clrscr;   {clear screen}
  writeln;
  SPARSE(B,N,X,RSQ);

  writeln;
  writeln(' Solution vector:');
  writeln;
  for i:=1 to 5 do write(X[I]:8:4); writeln;
  for i:=6 to N do write(X[I]:8:4); writeln;
  writeln;
  writeln(' RSQ = ', RSQ);
  writeln;

{verification A.x = B}
  ASUB(X,B1);
  writeln;
  writeln(' Product A.X:');
  writeln;
  for i:=1 to 5 do write(B1[I]:8:4); writeln;
  for i:=6 to N do write(B1[I]:8:4); writeln;
  writeln;

  Readkey; DoneWinCrt

END. 

{end of file tsparse.pas}