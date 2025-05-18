{***************************************************************************
*        Solve a symmetric linear system using subroutine SYMSOL           *
* ------------------------------------------------------------------------ *
* SAMPLE RUN:                                                              *
*  ( Solve symmetric linear system:                                        *
*      X1 -2X2 +3X3 -  4X4 +     5X5 -     6X6 = -3                        *
*    -2X1 +3X2 -5X3 +  7X4 -     4X5 +     8X6 =  7                        *
*     3X1 -5X2 +5X3 +  6X4 -     2X5 +      X6 =  8                        *
*    -4X1 +7X2 +6X3 -  9X4 +    20X5 +   0.5X6 = 20.5                      *
*     5X1 -4X2 -2X3 + 20X4 +    11X5 +0.3333X6 = 30.3333                   *
*    -6X1 +8X2 + X3 +0.5X4 +0.3333X5 +    13X6 = 16.8333 )                 *
*                                                                          *
*  SYMMETRIC SYSTEM TO SOLVE:                                              *
*  1.000000 -2.000000  3.000000 -4.000000  5.000000  -6.000000   -3.000000 *
* -2.000000  3.000000 -5.000000  7.000000 -4.000000   8.000000    7.000000 *
*  3.000000 -5.000000  5.000000  6.000000 -2.000000   1.000000    8.000000 *
* -4.000000  7.000000  6.000000 -9.000000 20.000000   0.500000   20.500000 *
*  5.000000 -4.000000 -2.000000 20.000000 11.000000   0.333300   30.333300 *
* -6.000000  8.000000  1.000000  0.500000  0.333300  13.000000   16.833300 *
*                                                                          *
*  SOLUTIONS:                                                              *
*  9.99999999999996E-0001                                                  *
*  1.00000000000002E+0000                                                  *
*  1.00000000000001E+0000                                                  *
*  1.00000000000001E+0000                                                  *
*  9.99999999999995E-0001                                                  *
*  9.99999999999989E-0001                                                  *
*                                                                          *
* ------------------------------------------------------------------------ *
* Reference:  From Numath Library By Tuan Dang Trong in Fortran 77         *
*             [BIBLI 18].                                                  *
*                                                                          *
*                                  TPW Release 1.0 By J-P Moreau, Paris    *
*                                           (www.jpmoreau.fr)              *
***************************************************************************}
PROGRAM TEST_SYMSOL;

Uses WinCrt;

Const NMAX = 75;

Type
      pMAT = ^MAT;
      MAT = Array[1..NMAX,1..NMAX] of Double;
      pVEC = ^VEC;
      VEC = Array[1..NMAX] of Double;

Var
      A    : pMAT;
      B,X,D: pVEC;
      I,J,N: Integer;


    Procedure SYMSOL(N:Integer; A:pMAT; B:pVEC; Var X:pVEC; D:pVEC);
{-----------------------------------------------------------------------
!     SOLVE A SYMMETRIC LINEAR SYSTEM A.X = B
!     A IS SUPPOSED WELL BALANCED AND CAN BE NOT POSITIVE DEFINITE,
!     THE COEFFICIENTS OF A ARE STORED IN THE UPPER TRIANGULAR HALF
!     OF THE MATRIX AND ARE LOST DURING THE PROCESS.
!     THE USED METHOD IS BASED UPON THE DECOMPOSITION OF THE INITIAL 
!     MATRIX INTO A PRODUCT OF 3 MATRICES: A = (U TRANSPOSE).D.U
!     INPUTS:
!     N       NUMBER OF LINES OF THE LINEAR SYSTEM
!     A       TABLE OF DIMENSION N X N STORING THE COEFFICIENTS
!             OF THE SYSTEM
!     B       VECTOR OF DIMENSION N (SECOND MEMBER)
!     OUTPUT:
!     X       SYSTEM SOLUTION
!     WORKING ZONE:
!     D       DIAGONAL D OF DECOMPOSITION
!----------------------------------------------------------------------}
    Label  2, 4, 6, 7, 8, 10;
    Var
        S: Double;
        I,K,L,M: Integer;
    Begin
      For I := 1 to N do
      begin
        L := I-1;
        M := I+1;
        S := A^[I,I];
        IF I=1 THEN GOTO 2;
        For K := 1 to L do
          S := S-A^[K,I]*D^[K]*A^[K,I];
2:      D^[I] := S;
        IF M > N THEN GOTO 6;
        For J := M to N do
        begin
          S := A^[I,J];
          IF I=1 THEN GOTO 4;
          For K := 1 to L do
            S := S-A^[K,I]*D^[K]*A^[K,J];
4:        A^[I,J] := S/D^[I];
        end
      end;
6:    for I := 1 to N do
      begin
        S := B^[I];
        L := I-1;
        IF I=1 THEN GOTO 8;
        For K := 1 to L do
          S := S-A^[K,I]*X^[K];
8:      X^[I] := S;
      end;
      For L := 1 to N do
      begin
        I := N-L+1;
        S := X^[I]/D^[I];
        M := I+1;
        IF M>N THEN GOTO 10;
        For K := M to N do
          S := S-A^[I,K]*X^[K];
10:     X^[I] := S;
      end
    End; {SYMSOL}

{main program}
BEGIN

  N := 6;    {size of linear system}

  {allocate memory for matrix A and vectors B, X, D}
  New(A); New(B); New(X); New(D);

  {define upper half part of matrix A}
  A^[1,1]:= 1.0; A^[1,2]:=-2.0; A^[1,3]:=3.0; A^[1,4]:=-4.0; A^[1,5]:=5.0; A^[1,6]:=-6.0;
  A^[2,2]:= 3.0; A^[2,3]:=-5.0; A^[2,4]:=7.0; A^[2,5]:=-4.0; A^[2,6]:=8.0;                
  A^[3,3]:= 5.0; A^[3,4]:= 6.0; A^[3,5]:=-2.0; A^[3,6]:=1.0;
  A^[4,4]:=-9.0; A^[4,5]:=20.0; A^[4,6]:=0.50;
  A^[5,5]:=11.0; A^[5,6]:= 0.3333;
  A^[6,6]:=13.0;

  {define 2nd member}
  B^[1]:=-3.0;
  B^[2]:=7.0;
  B^[3]:=8.0;
  B^[4]:=20.5;
  B^[5]:=30.3333;
  B^[6]:=16.8333;

  {optional section to print full data}
  For I:=2 to N do
    For J:=1 to I-1 do
      A^[I,J] := A^[J,I];

  writeln;
  writeln(' SYMMETRIC SYSTEM TO SOLVE:');
  For I:=1 to N do
  begin
    For J:=1 to N do write(A^[I,J]:10:6);
    writeln('  ',B^[I]:10:6)
  end;
  {end optional section}

  {call solver for symmetric linear systems}
  SYMSOL(N,A,B,X,D);

  {print solutions}
  writeln(' SOLUTIONS:');
  For I:=1 to N do writeln(X^[I]);

  {ending section}
  ReadKey;
  {free memory}
  Dispose(A); Dispose(B); Dispose(X); Dispose(D);
  DoneWinCrt

END.

{end of file tsymsol.pas}