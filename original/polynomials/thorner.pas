{************************************************************************
*     EVALUATE A POLYNOMIAL AND ITS DERIVATIVES BY HORNER'S METHOD      *
* --------------------------------------------------------------------- *
* SAMPLE RUN:                                                           *
*                                                                       *
*    EVALUATE A POLYNOMIAL AND ITS DERIVATIVES BY HORNER'S METHOD       *
*                                                                       *
*            Example:  P(X) = X^4 + 2X^3 +3X^2 + 4X + 5                 *
*                                                                       *
*      X          P(X)          P'(X)        P"(X)         P"'(X)       *
* ------------------------------------------------------------------    *
*  0.00000       5.00000      4.00000       6.00000      12.00000       *
*  0.20000       5.93760      5.47200       8.88000      16.80000       *
*  0.40000       7.23360      7.61600      12.72000      21.60000       *
*  0.60000       9.04160     10.62400      17.52000      26.40000       *
*  0.80000      11.55360     14.68800      23.28000      31.20000       *
*  1.00000      15.00000     20.00000      30.00000      36.00000       *
*  1.20000      19.64960     26.75200      37.68000      40.80000       *
*  1.40000      25.80960     35.13600      46.32000      45.60000       *
*  1.60000      33.82560     45.34400      55.92000      50.40000       *
*  1.80000      44.08160     57.56800      66.48000      55.20000       *
*  2.00000      57.00000     72.00000      78.00000      60.00000       *
* -------------------------------------------------------------------   *
*                                                                       *
* Press any key to continue                                             *
*                                                                       *
* --------------------------------------------------------------------- *
* Reference: From Numath Library By Tuan Dang Trong in Fortran 77       *
*            [BIBLI 18].                                                *  
*                                                                       *
*                                TPW Release 1.0 By J-P Moreau, Paris   *
*                                         (www.jpmoreau.fr)             *
************************************************************************}  
PROGRAM TEST_HORNER;

Uses WinCrt;

Const  NMAX = 25;

Type
       pVEC = ^VEC;
       VEC  = Array[1..NMAX] of Double;

Var
       A, RES: pVEC;
       DX,X,X1,X2: Double;
       I,J,K,N,NP:Integer;
       FLG: Boolean;


    Procedure HORNER (N:Integer; A:pVEC; X0:Double; K:Integer; B:Boolean; Var R:pVEC);
{---------------------------------------------------------------------
!     EVALUATE A POLYNOMIAL AND ITS DERIVATIVES BY HORNER'S METHOD
!
!     INPUTS:
!     N       ORDER OF POLYNOMIAL
!     A       VECTOR OF SIZE N+1 STORING THE COEFFICIENTS OF
!             POLYNOMIAL IN DECREASING ORDERS OF POWERS
!             I.E. P(X) = A(1)*X**N+A(2)*X**(N-1)+...+A(N)*X+A(N+1)
!     X0      GIVEN ARGUMENT
!     K       MAXIMUM ORDER OF DERIVATIVES ASKED FOR
!     B       FLAG
!             = .TRUE.  EVALUATION ONLY
!             = .FALSE. DETERMINATION OF POLYNOMIAL COEFFICIENTS
!                NEAR X0
!     R       VECTOR OF SIZE K+1 CONTAINS:
!             -  VALUES P(X0), P'(X0), P''(X0),..PK(X0),
!                IF B IS TRUE AND IF K < N
!             -  THE N+1 COEFFICIENTS OF POLYNOMIAL
!                Q(X-X0) = R(1)+R(2)*(X-X0)+...+R(N+1)*(X-X0)**(N),
!                IF B IS FALSE AND IF K = N.
!
!     REFERENCE:
!     ALGORITHM 337, COLLECTED ALGORITHMS FROM CACM, W.PANKIEWICKZ
!--------------------------------------------------------------------}
    Var RR: Double;
        I,J,L,NMJ: Integer;
    Begin
      RR := A^[1];
      For I := 1 to K+1 do R^[I] := RR;
      For J := 2 to N+1 do
      begin
        R^[1] := R^[1]*X0+A^[J];
        NMJ := N-J+1;
        IF NMJ > K THEN
          L := K
        ELSE
          L := NMJ;
        For I := 2 to L+1 do R^[I] := R^[I]*X0+R^[I-1]
      end;
      IF (B) THEN
      begin
        L := 1;
        For I := 2 to K+1 do
        begin
          L := L*(I-1);
          R^[I] := R^[I]*L
        end
      end
    End;


{main program}
BEGIN

  N := 4;      {ORDER OF POLYNOMIAL}
  NP := 11;    {NUMBER OF POINTS}

  New(A); New(RES);

  {define coefficients}
  A^[1] := 1.0;
  A^[2] := 2.0;
  A^[3] := 3.0;
  A^[4] := 4.0;
  A^[5] := 5.0;

  X1 := 0.0;   {X BEGIN}
  X2 := 2.0;   {X END}
  K:=3;        {NUMBER OF DERIVATIVES}
  FLG:=TRUE;   {EVALUATE ONLY}

  DX := (X2-X1) / (1.0*(NP-1));

  X := X1 - DX;

  writeln;
  writeln('     EVALUATE A POLYNOMIAL AND ITS DERIVATIVES BY HORNER''S METHOD');
  writeln;
  writeln('             Example:  P(X) := X^4 + 2X^3 +3X^2 + 4X + 5');
  writeln;
  writeln('       X          P(X)          P''(X)         P"(X)         P"''(X)');
  writeln(' -------------------------------------------------------------------');

  For I:=1 to NP do
  begin
    X := X + DX;
    HORNER(N,A,X,K,FLG,RES);
    write(X:10:5);
    For J:=1 to K+1 do write(RES^[J]:14:5);
    writeln
  end;
  writeln(' -------------------------------------------------------------------');

  Dispose(A); Dispose(RES);
  ReadKey;
  DoneWinCrt

END.

{end of file thorner.pas}