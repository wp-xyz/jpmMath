{*******************************************************
*  CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION   *
*  F(X), GIVEN BY N POINTS X(I), Y(I)                  *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
*                                                      *
* (Interpolate F(X), defined by 11 points:             *
*      0.0           0.0                               *
*      1.0           1.01                              *
*      2.0           3.99                              *
*      3.0           8.85                              *
*      4.0          15.0                               *
*      5.0          25.1                               *
*      6.0          37.0                               *
*      7.0          50.0                               *
*      8.0          63.9                               *
*      9.0          81.2                               *
*     10.0         100.5  for X = 8.7 )                *
*                                                      *
*  CUBIC SPLINE INTERPOLATION:                         *
*                                                      *
*  For X=  8.70000000   Y= 75.69211348                 *
*                                                      *
* ---------------------------------------------------- *
* Ref.: From Numath Library By Tuan Dang Trong in      *
*       Fortran 77.                                    *
*                                                      *
*                   TPW Release By J-P Moreau, Paris.  *
*                          (www.jpmoreau.fr)           *
*******************************************************}
PROGRAM TSEVAL;

Uses WinCrt;

Const NMAX = 25;

Type
      pVEC = ^VEC;
      VEC = Array[1..NMAX] of Double;

Var
      N: Integer;    {Number of given points}
      X,Y: pVEC;     {Coordinates of N points}
      B,C,D:pVEC;    {Coefficients of cubic spline}
      U,V:Double;    {Coordinates of interpolated point}


    FUNCTION SEVAL(N:Integer;U:Double;X,Y,B,C,D:pVEC): Double;
{------------------------------------------------------------------------
!     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
!     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
!     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
!     BY THE SPLINE SUBROUTINE.
!
!     INPUTS:
!     N       NUMBER OF POINTS OF CURVE Y = F(X)
!     U       ABSCISSA OF POINT TO BE INTERPOLATED
!     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
!             OF CURVE F(X)
!     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
!             CUBIC SPLINE
!
!     OUTPUTS:
!     SEVAL   INTERPOLATED VALUE
!             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
!
!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!-----------------------------------------------------------------------}
    Label 10, 20, 30;
    Var
        DX: Double; I,J,K: Integer;
    Begin

      I:=1;

{     BINARY SEARCH }

      IF I >= N Then I := 1;
      IF U < X^[I] then GOTO 10;
      IF U <= X^[I+1] then GOTO 30;
10:   I := 1;
      J := N+1;
20:   K := (I+J) Div 2;
      IF U < X^[K] then J := K;
      IF U >= X^[K] then I := K;
      IF J > I+1 then GOTO 20;

{     SPLINE EVALUATION }

30:   DX := U-X^[I];
      SEVAL := Y^[I]+DX*(B^[I]+DX*(C^[I]+DX*D^[I]))
    End;


    PROCEDURE SPLINE (N:Integer; X,Y:pVEC; Var B,C,D:pVEC);
{---------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
!     SPLINE TO BEST APPROXIMATE A DISCRETE FUNCTION GIVEN BY N POINTS
!
!     INPUTS:
!     N       NUMBER OF GIVEN POINTS
!     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
!             OF FUNCTION F(X)
!
!     OUTPUTS:
!     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
!             OF THE CUBIC SPLINE
!
!     REFERENCE:
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!--------------------------------------------------------------------}
    Label  15,50,60;
    Var
        I,L,NM1: Integer;
        T:Double;
    Begin
      NM1 := N-1;
      IF N < 2 then goto 60;
      IF N < 3 then GOTO 50;

{     BUILD THE TRIDIAGONAL SYSTEM
      B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)  }

      D^[1] := X^[2]-X^[1];
      C^[2] := (Y^[2]-Y^[1])/D^[1];
      For I := 2 to NM1 do
      begin
        D^[I] := X^[I+1]-X^[I];
        B^[I] := 2.0*(D^[I-1]+D^[I]);
        C^[I+1] := (Y^[I+1]-Y^[I])/D^[I];
        C^[I] := C^[I+1]-C^[I]
      end;

{     CONDITIONS AT LIMITS
      THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES }

      B^[1] := -D^[1];
      B^[N] := -D^[N-1];
      C^[1] := 0.0;
      C^[N] := 0.0;

      IF N = 3 then GOTO 15;
      C^[1] := C^[3]/(X^[4]-X^[2])-C^[2]/(X^[3]-X^[1]);
      C^[N] := C^[N-1]/(X^[N]-X^[N-2])-C^[N-2]/(X^[N-1]-X^[N-3]);
      C^[1] := C^[1]*D^[1]*D^[1]/(X^[4]-X^[1]);
      C^[N] := -C^[N]*SQR(D^[N-1])/(X^[N]-X^[N-3]);

{     FORWARD ELIMINATION }

15:   For I := 2 to N do
      begin
        T := D^[I-1]/B^[I-1];
        B^[I] := B^[I]-T*D^[I-1];
        C^[I] := C^[I]-T*C^[I-1]
      end;

{     BACK SUBSTITUTION }

      C^[N] := C^[N]/B^[N];
      For  L := 1 to NM1 do
      begin
        I := N-L;
        C^[I] := (C^[I]-D^[I]*C^[I+1])/B^[I]
      end;

{     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL }

      B^[N] := (Y^[N]-Y^[NM1])/D^[NM1]+D^[NM1]*(C^[NM1]+2.0*C^[N]);
      For I := 1 to NM1 do
      begin
        B^[I] := (Y^[I+1]-Y^[I])/D^[I]-D^[I]*(C^[I+1]+2.0*C^[I]);
        D^[I] := (C^[I+1]-C^[I])/D^[I];
        C^[I] := 3.0*C^[I]
      end;
      C^[N] := 3.0*C^[N];
      D^[N] := D^[NM1];
      goto 60;

{     CAS N = 2 }

50:   B^[1] := (Y^[2]-Y^[1])/(X^[2]-X^[1]);
      C^[1] := 0.0;
      D^[1] := 0.0;
      B^[2] := B^[1];
      C^[2] := 0.0;
      D^[2] := 0.0;
60: End;

{main program}
BEGIN

  N:=11;

  New(X); New(Y); New(B); New(C); New(D);

  X^[1]:=0.0; X^[2]:=1.0; X^[3]:=2.0; X^[4]:=3.0; X^[5]:=4.0; X^[6]:=5.0;
  X^[7]:=6.0; X^[8]:=7.0; X^[9]:=8.0; X^[10]:=9.0; X^[11]:=10.0;

  Y^[1]:=0.01; Y^[2]:=1.01; Y^[3]:=3.99; Y^[4]:=8.85; Y^[5]:=15.0; Y^[6]:=25.1;
  Y^[7]:=37.0; Y^[8]:=50.0; Y^[9]:=63.9; Y^[10]:=81.2; Y^[11]:=100.5;

  SPLINE (N,X,Y,B,C,D);

  U := 9.2;

  V := SEVAL(N,U,X,Y,B,C,D);

  writeln;
  writeln(' CUBIC SPLINE INTERPOLATION:');
  writeln;
  writeln(' For X=',U:12:8,'   Y=', V:12:8);
  writeln;

  ReadKey;
  Dispose(X); Dispose(Y); Dispose(B); Dispose(C); Dispose(D);
  DoneWinCrt

END.

{ end of file tseval.pas }