{****************************************************
*  Program to demonstrate Chebyshev economization   *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*            Pascal Version by J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)             *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* What is the degree of the input polynomial: ? 15  *
*                                                   *
* What is the degree of the desired economized      *
* polynomial: ? 9                                   *
*                                                   *
* What is the range of the input polynomial: ? 1.57 *
*                                                   *
* Input the coefficients:                           *
* C( 0) = ? 0                                       *
* C( 1) = ? 1                                       *
* C( 2) = ? 0                                       *
* C( 3) = ? -0.166666666                            *
* C( 4) = ? 0                                       *
* C( 5) = ? 0.00833333333                           *
* C( 6) = ? 0                                       *
* C( 7) = ? -0.0001984127                           *
* C( 8) = ? 0                                       *
* C( 9) = ? 0.000002755732                          *
* C( 10) = ? 0                                      *
* C( 11) = ? -0.000000025052109                     *
* C( 12) = ? 0                                      *
* C( 13) = ? 0.00000000016059045                    *
* C( 14) = ? 0                                      *
* C( 15) = ? -0.00000000000076471635                *
*                                                   *
* The Chebyshev series coefficients are:            *
*                                                   *
* A( 0) =  0.0000000000                             *
* A( 1) =  1.1334708982                             *
* A( 2) =  0.0000000000                             *
* A( 3) = -0.1378841454                             *
* A( 4) =  0.0000000000                             *
* A( 5) =  0.0044798168                             *
* A( 6) =  0.0000000000                             *
* A( 7) = -0.0000674667                             *
* A( 8) =  0.0000000000                             *
* A( 9) =  0.0000005865                             *
* A(10) =  0.0000000000                             *
* A(11) = -0.0000000033                             *
* A(12) =  0.0000000000                             *
* A(13) =  0.0000000000                             *
* A(14) =  0.0000000000                             *
* A(15) =  0.0000000000                             *
*                                                   *
* The economized polynomial coefficients are:       *
*                                                   *
* C( 0) =  0.0000000000                             *
* C( 1) =  0.9999999767                             *
* C( 2) =  0.0000000000                             *
* C( 3) = -1.6666647620                             *
* C( 4) =  0.0000000000                             *
* C( 5) =  0.0083329009                             *
* C( 6) =  0.0000000000                             *
* C( 7) = -0.0001980098                             *
* C( 8) =  0.0000000000                             *
* C( 9) =  0.0000025906                             *
*                                                   *
****************************************************}
PROGRAM CHEBECON;
Uses WinCrt;

CONST   SIZE = 40;

VAR

        A, C : Array[0..SIZE] of DOUBLE;
        B    : Array[0..SIZE,0..SIZE] of DOUBLE;

        i,m,m1 : INTEGER;
        x0     : DOUBLE;


{*******************************************************
* Chebyshev series coefficients evaluation subroutine  *
* The order of the polynomial is n. The coefficients   *
* are returned in the array B[i,j), i is the degree of *
* the polynomial, j is the coefficient order.          *
*******************************************************}
PROCEDURE Cheby_Ser(n:INTEGER);
Var i, j : integer;
Begin
  {Establish t0 and t1 coefficients}
  B[0, 0] := 1; B[1, 0] := 0; B[1, 1] := 1;
  {Return if order is less than two}
  IF n < 2 THEN exit;
  FOR i := 2 TO n DO
  begin
    FOR j := 1 TO i DO
      {Basic recursion relation}
      B[i, j] := 2 * B[i - 1, j - 1] - B[i - 2, j];
    B[i, 0] := -B[i - 2, 0]
  end
End;


{***********************************************************
* Chebyshev economization subroutine. The program takes    *
* the input polynomial coefficients, C(i), and returns the *
* Chebyshev series coefficients, A(i). The degree of the   *
* series passed to the routine is m. The degree of the     *
* series returned is m1. The maximum range of x is x0 used *
* for scaling. Note that the input series coefficients are *
* nulled during the process, and then set equal to the     *
* economized series coefficients.                          *
***********************************************************}
PROCEDURE Cheby_Econ;
Var bb : DOUBLE;  i,j,l,n : integer;
Begin
  {Start by scaling the input coefficients according to C(i) }
  bb := x0;
  FOR i := 1 TO m DO
  begin
    C[i] := C[i] * bb; bb := bb * x0
  end;
  {Call Chebyshev series coefficients subroutine}
  FOR n := m DOWNTO 0 DO
  begin
    Cheby_Ser(n);
    A[n] := C[n] / B[n, n];
    FOR l := 0 TO n DO
      {Chebyshev series of order l is substracted out of the polynomial}
      C[l] := C[l] - A[n] * B[n, l]
  end;
  {Perform truncation}
  FOR i := 0 TO m1 DO
    FOR j := 0 TO i DO
      C[j] := C[j] + A[i] * B[i, j];
  {Convert back to the interval X0}
  bb := 1.0 / x0;
  FOR i := 1 TO m1 DO
  begin
    C[i] := C[i] * bb;
    bb := bb / x0
  end
End;


{main program}
BEGIN
  ClrScr;
  Writeln;
  Write(' What is the degree of the input polynomial: '); read(m);
  Writeln;
  Write(' What is the degree of the desired economized polynomial: '); read(m1);
  Writeln;
  Write(' What is the range of input polynomial: '); read(x0);
  Writeln;
  Writeln(' Input the coefficients:');
  Writeln;
  FOR i := 0 TO m DO
  begin
    Write(' C(',i:2,') = '); readln(C[i])
  end;
  Writeln;

  Cheby_Econ;

  Writeln(' The Chebyshev series coefficients are:');
  Writeln;
  FOR i := 0 TO m DO
    Writeln(' A(',i:2,') = ',A[i]:13:10); 
    Writeln;
  Readkey;
  Writeln(' The economized polynomial coefficients are:');
  Writeln;
  FOR i := 0 TO m1 DO
    Writeln(' C(',i:2,') = ',C[i]:13:10);
  Writeln;
  Readkey; DoneWinCrt
END.

{End of file chebecon.pas}