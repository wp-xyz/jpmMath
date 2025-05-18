{***************************************************
*    Program to demonstrate Cheby_Ser Procedure    *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*          Pascal Version by J.-P. Moreau, Paris.  *
*                    (www.jpmoreau.fr)             *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* Chebyshev polynomial coefficients for degree 2   *
*    A(0) = -1                                     *
*    A(1) = 0                                      *
*    A(2) = 2                                      *
*                                                  *
* Chebyshev polynomial coefficients for degree 3   *
*    A(0) = 0                                      *
*    A(1) = -3                                     *
*    A(2) = 0                                      *
*    A(3) = 4                                      *
*                                                  *
* Chebyshev polynomial coefficients for degree 4   *
*    A(0) = 1                                      *
*    A(1) = 0                                      *
*    A(2) = -8                                     *
*    A(3) = 0                                      *
*    A(4) = 8                                      *
*                                                  *
* Chebyshev polynomial coefficients for degree 5   *
*    A(0) = 0                                      *
*    A(1) = 5                                      *
*    A(2) = 0                                      *
*    A(3) = -20                                    *
*    A(4) = 0                                      *
*    A(5) = 16                                     *
*                                                  *
* Chebyshev polynomial coefficients for degree 6   *
*    A(0) = -1                                     *
*    A(1) = 0                                      *
*    A(2) = 18                                     *
*    A(3) = 0                                      *
*    A(4) = -48                                    *
*    A(5) = 0                                      *
*    A(6) = 32                                     *
*                                                  *
* Chebyshev polynomial coefficients for degree 7   *
*    A(0) = 0                                      *
*    A(1) = -7                                     *
*    A(2) = 0                                      *
*    A(3) = 56                                     *
*    A(4) = 0                                      *
*    A(5) = -112                                   *
*    A(6) = 0                                      *
*    A(7) = 64                                     *
*                                                  *
* Chebyshev polynomial coefficients for degree 8   *
*    A(0) = 1                                      *
*    A(1) = 0                                      *
*    A(2) = -32                                    *
*    A(3) = 0                                      *
*    A(4) = 160                                    *
*    A(5) = 0                                      *
*    A(6) = -256                                   *
*    A(7) = 0                                      *
*    A(8) = 128                                    *
*                                                  *
* Chebyshev polynomial coefficients for degree 9   *
*    A(0) = 0                                      *
*    A(1) = 9                                      *
*    A(2) = 0                                      *
*    A(3) = -120                                   *
*    A(4) = 0                                      *
*    A(5) = 432                                    *
*    A(6) = 0                                      *
*    A(7) = -576                                   *
*    A(8) = 0                                      *
*    A(9) = 256                                    *
*                                                  *
* Chebyshev polynomial coefficients for degree 10  *
*    A(0) = -1                                     *
*    A(1) = 0                                      *
*    A(2) = 50                                     *
*    A(3) = 0                                      *
*    A(4) = -400                                   *
*    A(5) = 0                                      *
*    A(6) = 1120                                   *
*    A(7) = 0                                      *
*    A(8) = -1280                                  *
*    A(9) = 0                                      *
*    A(10) = 512                                   *
*                                                  *
***************************************************}
PROGRAM CHEBYSER;
Uses WinCrt;

CONST   SIZE = 10;

VAR
        B : Array[0..SIZE,0..SIZE] of DOUBLE;

        i, n : INTEGER;


{*******************************************************
* Chebyshev series coefficients evaluation subroutine  *
* ---------------------------------------------------- *
* The order of the polynomial is n. The coefficients   *
* are returned in the array B[i,j), i is the degree of *
* the polynomial, j is the coefficient order.          *
*******************************************************}
PROCEDURE Cheby_Ser;
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


{main program}
BEGIN
  Clrscr;
  Writeln;
  FOR n := 2 TO SIZE DO
  begin
    Cheby_Ser;
    Writeln(' Chebyshev polynomial coefficients for degree ',n);
    Writeln;
    FOR i := 0 TO n DO Writeln('  A(',i,') = ',Round(B[n, i]));
    Writeln;
    IF n < SIZE THEN Readkey
  end;
  Writeln;
  Readkey; DoneWinCrt
END.

{End of file chebyser.pas }