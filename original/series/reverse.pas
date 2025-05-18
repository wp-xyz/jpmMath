***************************************************
*        Program to demonstrate the series         *
*              reversion procedure                 *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*         Pascal Version By J.-P. Moreau, Paris.   *
*                    (www.jpmoreau.fr)             *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* What is the degree of the input polynomial: 3    *
*                                                  *
* Input the coefficients as prompted:              *
*                                                  *
*    A( 0) = ? 1                                   *
*    A( 1) = ? 1                                   *
*    A( 2) = ? 1                                   *
*    A( 3) = ? 1                                   *
*                                                  *
* The reversed polynomial coefficients are:        *
*                                                  *
*    B( 0) =  9                                    *
*    B( 1) =  1                                    *
*    B( 2) =  0                                    *
*    B( 3) = -1                                    *
*    B( 4) =  0                                    *
*    B( 5) =  3                                    *
*    B( 6) =  0                                    *
*    B( 7) = -12                                   *
*                                                  *
***************************************************}
PROGRAM Reverse;
Uses WinCrt;


VAR
        A, B  : Array[0..10] of DOUBLE;
        i, n : INTEGER;
 

{*****************************************************
*           Series reversion subroutine              *
* -------------------------------------------------- *
* This routine takes a polynomial                    *
* Y = A(0) + A(1) * X + ... and returns a polynomial *
* X = B(0) + B(1) * Y + ... A(1) must be <> 0.       *
* The degree of reversion is limited to seven.       *
 --------------------------------------------------- *
* Reference: CRC Standard Mathematical Tables,       *
* 24th edition.                                      *
*****************************************************}
PROCEDURE Reverse_coeff;
Var
        a1,a2,a3,a4,a5,a6,a7,aa,bb : DOUBLE;
Begin
  a1 := A[1];
  IF a1 = 0.0 THEN
  begin
    Writeln(' Divide zero error (A(1) must be nonzero).');
    exit
  end;
  B[1] := 1.0 / a1;
  aa := 1.0 / a1;
  bb := aa * aa; aa := aa * bb;
  B[2] := -a2 / aa;
  a3 := A[3]; aa := aa * bb;
  B[3] := aa * (2.0 * a2 * a2 - a1 * a3);
  a4 := A[4]; aa := aa * bb;
  B[4] := aa * (5.0 * a1 * a2 * a3 - a1 * a1 * a4 - 5.0 * a2 * a2 * a2);
  a5 := A[5]; aa := aa * bb;
  B[5] := 6.0 * a1 * a1 * a2 * a4 + 3.0 * a1 * a1 * a3 * a3 + 14.0 * a2 * a2 * a2 * a2;
  B[5] := B[5] - a1 * a1 * a1 * a5 - 21.0 * a1 * a2 * a2 * a3;
  B[5] := aa * B[5];
  a6 := A[6]; aa := aa * bb;
  B[6] := 7.0 * a1 * a1 * a1 * a2 * a5 + 7.0 * a1 * a1 * a1 * a3 * a4 + 84 * a1 * a2 * a2 * a2 * a3;
  B[6] := B[6] - a1 * a1 * a1 * a1 * a6 - 28.0 * a1 * a1 * a2 * a2 * a4;
  B[6] := B[6] - 28.0 * a1 * a1 * a2 * a3 * a3 - 42.0 * a2 * a2 * a2 * a2 * a2;
  B[6] := aa * B[6];
  a7 := A[7]; aa := aa * bb;
  B[7] := 8.0 * a1 * a1 * a1 * a1 * a2 * a6 + 8.0 * a1 * a1 * a1 * a1 * a3 * a5;
  B[7] := B[7] + 4.0 * a1 * a1 * a1 * a1 * a4 * a4 + 120.0 * a1 * a1 * a2 * a2 * a2 * a4;
  B[7] := B[7] + 180.0 * a1 * a1 * a2 * a2 * a3 * a3 + 132.0 * a2 * a2 * a2 * a2 * a2 * a2;
  B[7] := B[7] - a1 * a1 * a1 * a1 * a1 * a7 - 36.0 * a1 * a1 * a1 * a2 * a2 * a5;
  B[7] := B[7] - 72.0 * a1 * a1 * a1 * a2 * a3 * a4 - 12.0 * a1 * a1 * a1 * a3 * a3 * a3;
  B[7] := B[7] - 330.0 * a1 * a2 * a2 * a2 * a2 * a3;
  B[7] := aa * B[7];
  B[0] := 0.0;
  aa := A[0];
  FOR i := 1 TO 7 DO
  begin
    B[0] := B[0] - B[i] * aa;
    aa := aa * A[0]
  end
End;


{main program}
BEGIN
  ClrScr;
  Writeln;
  Write(' What is the degree of the input polynomial: '); read(n);
  Writeln;
  Writeln(' Input the coefficients as prompted:');
  Writeln;
  FOR i := 0 TO n DO
  begin
    Write('  A(',i,') = '); read(A[i])
  end;
  Writeln;

  Reverse_coeff;

  IF A[1] <> 0 THEN
  begin
    Writeln(' The reversed polynomial coefficients are:');
    Writeln;
    FOR i := 0 TO 7 DO
    begin
      if abs(B[i]) < 1e-20 then B[i] := 0.0;
      Writeln('  B(',i,') = ',B[i]:3:0);
    end
  end;
  Writeln;
  ReadKey; DoneWinCrt
END.

{End of file reverse.pas}