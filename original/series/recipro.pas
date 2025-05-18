{***************************************************
*        Program to demonstrate the series         *
*              inversion procedure                 *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version By J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* What is the degree of the input polynomial: 2    *
*                                                  *
* What is the degree of the inverted polynomial: 6 *
*                                                  *
* Input the polynomial coefficients:               *
*                                                  *
*    A( 0) = ? 1                                   *
*    A( 1) = ? .1                                  *
*    A( 2) = ? .01                                 *
*                                                  *
* The inverted polynomial coefficients are:        *
*                                                  *
*    B( 0) =  1.000000                             *
*    B( 1) = -0.100000                             *
*    B( 2) =  0.000000                             *
*    B( 3) =  0.001000                             *
*    B( 4) = -0.000100                             *
*    B( 5) =  0.000000                             *
*    B( 6) =  0.000001                             *
*                                                  *
***************************************************}
PROGRAM Recipro;
Uses WinCrtMy, Typedef;

VAR
        A, B : Array[0..10] of DOUBLE;
        i, m, n : INTEGER;


{*****************************************************
*        Reciprocal power series procedure           *
* -------------------------------------------------- *
* The input series coefficients are A(i). The output *
* series coefficients are B(i). The degree of the    *
* input polynomial is n. The degree of the inverted  *
* polynomial is m. The routine takes care of norma-  *
* lization using l = A(0).                           *
* -------------------------------------------------- *
* Reference; Computational analysis by Henrici       *
*****************************************************}
PROCEDURE Reciprocal;
Label 100;
Var l   : DOUBLE;
    i,j : integer;
Begin  
  l := A[0];
  FOR i := 0 TO n DO
  begin
    A[i] := A[i] / l; B[i] := 0.0
  end;
  {Clear arrays}
  FOR i := n + 1 TO m DO
  begin
    A[i] := 0; B[i] := 0
  end;
  {Calculate the B(i) coefficients}
  B[0] := 1;
  FOR i := 1 TO m DO
  begin
    j := 1;
100: B[i] := B[i] - A[j] * B[i - j];
    j := j + 1;
    IF j <= i THEN GOTO 100
  end;
  {Un-normalize the A(i) and B(i)}
  FOR i := 0 TO m DO
  begin
    A[i] := A[i] * l;
    B[i] := B[i] / l
  end
End;


{main program}
BEGIN
  ClrScr;
  Writeln;
  Write(' What is the degree of the input polynomial; '); read(n);
  Writeln;
  Write(' What is the degree of the inverted polynomial; '); read(m);
  Writeln;
  Writeln(' Input the polynomial coefficients;');
  Writeln;
  FOR i := 0 TO n DO
  begin
    Write('  A(',i,') = '); read(A[i])
  end;

  Reciprocal;

  Writeln;
  Writeln(' The inverted polynomial coefficients are;');
  Writeln;
  FOR i := 0 TO m DO
    Writeln('  B(',i,') = ',B[i]:9:6);
  Writeln;
  ReadKey; DoneWinCrt
END.

{End of file recipro.pas}