{***************************************************
*    Program to demonstrate Bessel Coefficients    *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*          Pascal Version by J.-P. Moreau, Paris.  *
*                    (www.jpmoreau.fr)             *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* BESSEL COEFFICIENTS                              *
*                                                  *
* What is the order of the Bessel function ?  0    *
* What degree is desired ? 41                      *
*                                                  *
* The coefficients are:                            ********************
*                                                                     *       
* A( 0) =  0.100000E+01  A( 1) =  0.000000E+00  A( 2) = -0.250000E+00 *
* A( 3) =  0.000000E+00  A( 4) =  0.156250E-01  A( 5) =  0.000000E+00 *
* A( 6) = -0.434028E-03  A( 7) =  0.000000E+00  A( 8) =  0.678168E-05 *
* A( 9) =  0.000000E+00  A(10) = -0.678168E-07  A(11) =  0.000000E+00 *
* A(12) =  0.470950E-09  A(13) =  0.000000E+00  A(14) = -0.240281E-11 *
* A(15) =  0.000000E+00  A(16) =  0.938597E-14  A(17) =  0.000000E+00 *
* A(18) = -0.289690E-16  A(19) =  0.000000E+00  A(20) =  0.724226E-19 *
* A(21) =  0.000000E+00  A(22) = -0.149633E-21  A(23) =  0.000000E+00 *
* A(24) =  0.259780E-24  A(25) =  0.000000E+00  A(26) = -0.384290E-27 *
* A(27) =  0.000000E+00  A(28) =  0.490166E-30  A(29) =  0.000000E+00 *
* A(30) = -0.544629E-33  A(31) =  0.000000E+00  A(32) =  0.531864E-36 *
* A(33) =  0.000000E+00  A(34) = -0.460090E-39  A(35) =  0.000000E+00 *
* A(36) =  0.355008E-42  A(37) =  0.000000E+00  A(38) = -0.245850E-45 *
* A(39) =  0.000000E+00  A(40) =  0.153657E-48  A(41) =  0.000000E+00 *
*                                                                     *
* Argument ?  1                                                       *
* Y=  0.765198                                                        *
**********************************************************************}  
PROGRAM Bessel;
Uses WinCrtMy, Typedef, Math, Utilit;

Var
        A, B : Array[0..50] of ar_reel;
        i, m, n : integer;
        x, y    : ar_reel;


{************************************************************
* Bessel function series coefficient evaluation subroutine  *
* m+1 is the number of coefficients desired, n is the order *
* of the Bessel function. The coefficients are returned in  *
* A(i).                                                     *
************************************************************}
PROCEDURE Bessel_coeff;
Var a1,b1 : ar_reel;
    i : integer;
Begin
  a1 := 1.0; b1 := 1.0;
  FOR i := 1 TO n DO
  begin
    B[i - 1] := 0; b1 := b1 * i; a1 := a1 / 2.0
  end;
  b1 := a1 / b1; a1 := 1.0;
  i := 0;
  while i <= m do
  begin
    A[i] := a1 * b1; A[i + 1] := 0;
    a1 := -a1 / ((i + 2) * (n + n + i + 2));
    Inc(i,2)
  end;
  a1 := a1 / 2.0;
  FOR i := 0 TO m  DO  B[i + n] := A[i];
  FOR i := 0 TO n + m DO  A[i] := B[i]
End;


{main program}
BEGIN
  CLRSCR;
  Writeln(' BESSEL COEFFICIENTS');
  Writeln;
  Write(' What is the order of the Bessel function ? '); read(n);
  Write(' What degree is desired ? '); read(m);
  Writeln;
  
  Bessel_coeff;

  Writeln(' The coefficients are:');
  Writeln;
  FOR i := 0 TO m DO
  begin
    Write(' A(',i:2,') = '); Aff_reel(A[i]); write('  ');
    if i=2 then writeln;
    if i > 3 then
       if ((i+1) MOD 3) = 0  then writeln
  end;
  Writeln;
  Write(' Argument ? '); read(x);
  y := A[0] + A[1] * x;
  IF m > 1 THEN y := y + A[2] * x * x;
  IF m > 2 THEN y := y + A[3] * x * x * x;
  IF m > 3 THEN y := y + A[4] * Puissance(x,4);
  FOR i := 4 TO m - 1 DO
    IF m > i THEN y := y + A[i + 1] * Puissance(x,i+1);
  Write(' Y='); Aff_reel(y); writeln;
  ReadKey; DoneWinCrt
END.


{End file bessel.pas}



