{***************************************************
*       Program to demonstrate Bessel Series       *
*               Summation Subroutine               *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version by J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
*  BESSEL SERIES SUMMATION                         *
*                                                  *
*  What is the order of the Bessel function ? 2    *
*  Argument ? 1                                    *
*  Convergence criterion ? 1e-6                    *
*                                                  *
*                                                  *
*  J( 1.000) of order  2 = 0.114903                *
*                                                  *
*  Number of terms used:  4                        *
*                                                  *
***************************************************}
PROGRAM BESSLSER;
Uses WinCrt;

VAR  m,n   : integer;
     e,x,y : DOUBLE;

{************************************
* Bessel function series subroutine *
* The order is n, the argument x.   *
* The returned value is in y.       *
* The number of terms used is in m. *
* e is the convergence criterion    *
* Example: e = 1e-6.                *
************************************}
PROCEDURE Bessel_Series;
Label 100, 200, 210;
Var a,b0,b1,b2 : DOUBLE;
    i : integer;
Begin
  a := 1;
  IF n <= 1 THEN GOTO 100;
  {Calculate N! }
  FOR i := 1 TO n DO a := a * i;
100: a := 1 / a;
  IF n = 0 THEN GOTO 200;
  {Calculate multiplying term}
  FOR i := 1 TO n DO a := a * x / 2.0;
200: b0 := 1.0;
  b2 := 1.0; m := 0;
  {Assemble series sum}
210: m := m + 1;
  b1 := -(x * x * b0) / (m * (m + n) * 4);
  b2 := b2 + b1;
  b0 := b1;
  {Test for convergence}
  IF ABS(b1) > e THEN GOTO 210;
  {Form final answer}
  y := a * b2
End;


{main program}
BEGIN
  Writeln;
  Writeln(' BESSEL SERIES SUMMATION');
  Writeln;
  Write(' What is the order of the Bessel function ? '); read(n);
  Write(' Argument ? '); read(x);
  Write(' Convergence criterion ? '); read(e);
  Writeln;
  Writeln;
  Bessel_Series;
  Writeln(' J(',x:7:4,') of order ',n:2,' = ',y:9:6);
  Writeln;
  Writeln(' Number of terms used: ',m);
  Writeln;
  Readkey; DoneWinCrt
END.

{End file besslser.pas}



