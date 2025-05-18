{***************************************************
*         Program to demonstrate ASYMERF           *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version by J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* Find the value of ERF(X)=2*Exp(-X*X)/SQRT(PI)    *
*                                                  *
* Input X ? 3                                      *
* ERF(X)= .9999779 with error estimate=-0.00000000 *
* Number of terms evaluated was 10                 *
*                                                  *
* Input X ? 4                                      *
* ERF(X)= 1.0000000 with error estimate= 0.0000000 *
* Number of terms evaluated was 17                 *
*                                                  *
***************************************************}
PROGRAM ASYMERF;
Uses WinCrt;

Var
        e, x, y : DOUBLE;
        n       : integer;


{**********************************************************
* Asymptotic series expansion of the integral of          *
* 2 EXP(-X*X)/(X*SQRT(PI)), the normalized error function *
* (ASYMERF). This program determines the values of the    *
* above integrand using an asymptotic series which is     *
* evaluated to the level of maximum accuracy.             *
* The integral is from 0 to X. The input parameter, X     *
* must be > 0. The results are returned in Y and Y1,      *
* with the error measure in E. The number of terms used   *
* is returned in N. The error is roughly equal to first   *
* term neglected in the series summation.                 *
* ------------------------------------------------------- *
* Reference: A short table of integrals by B.O. Peirce,   *
* Ginn and Company, 1957.                                 *
**********************************************************}
PROCEDURE Asym_erf;
Label 100, 200;
Var c1, c2, y1 : DOUBLE;
Begin
  n := 1; y := 1;
  c2 := 1 / (2 * x * x);
  100: y := y - c2;
  n := n + 2; c1 := c2;
  c2 := -c1 * n / (2 * x * x);
  {Test for divergence - The break point is roughly N:=X*X }
  IF ABS(c2) > ABS(c1) THEN GOTO 200;
  {Continue summation}
  GOTO 100;
  200: n := (n + 1) DIV 2;
  e := EXP(-x * x) / (x * SQRT(PI));
  y1 := y * e;
  y := 1.0 - y1;
  e := e * c2
End;

{main program}
BEGIN
  ClrScr;
  Writeln;
  Write(' Input X: '); read(x);

  Asym_erf;

  Writeln;
  Writeln(' ERF(x)= ',y:9:7,'  with error estimate= ',e:9:7);
  Writeln;
  Writeln;
  Writeln(' Number of terms evaluated was ', n:2);
  Writeln;
  ReadKey; DoneWinCrt
END.

{End of file asymerf.pas}