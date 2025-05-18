{***************************************************
* Program to demonstrate inverse normal subroutine *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*          Pascal Version by J.-P. Moreau, Paris.  *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
*  P(Z>X)     X                                    *
* ----------------                                 *
*  0.50    0.0000                                  *
*  0.48    0.0500                                  *
*  0.46    0.1002                                  *
*  0.44    0.1507                                  *
*  0.42    0.2015                                  *
*  0.40    0.2529                                  *
*  0.38    0.3050                                  *
*  0.36    0.3580                                  *
*  0.34    0.4120                                  *
*  0.32    0.4673                                  *
*  0.30    0.5240                                  *
*  0.28    0.5825                                  *
*  0.26    0.6430                                  *
*  0.24    0.7060                                  *
*  0.22    0.7719                                  *
*  0.20    0.8414                                  *
*  0.18    0.9152                                  *
*  0.16    0.9944                                  *
*  0.14    1.0804                                  *
*  0.12    1.1751                                  *
*  0.10    1.2817                                  *
*  0.08    1.4053                                  *
*  0.06    1.5551                                  *
*  0.04    1.7511                                  *
*  0.02    2.0542                                  *
* -0.00    INF.                                    *
*                                                  *
***************************************************}
PROGRAM INVNORM;
Uses WinCrt;

Var  i   : integer;
     x,y : DOUBLE;


{**********************************************
*   Inverse normal distribution subroutine    *
* ------------------------------------------- *
* This program calculates an approximation to *
* the integral of the normal distribution     *
* function from x to infinity (the tail).     *
* A rational polynomial is used. The input is *
* in y, with the result returned in x. The    *
* accuracy is better then 0.0005 in the range *
* 0 < y < 0.5.                                *
* ------------------------------------------- *
* Reference: Abramowitz and Stegun.           *
**********************************************}
PROCEDURE Inverse_Normal;
Label fin;
Var  c0,c1,c2,d1,d2,d3,z : DOUBLE;
Begin
  {Define coefficients}
  c0 := 2.515517;
  c1 := 0.802853;
  c2 := 0.010328;
  d1 := 1.432788;
  d2 := 0.189269;
  d3 := 0.001308;
  IF y <= 0 THEN
  begin
    x := 1e13;
    goto fin
  end;
  z := SQRT(-LN(y * y));
  x := 1.0 + d1 * z + d2 * z * z + d3 * z * z * z;
  x := (c0 + c1 * z + c2 * z * z) / x;
  x := z - x;
fin:End;

{main program}
BEGIN
  ClrScr;
  Writeln('  P(Z>X)     X   ');
  Writeln(' ----------------');
  i := 1; y := 0.5;
  While y > -0.01 do
  begin
    Inverse_Normal;
    IF x < 0.000001 THEN x := 0;
    IF x < 1e10 then
      Writeln(' ',y:5:2,'    ',x:6:4)
    else
      Writeln(' ',y:5:2,'    INF.');
    IF i = 20 THEN
    begin
      Readkey; i := 1
    end;
    i := i + 1;
    y := y - 0.02;
  end;
  Writeln;
  Readkey; DoneWinCrt
END.

{End of file invnorm.pas}