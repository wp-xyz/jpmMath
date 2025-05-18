{***************************************************
* Program to demonstrate the chi-square procedure  *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version by J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
* How many degrees of freedom: 100                 *
* What is the range (X1,X2):                       *
*   X1: 50                                         *
*   X2: 150                                        *
* What is the table step size: 5                   *
*                                                  *
*    X       Chi-Square PDF                        *
*   ------------------------                       *
*    50        0.00000                             *
*    55        0.00003                             *
*    60        0.00018                             *
*    65        0.00076                             *
*    70        0.00237                             *
*    75        0.00571                             *
*    80        0.01107                             *
*    85        0.01772                             *
*    90        0.02393                             *
*    95        0.02779                             *
*   100        0.02816                             *
*   105        0.02525                             *
*   110        0.02025                             *
*   115        0.01468                             *
*   120        0.00970                             *
*   125        0.00588                             *
*   130        0.00330                             *
*   135        0.00172                             *
*   140        0.00084                             *
*   145        0.00038                             *
*   150        0.00017                             *
*                                                  *
***************************************************}
PROGRAM CHISQR;
Uses WinCrt;

VAR
        m1,x,x1,x2,x3,y : DOUBLE;
        m : integer;


{**************************************************
* Series approximation subroutine LN(X!)          *
* Accuracy better then 6 places for x>=3          *
* Accuracy better than 12 places for x>10         *
* Advantage is that very large values of the      *
* argument can be used without fear of over flow. *
* Reference: CRC Math Tables.                     *
* x is the Write, y is the output.                *
* ************************************************}
PROCEDURE LN_FACTX;
Var x1 : DOUBLE;
begin
  x1 := 1 / (x * x);
  y := (x + 0.5) * LN(x) - x * (1 - x1 / 12 + x1 * x1 / 360 - x1 * x1 * x1 / 1260 + x1 * x1 * x1 * x1 / 1680);
  y := y + 0.918938533205
end;

{*****************************************************
* Chi-square function subroutine. This program takes *
* a given degree of freedom, m and value, x, and     *
* calculates the chi-square density distribution     *
* function value, y.                                 *
* -------------------------------------------------- *
* Reference: Texas Instruments SR-51 owners Manual,  *
* 1974.                                              *
* -------------------------------------------------- *
* Subroutine used: LN(X!).                           *
*****************************************************}
PROCEDURE Chi_Square;
var c : DOUBLE;
begin
  {Save X}
  m1 := x;
  {Perform calculation}
  x := m / 2 - 1;
  {Call LN(X!) subroutine}
  LN_FACTX;
  x := m1;
  c := -x / 2 + (m / 2 - 1) * LN(x) - (m / 2) * LN(2) - y;
  y := EXP(c)
end;

{main program}
BEGIN
  CLRSCR;
  Writeln;
  Write(' How many degrees of freedom: '); read(m);
  Writeln(' What is the range (X1,X2):');
  Write('   X1: '); read(x1);
  Write('   X2: '); read(x2);
  Write(' What is the table step size: '); read(x3);
  CLRSCR;
  Writeln;
  Writeln('   X       Chi-Square PDF ');
  Writeln(' -------------------------');
  x:=x1;
  while x <= x2 do
  begin
    Chi_Square;
    Writeln('  ',x:3:0,'        ',y:7:5);
    x := x + x3
  end;
  ReadKey; DoneWinCrt
END.

{End of file chi-sqr.pas}