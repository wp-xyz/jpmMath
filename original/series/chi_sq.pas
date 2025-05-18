{***************************************************
* Program to demonstrate the chi-square procedure  *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*          Pascal Version by J.-P. Moreau, Paris.  *
*                    (www.jpmoreau.fr)             *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
* How many degrees of freedom: 100                 *
* What is the range (X1,X2):                       *
*   X1: 50                                         *
*   X2: 150                                        *
* What is the table step size: 5                   *
* Summation truncation error bound: 1e-6           *
*                                                  *
*   X      Chi-Square CDF                          *
* ------------------------                         *
*   50       0.00001                               *
*   55       0.00007                               *
*   60       0.00052                               *
*   65       0.00261                               *
*   70       0.00985                               *
*   75       0.02918                               *
*   80       0.07034                               *
*   85       0.14206                               *
*   90       0.24680                               *
*   95       0.37742                               *
*  100       0.51881                               *
*  105       0.65350                               *
*  110       0.76780                               *
*  115       0.85507                               *
*  120       0.91559                               *
*  125       0.95401                               *
*  130       0.97649                               *
*  135       0.98869                               *
*  140       0.99486                               *
*  145       0.99779                               *
*  150       0.99910                               *
*                                                  *
***************************************************}
PROGRAM CHISQR;
Uses WinCrt;

VAR
        ee,m1,x,x1,x2,xx2,x3,y : DOUBLE;
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
* Reference: Texas Instruments SR-51 owners Manual,  *
* 1974.                                              *
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

{****************************************************************
*      Chi-square cummulative distribution subroutine           *
* ------------------------------------------------------------- *
* The program is fairly accurate and calls upon the chi-square  *
* probability density function subroutine. The input parameter  *
* is m, the number of degrees of freedom. Also required is the  *
* ordinate value, x. The subroutine returns y, the cummulative  *
* distribution integral from 0 to x. This program also requires *
* an accuracy parameter, e, to determine the level of summation.*
* Reference: Hewlett-Packard statistics programs, 1974.         *
****************************************************************}
Procedure Chi_Square_Cumul;
Label 100, 200;
Var m2:integer; y1:DOUBLE;
begin
  y1 := 1; x2 := x; m2 := m + 2;
  x2 := x2 / m2;
  100: y1 := y1 + x2;
  IF x2 < ee THEN GOTO 200;
  m2 := m2 + 2;
  {This form is used to avoid overflow}
  x2 := x2 * (x / m2);
  {Loop to continue sum}
  GOTO 100;
  {Obtain y, the probability density function}
  200: Chi_Square;
  y := (y1 * y * 2) * (x / m)
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
  Write(' Summation truncation error bound: '); read(ee);
  CLRSCR;
  writeln(' x1=',x1:3:0,' x2=',x2:3:0,' x3=',x3:3:0);
  Writeln;
  Writeln('   X       Chi-Square CDF ');
  Writeln(' -------------------------');
  x:=x1; xx2:=x2;
  while x <= xx2 do
  begin
    Chi_Square_Cumul;
    Writeln('  ',x:3:0,'        ',y:7:5);
    x := x + x3
  end;
  ReadKey; DoneWinCrt
END.

{End of file chi-sqr.pas}