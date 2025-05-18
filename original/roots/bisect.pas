{***************************************************
*   Program to demonstrate Bisection subroutine    *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*             Pascal Release By J-P Moreau, Paris  *
*                      (www.jpmoreau.fr)           *
* ------------------------------------------------ *
* Example:  Find a real root of f(x)=(x+1)^5       *
*                                                  *
* Sample run:                                      *
*                                                  *
* What is the initial range (X0,X1):               *
*    X0 = -5                                       *
*    X1 = 0                                        *
* Convergence criterion: 1e-6                      *
*                                                  *
* The calculated zero is X = -9.99755859375E-001   *
* The associated Y value is Y = 0.000000000E-000   *
* The number of steps was: 12                      *
*                                                  *
***************************************************}
PROGRAM Bisection;
Uses WinCrt;

VAR
        e,x,x0,x1 : DOUBLE;
        m : INTEGER;


{**********************************************}
Function Y(x:DOUBLE): DOUBLE;
Begin
  Y := 1+5*x+10*x*x+10*x*x*x+5*x*x*x*x+x*x*x*x*x
End;
{**********************************************}


{**********************************************
*        Bisection method subroutine          *
* ------------------------------------------- *
* This routine iteratively seeks the zero of  *
* function Y(x) using the method of interval  *
* halving until the interval is less than e   *
* in width. It is assumed that the function   *
* Y(x) is available from a function routine.  *
* ------------------------------------------- *
* Input values: range (x0,x1), and e.         *
* Output values: root x, Y(x) and number of   *
* steps m.                                    *
**********************************************}
PROCEDURE Bisect;
Label 100, fin;
Var   y0,yy:double;
Begin
  m:=0;
100: y0:=Y(x0);
  x:=(x0+x1)/2; yy:=Y(x);
  m:=m+1;
  if yy*y0=0 then goto fin;
  if yy*y0<0 then x1:=x;
  if yy*y0>0 then x0:=x;
  if ABS(x1-x0)>e then goto 100;
fin: End;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' What is the initial range (X0,X1):');
  writeln;
  write('    X0 = '); read(x0);
  write('    X1 = '); read(x1);
  writeln;
  write(' Convergence criterion: '); read(e);

  Bisect;   {Call bisection routine}

  writeln;
  writeln;
  writeln(' The calculated zero is X = ',x);
  writeln;
  writeln(' The associated Y value is Y = ',Y(x));
  writeln;
  writeln(' The number of steps was: ',m);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file bisect.pas}