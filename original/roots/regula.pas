{***********************************************************
*    Program to demonstrate the modified false position    *
*               method subroutine                          *
* -------------------------------------------------------- *
*     Reference: BASIC Scientific Subroutines, Vol. II     *
*     By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981         *
*     [BIBLI 01].                                          *
*                                                          *
*                     Pascal Version By J-P Moreau, Paris. *
*                             (www.jpmoreau.fr)            *
* -------------------------------------------------------- *
* SAMPLE RUN:                                              *
* (Example: Find a real root of (x+1)^5)                   *
*                                                          *
* Input the initial guesses (root must be between x0, x1): *
*                                                          *
*     X0 = -2                                              *
*     X1 = 2                                               *
*                                                          *
* Convergence factor: .0001                                *
*                                                          *
* Maximum number of iterations: 50                         *
*                                                          *
*                                                          *
* The calculated zero is X = -1.000646                     *
*                                                          *
* The associated Y value is Y =  0.00000000                *
*                                                          *
* The number of steps was: 46                              *
*                                                          *
***********************************************************}
PROGRAM Zero_Regula;
Uses WinCrt;

Var
        e,x,x0,x1 : DOUBLE;
        m, n      : integer;

{**************************************************
  Function subroutine                             }
  FUNCTION Y(x:DOUBLE):DOUBLE;
  Begin
    Y := 1+5*x+10*x*x+10*x*x*x+5*x*x*x*x+x*x*x*x*x
  end;
{*************************************************}


{**********************************************
*     Modified false position subroutine      *
* ------------------------------------------- *
* This subroutine calculates the zeroes of a  *
* function Y(x) uses Hamming's modification   *
* to speed convergence.                       *
* Two initial guess are required, x0 and x1,  *
* bracketting the root, and a convergence     *
* criterion, e. Also required is the maximum  *
* number of iterations, m. The root is retur- *
* ned in x, the actual number of iterations   *
* used in n.                                  *
**********************************************}
PROCEDURE False_position;
Label 100,200,300,fin;
Var yy,y0,y1 : DOUBLE;
Begin
  n:=0;
  {Make x0<x1}
  if x0<x1 then goto 100;
  x:=x0; x1:=x;
  {Get y0 and y1}
100: y0:=Y(x0); y1:=Y(x1);
  {Calculate a new estimate, x}
200: x:=(x0*y1-x1*y0)/(y1-y0);
  {Test for convergence}
  n:=n+1;
  if n>=m then goto fin;
  if ABS(x1-x)<e then goto fin;
  {Get anew Y(x) value}
  yy:=Y(x);
  {Apply Hamming's modification}
  if y1*yy=0 then goto fin;
  if y0*yy>0 then goto 300;
  x1:=x; y1:=yy; y0:=y0/2.0;
300: x0:=x; y0:=yy; y1:=y1/2.0;
  goto 200;
fin: End;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Input the initial guesses:');
  writeln;
  write('    X0 = '); read(x0);
  write('    X1 = '); read(x1);
  writeln;
  write(' Convergence factor: '); read(e);
  writeln;
  write(' Maximum number of iterations: '); read(m);

  False_position;   {Call Secant routine}

  writeln;
  writeln;
  writeln(' The calculated zero is X = ',x:9:6);
  writeln;
  writeln(' The associated Y value is Y = ',Y(x):11:8);
  writeln;
  writeln(' The number of steps was: ',n);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file Regula.pas}

