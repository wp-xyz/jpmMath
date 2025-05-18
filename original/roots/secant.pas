{****************************************************
*  Program to demonstrate secant method subroutine  *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*            Pascal Version By J-P Moreau, Paris.   *
*                     (www.jpmoreau.fr)             *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
* ( Example: Find a real root of f(x)=(x+1)^5 )     *
*                                                   *
* Input the initial guesses:                        *
*                                                   *
*    X0 = -2                                        *
*    X1 = 0                                         *
*                                                   *
* Convergence factor: 1e-6                          *
* Maximum number of iterations: 50                  *
*                                                   *
*                                                   *
* The calculated zero is X = -1.000000              *
*                                                   *
* The associated Y value is Y =  0.000000           *
*                                                   *
* The number of steps was: 3                        *
*                                                   *
****************************************************}
PROGRAM Zero_Secant;
Uses WinCrt;

Var
        e,x0,x1 : double;
        m, n    : integer;

{**************************************************
  Function subroutine                             }
  FUNCTION Y(x:double):double;
  Begin
    Y := 1+5*x+10*x*x+10*x*x*x+5*x*x*x*x+x*x*x*x*x
  end;
{*************************************************}

{**********************************************
*          Secant method subroutine           *
* ------------------------------------------- *
* This subroutine calculates the zeroes of a  *
* function Y(x) using the secant method.      *
* Two initial guess are required, x0 and x1,  *
* and a convergence criterion, e. Also requi- *
* red is the maximum number of iterations, m. *
* The root is returned in x, the number of    *
* iterations performed in n.                  *
**********************************************}
PROCEDURE Secant;
Label 100,fin;
Var x,y0,y1:double;
Begin
  n:=0;
  {Start iteration}
100: y0:=Y(x0);
  y1:=Y(x1);
  {Calculate new estimate
   Guard against y1-y0 too small }
  if ABS(y1-y0)<0.001 then y1:=y1+0.001;
  x:=(x0*y1-x1*y0)/(y1-y0);
  n:=n+1;
  {Test for convergence}
  if n>=m then goto fin;
  if ABS(x1-x0)<e then goto fin;
  {Update positions}
  x0:=x1; x1:=x;
  goto 100;

fin: end;


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
  write(' Maximum number of iterations: '); read(m);

  Secant;   {Call Secant routine}

  writeln;
  writeln;
  writeln(' The calculated zero is X = ',x0:9:6);
  writeln;
  writeln(' The associated Y value is Y = ',Y(x0):9:6);
  writeln;
  writeln(' The number of steps was: ',n);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file secant.pas}