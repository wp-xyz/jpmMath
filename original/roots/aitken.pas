{****************************************************
*        Program to demonstrate the Aitken          *
*         acceleration method subroutine            *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*             Pascal version by J-P Moreau, Paris.  *
*                     (www.jpmoreau.fr)             *
* ------------------------------------------------- *
* Example: Find a real root of f(x)=(x+1)^5         *
*                                                   *
* Sample run:                                       *
*                                                   *
*  Input the initial guess:                         *
*     X0 = 0                                        *
*  Convergence criterion: 0.000001                  *
*  Convergence factor: -1                           *
*  Maximum number of iterations: 100                *
*                                                   *
*  The calculated zero is X = -1                    *
*  The associated Y value is Y =  0                 *
*  The number of iterations was:  2                 *
*                                                   *
****************************************************}
PROGRAM Demo_Aitken;
Uses WinCrt;

Var
        c,e,x,x0  : double;
        m, n      : integer;

{**************************************************
  Function subroutine                             }
  FUNCTION Y(x:double):double;
  Begin
    Y := 1+5*x+10*x*x+10*x*x*x+5*x*x*x*x+x*x*x*x*x
  end;
{*************************************************}


{**********************************************
*          Aitken Method Subroutine           *
* ------------------------------------------- *
* This subroutine calculates the zeroes of a  *
* function Y(x) by iterations, and employs    *
* Aitken acceleration to speed up convergence.*
* An initial guess is required, x0, and two   *
* convergence factors, c and e. e relates to  *
* the accuracy of the estimate, and c is used *
* to aid convergence. Also required is the    *
* maximum number of iterations, m. c=-1 is a  *
* normal value, if divergence occurs, smaller *
* and/or positive values should be tried.     *
* The root is returned in x, the number of    *
* iterations in n.                            *
**********************************************}
PROCEDURE Aitken;
Label 100,200,fin;
Var x1,x2,xk,yy : double;
Begin
  n:=0;
  x:=x0;
  {Get y }
100: yy:=Y(x);
  yy:=x+c*yy;
  {Enough points for acceleration ? }
  if n>0 then goto 200;
  x1:=yy; x:=x1; n:=n+1;
  goto 100;
200: x2:=yy; n:=n+1;
  {Guard against a zero denominator}
  if x2-2.0*x1+x0=0 then x0:=x0+0.001;
  {Perform acceleration}
  xk:=(x2-x1)*(x2-x1)/(x2-2.0*x1+x0);
  x2:=x2-xk;
  {Test for convergence}
  if n>=m then goto fin;
  if ABS(xk)<e then goto fin;
  x0:=x1; x1:=x2; x:=x1;
  goto 100;
fin: End;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Input the initial guess:');
  writeln;
  write('    X0 = '); read(x0);
  writeln;
  write(' Convergence criterion: '); read(e);
  writeln;
  write(' Convergence factor: '); read(c);
  writeln;
  write(' Maximum number of iterations: '); read(m);

  Aitken;   {Call Aitken routine}

  writeln;
  writeln;
  writeln(' The calculated zero is X := ',x);
  writeln;
  writeln(' The associated Y value is Y := ',Y(x));
  writeln;
  writeln(' The number of steps was: ',n);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file aitken.pas}