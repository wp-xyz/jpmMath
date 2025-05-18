{****************************************************
*   Program to demonstrate the Aitken Steffenson    *
*              iteration subroutine                 *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*             Pascal Version By J-P Moreau, Paris.  *
*                       (www.jpmoreau.fr)           *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
* (Example: find a real root of f(x) = x-2*SIN(x))  *
*                                                   *
* Input the initial guess:                          *
*                                                   *
*    X0 = -2                                        *
*                                                   *
* Convergence criterion: 1e-6                       *
*                                                   *
* Convergence factor: -1                            *
*                                                   *
* Maximum number of iterations: 25                  *
*                                                   *
*                                                   *
* The calculated zero is X = -1.895494              *
*                                                   *
* The associated Y value is Y = -0.000000           *
*                                                   *
* The number of iterations was: 3                   *
*                                                   *
* ------------------------------------------------- *
* Note: this algorithm fails with function (x+1)^5. *
****************************************************}
PROGRAM Zero_Aitken_Steffenson;
Uses WinCrt;

Var
        c,e,x,x0  : double;
        m, n      : integer;

{**************************************
  Function subroutine                 }
  FUNCTION Y(x:double):double;
  Begin
    Y := x-2.0*SIN(x);
    c := 1.0-2.0*COS(x);
    c := -1.0/c
  End;
{*************************************}

{**********************************************
*   Aitken Steffenson iteration subroutine    *
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
PROCEDURE Steffenson;
Label 50,100,200,fin;
Var x1,x2,xk,yy:double;
    m1:integer;
Begin
  n:=0;
50: m1:=0;
  x:=x0;
  {Get yy}
100: yy := Y(x);
  yy:=x+c*yy;
  {Enough points for acceleration ? }
  if m1>0 then goto 200;
  x1:=yy; x:=x1; n:=n+1; m1:=m1+1;
  goto 100;
200: x2:=yy;
  {Perform acceleration
   Guard against a zero denominator}
  xk:=x2-2.0*x1+x0;
  if xk=0 then xk:=xk+0.001;
  xk:=(x1-x0)*(x1-x0)/xk;
  x0:=x0-xk;
  {Test for convergence}
  if n>=m then goto fin;
  if ABS(x-x0)<e then goto fin;
  x0:=x1; x1:=x2; x:=x1;
  {Repeat process}
  goto 50;
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

  Steffenson;   {Call Aitken_Steffenson routine}

  writeln;
  writeln;
  writeln(' The calculated zero is X = ',x:9:6);
  writeln;
  writeln(' The associated Y value is Y = ',Y(x):9:6);
  writeln;
  writeln(' The number of steps was: ',n);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file steffen.pas}