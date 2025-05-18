{******************************************************
*  Program to demonstrate Newton's method subroutine  *
* --------------------------------------------------- *
*  Reference: BASIC Scientific Subroutines, Vol. II   *
*  By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                     *
* SAMPLE RUN:                                         *
*                                                     *
* (Example: find a real root of (x+1)^5)              *
*                                                     *
* Input the initial guess:                            *
*                                                     *
*    X0 = 0                                           *
*                                                     *
* Convergence factor: 1e-6                            *
*                                                     *
* Maximum number of iterations: 30                    *
*                                                     *
*                                                     *
* The calculated zero is X = -9.98744923787032E-0001  *
*                                                     *
* The associated Y value is Y = 8.32667268468867E-0015*
*                                                     *
* The number of steps was: 30                         *
*                                                     *
*               Pascal version By J-P Moreau, Paris.  *
*                        (www.jpmoreau.fr)            *
******************************************************}
PROGRAM Zero_Newton;
Uses WinCrt;

Var
        e,x0,yy : double;
        m,n : integer;

{**************************************************
  Function subroutine                             }
  FUNCTION Y(x:double;var y1:double):double;
  Begin
    Y := 1+5*x+10*x*x+10*x*x*x+5*x*x*x*x+x*x*x*x*x;
    {derivative}
    y1 := 5+20*x+30*x*x+20*x*x*x+5*x*x*x*x
  end;
{**************************************************}

{**********************************************
*         Newton's method subroutine          *
* ------------------------------------------- *
* This routine calculates the zeros of a      *
* function Y(x) by Newton's method.           *
* The routine requires an initial guess, x0,  *
* and a convergence factor, e. Also required  *
* is a limit on the number of iterations, m.  *
**********************************************}
PROCEDURE Newton;
Label 100,fin;
Var y1:double;
Begin
  n:=0;
  {Get y and y1}
100: yy:=Y(x0,y1);
  {Update estimate}
  x0:=x0-(yy/y1);
  n:=n+1;
  if n>=m then goto fin;
  if ABS(yy/y1)>e then goto 100;
fin: end;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Input the initial guess:');
  writeln;
  write('    X0 = '); read(x0);
  writeln;
  write(' Convergence factor: '); read(e);
  writeln;
  write(' Maximum number of iterations: '); read(m);

  Newton;   {Call Newton routine}

  writeln;
  writeln;
  writeln(' The calculated zero is X = ',x0);
  writeln;
  writeln(' The associated Y value is Y = ',yy);
  writeln;
  writeln(' The number of steps was: ',n);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file newton1.pas}