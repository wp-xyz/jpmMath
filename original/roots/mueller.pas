{****************************************************
*     Program to demonstrate Mueller's method       *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*             Pascal version By J-P Moreau, Paris.  *
*                      (www.jpmoreau.fr)            *
* ------------------------------------------------- *
* Example:  Find a real root of f(x)=(x+1)^5        *
*                                                   *
* Sample run:                                       *
*                                                   *
* Input the initial guess:                          *
*    X0 = 0                                         *
* What is the bond of this guess: 3                 *
* Error criterion: 1e-6                             *
* Maximum number of iterations: 100                 *
*                                                   *
* The calculated zero is X = -9.99525660639951E-001 *
* The associated value is Y = -2.2204460492503E-016 *
* The number of steps was: 42                       *
*                                                   *
****************************************************}
PROGRAM DEMO_MUELLER;
Uses WinCrt;

Var     {global variables}
        k,n : INTEGER;
        d,e,x,x0 : DOUBLE;


{***********************************************
  Function subroutine                          }
FUNCTION Y(x:DOUBLE): DOUBLE;
Begin   
  Y := 1+5*x+10*x*x+10*x*x*x+5*x*x*x*x+x*x*x*x*x
End;
{**********************************************}


{**********************************************
*         Mueller's method subroutine         *
* ------------------------------------------- *
* This routine iteratively seeks the root of  *
* a function Y(x) by fitting a parabola to    *
* three points and calculating the nearest    *
* root as described in Becket and Hurt, nume- *
* rical calculations and algorithms.          *
* writeS:                                     *
* The routine requires an initial guess, x0,  *
* a bound on the error on this guess, d and a *
* convergence criteria, e. Also required is a *
* limit on the number of iterations, n.       *
* OUTPUTS:                                    *
* The routine returns the value of the root   *
* found, x and the number of iterations per-  *
* formed, k.                                  *
**********************************************}
PROCEDURE Mueller;
Label 100,200,fin;
Var a1,b,c1,d1,e1,e2,e3,x1,x2,x3,xl,xl1 : DOUBLE;
Begin
  {Set up the three evaluation points}
  k:=1; x3:=x0; x1:=x3-d; x2:=x3+d;
  {Calculate Mueller parameters and guard against divide by zero}
  if x2-x1=0 then x2:=x2*1.0000001;
100: if x2-x1=0 then x2:=x2+0.0000001;
  xl1:=(x3-x2)/(x2-x1);
  d1:=(x3-x1)/(x2-x1);
  if k>1 then goto 200;
  {Get values of function}
  e1:=Y(x1); e2:=Y(x2);
200: e3:=Y(x3);
  a1:=xl1*xl1*e1-d1*d1*e2+(xl1+d1)*e3;
  c1:=xl1*(xl1*e1-d1*e2+e3);
  b:=a1*a1-4.0*d1*c1*e3;
  {Test for complex root, meaning the parabola is inverted}
  if b<0 then b:=0.0;
  {Choose closest root}
  if a1<0 then a1:=a1-SQRT(b);
  if a1>0 then a1:=a1+SQRT(b);
  {Guard against divide by zero}
  if a1=0 then a1:=1e-7;
  xl:=-2.0*d1*e3/a1;
  {Calculate next estimate}
  x:=x3+xl*(x3-x2);
  {Test for convergence}
  if ABS(x-x3)<e then goto fin;
  {Test for number of iterations}
  if k>=n then goto fin;
  {Otherwise, make another pass}
  k:=k+1;
  {Some saving}
  x1:=x2; x2:=x3; x3:=x; e1:=e2; e2:=e3;
  goto 100;
fin: End; {Mueller}


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Input the initial guess:');
  writeln;
  write('    X0 = '); read(x0);
  writeln;
  write(' What is the bond of this guess: '); read(d);
  writeln;
  write(' Error criterion: '); read(e);
  writeln;
  write(' Maximum number of iterations: '); read(n);

  Mueller;  {Call Mueller's routine}

  writeln;
  writeln;
  writeln(' The calculated zero is X = ', x);
  writeln;
  writeln(' The associated Y value is Y = ', Y(x));
  writeln;
  writeln(' The number of steps was: ', k);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file mueller.pas}