{******************************************************
*       Program to demonstrate two dimensional        *
*               Mueller's method                      *
* --------------------------------------------------- *
*  Reference: BASIC Scientific Subroutines, Vol. II   *
*  By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                     *
*                Pascal version by J-P Moreau, Paris. *
*                         (www.jpmoreau.fr)           *
* --------------------------------------------------- *
* Example: Find a real root of f(x,y)=(x+1)^5*(y-1)^5 *
*                                                     *
* Sample run:                                         *
*                                                     *
* Input the initial guesses and their bonds:          *
*    X0 = 0                                           *
*    Bond on X0 = 3                                   *
*    Y0 = 0                                           *
*    Bond on Y0 = 3                                   *
* Error criterion: 1e-6                               *
* Maximum number of iterations: 100                   *
*                                                     *
* The calculated zero is X= -.99999812  Y= 0.99999812 *
* The associated W value is W = 0.00000000            *
* The number of steps was: 43                         *
*                                                     *
******************************************************}
PROGRAM DEMO_MUELLER2;
Uses WinCrt;

Var     {global variables}
        k,n : INTEGER;
        b1,b2,e,x,y,x0,y0 : DOUBLE;


{***********************************************
  Function subroutine                          }
FUNCTION W(x,y:DOUBLE): DOUBLE;
Var temp : double;
Begin   
  temp := (x+1)*(x+1)*(x+1)*(x+1)*(x+1);
  W := temp * (y-1)*(y-1)*(y-1)*(y-1)*(y-1)
End;
{**********************************************}

{************************************************
*  Two Dimensional Mueller's method subroutine  *
* --------------------------------------------- *
* This routine iteratively seeks the root of a  *
* function W(x,y) by fitting a parabola to      *
* three points and calculating the nearest root *
* as described in Becket and Hurt, numerical    *
* calculations and algorithms.                  *
* INPUTS:                                       *
*   x0, y0 - The initial guess                  *
*   b1, b2 - A bound on the error in this guess *
*   e - The convergence criteria                *
*   n - The maximum number of iterations        *
* OUTPUTS:                                      *
*   x,y - The value of the root found           *
*   k - The number of iterations performed,     *
************************************************}
PROCEDURE Mueller2D;
Label 100,fin;
Var   a1,b,c1,d1,e1,e2,e3,xl,xl1 : DOUBLE;
      x1,x2,x3,y1,y2,y3 : DOUBLE;

  Procedure Utility;
  Begin
    a1:=xl1*xl1*e1-d1*d1*e2+(xl1+d1)*e3;
    c1:=xl1*(xl1*e1-d1*e2+e3);
    b:=a1*a1-4.0*d1*c1*e3;
    {Test for complex root, meaning the parabola is inverted}
    if b<0 then b:=0;
    {Choose closest root}
    if a1<0 then a1:=a1-SQRT(b);
    if a1>0 then a1:=a1+SQRT(b);
    {Guard against a divide by zero}
    if ABS(a1)+ABS(b)=0 then a1:=4.0*d1*e3;
    {Calculate a relative distance of next guess
     and guard against a divide by zero }
    if a1=0 then a1:=1e-7;
    xl:=-2.0*d1*e3/a1
  End;

Begin {Mueller2D}
  {Set up the three evaluation points}
  k:=1;
100: x3:=x0; x1:=x3-b1; x2:=x3+b1;
  {Calculate Mueller parameters and guard against divide by zero}
  if x2-x1=0 then x2:=x2*1.0000001;
  if x2-x1=0 then x2:=x2+0.0000001;
  xl1:=(x3-x2)/(x2-x1);
  d1:=(x3-x1)/(x2-x1);
  {Get values of function}
  e1:=W(x1,y0);
  e2:=W(x2,y0);
  e3:=W(x3,y0);
  Utility;
  {Calculate new x estimate}
  b1:=xl*(x3-x2);
  x:=x3+b1;
  {Test for convergence}
  if ABS(b1)+ABS(b2)<e then goto fin;
  x0:=x;
  {Repeat for the y direction}
  y3:=y0; y1:=y3-b2; y2:=y3+b2;
  {Calculate Mueller parameters and guard against divide by zero}
  if y2-y1=0 then y2:=y2*1.0000001;
  if y2-y1=0 then y2:=y2+0.0000001;
  xl1:=(y3-y2)/(y2-y1);
  d1:=(y3-y1)/(y2-y1);
  {Get values of function}
  e1:=W(x0,y1);
  e2:=W(x0,y2);
  e3:=W(x0,y3);
  Utility;
  {Calculate new y estimate}
  b2:=xl*(y3-y2);
  y:=y3+b2;
  {Test for convergence}
  if ABS(b1)+ABS(b2)<e then goto fin;
  {Test for number of iterations}
  if k>=n then goto fin;
  y0:=y; k:=k+1;
  {Start another pass}
  goto 100;
fin: End; {Mueller2D}



{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Input the initial guess:');
  writeln;
  write('    X0 = '); read(x0);
  write('    Bond on X0 = '); read(b1);
  writeln;
  write('    Y0 = '); read(y0);
  write('    Bond on Y0 = '); read(b2);
  writeln;
  write(' Error criterion: '); read(e);
  writeln;
  write(' Maximum number of iterations: '); read(n);

  Mueller2D;  {Call 2D Mueller's routine}

  writeln;
  writeln;
  writeln(' The calculated zero is (X,Y) = ',x:10:8,'  ',y:10:8);
  writeln;
  writeln(' The associated Y value is Y = ', W(x,y):10:8);
  writeln;
  writeln(' The number of steps was: ', k);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file mueller2.pas}