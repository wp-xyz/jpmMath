{****************************************************
*    Program to demonstrate the complex domain      *
*            Mueller's subroutine                   *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*              Pascal Version By J-P Moreau, Paris. *
*                       (www.jpmoreau.fr)           *
* ------------------------------------------------- *
* Example:   Find a complex root of z² + 1 = 0      *
*                                                   *
* SAMPLE RUN:                                       *
*                                                   *
* Write the initial guesses and their bounds:       *
*                                                   *
*    X0 = 0                                         *
*    Bond on X0 = 3                                 *
*    Y0 = 0                                         *
*    Bond on Y0 = 3                                 *
*                                                   *
* Convergence criterion: 0                          *
*                                                   *
* Maximum number of iterations: 100                 *
*                                                   *
* The calculated zero is (X,Y) =                    *
* -1.43254447933788E-0005 -1.00003709578020E+0000   *
*                                                   *
* The associated Z value is (U,V) =                 *
* -7.41927312872814E-0005  2.86519524138603E-0005   *
*                                                   *
* The number of steps was: 100                      *
****************************************************}
PROGRAM DEMO_ZMUELLER;
Uses WinCrt;

Var  k,n : INTEGER;
     b1,b2,e,u,v,x,y,x0,y0 : DOUBLE;


{*************************************
  Function subroutine                }
PROCEDURE Z(x,y:DOUBLE; VAR u,v:DOUBLE);
Begin 
  u := x*x-y*y+1;
  v := 2.0*x*y
End;
{************************************}


{************************************************
* Mueller's method for complex roots subroutine *
* --------------------------------------------- *
* This routine uses the parabolic fitting tech- *
* nique associated with Mueller's method, but   *
* does it in the complex domain.                *
* --------------------------------------------- *
* writeS:                                       *
*   x0, y0 - The initial guess                  *
*   b1, b2 - A bound on the error in this guess *
*   e - The convergence criteria                *
*   n - The maximum number of iterations        *
* OUTPUTS:                                      *
*   x,y - The value of the complex root found   *
*   k - The number of iterations performed,     *                                  *
************************************************}
PROCEDURE ZMueller;
Label 100,200,300,fin;
Var d,d1,d2,x1,x2,x3,y1,y2,y3 : DOUBLE;
    u1,u2,u3,v1,v2,v3,xl1,xl2 : DOUBLE;
    a,a1,a2,b,c1,c2,e1,e2,e3  : DOUBLE;
Begin
  {Start calculations}
  k:=1;
  x3:=x0; y3:=y0; x1:=x3-b1; y1:=y3-b2; x2:=x3+b1; y2:=y3+b2;
100: d:=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
  {Avoid divide by zero}
  if d=0 then d:=1e-7;
  xl1:=(x3-x2)*(x2-x1)+(y3-y2)*(y2-y1);
  xl1:=xl1/d;
  xl2:=(x2-x1)*(y3-y2)+(x3-x2)*(y2-y1);
  xl2:=xl2/d;
  d1:=(x3-x1)*(x2-x1)+(y3-y1)*(y2-y1);
  d1:=d1/d;
  d2:=(x2-x1)*(y3-y1)+(x3-x1)*(y2-y1);
  d2:=d2/d;
  {Get function values}
  Z(x1,y1,u1,v1);
  Z(x2,y2,u2,v2);
  Z(x3,y3,u3,v3);
  {Calculate Mueller parameters}
  e1:=u1*(xl1*xl1-xl2*xl2)-2.0*v1*xl1*xl2-u2*(d1*d1-d2*d2);
  e1:=e1+2.0*v2*d1*d2+u3*(xl1+d1)-v3*(xl2+d2);
  e2:=2.0*xl1*xl2*u1+v1*(xl1*xl1-xl2*xl2)-2.0*d1*d2*u2-v2*(d1*d1-d2*d2);
  e2:=e2+u3*(xl2+d2)+v3*(xl1+d1);
  c1:=xl1*xl1*u1-xl1*xl2*v1-d1*xl1*u2+xl1*d2*v2+u3*xl1;
  c1:=c1-u1*xl2*xl2-v1*xl1*xl2+u2*xl2*d2+v2*d1*xl2-v3*xl2;
  c2:=xl1*xl2*u1+xl1*xl1*v1-d2*xl1*u2-xl1*d1*v2+v3*xl1;
  c2:=c2+u1*xl1*xl2-v1*xl2*xl2-u2*xl2*d1+v2*d2*xl2+u3*xl2;
  b1:=e1*e1-e2*e2-4.0*(u3*d1*c1-u3*d2*c2-v3*d2*c1-v3*d1*c2);
  b2:=2.0*e1*e2-4.0*(u3*d2*c1+u3*d1*c2+v3*d1*c1-v3*d2*c2);
  {Guard against divide by zero}
  if b1=0 then b1:=1e-7;
  a:=ARCTAN(b2/b1); a:=a/2.0;
  b:=SQRT(SQRT(b1*b1+b2*b2));
  b1:=b*COS(a); b2:=b*SIN(a);
  a1:=(e1+b1)*(e1+b1)+(e2+b2)*(e2+b2);
  a2:=(e1-b1)*(e1-b1)+(e2-b2)*(e2-b2);
  if a1>a2 then goto 200;
  a1:=e1-b1; a2:=e2-b2; goto 300;
200: a1:=e1+b1; a2:=e2+b2;
300: a:=a1*a1+a2*a2;
  xl1:=a1*d1*u3-a1*d2*v3+a2*u3*d2+a2*v3*d1;
  {Guard against divide by zero}
  if a=0 then a:=1e-7;
  xl1:=-2.0*xl1/a;
  xl2:=-d1*u3*a2+d2*v3*a2+a1*u3*d2+a1*v3*d1;
  xl2:=-2.0*xl2/a;
  {Calculate new estimate}
  x:=x3+xl1*(x3-x2)-xl2*(y3-y2);
  y:=y3+xl2*(x3-x2)+xl1*(y3-y2);
  {Test for convergence}
  if ABS(x-x0)+ABS(y-y0)<e then goto fin;
  {Test for number of iterations}
  if k>=n then goto fin;
  {Continue}
  k:=k+1; x0:=x; y0:=y; x1:=x2; y1:=y2; x2:=x3; y2:=y3; x3:=x; y3:=y;
  goto 100;
fin: End; {of ZMueller}


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Write the initial guesses and their bounds:');
  writeln;
  write('    X0 = '); read(x0);
  write('    Bond on X0 = '); read(b1);
  writeln;
  write('    Y0 = '); read(y0);
  write('    Bond on Y0 = '); read(b2);
  writeln;
  write(' Convergence criterion: '); read(e);
  writeln;
  write(' Maximum number of iterations: '); read(n);

  ZMueller;   {Call complex domain Mueller subroutine}

  writeln; writeln;
  writeln(' The calculated zero is (X,Y) ='); writeln(' ',x,' ',y);
  writeln;
  Z(x,y,u,v);
  writeln(' The associated Z value is (U,V) ='); writeln(' ',u,' ',v);
  writeln;
  writeln(' The number of steps was: ',k);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file zmueller.pas}