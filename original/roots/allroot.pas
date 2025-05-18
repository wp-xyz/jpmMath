{****************************************************
*    Program to demonstrate the complex domain      *
*             Allroot subroutine                    *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*           Pascal version by J-P Moreau, Paris.    *
*                    (www.jpmoreau.fr)              *
* ------------------------------------------------- *
* Example: Find the two complex roots of z²+1 = 0   *
*                                                   *
* SAMPLE RUN:                                       *
*                                                   *
* write the initial guesses and their bounds:       *
*                                                   *
*    X0 = 4                                         *
*    Bond on X0 = 1                                 *
*    Y0 = 4                                         *
*    Bond on Y0 = 1                                 *
*                                                   *
* Convergence criterion: 1e-9                       *
* Maximum number of iterations: 30                  *
* How many roots are to be sought: 2                *
* Is the function defined in 1000 (1)               *
* or is it a series (2): 1                          *
*                                                   *
* The estimated roots are:                          *
*                                                   *
*      X = 5E-22                                    *
*      Y = 1                                        *
*                                                   *
*      X = -2E-22                                   *
*      Y = -.99999999                               *
*                                                   *
* The last number of iterations was: 6              *
*                                                   *
****************************************************}
PROGRAM DEMO_ALLROOT;
Uses WinCrt;

VAR
        A, X, Y : ARRAY[0..10] of DOUBLE;
        b1,b2,e,u,v,u1,v1,u2,v2,u3,v3,x0,y0,xx,yy,x4,y4,z1,z2 : DOUBLE;
        i,j1,k,m,n,n2,n3 : INTEGER;



{Rectangular to polar conversion}
Procedure RectPol(x,y:DOUBLE; VAR u,v:DOUBLE);
Begin
  u:=SQRT(x*x+y*y);
  {Guard against ambiguous vector}
  if y=0 then y:=1e-16;
  {Guard against divide by zero}
  if x=0 then x:=1e-16;
  v:=ARCTAN(y/x);
  {Check quadrant and adjust}
  if x<0 then v:=v+PI;
  if v<0 then v:=v+2.0*PI
End;

{Polar to rectangular conversion}
Procedure PolRect(u,v:DOUBLE;VAR x,y:DOUBLE);
Begin
  x:=u*COS(v);
  y:=u*SIN(v)
End;

{Polar power}
Procedure PolPower(n:INTEGER;u,v:DOUBLE;VAR u1,v1:DOUBLE);
Var i : INTEGER;
Begin
  u1:=u;
  for i:=2 to n do u1:=u1*u;
  v1:=n*v;
  v1:=v1-(2.0*PI)*INT(v1/2.0/PI)
End;

{Rectangular complex number power}
Procedure ComplexPower(n:INTEGER;x,y:DOUBLE;VAR x1,y1:DOUBLE);
Var u,v,u1,v1:DOUBLE;
Begin
  {Rectangular to polar conversion}
  RectPol(x,y,u,v);
  {Polar power}
  PolPower(n,u,v,u1,v1);
  {Polar to rectangular conversion}
  PolRect(u1,v1,x1,y1)
End;

{*****  Function subroutine  *****}
Procedure Eval(x,y:DOUBLE;VAR u,v:DOUBLE);
Begin
  u:=x*x-y*y+1.0;
  v:=2.0*x*y
End;

{**********************************************
*    Complex series evaluation subroutine     *
* ------------------------------------------- *
* The series coefficients are A(i), assumed   *
* real. The order of the polynomial is m. The *
* routine uses repeated calls to the nth      *
* power (Z^N) complex number subroutine.      *
* writeS:    x,y,m and A(i)                   *
* OUTPUTS:   z1(real) and z2(imaginary)       *
**********************************************}
Procedure Series;
Var a1,a2,x1,y1 : DOUBLE;
Begin
  z1:=A[0]; z2:=0.0;
  {Store x and y}
  a1:=xx; a2:=yy;
  for n:=1 to m do
  begin
    {Recall original x and y}
    xx:=a1; yy:=a2;
    {Call Z^N routine}
    ComplexPower(n,xx,yy,x1,y1);
    {Form partial sum}
    z1:=z1+A[n]*x1; z2:=z2+A[n]*y1
  end;
  {Restore x and y}
  xx:=a1; yy:=a2
End;


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
*   k - The number of iterations performed,     *
************************************************}
Procedure ZMueller;
Label 100,200,300,fin;
Var aa,x1,y1,x2,y2,x3,y3 : DOUBLE;
    a1,a2,b,b1,c1,c2,d,d1,d2,e1,e2,xl1,xl2 : DOUBLE;

  {***** Supervisor subroutine *****}
  Procedure Supervisor;
  Label 100,fin;
  Var j2,n5:INTEGER; a4,u5,v5:DOUBLE;
  Begin
    n5:=n; u5:=u1; v5:=v1;
    {Function or series ? }
    if n3=1 then Eval(xx,yy,u,v);
    if n3=2 then Series;
    if n3=1 then goto 100;
    u:=z1; v:=z2;
    {Restore parameters}
    n:=n5; u1:=u5; v1:=v5;
  100: if j1=0 then goto fin;
    {Divide by the j1 roots already found}
    for j2:=1 to j1 do
    begin
      u5:=u;
      u:=(xx-X[j2])*u+(yy-Y[j2])*v;
      v:=(xx-X[j2])*v-(yy-Y[j2])*u5;
      a4:=(xx-X[j2])*(xx-X[j2])+(yy-Y[j2])*(yy-Y[j2]);
      {Gard against divide by zero}
      if a4=0 then a4:=1e-7;
      v:=v/a4; u:=u/a4;
    end;
  fin: End;  {Supervisor}

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
  xx:=x1; yy:=y1; Supervisor; u1:=u; v1:=v;
  xx:=x2; yy:=y2; Supervisor; u2:=u; v2:=v;
  xx:=x3; yy:=y3; Supervisor; u3:=u; v3:=v;
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
  aa:=ARCTAN(b2/b1); aa:=aa/2.0;
  b:=SQRT(SQRT(b1*b1+b2*b2));
  b1:=b*COS(aa); b2:=b*SIN(aa);
  a1:=(e1+b1)*(e1+b1)+(e2+b2)*(e2+b2);
  a2:=(e1-b1)*(e1-b1)+(e2-b2)*(e2-b2);
  if a1>a2 then goto 200;
  a1:=e1-b1; a2:=e2-b2; goto 300;
200: a1:=e1+b1; a2:=e2+b2;
300: aa:=a1*a1+a2*a2;
  xl1:=a1*d1*u3-a1*d2*v3+a2*u3*d2+a2*v3*d1;
  {Guard against divide by zero}
  if aa=0 then aa:=1e-7;
  xl1:=-2.0*xl1/aa;
  xl2:=-d1*u3*a2+d2*v3*a2+a1*u3*d2+a1*v3*d1;
  xl2:=-2.0*xl2/aa;
  {Calculate new estimate}
  xx:=x3+xl1*(x3-x2)-xl2*(y3-y2);
  yy:=y3+xl2*(x3-x2)+xl1*(y3-y2);
  {Test for convergence}
  if ABS(xx-x0)+ABS(yy-y0)<e then goto fin;
  {Test for number of iterations}
  if k>=n then goto fin;
  {Continue}
  k:=k+1; x0:=xx; y0:=yy; x1:=x2; y1:=y2; x2:=x3; y2:=y3; x3:=xx; y3:=yy;
  goto 100;
fin: End;  {Zmueller}


{***************************************************
*     General root determination subroutine        *
* ------------------------------------------------ *
* This routine attempts to calculate the several   *
* roots of a given series or function by repea-    *
* tedly using the Zmueller subroutine and removing *
* the roots already found by division.             *
* ------------------------------------------------ *
* writeS:     X0, Y0 - The initial guess           *
*	      b1, b2 - The bounds on this guess    *
*	      e - The convergence criteria         *
*	      n - The maximum number of iterations *
*                 for each root                    *
*	      n2 - the number of roots to seek     *
*      	      n3 - A flag (1) for a function F(Z)  *
*                         (2) for a polynomial.    *
* ------------------------------------------------ *
* OUTPUTS:    The n2 roots found X(i), Y(i)        *
* 	      The last number of iterations, k.    *
***************************************************}
Procedure Allroot;
Label 100,200,fin;

  {***** Coefficients subroutine *****}
  Procedure Coefficients;
  Begin
    m:=5; A[0]:=0; A[1]:=24; A[2]:=-50; A[3]:=35; A[4]:=-10; A[5]:=1
  End;

Begin
  k:=0;
  if n3=1 then goto 100;
  if n3<>2 then goto fin;  {error in n3}
100: j1:=0;
  {Save the initial guess}
  x4:=x0; y4:=y0;
  {If n3=2 get the series coefficients}
  if n3=2 then Coefficients;
  {Test for completion}
200: if j1=n2 then goto fin;
  {Call Zmueller subroutine}
  Zmueller;
  j1:=j1+1;
  X[j1]:=xx; Y[j1]:=yy; x0:=x4; y0:=y4;
  {Try another pass}
  goto 200;
fin: End;

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
  writeln;
  write(' How many roots are to be sought: '); read(n2);
  writeln;
  writeln(' Is the function defined in Eval (1)');
  write(' or is it a series (2) '); read(n3);

  AllRoot;   {Call complex domain Allroot subroutine}

  writeln; writeln;
  writeln(' The estimated roots are:');
  for i:=1 to n2 do
  begin
    writeln;
    writeln('   X = ',X[i]);
    writeln('   Y = ',Y[i])
  end;
  writeln;
  writeln(' The last number of iterations was: ',k);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file allroot.pas}