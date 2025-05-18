{*********************************************
*        Solving algebraic equations         *
*           by Bairstow's Method             *
* ------------------------------------------ *
* Ref.: "Mathematiques par l'informatique    *
*        individuelle - Programmes en basic  *
*        By H. Lehning and D. Jakubowicz,    *
*        MASSON Paris, 1982" [BIBLI 06].     *
*                                            *
* SAMPLE RUN:                                *
* (Find real and complex roots of equation:  *
*  x^6 -127x^5 +215x^4 +28x^3 -39x^2 +20x    *
*  -15 = 0)                                  *
*                                            *
* Input degree of equation: 6                *
*   A(6) = 1                                 *
*   A(5) = -127                              *
*   A(4) = 215                               *
*   A(3) = 28                                *
*   A(2) = -39                               *
*   A(1) = 20                                *
*   A(0) = -15                               *
*                                            *
*             ROOTS                ERROR     *
*   0.03989619     0.44667179  1.874778E-011 *
*   0.03989619    -0.44667179  1.874778E-011 *
*   0.52383379     0.00000000  1.297698E-006 *
*  -0.64575255     0.00000000  3.497111E-006 *
* 125.28210889     0.00000000  3.889075E-009 *
*   1.76001412     0.00000000  1.438088E-006 *
*                                            *
* NOTE: if a better precision is desired,    *
*       you can set e to a smaller value or  *
*       use Newton's method for a real root. *
*                                            *
*             Pascal version by J-P Moreau.  *
*                   (www.jpmoreau.fr)        *
**********************************************
See explanation file Bairstow.txt
---------------------------------    }
Program Bairstow;
Uses WinCrt;

Label 100,200,500,550,1000,fin;

Var
    A,B,C,P: Array[0..20] of double;
    a1,b1,d,e,f,p1,q,r,t,u,v,x,y,z: double;
    i,j,k,m,n: integer;

Procedure Print_results;
Label 10;
Var a1:double;
Begin
  write(x:13:8,'  ',y:13:8);
{calculate error estimation (optional) }
  u:=P[0]; v:=0.0;
  for I:=1 to M do
  begin
    r:=u*x-v*y; v:=u*y+v*x; u:=r+P[I];
  end;
  a1:=SQRT(u*u+v*v);
  u:=0.0; v:=0.0;
  for I:=1 to M do
  begin
    r:=u*x-v*y; v:=u*y+v*x; u:=r+(M-I+1)*P[I-1]
  end;
  if A1=0 then goto 10;
  a1:=a1/SQRT(u*u+v*v);
10: writeln('  ', a1)
End;

{Solve x^2+p1x+q:=0}
Procedure Solve2;
Label 10,fin;
Var d:double;
Begin
  d:=p1*p1-4*q;
  if d<0 then goto 10;
  d:=SQRT(d); y:=0.0;
  x:=(-p1+d)/2.0; Print_results;
  x:=(-p1-d)/2.0; Print_results;
  goto fin;
10: d:=SQRT(-d)/2.0; x:=-p1/2.0;
  y:=d; Print_results;
  y:=-d; Print_results;
fin:End;


{main program}
BEGIN
{Enter polynomial A(i) and copy in P(i) }
  clrscr;
  writeln;
  write(' Input degree of equation: '); readln(N);
  M:=N;  {save N for error estimation}
  for I:=0 to N do
  begin
    write('   A(',N-i,') = '); readln(A[i]);
    P[I]:=A[I]
  end;
  writeln;
  writeln('             ROOTS                 ERROR');
{Init section}
  p1:=0.0; q:=0.0; k:=100; e:=0.001;
{Factoring main loop}
100:if N<=2 then goto 500;
  J:=0;
200: if J>K then goto 1000;
J:=J+1;
{calculate B(i) and C(i) }
  B[0]:=0.0; B[1]:=0.0; C[0]:=0.0; C[1]:=0.0;
  for I:=2 to N+2 do
  begin
    B[I]:=A[I-2]-p1*B[I-1]-q*B[I-2];
    C[I]:=-B[I-1]-p1*C[I-1]-q*C[I-2]
  end;
{calculate dp=a1 and dq=b1}
  x:=B[N+1]; y:=B[N+2]; z:=C[N];
  t:=C[N+1]; u:=C[N+2];
  d:=t*t-z*(u+x);
  if d=0 then goto 1000;
  a1:=(z*y-x*t)/d;
  b1:=(-x*(q*z+p1*t)-y*t)/d;
{New p1 and q}
  p1:=p1+a1; q:=q+b1;
  f:=(ABS(a1)+ABS(b1))/(ABS(p1)+ABS(q));
  if f>e then goto 200;
{A factor has been found}
  Solve2;
{Update polynomial}
  N:=N-2;
  for I:=0 to N do A[I]:=B[I+2];
  goto 100;
500: {Last factor, first or second degree}
  if N=2 then goto 550;
  x:=-A[1]/A[0]; Y:=0.0; {first degree}
  Print_results;
  goto fin;
550: p1:=A[1]/A[0]; q:=A[2]/A[0];
  Solve2;
  goto fin;
1000: writeln(' Process not convergent !');
fin:ReadKey; DoneWinCrt
END.

{end of file bairsto1.pas}