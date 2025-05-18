{**************************************************
*          Statistical distributions              *
* ----------------------------------------------- *
*   This program allows computing several distri- *
*   butions:                                      *
*     1.  binomial distribution                   *
*     2.  Poisson distribution                    *
*     3.  normal distribution                     *
*     4.  normal distribution (2 variables)       *
*     5.  chi-square distribution                 *
*     6.  Student T distribution                  *
* ----------------------------------------------- *
* REFERENCE: "Mathematiques et statistiques By H. *
*             Haut, PSI Editions, France, 1981"   *
*             [BIBLI 13].                         *
*                                                 *
*                 Pascal version By J-P Moreau.   *
*                      (www.jpmoreau.fr)          *
* ----------------------------------------------- *
* SAMPLE RUN:                                     *
*                                                 *
*       Statistical distributions                 *
*                                                 *
* Tutorial                                        *
*                                                 *
* 1. Input distribution number:                   *
*                                                 *
*    1: binomial distribution                     *
*    2: Poisson distribution                      *
*    3: normal distribution                       *
*    4: normal distribution (2 variables)         *
*    5: chi-square distribution                   *
*    6: Student T distribution                    *
*                                                 *
* 2. Define the parameters of chosen distribution *
*                                                 *
* 3. Input value of random variable               *
*                                                 *
*                                                 *
* Input distribution number (1 to 6): 3           *
*                                                 *
*   Normal distribution                           *
*                                                 *
*   MU=mean                                       *
*   S =standard deviation                         *
*   X =value of random variable                   *
*                                                 *
*   MU = 2                                        *
*   S  = 3                                        *
*   X  = 5                                        *
*                                                 *
*   Probability of random variable  = X: .0806569 *
*   Probability of random variable <= X: .8413447 *
*   Probability of random variable >= X: .1586552 *
*                                                 *
**************************************************}
Program stat_distributions;
Uses WinCrt;

Label   100,200,300,400,500,600,700,800,1000;

Var
        B       : Array[0..4] of double;    {used by normal distribution}
        ch,n,nu : integer;
        as,bs,cs,ds : string;

        bt,fx,fxy,p,px,qx,ro,s,sx,sy,x,xm,y,ym : double;

{ Calculates x power n}
Function PowerI(x:double; n:integer): double;
var i,m : integer;
    result :double;
begin
  result := 1.0;
  if n=0 then
  begin
    PowerI:=result;
    exit;
  end;
  m:=  n;
  if n<0 then m:=-n;
  for i:=1 to m do result :=x*result;
  PowerI :=result;
  if n<0 then PowerI:=1.0/result
end;

{y power x}
FUNCTION Power(y,x:DOUBLE):DOUBLE;
BEGIN
  IF x<0 THEN EXIT;
  Power:=Exp(x*Ln(y))
END;

{***************************************************
*            BINOMIAL  DISTRIBUTION                *
* ------------------------------------------------ *
* INPUTS:                                          *
*          p:  probability of success              *
*          n:  number of trials                    *
*          x:  value of random variable            *
* OUTPUTS:                                         *
*          fx: probability of random variable  = X *
*          px: probability of random variable <= X *
*          qx: probability of random variable >= X *
***************************************************}
Procedure Binomial(p:double; n:integer; x:double;
                   var fx:double; var px:double; var qx:double);
Label 100;
Var q,t : double; i,ix,n1:integer;
Begin
  q:=1.0-p; t:=p/q;
  n1:=n+1;
  fx:=Power(q,n);   {fx:=prob(0) }
  px:=fx;
  if x=0 then goto 100;
  ix:=Round(x);
  for i:=1 to ix do
  begin
    fx := (n1-i)*t*fx/i;
    px := px + fx
  end;
100: qx:=1.0+fx-px
End;


{***************************************************
*              POISSON  DISTRIBUTION               *
* ------------------------------------------------ *
* INPUTS:                                          *
*          xm:  mean                               *
*           x:  value of random variable           *
* OUTPUTS:                                         *
*          fx: probability of random variable  = X *
*          px: probability of random variable <= X *
*          qx: probability of random variable >= X *
***************************************************}
Procedure Poisson(xm,x:double; var fx,px,qx:double);
Label 100;
Var i,ix:integer;
Begin
  fx:=EXP(-xm);
  px:=fx;
  if x=0 then goto 100;
  ix:=Round(x);
  for i:=1 to ix do
  begin
    fx:=fx*xm/i;
    px:=px+fx
  end;
100: qx:=1.0+fx-px
End;

{***************************************************
* INPUTS:                                          *
*          xm: mean                                *
*           s: standard deviation                  *
*           x: value of random variable            *
* OUTPUTS:                                         *
*          fx: probability of random variable  = X *
*          px: probability of random variable <= X *
*          qx: probability of random variable >= X *
***************************************************}
Procedure normal(xm,s,x:double; var fx,px,qx:double);
Var po,pp,t,y,zy : double; i:integer;
Begin
  {define coefficients B[i) }
  B[0] :=  1.330274429;
  B[1] := -1.821255978;
  B[2] :=  1.781477937;
  B[3] := -0.356563782;
  B[4] :=  0.319381530;
  pp:=0.2316419;
  y := (x-xm)/s;   {reduced variable}
  zy:=EXP(-y*y/2.0)/2.506628275;
  fx:=zy/s;
  {calculate qx by Horner's method}
  t:=1.0/(1.0+pp*y);
  po:=0.0;
  for i:=0 to 4 do
    po:=po*t + B[i];
  po:=po*t;
  qx:=zy*po;
  px:=1.0-qx
End;

{***************************************************
*        Normal distribution (2 variables)         *
* ------------------------------------------------ *
* INPUTS:                                          *
*          xm: mean of x                           *
*          ym: mean of y                           *
*          sx: standard deviation of x             *
*          sy: standard deviation of y             *
*          ro: correlation coefficient             *
*           x: value of 1st random variable        *
*           y: value of 2nd random variable        *
* OUTPUTS:                                         *
*         fxy: probability of random variables     *
*              = (X,Y)                             *
***************************************************}
Procedure Normal2(xm,ym,sx,sy,ro,x,y:double; var fxy:double);
Var r,s,t,vt:double;
Begin
  s:=(x-xm)/sx;  {1st reduced variable}
  t:=(y-ym)/sy;  {2nd reduced variable}
  r:=2.0*(1.0-ro*ro);
  vt:=s*s+t*t-2.0*ro*s*t;
  vt:=EXP(-vt/r);
  s:=sx*sy*sqrt(2.0*r)*pi;
  fxy:=vt/s
End;

{***************************************************
*             Chi-square distribution              *
* ------------------------------------------------ *
* INPUTS:                                          *
*          nu: number of degrees of freedom        *
*           x: value of random variable            *
* OUTPUTS:                                         *
*          fx: probability of random variable = X  *
*          px: probability of random variable <= X *
*          qx: probability of random variable >= X *
***************************************************}
Procedure Chisqr(nu:integer;x:double;var fx,px,qx:double);
Label 100,200,300,400;
Var aa,bb,fa,fm,gn,xn2:DOUBLE; i:INTEGER;
Begin
  {Calculate gn=Gamma(nu/2) }
  xn2:=nu/2;
  gn:=1.0;
  n:=nu DIV 2;
  if nu=2*n then goto 200;
  {Calculate gn for nu odd}
  if nu=1 then goto 100;
  i:=1;
  while i<=2*n-1 do
  begin
    gn := gn * i / 2.0;
    Inc(i,2)
  end;
100: gn := gn * 1.7724528509;
  goto 300;
200: {Calculate gn for nu even}
  if nu=2 then goto 300;
  for i:=1 to n-1 do gn := gn * i;
300: {Calculate fx}
  fx := Power(x,xn2-1);
  aa := EXP(-x/2.0);
  bb := gn*Power(2.0,xn2);
  fx := fx * aa / bb;
  {Calculate px}
  fa := Power((x/2.0),xn2);
  fa := fa * aa / (xn2*gn);;
  fm := 1.0; px := 1.0; i := 2;
400:  fm := fm * x / (nu+i);
  px := px + fm;
  if fm*fa > 1e-8 then
  begin
    i:=i+2;
    goto 400
  end;
  {Accuracy is obtained}
  px := px * fa;
  qx := 1.0 - px
End;

{***************************************************
*            Student's T distribution              *
* ------------------------------------------------ *
* INPUTS:                                          *
*          nu: number of degrees of freedom        *
*           x: value of random variable            *
* OUTPUTS:                                         *
*          bt: probability of random variable      *
*              between -X and +X                   *
*          px: probability of random variable <= X *
*          qx: probability of random variable >= X *
***************************************************}
Procedure Student(nu:integer;x:double;var bt,px,qx:double);
Label 100,200,300,500,600;
Var c,t,tt,vt : double;
    i,nt : integer;
Begin
  c:=0.6366197724;
  x:=x/sqrt(nu);
  t:=ARCTAN(x);    {t:=theta}
  if nu=1 then
  begin
    bt:=t*c;
    goto 500
  end;
  nt:=2*(nu DIV 2);
  if nt=nu then goto 200;
  {calculate A(X/NU) for nu odd}
  bt:=COS(t);
  tt:=bt*bt;
  vt:=bt;
  if nu=3 then goto 100;
  i:=2;
  while i<=nu-3 do
  begin
    vt:=vt*i*tt/(i+1);
    bt:=bt+vt;
    Inc(i,2)
  end;
100: bt:=bt*SIN(t);
  bt:=(bt+t)*c;
  goto 500;
  {calculate A(X/NU) for nu even}
200: tt:=Power(COS(t),2);
  bt:=1.0; vt:=1.0;
  if nu=2 then goto 300;
  i:=1;
  while i<=nu-3 do
  begin
    vt:=vt*i*tt/(i+1);
    bt:=bt+vt;
    Inc(i,2)
  end;
300: bt:=bt*SIN(t);
500: px := (1.0+bt)/2.0;
  qx:=1.0-px
End;

{main program}
BEGIN
  clrscr;
  writeln('       Statistical distributions');
  writeln;
  writeln(' Tutorial');
  writeln;
  writeln(' 1. Input distribution number:');
  writeln;
  writeln('    1: binomial distribution');
  writeln('    2: Poisson distribution');
  writeln('    3: normal distribution');
  writeln('    4: normal distribution (2 variables)');
  writeln('    5: chi-square distribution');
  writeln('    6: Student T distribution');
  writeln;
  writeln(' 2. Define the parameters of chosen distribution');
  writeln;
  writeln(' 3. Input value of random variable');
  writeln;
  writeln;
  write(' Input distribution number (1 to 6): '); read(ch);

  clrscr;
  if (ch<1) or (ch>6) then
  begin
    writeln(' Error: Invalid choice!');
    Halt
  end;

  {call appropriate input section}
  case ch of
    1: goto 100;
    2: goto 200;
    3: goto 300;
    4: goto 400;
    5: goto 500;
    6: goto 600
  end;

100: writeln;
  writeln('   Binomial distribution:');
  writeln;
  writeln('   P=probability of success');
  writeln('   N=number of trials');
  writeln('   X=value of random variable');
  writeln;
  write('   P = '); read(P);
  write('   N = '); read(N);
  write('   X = '); read(X);
  goto 700;

200: writeln;
  writeln('   Poisson distribution');
  writeln;
  writeln('   MU=mean');
  writeln('   X =value of random variable');
  writeln;
  write('   MU = '); read(xm);
  write('   X  = '); read(x);
  goto 700;

300: writeln;
  writeln('   Normal distribution');
  writeln;
  writeln('   MU=mean');
  writeln('   S =standard deviation');
  writeln('   X =value of random variable');
  writeln;
  write('   MU = '); read(xm);
  write('   S  = '); read(s);
  write('   X  = '); read(x);
  goto 700;

400: writeln;
  writeln('   Normal distribution (2 variables)');
  writeln;
  writeln('   MX=mean of X');
  writeln('   MY=mean of Y');
  writeln('   SX=standard deviation of X');
  writeln('   SY=standard deviation of Y');
  writeln('   RO=correlation coefficient');
  writeln('   X =value of 1st random variable');
  writeln('   Y =value of 2nd random variable');
  writeln;
  write('   MX = '); read(xm);
  write('   SX = '); read(sx);
  write('   MY = '); read(ym);
  write('   SY = '); read(sy);
  write('   RO = '); read(ro);
  write('   X  = '); read(x);
  write('   Y  = '); read(y);
  goto 700;

500: writeln;
  writeln('   Chi-square distribution');
  writeln;
  writeln('   NU=number of degrees of freedom');
  writeln('   X =value of random variable');
  writeln;
  write('   NU = '); read(nu);
  write('   X  = '); read(x);
  goto 700;

600: writeln;
  writeln('   Student''s T distribution');
  writeln;
  writeln('   NU=number of degrees of freedom');
  writeln('   X =value of random variable');
  writeln;
  write('   NU = '); read(nu);
  write('   X  = '); read(x);

700: {call appropriate subroutine}
  case ch of
    1: Binomial(p,n,x,fx,px,qx);
    2: Poisson(xm,x,fx,px,qx);
    3: Normal(xm,s,x,fx,px,qx);
    4: Normal2(xm,ym,sx,sy,ro,x,y,fxy);
    5: Chisqr(nu,x,fx,px,qx);
    6: Student(nu,x,bt,px,qx)
  end;

  {print results}
  writeln;
  as:='   Probability of random variable  = X: ';
  bs:='   Probability of random variable <= X: ';
  cs:='   Probability of random variable >= X: ';
  ds:='   Probability of random variables = X,Y: ';
  if ch=4 then begin write(ds,fxy); goto 1000 end;
  if ch<>6 then goto 800;
  writeln('   Prob(-X<=random variable<=X) = ',bt);
  writeln(bs, px);
  writeln(cs, qx);
  goto 1000;
800: writeln(as, fx);
  writeln(bs, px);
  writeln(cs, qx);

1000: ReadKey; DoneWinCrt

END.
    
{End of file distri.pas}