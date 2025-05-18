{********************************************************** 
* This program calculates the F Distribution Q(F|nu1,nu2) *
* for F, nu1, nu2 given.                                  *
* ------------------------------------------------------- *
* Ref.: "Numerical Recipes, By W.H. Press, B.P. Flannery, *
* S.A. Teukolsky and T. Vetterling, Cambridge University  *
* Press, 1986, page 166" [BIBLI 08].                      *
* ------------------------------------------------------- *
* SAMPLE RUN:                                             *
*                                                         *
* Calculate F Distribution Function Q(F|nu1,nu2) for      * 
* F, nu1 and nu2 given.                                   *
*                                                         *
* Ratio F of dispersion of 1st sample / dispersion of     *
* 2nd sample: 1.25                                        *
* Degree of freedom nu1 of first  sample: 25              *
* Degree of freedom nu2 of second sample: 15              *
*                                                         *
* Number of iterations: 6                                 *
*                                                         *
*                                                         *
* Q = 0.332373                                            *                                                   *
*                                                         *
* ------------------------------------------------------- *
*                    Pascal version by J-P Moreau, Paris. *
*                             (www.jpmoreau.fr)           *
**********************************************************}
PROGRAM FDISTRI;
Uses WinCrt;

CONST
        MAXIT=  100;
        EPS  = 3e-7;
        FPMIN=1e-30;

VAR
        f,nu1,nu2,q,tmp: REAL;


{*********************************************************************
  Returns the value ln(Gamma(xx)) for xx>0. Full accuracy is obtained
  for xx > 1. For 0<xx<1, the reflection formula can be used first:

     Gamma(1-z) := pi/Gamma(z)/sin(pi*z) := pi*z/Gamma(1+z)/sin(pi*z)
*********************************************************************} 
FUNCTION gammln(xx:REAL): REAL;
VAR
cof: array[1..6] of REAL;
stp,half,one,fpf,x,tmp,ser: REAL; {internal arithmetic in double precision}
j: INTEGER;
Begin
  cof[1]:=76.18009173; cof[2]:=-86.50532033; cof[3]:=24.01409822;
  cof[4]:=-1.231739516; cof[5]:=0.120858003e-2; cof[6]:=-0.536382e-5;
  stp:=2.50662827465;
  half:=0.5; one:=1.0; fpf:=5.5;
  x:=xx-one;
  tmp:=x+fpf;
  tmp:=(x+half)*Ln(tmp)-tmp;
  ser:=one;
  for j:=1 to 6 do
  begin
    x:=x+one;
    ser:=ser+cof[j]/x
  end;
  gammln:=tmp+Ln(stp*ser)
End;

FUNCTION betacf(a,b,x:REAL):REAL; Forward;

{**********************************************
  Returns the incomplete beta function Ix(a,b)
**********************************************}
FUNCTION betai(a, b, x: REAL): REAL;
VAR  bt: REAL;
Begin
  if (x < 0) or (x > 1) then
  begin
    writeln(' Bad argument x in function BetaI.');
    betai:=0
  end;

  if (x=0) or (x=1) then
    bt:=0
  else 
    bt:=exp(gammln(a+b)-gammln(a)-gammln(b)+a*Ln(x)+b*Ln(1-x));

  if x < (a+1)/(a+b+2) then            {use continued fraction directly}
    betai:=bt*betacf(a,b,x)/a
  else
    betai:=1-bt*betacf(b,a,1.-x)/b     {use continued fraction after
	 			        making the symmetry transformation}
End;

{*****************************************************************
  Continued fraction for incomplete beta function (used by BETAI)
*****************************************************************}
FUNCTION betacf(a, b, x: REAL): REAL;
LABEL 10;
VAR
  m,m2: INTEGER;
  aa,c,d,del,h,qab,qam,qap: REAL;
Begin
  qab:=a+b; qap:=a+1.0; qam:=a-1.0;
  c:=1.0;
  d:=1.0-qab*x/qap;
  if abs(d) < FPMIN then d:=FPMIN;
  d:=1.0/d; h:=d;
  for m:=1 to MAXIT do
  begin
    m2:=2*m;
    aa:=m*(b-m)*x/((qam+m2)*(a+m2));
    d:=1.0+aa*d;
    if abs(d) < FPMIN then d:=FPMIN;
    c:=1.0+aa/c;
    if abs(c) < FPMIN then c:=FPMIN;
    d:=1.0/d;
    h:=h*d*c;
    aa:=-(a+m)*(qab+m)*x/((qap+m2)*(a+m2));
    d:=1.0+aa*d;
    if abs(d) < FPMIN then d:=FPMIN;
    c:=1.0+aa/c;
    if abs(c) < FPMIN then c:=FPMIN;
    d:=1.0/d;
    del:=d*c;
    h:=h*del;
    if abs(del-1.0) < EPS then goto 10
  end; {of m loop}
10:if m > MAXIT then 
    writeln(' a or b too big, or MAXIT too small.')
  else 
    writeln(' Number of iterations: ', m); 
  betacf:=h
End;


{main program}
BEGIN
{we use the formula: Q(F|nu1,nu2) = BetaI(nu2/2,nu1/2,nu2/(nu2+nu1*F)) }
  writeln;
  writeln(' Calculate F Distribution Function Q(F|nu1,nu2) for');
  writeln(' F, nu1 and nu2 given.');
  writeln;
  write(' Ratio F of dispersion of 1st sample / dispersion of 2nd sample: ');
  readln(f);
  write(' Degree of freedom nu1 of first  sample: '); readln(nu1);
  write(' Degree of freedom nu2 of second sample: '); readln(nu2);
  writeln;

  tmp:=nu2/(nu2+nu1*f);

  q:=betai(nu2/2,nu1/2,tmp);

  writeln;
  writeln(' Q = ', q:8:6);

  Readln; DoneWinCrt
END.

{end of file fdistri.pas}