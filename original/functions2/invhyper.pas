{***************************************************
*  Program to demonstrate the inverse hyperbolic   *
*  functions subroutines                           *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*                   Pascal version by J-P Moreau   *
*                         (www.jpmoreau.fr)        *
***************************************************}
PROGRAM InvHyper;
Uses WinCrt;


Var
    E:double; ligne,N:integer;
    A,W: Array[1..15] of double;
    x,y,Z:double;

Procedure Coeff(VAR N:integer); Forward;
Procedure WCoeff(X:double); Forward;

{*********************************************************
*        Modified cordic logarithm subroutine
* -------------------------------------------------------
* This subroutine takes an input value and returns Y=LN(X)
* X may be any positive real value.
*********************************************************}
Function XLn(X:double): double;
Label 10,20,30,fin;
Var I,K:integer;
    aa,X1,Y:double;
Begin
{ Get coefficients }
  Coeff(N);
{ If X<=0 then an error exists, return }
  if X<=0 then goto fin;
  K:=0;
{ Save X }
  X1:=X;
{ Reduce the range of X }
10: if X<E then goto 20;
{ Divide out a power of E }
    K:=K+1;
    X:=X/E;
  goto 10;
{ Test if X>=1, if so go to next step 
  Otherwise, bring X to >1 }
20: if X>=1.0 then goto 30;
  K:=K-1;
  X:=X*E;
  goto 20;
{ Determine the weighting coefficients, W(I) }
30: WCoeff(X);
{ Calculate residual factor based on Z
  Want LN(Z), where Z is near unity  }
  Z:=Z-1;
  Z:=Z*(1.0-(Z/2.0)*(1.0+(Z/3.0)*(1.0-Z/4.0)));
{ Assemble results }
  aa:=0.5;
  for I:=1 to N do
  begin
    Z:=Z+W[I]*aa;
    aa:=aa/2.0
  end;
{ Z is now the mantissa, K the characteristic }
  Y:=K+Z;
{ Restore X }
  X:=X1;
  XLn:=Y;
fin:End;

{ Weight determination subroutine }
Procedure WCoeff(X:double);
Var I:integer;
Begin
  Z:=X;
  for I:=1 to N do
  begin
    W[I]:=0.0;
    if Z>A[I] then W[I]:=1.0;
    if W[I]=1.0 then Z:=Z/A[I]
  end
End;

{ Exponential coefficients subroutine }
Procedure Coeff(VAR N:Integer);
Begin
  N:=15;
  E:=2.718281828459045;
  A[1]:=1.648721270700128;
  A[2]:=1.284025416687742;
  A[3]:=1.133148453066826;
  A[4]:=1.064494458917859;
  A[5]:=1.031743407499103;
  A[6]:=1.015747708586686;
  A[7]:=1.007843097206448;
  A[8]:=1.003913889338348;
  A[9]:=1.001955033591003;
  A[10]:=1.000977039492417;
  A[11]:=1.000488400478694;
  A[12]:=1.000244170429748;
  A[13]:=1.000122077763384;
  A[14]:=1.000061037018933;
  A[15]:=1.000030518043791
End;

{----------------------------------------------
*    Inverse hyperbolic sine subroutine
* -------------------------------------------
* This subroutine calculates the inverse of
* the hyperbolic sine using the modified cordic
* natural logarithm subroutine 500.
* Input: X - Output: Y := ARCSINH(X)
* Formula: ARCSINH(X]:=LN(X+SQRT(X*X+1))
----------------------------------------------}
Function ArcSinH(X:double):double;
Label 10,fin;
Var X2,Y:double;
Begin
{ Test for zero argument }
  if X<>0.0 then goto 10;
  Y:=0.0;
  goto fin;
10: {Save X}
  X2:=X;
  X:=ABS(X);
  X:=X+SQRT(X*X+1.0);
{ Get LN(X) }
  Y:=XLn(X);
{ Insert sign }
  Y:=(X2/ABS(X2))*Y;
{ Restore X }
  X:=X2;
fin:ArcSinH:=Y
End;

{Inverse hyperbolic cosine subroutine
    ARCCOSH(X]:=LN(X+SQRT(X*X-1)) }
Function ArcCosH(X:double):double;
Label 10,fin;
Var X2,Y:double;
Begin
{ Test for argument less than or equal to unity }
  if X>1.0 then goto 10;
  Y:=0.0;
  goto fin;
10: {Save X}
  X2:=X;
  X:=ABS(X);
  X:=X+SQRT(X-1.0)*SQRT(X+1.0);
  Y:=XLn(X);
{ Restore X }
  X:=X2;
fin:ArcCosH:=Y
End;

{Inverse hyperbolic tangent subroutine
   ARCTANH(X]:=1/2*LN(1+x/1-x)        }
Function ArcTanH(X:double):double;
Label 10,20,fin;
Var X2,Y:double;
Begin
{ Test for X>:= +/- 1 }
  if ABS(X)<=0.999999 then goto 10;
{ Y is BIG! (here +/- 1E18) }
  Y:=(X/ABS(X))*1e18;
  goto fin;
10: {Test for zero argument }
  if X<>0.0 then goto 20;
  Y:=0.0;
  goto fin;
20: {Save X}
  X2:=X;
  X:=(1.0+X)/(1.0-X);
  Y:=XLn(X);
{ Restore X }
  X:=X2;
fin:ArcTanH:=Y
End;


BEGIN
  clrscr;
  ligne:=1;
  writeln;
  writeln('   X      ARCSINH(X)      ARCCOSH(X)      ARCTANH(X) '); 
  writeln(' ----------------------------------------------------');
  x:=-3.0;
  while x<=3.2 do
  begin
    write(X:6:1,'  ');
    y:=arcsinh(x);
    write(Y:14:10,'  ');
    y:=arccosh(x);
    write(Y:14:10,'  ');
    y:=arctanh(x);
    if ABS(y)<1e12 then
      writeln(Y:14:10)
    else
      writeln(' ',Y:6);
    if ligne=20 then
    begin
      ligne:=0;
      Readkey;
      clrscr;
      writeln;writeln
    end;
    Inc(ligne);
    x:=x+0.2
  end;
  writeln;
  Readkey; DoneWinCrt
END.

{End of file invhyper.pas}