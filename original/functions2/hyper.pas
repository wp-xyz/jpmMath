{****************************************************
* Program to demonstrate the hyperbolic subroutines *
* ------------------------------------------------- *
* Reference: BASIC Scientific Subroutines, Vol. II  *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                   *
*                   Pascal version by J-P Moreau.   *
*                        (www.jpmoreau.fr)          *
****************************************************}
PROGRAM Hyper;
Uses WinCrt;

VAR  A, W: Array[1..9] of DOUBLE;
     N,ligne:Integer;
     aa,E,K,XX,Y,Z:DOUBLE;

Procedure Coeffs; Forward;
Procedure WCoeffs(X:DOUBLE); Forward;

{Modified cordic exponential subroutine}
Function XExp(X:DOUBLE): DOUBLE;
{ This subroutine takes an input value and returns Y:=EXP(X)
  X may be any positive or negative real value
  Get coefficients }
Label fin;
Var Y:DOUBLE;
    I:Integer;
Begin
  Coeffs;
{ Reduce the range of X }
  K:=INT(X);
  X:=X-K;
{ Determine the weighting coeffs, W(I) }
  WCoeffs(X);
{ Calculate products }
  Y:=1.0;
  for I:=1 to N do
    if W[I]>0.0 then Y:=Y*A[I];
{ Perform residual multiplication }
  Y:=Y*(1.0+Z*(1.0+Z/2.0*(1.0+Z/3.0*(1.0+Z/4.0))));
{ Account for factor EXP(K) }
  if K<0 then E:=1.0/E;
  if ABS(K)<1 then goto fin;
  for I:=1 to ABS(Round(K)) do Y:=Y*E;
{ Restore X }
  X:=X+K;
fin:XExp:=Y
End;

Procedure WCoeffs(X:DOUBLE);
Var I:Integer;
Begin
  aa:=0.5;
  Z:=X;
  for I:=1 to N do
  begin
    W[I]:=0.0;
    if Z>aa then W[I]:=1.0;
    Z:=Z-W[I]*aa;
    aa:=aa/2.0
  end
End;

Procedure Coeffs;
Begin
  N:=9;
  E:=2.718281828459045;
  A[1]:=1.648721270700128;
  A[2]:=1.284025416687742;
  A[3]:=1.133148453066826;
  A[4]:=1.064494458917859;
  A[5]:=1.031743407499103;
  A[6]:=1.015747708586686;
  A[7]:=1.007843097206448;
  A[8]:=1.003913889338348;
  A[9]:=1.001955033591003
End;

{---------------------------------------------*
*          Hyperbolic sine Function           *
* ------------------------------------------- *
* This Procedure uses the definition of the   *
* hyperbolic sine and the modified cordic     *
* exponential Function XExp(X) to approximate *
* SINH(X) over the entire range of real X.    *
* ------------------------------------------- }
Function SinH(X:DOUBLE): DOUBLE;
Label 10,fin;
Var Y:DOUBLE;
    I:Integer;
Begin
{ Is X small enough to cause round off erroe ? }
  if ABS(X)<0.35 then goto 10;
{ Calculate SINH(X) using exponential definition
  Get Y:=EXP(X) }
  Y:=XExp(X);
  Y:=(Y-(1.0/Y))/2.0;
  goto fin;
10: {series approximation (for X small) }
  Z:=1.0; Y:=1.0;
  for I:=1 to 8 do
  begin
    Z:=Z*X*X/((2*I)*(2*I+1));
    Y:=Y+Z
  end;
  Y:=X*Y;
fin:SinH:=Y
End;

{ hyperbolic cosine Function }
Function CosH(X:DOUBLE): DOUBLE;
Var Y:DOUBLE;
Begin
  Y:=XExp(X);
  Y:=(Y+(1.0/Y))/2.0;
  CosH:=Y
End;

{ hyperbolic tangent Function
    TANH(X]:=SINH(X)/COSH(X)  }
Function TanH(X:DOUBLE): DOUBLE;
Var V,Y: DOUBLE;
Begin
  V:=SinH(X);
  Y:=CosH(X);
  TanH:=V/Y
End;


{main program}
BEGIN
ligne:=1;
writeln;
writeln('   X        SINH(X)         COSH(X)         TANH(X) ');
writeln(' ---------------------------------------------------');
XX:=-5.0;
Repeat
  write(XX:6:1,'  ');
  Y:=SinH(XX);
  write(Y:14:10,'  ');
  Y:=CosH(XX);
  write(Y:14:10,'  ');
  Y:=TanH(XX);
  writeln(Y:14:10);
  if ligne=20 then
  begin
    ligne:=0;
    Readkey;
    clrscr;
    writeln;writeln
  end;
  Inc(ligne);
  XX:=XX+0.2
Until XX>=5.2;
writeln;
Readkey; DoneWinCrt
END.

{End of file hyper.pas}