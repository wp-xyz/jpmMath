{***************************************************************
* Integration of a discrete function F(x) by Simpson's method  *
* ------------------------------------------------------------ *
* SAMPLE RUN                                                   *
*                                                              *
* (Find integral of digitalized function sin(x) from x=0 to    *
*  x=1).                                                       *
*                                                              *
* Integration of a discrete real function F(x)                 *
*             by Simpson's Method                              *
*                                                              *
* Number of points: 25                                         *
* Begin x         : 0                                          *
* End x           : 1                                          *
*                                                              *
* Value of integral: 4.59697701833193E-0001                    *
*                                                              *
*                                                              *
*                               TPW version by J-P Moreau.     *
*                                   (www.jpmoreau.fr)          *
***************************************************************}
PROGRAM discrete_Simpson;
Uses WinCrt;

Const  SIZE = 100;

Type   Tab = Array[0..SIZE] of REAL;

Var
       Xi,Yi: Tab;
       i,npoints: Integer;
       dx,result,x,x1,x2: REAL;

{**********************************************************
* This function returns the integral value of a discrete  *
* real function F(x) defined by n points from x1 to xn.   *
* ------------------------------------------------------- *
* INPUTS:                                                 *
*          n      Number of points (xi, yi)               *
*          X      Table of n xi values in ascending order *
*          Y      Table of n yi values                    *
*                                                         *
* OUTPUT          The function returns the integral of    *
*                 F(x) from x1 to xn.                     *
* ------------------------------------------------------- *
* Ref.: "Mathematiques en Turbo Pascal By Alain Reverchon *
*        and Marc Ducamp, Editions Eyrolles, Paris, 1991  *
*        [BIBLI 03].                                      *
**********************************************************}
Function discreteIntegral(n:Integer; X,Y:Tab): REAL;
Var i,i1,i2: Integer; h,k,j,r: REAL;
Begin
  r:=0;
  i:=0;
  Repeat
    i1:=i+1; i2:=i+2; h:=X[i1]-X[i];
    k:=X[i2]-X[i1]; j:=h+k;
    if h=0 then exit;
    r:=r+((h+h-k)*Y[i]+(j*j*Y[i1]+h*(k+k-h)*Y[i2])/k)*j/6/h;
    i:=i2
  Until i>n-3;
  if Not Odd(n) then
    r:=r+(X[n-1]-X[n-2])*(Y[n-1]+Y[n-2])/2;
  discreteIntegral:=r
End;

{main program}
BEGIN
  writeln('    Integration of a discrete real function F(x)');
  writeln('                by Simpson''s Method');
  writeln;
  write(' Number of points: '); readln(npoints);
  write(' Begin x         : '); readln(x1);
  write(' End x           : '); readln(x2);
  writeln;
  {build up test discrete curve F(x)=sin(x) }
  dx:=(x2-x1)/(npoints-1);
  x:=x1-dx;
  for i:=0 to npoints-1 do
  begin
    x:=x+dx;
    Xi[i]:=x; Yi[i]:=sin(x)
  end;

  {call integration function}
  result := discreteIntegral(npoints,Xi,Yi);

  {print value of integral}
  writeln;
  writeln(' Value of integral: ', result);
  writeln;
  ReadKey;
  DoneWinCrt
    
END.

{end of file disinteg.pas}