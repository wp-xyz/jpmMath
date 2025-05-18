{********************************************************
*  Calculate a limited development of a real function   *
*  f(x)*g(x) at x=0 up to order 25, knowing the limited *
*  developments of f(x) and g(x).                       *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
*                                                       *
* Limited development of a real function f(x)*g(x):     *
*                                                       *
* Function to develop: exp(x)*(1/1+x)                   *
*                                                       *
* Limited development of f(x)=exp(x) is:                *
*   1 + x + x^2/2! + x^3/3! + ... + x^n/n! + ...        *
*                                                       *
* Limited development of g(x)=1/(1+x) is:               *
*   1 -x + x^2 + ... + (-1)^n*x^n + ...                 *
*                                                       *
* Order of development (max=25): 5                      *
*                                                       *
* Input coefficients of limited dev. of f(x):           *
*  a0 = 1                                               * 
*  a1 = 1                                               *
*  a2 = 0.5                                             *
*  a3 = 0.1666667                                       *
*  a4 = 0.0416667                                       *
*  a5 = 0.0083333                                       *
*                                                       *
* Input coefficients of limited dev. of g(x):           *
*  b0 = 1                                               * 
*  b1 = -1                                              *
*  b2 = 1                                               *
*  b3 = -1                                              *
*  b4 = 1                                               *
*  b5 = -1                                              *
*                                                       *
* The coefficients of limited dev. of f(x)*g(x) are:    *
*  c0 =  1.0000000                                      *
*  c1 =  0.0000000                                      *
*  c2 =  0.5000000                                      *
*  c3 = -0.3333333                                      *
*  c4 =  0.3750000                                      *
*  c5 = -0.3666667                                      *
*                                                       *
* ----------------------------------------------------- *
* Ref.: "Mathematiques en Turbo Pascal By Alain         *
*        Reverchon and Marc Ducamp, Editions Eyrolles,  *
*        Paris, 1991" [BIBLI 03].                       *
*                                                       *
*                            TPW version by J-P Moreau. *
*                                (www.jpmoreau.fr)      *
********************************************************}
Program Ltddev1;
Uses WinCrt;

Const
     SIZE = 25;

Type
     Table = Array[0..SIZE] of Double;

Var
     i,m: Integer;
     T1,T2,R: Table;

Procedure ar_dlmult(n:Integer; var t1:Table; var t2:Table; var res:Table);
Var i,j:Integer; x: Double;
Begin
  if n > SIZE then exit;
  for i:=0 to n do
  begin
    x:=0;
    for j:=0 to i do x:=x+t1[j]*t2[i-j];
    res[i]:=x
  end
End;

{main program}
BEGIN
  writeln;
  writeln(' Limited development of a real function f(x)*g(x):');
  writeln;
  writeln(' Function to develop: exp(x)*(1/1+x)');
  writeln;
  writeln(' Limited development of f(x)=exp(x) is:');
  writeln('  1 + x + x^2/2! + x^3/3! + ... + x^n/n! + ...');
  writeln;
  writeln(' Limited development of g(x)=1/(1+x) is:');
  writeln('  1 -x + x^2 + ... + (-1)^n*x^n + ...');
  writeln;
  write(' Order of development (max=25): '); readln(m);
  writeln;
  writeln(' Input coefficients of limited dev. of f(x):');
  for i:=0 to m do
  begin
    write('  a',i,'= '); readln(T1[i])
  end;
  writeln;
  writeln(' Input coefficients of limited dev. of g(x):');
  for i:=0 to m do
  begin
    write('  b',i,'= '); readln(T2[i])
  end;
  writeln;

  ar_dlmult(m,T1,T2,R);

  writeln(' The coefficients of limited dev. of f(x)*g(x) are:');
  for i:=0 to m do writeln('  c',i,'= ',R[i]:10:7);
  writeln;
  Readkey; DoneWinCrt
END.

{end of file ltddev1.pas}