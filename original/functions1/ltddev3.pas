{********************************************************
*  Calculate a limited development of a real function   *
*  f(x)og(x) at x=0 up to order 25, knowing the limited *
*  developments of f(x) and g(x).                       *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
*                                                       *
* Limited development of a real function f(x)og(x):     *
*                                                       *
* Function to develop: exp(sin(x))                      *
*                                                       *
* Limited development of f(x)=exp(x) is:                *
* 1 + x + x^2/2! + x^3/3! + ... + x^n/n! + ...          *
*                                                       *
* Limited development of g(x)=sin(x) is:                *
* x - x^3/3! +x^5/5! - ... + (-1)^n x^2n+1/(2n+1)! + ...*
*                                                       *
* Order of development (max=25): 10                     *
*                                                       *
* Input coefficients of limited dev. of f(x):           *
*  a0 = 1                                               * 
*  a1 = 1                                               *
*  a2 = 0.5                                             *
*  a3 = 0.16666667                                      *
*  a4 = 0.04166667                                      *
*  a5 = 0.00833333                                      *
*  a6 = 0.00138889                                      *
*  a7 = 0.00019841                                      *
*  a8 = 0.00002480                                      *
*  a9 = 0.00000276                                      *
*  a10 = 0.000000276                                    *
*                                                       *
* Input coefficients of limited dev. of g(x):           *
*  b0 = 0                                               * 
*  b1 = 1                                               *
*  b2 = 0                                               *
*  b3 = -0.16666667                                     *
*  b4 = 0                                               *
*  b5 = 0.00833333                                      *
*  b6 = 0                                               *
*  b7 = -0.00019841                                     *
*  b8 = 0                                               *
*  b9 = 0.00000276                                      *
*  b10 = 0                                              *                          
*                                                       *
* The coefficients of limited dev. of f(x)og(x) are:    *
*  c0 =  1.0000000                                      *
*  c1 =  1.0000000                                      *
*  c2 =  0.5000000                                      *
*  c3 =  0.0000000                                      *
*  c4 = -0.1250000                                      *
*  c5 = -0.0666667                                      *
*  c6 = -0.0041667                                      *
*  c7 =  0.0111111                                      *
*  c8 =  0.0053819                                      *
*  c9 =  0.0001764                                      *
*  c10 = -0.0008132                                     *
*                                                       *
* ----------------------------------------------------- *
* Ref.: "Mathematiques en Turbo Pascal By Alain         *
*        Reverchon and Marc Ducamp, Editions Eyrolles,  *
*        Paris, 1991" [BIBLI 03].                       *
*                                                       *
*                            TPW version by J-P Moreau. *
*                                (www.jpmoreau.fr)      *
********************************************************}
Program Ltddev3;
Uses WinCrt;

Const
     SIZE = 25;

Type
     Table = Array[0..SIZE] of Double;

Var
     i,m: Integer;
     T1,T2,R: Table;

Procedure ar_dlcomp(n:Integer; var t1:Table; var t2:Table; var res:Table);
Var i,j,m,pos:Integer; x,s,h: Double;
    binome: array[0..SIZE,0..SIZE] of Double;
    id: array[0..SIZE] of Integer;

    Function ProductBin: Double;
    Var k,p,q:Integer; u:Double;
    Begin
      u:=t2[id[1]];
      p:=i; q:=1;
      for k:=2 to i do
      begin
        u:=u*t2[id[k]];
        if id[k]=id[k-1] then
          Inc(q)
        else
        begin
          u := u*binome[p,q];
          p:=p-q;
          q:=1
        end
      end;
      ProductBin:=u
    End;

Begin
  if n > SIZE then exit;
  if t2[0] <> 0 then exit;
  res[0]:=t1[0];
  {calculate coefficients of binome}
  binome[1,0]:=1;
  binome[1,1]:=1; 
  for i:=2 to n do
  begin
    binome[i,0]:=1;
    binome[i,i]:=1;
    for j:=1 to i-1 do
      binome[i,j]:=binome[i-1,j-1]+binome[i-1,j]
  end;
  {calculate coefficients of gof}
  for m:=1 to n do
  begin
    x:=t1[1]*t2[m];
    for i:= 2 to m do
    begin
      {calculate s=gamma(i,m) }
      s:=0;
      for j:=1 to i-1 do id[j]:=1;
      pos:=1;
      Repeat
        for j:=i-pos+1 to i-1 do id[j]:=id[i-pos];
        id[i]:=m;
        for j:=1 to i-1 do id[i]:=id[i]-id[j];
        while id[i] >= id[i-1] do
        begin
          pos:=1;
          s:=s+ProductBin;
          id[i]:=id[i]-1;
          id[i-1]:=id[i-1]+1
        end;
        Inc(pos);
        if pos < i then  id[i-pos]:=id[i-pos]+1
      Until pos >= i;
      x:=x+s*t1[i]
    end;
    res[m]:=x
  end
End;

{main program}
BEGIN
  writeln;
  writeln(' Limited development of a real function f(x)og(x):');
  writeln;
  writeln(' Function to develop: exp(sin(x))');
  writeln;
  writeln(' Limited development of f(x)=exp(x) is:');
  writeln('  1 + x + x^2/2! + x^3/3! + ... + x^n/n! + ...');
  writeln;
  writeln(' Limited development of g(x)=sin(x) is:');
  writeln('  x - x^3/3! +x^5/5! - ... + (-1)^n x^2n+1/(2n+1)! + ...');
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

  ar_dlcomp(m,T1,T2,R);

  writeln(' The coefficients of limited dev. of f(x)og(x) are:');
  for i:=0 to m do writeln('  c',i,'= ',R[i]:10:7);
  writeln;
  Readkey; DoneWinCrt
END.

{end of file ltddev3.pas}