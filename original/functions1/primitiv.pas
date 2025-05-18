{*******************************************************
* This program calculates the discrete primitive of a  *
* user defined real function F(x) in n points by using *
* a classical Runge-Kutta method of order 4 (better    *
* than a Simpson method of order 2).                   *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
* (Calculate 20 points of the primitive of F(x)=x      *
*  from x=-5 to x=5 with a finesse of 2. the value of  *
*  the primitive at x=-5 will be given as 25/2).       *
*                                                      *
* Calculate the discrete primitive of a function F(x)  *
*                                                      *
*   Value of x begin: -5                               *
*   Value of x end  : 5                                *
*   Value of primitive in x begin: 12.5                *
*   Number of points to calculate: 20                  *
*   Number of intermediate points: 2                   *
*                                                      *
*  x1 =  -4.500000  y1 =  10.125000                    *
*  x2 =  -4.000000  y2 =   8.000000                    *
*  x3 =  -3.500000  y3 =   6.125000                    *
*  x4 =  -3.000000  y4 =   4.500000                    *
*  x5 =  -2.500000  y5 =   3.125000                    *
*  x6 =  -2.000000  y6 =   2.000000                    *
*  x7 =  -1.500000  y7 =   1.125000                    *
*  x8 =  -1.000000  y8 =   0.500000                    *
* -----------------------------------                  *
*  x16 =   3.000000  y16 =   4.500000                  *
*  x17 =   3.500000  y17 =   6.125000                  *
*  x18 =   4.000000  y18 =   8.000000                  *
*  x19 =   4.500000  y19 =  10.125000                  *
*  x20 =   5.000000  y20 =  12.500000                  *
*                                                      *
* (Exact values are yi = xi*xi/2).                     *
*                                                      *
*                  TPW version by J-P Moreau, Paris.   *
*                          (www.jpmoreau.fr)           *
*******************************************************}
Program discrete_Primitive;
Uses WinCrt;

VAR
    fa,xa,xb: real;
    finesse,n: integer;

{given user defined function F(x) }
Function F(x:real): real;
begin
  F:=x
end;

{******************************************************
* Calculate the discrete primitive of a real function *
* F(x) by using a Runge-Kutta method of order 4.      *
* --------------------------------------------------- *
* INPUTS:                                             *
*           a        begin x                          *
*           b        end x                            *
*           fa       value of primitive at x=a        *
*           n        number of calculated points      *
*           fi       number of intermediate points    *
*                    (finesse)                        *
*                                                     *
* OUTPUTS:                                            *
*           The procedure displays the n calculated   *
*           points (xi,yi) of the primitive of the    *
*           user defined function F(x).               *
* --------------------------------------------------- *
* Ref.: "Mathematiques en Turbo Pascal By Alain Rever-*
*        chon and Marc Ducamp, Editions Eyrolles,     *
*        Paris, 1991" [BIBLI 03].                     *
******************************************************}
Procedure Primitive(a,b,fa:real; n,fi:Integer);
Var
    h,h2,h6,ly,lr,x,y,u: real;
    i,j: integer; ni: longint;
Begin
  h:=(b-a)/fi/n; h2:=h/2; h6:=h/6;
  lr:=fa; ly:=F(a);
  for i:=1 to n do
  begin
    ni:=(i-1)*fi-1;
    for j:=1 to fi do
    begin
      x:=a+h*(ni+j);
      u:=4*F(x+h2);
      x:=x+h;
      y:=F(x);
      lr:=lr+h6*(ly+u+y); ly:=y
    end;
    writeln('  x',i,' = ',x:10:6,'  y',i,' = ',lr:10:6)
  end
End;


BEGIN

  writeln;
  writeln(' Calculate the discrete primitive of a function F(x)');
  writeln;
  write('   Value of x begin: '); readln(xa);
  write('   Value of x end  : '); readln(xb);
  write('   Value of primitive in x begin: '); readln(fa);
  write('   Number of points to calculate: '); readln(n);
  write('   Number of intermediate points: '); readln(finesse);
  writeln;

  Primitive(xa,xb,fa,n,finesse);

  writeln; Readkey;
  DoneWinCrt

END.

{end of file primitiv.pas}