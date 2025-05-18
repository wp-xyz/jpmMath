{*******************************************************
* This program calculates the discreet primitive of a  *
* user defined real function F(x) in n points by using *
* a classical Runge-Kutta method of order 4 (better    *
* than a Simpson method of order 2) and graphically    *
* displays the result.                                 *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
* (Calculate 20 points of the primitive of F(x)=x      *
*  from x=-5 to x=5 with a finesse of 2. the value of  *
*  the primitive at x=-5 will be given as 25/2).       *
*                                                      *
* Calculate the discreet primitive of a function F(x)  *
*                                                      *
*   Value of x begin: -5                               *
*   Value of x end  : 5                                *
*   Value of primitive in x begin: 12.5                *
*   Number of points to calculate: 100                 *
*   Number of intermediate points: 5                   *
*                                                      *
* A curve P(x) is displayed from x=-5 to x=5.          * 
*                                                      *
*                  TPW version by J-P Moreau, Paris.   *
*                          (www.jpmoreau.fr)           *
*******************************************************}
Program Discreet_Primitive;
Uses WinCrtMy, Type_def, Graph_2d;

VAR
    fa,xa,xb: real;
    finesse,n: integer;
    Yi: RV;              {pointer to a real vector, see Type_def.pas}

{given user defined function F(x) }
Function F(x:real): real;
begin
  F:=x
end;

{******************************************************
* Calculate the discreet primitive of a real function *
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
*           The n calculated ordinates yi of the      *
*           primitive of the user defined function    *
*           F(x) are stored in Table Y(i).            * 
* --------------------------------------------------- *
* Ref.: "Mathematiques en Turbo Pascal By Alain Rever-*
*        chon and Marc Ducamp, Editions Eyrolles,     *
*        Paris, 1990" [BIBLI 03].                     *
******************************************************}
Procedure Primitive(a,b,fa:real; n,fi:Integer);
Var
    h,h2,h6,ly,lr,x,y,u: real;
    i,j: integer; ni: longint;
Begin
  h:=(b-a)/fi/n; h2:=h/2; h6:=h/6;
  lr:=fa; ly:=F(a);
  New(Yi);
  Yi^[1]:=lr;
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
    Yi^[i+1]:=lr
  end
End;


{main program}
BEGIN

  WinCrtInit('GPRIMITIV');
  writeln;
  writeln(' Calculate the discreet primitive of a function F(x)');
  writeln;
  write('   Value of x begin: '); readln(xa);
  write('   Value of x end  : '); readln(xb);
  write('   Value of primitive in x begin: '); readln(fa);
  write('   Number of points to calculate: '); readln(n);
  write('   Number of intermediate points: '); readln(finesse);
  writeln;

  Primitive(xa,xb,fa,n,finesse);

  ClrScr;
  CourbeXY(CrtDc,n+1,10,Yi,xa,xb);
  Legendes(CrtDc,'Primitive of F(x)','X','Y');
  SortieGraphique;
  DoneWinCrt

END.

{end of file gprimiti.pas}