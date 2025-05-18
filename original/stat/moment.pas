{*****************************************************
*         Computing the means and moments            *
*           of a statistical variable                *
*                                                    *
* -------------------------------------------------- *
* REFERENCE:  "Mathematiques et statistiques By H.   *
*              Haut, PSI Editions, France, 1981"     *
*              [BIBLI 13].                           *
*                                                    *
*              Pascal version by J-P Moreau, Paris   *
*                       (www.jpmoreau.fr)            *
* -------------------------------------------------- *
* SAMPLE RUN:                                        *
*                                                    *
* TUTORIAL                                           *
*                                                    *
* Means and moments of a statistical variable X[i]   *
* with frequency F[i]                                *
*                                                    *
* 1. Input number n of data                          *
*                                                    *
* 2. Input successively the n values X[i] and F[i]   *
*                                                    *
* Number of data: 3                                  *
*                                                    *
*   1  1 6                                           *
*   2  3 4                                           *
*   3  5 8                                           *
*                                                    *
* Calculate generalized mean for t= 2                *
*                                                    *
*                                                    *
* Arithmetic mean:  3.22222222                       *
* Geometric mean :  2.61023904                       *
* Harmonic mean  :  2.01492537                       *
*                                                    *
* Moments:                                           *
*                                                    *
*   M1=    3.22222222                                *
*   M2=    3.06172840                                *
*   M3=   -1.16323733                                *
*   M4=   12.56881570                                *
*                                                    *
* Flatness coefficient....:   1.34079084             *
* Coefficient of asymmetry:  -0.21712925             *
*                                                    *
* Gen. mean M( 2): 3.66666667                        *
*                                                    *
*****************************************************}
PROGRAM Moments;
Uses WinCrt;

Const   SIZE = 100;

Var
        i,ite,n : INTEGER;
        F,X : Array[1..SIZE] of DOUBLE;
        a,ap,asy,g,h,t,xm1,xm2,xm3,xm4,xmt : DOUBLE;


{y power x}
FUNCTION Power(y,x:DOUBLE):DOUBLE;
BEGIN
  IF x<0 THEN EXIT;
  Power:=Exp(x*Ln(y))
END;

{****************************************************
* Procedure to calculate the means and moments of a *
*            statistical variable X[i]              *
* ------------------------------------------------- *
* INPUTS:                                           *
*          n: number of data X[i] and F[i]          *
*          X: n values X[i]                         *
*          F: n values F[i], frequency of X[i]      *
*          t: coefficient of generalized mean M(t)  *
*             to calculate                          *
* OUTPUTS:                                          *
*          a: arithmetic mean of X[i]               *
*          g: geometric mean of X[i]                *
*          h: harmonic mean of X[i]                 *
*        xm1: moment of first order                 *
*        xm2: moment of second order                *
*        xm3: moment of third order                 *
*        xm4: moment of fourth order                *
*         ap: flatness coefficient                  *
*        asy: coefficient of assymetry              *
*        xmt: generalized mean M(t)                 *
*        ite: test flag (if te=1, harmonic mean h   *
*             is not defined                        *
****************************************************}
Procedure Calculate;
Label 100;
Var
     v1,v2,v3,v4,v5,v6,v7,vt: DOUBLE;
     xm : DOUBLE;
     i: INTEGER;
Begin
  ite:=0;
  v1:=0.0; v2:=0.0; v3:=0.0; v4:=0.0;
  v5:=1.0; v6:=0.0; v7:=0.0; xm:=0.0;
  {Calculate necessary sums}
  for i:=1 to n do
  begin
    vt:=F[i]*X[i];
    v1:=v1+vt;
    v2:=v2+vt*X[i];
    v3:=v3+vt*X[i]*X[i];
    v4:=v4+vt*X[i]*X[i]*X[i];
    xm:=xm+F[i];
    v7:=v7+F[i]*Power(X[i],t);
    {test for a X[i]=0}
    if X[i]=0 then
    begin
      ite:=1;
      goto 100
    end;
    {if one X[i]:=0, the harmonic mean is undefined
     and the geometric mean= 0  }
    if ite=1 then goto 100;
    v5:=v5*Power(X[i],F[i]);
    v6:=v6+F[i]/X[i];
100: end;
  {prepare outputs}
  a:=v1/xm;
  xmt:=Power((v7/xm),(1.0/t));
  xm1:=a;
  xm2:=(v2-2.0*a*v1+a*a*xm)/xm;
  xm3:=(v3-3.0*a*v2+3.0*a*a*v1-xm*Power(a,3))/xm;
  xm4:=(v4-4.0*a*v3+6.0*a*a*v2-4.0*v1*Power(a,3)+xm*Power(a,4))/xm;
  ap:=xm4/xm2/xm2;
  asy:=xm3/Power(xm2,1.5);
  if ite=1 then {one X[i]:=0 }
  begin
    g:=0.0;
    exit
  end;
  g:=Power(v5,(1.0/xm));
  h:=xm/v6
End;


{main program}
BEGIN

  clrscr;
  writeln;
  writeln(' TUTORIAL');
  writeln;
  writeln(' Means and moments of a statistical variable X[i]');
  writeln(' with frequency F[i] ');
  writeln;
  writeln(' 1. Input number n of data');
  writeln;
  writeln(' 2. Input sucessively the n values X[i] and F[i]');
  writeln;
  write(' Number of data: '); readln(n);
  writeln;

  for i:=1 to n do
  begin
    write('  ',i,'  '); readln(X[i], F[i])
  end;
  writeln;
  write(' Calculate generalized mean for t= ');  readln(t);
  writeln;

  Calculate;

  {writeln results}
  clrscr;
  writeln;
  writeln(' Arithmetic mean: ', a:13:8);
  writeln;
  writeln(' Geometric mean : ', g:13:8);
  writeln;
  if ite=1 then
    writeln(' Harmonic mean undefined.')
  else
    writeln(' Harmonic mean  : ', h:13:8);
  writeln;
  writeln(' Moments:');
  writeln;
  writeln('   M1= ', xm1:13:8);
  writeln('   M2= ', xm2:13:8);
  writeln('   M3= ', xm3:13:8);
  writeln('   M4= ', xm4:13:8);
  writeln;
  writeln(' Flatness coefficient....: ', ap:13:8);
  writeln(' Coefficient of asymmetry: ', asy:13:8);
  writeln;
  writeln(' Gen. mean M(',t:2:0,'): ', xmt:13:8);
  writeln;
  ReadKey; DoneWinCrt

END.

{End of file moment.pas}