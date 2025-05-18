{********************************************
* Program to demonstrate the Gamma Function *
*                                           *
*             Pascal version by J-P Moreau  *
*                  (www.jpmoreau.fr)        *
* ----------------------------------------- *
* Reference:                                * 
* "Numerical Recipes, By W.H. Press, B.P.   *
*  Flannery, S.A. Teukolsky and T. Vetter-  *
*  ling, Cambridge University Press, 1986"  *
*  [BIBLI 08].                              *
* ----------------------------------------- *
* SAMPLE RUN:                               *
*                                           *
*       X        Gamma(X)                   *   
*  -------------------------                *
*    0.5000      1.772454                   *
*    1.0000      1.000000                   *
*    1.5000      0.886227                   *
*    2.0000      1.000000                   *
*    2.5000      1.329340                   *
*    3.0000      2.000000                   *
*    3.5000      3.323351                   *
*    4.0000      6.000000                   *
*    4.5000     11.631728                   *
*    5.0000     24.000000                   *
*                                           *
********************************************}
Program Test_Gamma;
Uses WinCrt;

Const
       half = 0.5;
       one  = 1.0;
       fpf  = 5.5;
       zero = 0.0;

Var
       x,y : double;
       i   : integer;


{******************************************
*           FUNCTION  GAMMA(X)            *
* --------------------------------------- *
* Returns the value of Gamma(x) in double *
* precision as EXP(LN(GAMMA(X))) for X>0. *
******************************************}
Function Gamma(xx:double):double;
Var
  cof:Array[1..6] of double;
  stp,x,tmp,ser:double;
  j:integer;
Begin
  cof[1]:=76.18009173;
  cof[2]:=-86.50532033;
  cof[3]:=24.01409822;
  cof[4]:=-1.231739516;
  cof[5]:=0.120858003e-2;
  cof[6]:=-0.536382e-5;
  stp:=2.50662827465;
  
  x:=xx-one;
  tmp:=x+fpf;
  tmp:=(x+half)*LN(tmp)-tmp;
  ser:=one;
  for j:=1 to 6 do
  begin
    x:=x+one;
    ser:=ser+cof[j]/x
  end;
  Gamma := EXP(tmp+LN(stp*ser))
End;


{main program}
BEGIN
  writeln;
  writeln('      X        Gamma(X)   ');
  writeln(' -------------------------');
  x:=zero;
  for i:=1 to 10 do
  begin
    x := x + half;
    y := Gamma(x);
    writeln('  ',x:7:4,'  ',y:12:6)
  end;
  writeln;
  Readkey; DoneWinCrt
END.

{end of file gamma.pas}