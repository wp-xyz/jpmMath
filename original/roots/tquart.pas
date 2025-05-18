{****************************************************************
*   Program to demonstrate Bisection & Quartile subroutines     *
* ------------------------------------------------------------- *
* References:                                                   *
*                                                               *
*  Bisection: BASIC Scientific Subroutines, Vol. II             *
*             By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981      *
*             [BIBLI 01].                                       *
*                                                               *
*  Quartile:  After an algorithm provided By Namir Shammas.     *
*                                                               *
*                          Pascal Release By J-P Moreau, Paris. *
*                                  (www.jpmoreau.fr)            *
* ------------------------------------------------------------- *
* Example:  Find a real root of f(x)=(x+1)^5                    *
*                                                               *
* Sample run:                                                   *
*                                                               *
* What is the initial range (X0,X1):                            *
*    X0 = -5                                                    *
*    X1 = 0                                                     *
* Convergence criterion: 1e-6                                   *
*                                                               *
*                                                               *
* Bisection: The calculated zero is X = -9.99999642372131E-0001 *
*                                                               *
* The associated Y value is Y =  5.85001220607623E-0033         *
*                                                               *
* The number of steps was: 23                                   *
*                                                               *
*                                                               *
* Quartile: The calculated zero is X = -1.00000008940697E+0000  *
*                                                               *
* The associated Y value is Y = -5.71290254499632E-0036         *
*                                                               *
* The number of steps was: 12                                   *
*                                                               *
*****************************************************************
Note: The quartile method is slightly faster and more accurate
      than the classical Bisection Method.
----------------------------------------------------------------}
PROGRAM TQuart;
Uses WinCrt;

VAR
        a,b,root, e,x,x0,x1 : DOUBLE;
        m : INTEGER;


{**********************************************}
Function Y(x:DOUBLE): DOUBLE;
Begin
  Y := (x+1)*(x+1)*(x+1)*(x+1)*(x+1)
End;
{**********************************************}


{**********************************************
*        Bisection method subroutine          *
* ------------------------------------------- *
* This routine iteratively seeks the zero of  *
* function Y(x) using the method of interval  *
* halving until the interval is less than e   *
* in width. It is assumed that the function   *
* Y(x) is available from a function routine.  *
* ------------------------------------------- *
* Input values: range (x0,x1), and e.         *
* Output values: root x, Y(x) and number of   *
* steps m.                                    *
**********************************************}
PROCEDURE Bisect;
Label 100, fin;
Var   y0,yy:double;
Begin
  m:=0;
  if Y(x0) * Y(x1) > 0.0 then
  begin
    writeln(' Bisection: No root found in given interval.');
    goto fin
  end
  else
  begin
100: y0:=Y(x0);
    x:=(x0+x1)/2; yy:=Y(x);
    m:=m+1;
    if yy*y0=0 then goto fin;
    if yy*y0<0 then x1:=x;
    if yy*y0>0 then x0:=x;
    if ABS(x1-x0)>e then goto 100
  end;
fin: End;

{***********************************************
*         Quartile method subroutine           *
* -------------------------------------------- *
* This routine iteratively seeks the zero of   *
* function Y(x) using the method of interval   *
* quarting until the interval is less than tol *
* in width. It is assumed that the function    *
* Y(x) is available from a function routine.   *
* -------------------------------------------- *
* Input values: range (a, b), and tol.         *
* Output values: found root and number of used *
* steps.                                       *
***********************************************}
PROCEDURE Quartile(a,b,tol:double; Var root:double; Var step:integer);
Label 10;
Var m, co: double;
Begin
  co:=0.25; {can be 0.3} step:=0;
  if Y(a) * Y(b) > 0.0 then
  begin
    writeln(' Quartile: No root found in given interval.');
    goto 10
  end
  else
  begin
    Repeat
      if abs(Y(a)) < abs(Y(b)) then
        m := a + co*(b-a)
      else
        m := a + (1.0-co)*(b-a);
      if Y(m)*Y(a) > 0.0 then a:= m else b:=m;
      Inc(step)
    Until abs(a-b) < tol;
    root := (a+b)/2.0;
  end;
10:End;

{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' What is the initial range (X0,X1):');
  writeln;
  write('    X0 = '); read(x0);
  write('    X1 = '); read(x1);
  writeln;
  write(' Convergence criterion: '); read(e);
  writeln;

  a:=x0; b:=x1;

  Bisect;   {Call bisection routine}

  if m<>0 then
  begin
    writeln;
    writeln(' Bisection: The calculated zero is X = ',x);
    writeln;
    writeln(' The associated Y value is Y = ',Y(x));
    writeln;
    writeln(' The number of steps was: ',m);
    writeln
  end;

  Quartile(a, b, e, root, m);

  if m<>0 then
  begin
    writeln;
    writeln(' Quartile: The calculated zero is X = ', root);
    writeln;
    writeln(' The associated Y value is Y = ',Y(root));
    writeln;
    writeln(' The number of steps was: ',m);
    writeln
  end;

  ReadKey; DoneWinCrt
END.

{End of file tquart.pas}