{******************************************************
*  Program to demonstrate the Newton root subroutine  *
*             in the complex domain                   *
* --------------------------------------------------- *
*  Reference; BASIC Scientific Subroutines, Vol. II   *
*  By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                     *
*               Pascal Version By J-P Moreau, Paris.  *
*                        (www.jpmoreau.fr)            *
* --------------------------------------------------- *
* SAMPLE RUN:                                         *
*                                                     *
* (Example: find a complex root of  z^2 + 1 = 0)      *
*                                                     *
* What is the initial guess:                          *
*                                                     *
*    X0 = .2                                          *
*    Y0 = 1.2                                         *
*                                                     *
* Convergence criterion: 1e-8                         *
* Maximum number of iterations: 10                    *
*                                                     *
* The root estimate is:                               *
*                                                     *
*    X0 = -0.000000                                   *
*    Y0 =  1.000000                                   *
*                                                     *
* The number of iterations performed was: 5           *
*                                                     *
******************************************************}
PROGRAM DEMO_ZNEWTON;
Uses WinCrt;

Var
        k,n : INTEGER;
        a,e,x0,y0,u,v,u1,v1,u2,v2,x,y : DOUBLE;


{*****************************************************
  Functions subroutine                               }
PROCEDURE Eval(x,y:double;VAR u,v,u1,v1,u2,v2:double);
Begin
  u:=x*x-y*y+1.0 ; v:=2.0*x*y;
  u1:=2.0*x ; u2:=-2.0*y;
  v1:=2.0*y ; v2:=2.0*x
End;
{****************************************************}

{************************************************
*  Complex root seeking using Newton's method   *
* --------------------------------------------- *
* This routine uses the complex domain form of  *
* Newton's method for iteratively searching     *
* for roots. The complex function and its first *
* partial derivatives must be available in the  *
* form; F(Z) = U(X,Y) + I V(X,Y). The required  *
* derivatives are DU/DX and DU/DY.              *
* --------------------------------------------- *
* INPUTS; initial guess, x0 and y0, convergence *
* criteria, e, maximum number of iterations, n. *
* OUTPUTS; approximation to the root, x and y,  *
* number of performed iterations, k.            *
************************************************}
PROCEDURE ZNewton;
Label 100,fin;
Const TINY = 1e-12;
Begin
  k:=0;
100: k:=k+1;
  {Get u,v and the derivatives u1,v1,u2,v2}
  Eval(x0,y0,u,v,u1,v1,u2,v2);
  a:=u1*u1+u2*u2;
  {Guard against a=0}
  if a < TINY then
  begin
    writeln(' ZERO DIVIDE ERROR - u1*u1+u2*u2 must be <> 0.');
    goto fin
  end;
  x:=x0+(v*u2-u*u1)/a;
  y:=y0-(v*u1+u*u2)/a;
  {Check for convergence in euclidean space}
  if (x0-x)*(x0-x)+(y0-y)*(y0-y)<=e*e then goto fin;
  if k>=n then goto fin;
  x0:=x ; y0:=y;
  goto 100;
fin: End;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' What is the initial guess:');
  writeln;
  write('    X0 = '); read(x0);
  write('    Y0 = '); read(y0);
  writeln;
  write(' Convergence criterion: '); read(e);
  write(' Maximum number of iterations: '); read(n);
  writeln;

  ZNewton;     {Call ZNewton subroutine}

  if a<>0 then
  begin
    writeln;
    writeln(' The root estimate is:');
    writeln;
    writeln('    X0 = ',x:9:6);
    writeln('    Y0 = ',y:9:6);
    writeln;
    writeln(' The number of iterations performed was: ',k);
    writeln
  end;
  ReadKey; DoneWinCrt
END.

{End of file znewton.pas}





