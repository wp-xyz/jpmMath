{**********************************************************
*        SOLVING DIFFERENTIAL EQUATIONS OF ORDER 1        *
*             (Prediction-correction method)              *
* ------------------------------------------------------- *
* Reference:                                              * 
*   "Analyse en Turbo Pascal versions 5.5 et 6.0 By Marc  *
*    DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991"    *
*    [BIBLI 03].                                          *
*                                                         *
*                       TPW Version By J-P Moreau, Paris  *
*                              (www.jpmoreau.fr)          *
* ------------------------------------------------------- *
* SAMPLE RUN:                                             *
* (Integrate y'=(4xy+4x*sqrt(y))/(1+x²) from x=0 to x=10  *
*  with initial condition y(0)=1)                         *
*                                                         *
*  Input x begin: 0                                       *                 
*  Input x end  : 10                                      *
*  Input y begin: 1                                       * 
*  Input number of points: 10                             *
*  Input finesse: 10                                      *
*                                                         *
*         X               Y                               *
*  -----------------------------                          *
*     1.000000         8.999984                           *
*     2.000000        80.999881                           *
*     3.000000       360.999496                           *
*     4.000000      1088.998512                           *
*     5.000000      2600.996484                           *
*     6.000000      5328.992838                           *
*     7.000000      9800.986875                           *
*     8.000000     16640.977767                           *
*     9.000000     26568.964559                           *
*    10.000000     40400.946171                           *
*  -----------------------------                          *
*                                                         *
**********************************************************}
PROGRAM TEQDIF_CP;
Uses WinCrt;

CONST   MAXDATA = 100;

TYPE
        pTab  = ^Table;
        Table = Array[0..MAXDATA] of REAL;

VAR
        h,x,xi,xf,yi:REAL;
        tx,ty : pTab;
        fi,i,n:INTEGER;

  { y'=(4xy+4x*sqrt(y))/(1+x²) }
  Function f(x, y:REAL):REAL;
  begin
    if y>=0.0 then
      f := (4.0*x*y+4.0*x*sqrt(y))/(1.0+x*x)
    else
      f := 0.0
  end; 

  {classical Runge-Kutta method of order 4}
  Procedure runge_kutta(VAR y:REAL);
  Var a,b,c,d:REAL;
  Begin
    a:=h*f(x,y);
    b:=h*f(x+h/2,y+a/2);
    c:=h*f(x+h/2,y+b/2);
    x:=x+h;
    d:=h*f(x,y+c);
    y := y + (a+b+b+c+c+d)/6
  End;

  {*********************************************************
  *             Prediction-correction method               *
  * ------------------------------------------------------ *
  * INPUTS:                                                *
  *	      xi      begin x value                        *
  *	      xf      end x value                          *
  *           y1      begin y value (at x=xi)              *
  *           m	      number of points to calculate        *
  *           fi      finesse (number of intermediate      *
  *                   points (for example 10)              *
  * OUTPUTS:                                               *
  *           tx      table of x values (m values)         *
  *           ty      table of y values (m values)         * 
  *                                                        *
  * DESCRIPTION:                                           *
  * The algorithm has the following steps:                 *
  *     1. calculate y2, y3, y4 using a classical Runge-   *
  *        Kutta method of order 4.                        *
  *     2. from point (x4, y4), first estimate y(n+1) by   *
  *        formula:                                        * 
  *        y(n+1)=y(n) + H/24(55y'(n)-59y'(n-1)+37y'(n-2)  *
  *               -9y'(n-3)                                *
  *        then continue with formula:                     *
  *        y(n+1)=y(n) + H/24(9y'(n+1)+19y'(n)-5y'(n-1)    *
  *               +y'(n-2),                                *
  *        noting that y'(n+1)=f(x(n+1),y(n+1)) with the   *
  *        estimated value of y(n+1), until convergence is *
  *        obtained.                                       *
  *********************************************************}
  Procedure equadiff_pc(VAR tx,ty:pTab; xi, xf, yi:REAL; m,fi:INTEGER);
  Var
    z,y,w:REAL; p:Array[0..4] of REAL;
    i,j,k:INTEGER; ni:LONGINT;
  Begin
    writeln(' Arrived in eqdifpc.');
    z:=yi;
    if (m>MAXDATA) or (fi<1) then exit;
    h := (xf-xi)/fi/m;
    p[3]:=f(xi,yi);
    tx^[0]:=xi; ty^[0]:=yi;
    k:=0;
    for i:=1 to m do
    begin
      ni:= (i-1)*fi-1;
      for j:=1 to fi do
      begin
	x:=xi+h*(ni+j);
        Inc(k);
	if k<4 then
        begin
	  runge_kutta(z);
	  p[3-k]:=f(x,z)
        end
	else
        begin
	  x:=x+h;
          w:=z+h/24.0*(55.0*p[0]-59.0*p[1]+37.0*p[2]-9.0*p[3]);
          Repeat
	    y:=w;
            w:=z+h/24.0*(9.0*f(x,y)+19.0*p[0]-5.0*p[1]+p[2])
          Until abs(y-w) < 1e-10;
	  z:=w; p[3]:=p[2]; p[2]:=p[1];
	  p[1]:=p[0]; p[0]:=f(x,z)
        end
      end; {j loop}
      tx^[i]:=x; ty^[i]:=z
    end {i loop}
  End;


  {main program}
  BEGIN

    New(tx); New(ty);
    writeln;
    write(' Input x begin: '); read(xi);
    write(' Input x end  : '); read(xf);
    write(' Input y begin: '); read(yi);
    write(' Input number of points: '); read(n);
    write(' Input finesse: '); read(fi);

    equadiff_pc(tx,ty,xi,xf,yi,n,fi);

    {print results}
    writeln;
    writeln('        X               Y     ');
    writeln(' -----------------------------');
    for i:=1 to n do 
      writeln(tx^[i]:12:6,'     ',ty^[i]:12:6);
    writeln(' -----------------------------');
    writeln;
    Dispose(tx); Dispose(ty);
    ReadKey; DoneWinCrt

  END.

{end of file teqdifpc.pas}