{*******************************************************************
*          Differential equations of order 1 to N                  *
*             by Runge-Kutta method of order 4                     *
* ---------------------------------------------------------------- *
* Reference: "Analyse en Turbo Pascal versions 5.5 et 6.0" by Marc *
*             DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991     *
*             [BIBLI 03].                                          *
*                                                                  *
*             See demonstration program teqdif1.pas.               *
*                                                                  *
* Release 1.1 (Sept. 2008): added procedure display.               *
*                                                                  *
*                                             J-P Moreau, Paris.   *
*                                             (www.jpmoreau.fr)    *
*******************************************************************}
UNIT eqdif1;

INTERFACE

USES WinCrtMy, Type_def;

CONST MAXSIZE = 9;

TYPE
     Table  = Array[0..MAXSIZE] of real_ar;
     FPrime = FUNCTION(X,Y:real_ar): real_ar;
     Fn     = FUNCTION(X:real_ar;Y:Table): real_ar;
     Fp1    = FUNCTION(K:INTEGER;X:real_ar;Y:Table): real_ar;

PROCEDURE Affiche(V:RV; ndata:INTEGER; xi,xf:real_ar);
PROCEDURE Display(V1,V2,V3,V4:RV; ndata,p:INTEGER; xi,xf:real_ar);
PROCEDURE Equadiff1( fp:FPrime; VAR t:RV; xi,xf,yi:real_ar; m,fi:INTEGER );
PROCEDURE Equadiffn( fp:Fn; VAR t:RV; xi,xf:real_ar; Yi:Table;m,n,fi:INTEGER );
PROCEDURE Equadiffp( fp:Fp1; t1,t2,t3,t4:RV; xi,xf:real_ar; Yi:Table;m,p,fi:INTEGER );

IMPLEMENTATION

PROCEDURE Aff_reel(x:real_ar);
Begin
  write(x:12:6)
End;

{Display results for size=1}
PROCEDURE Affiche;
Var i:INTEGER;
    h,x:real_ar;
Begin
  h:=(xf-xi)/(ndata-1); 
  x:=xi-h;
  writeln;
  writeln('      X           Y       ');
  writeln(' -------------------------');
  for i:=1 to ndata do
  begin
    x:=x+h;
    aff_reel(x); aff_reel(V^[i]); writeln
  end;
  writeln(' -------------------------')
End;

{Display results for size > 1}
PROCEDURE Display;
Var i:INTEGER;
    h,x:real_ar;
Begin
  h:=(xf-xi)/(ndata-1); 
  x:=xi-h;
  Clrscr;
  writeln;
  write('       X');
  For i:=1 to p do write('          Y',i:1);
  writeln;
  writeln('--------------------------------------------------------------');
  for i:=1 to ndata do
  begin
    x:=x+h;
    if p=2 then
    begin
      aff_reel(x); aff_reel(V1^[i]); aff_reel(V2^[i]); writeln
    end
    else if p=3 then
    begin
      aff_reel(x); aff_reel(V1^[i]); aff_reel(V2^[i]); aff_reel(V3^[i]); writeln
    end
    else if p>3 then
    begin
      aff_reel(x); aff_reel(V1^[i]); aff_reel(V2^[i]); aff_reel(V3^[i]); aff_reel(V4^[i]); writeln
    end
  end;
  writeln('--------------------------------------------------------------')
End;

{**************************************************************************
*        SOLVING DIFFERENTIAL EQUATIONS WITH 1 VARIABLE OF ORDER 1        *
*                       of type y' = f(x,y)                               *
* ----------------------------------------------------------------------- *
*  INPUTS:                                                                *
*    fp        Given equation to integrate (see test program)             *
*    xi, xf    Begin, end values of variable x                            *
*    yi        Begin Value of y at x=xi                                   *
*    m         Number of points to calculate                              *
*    fi        finesse (number of intermediate points)                    *
*                                                                         *
*  OUTPUTS:                                                               *
*    t       : pointer to real vector storing m results for function y    *
**************************************************************************}
PROCEDURE  Equadiff1;
Var
  h, a, b, c, d, x, y : real_ar;
  i, j, ni : INTEGER;
Begin
  if (m > SIZE) OR (fi < 1) then exit;
  h   := (xf - xi) / fi / (m-1);
  y   := yi;
  t^[1] := yi; 

  for i := 1 to m do
  begin
    ni := (i - 1) * fi - 1;
    for j := 1 to fi do
    begin
      x := xi + h * (ni + j);
      a := h * fp(x, y);                  
      b := h * fp(x + h / 2, y + a / 2);  
      c := h * fp(x + h / 2, y + b / 2);  
      x := x + h;
      d := h * fp(x, y + c);              
      y := y + (a + b + b + c + c + d) / 6
    end;
    t^[i+1] := y
  end;
  Affiche(t,m,xi,xf)
End;

{**************************************************************************
*        SOLVING DIFFERENTIAL EQUATIONS WITH 1 VARIABLE OF ORDER N        *
*                of type y(n) = f(y(n-1),y(n-2),...,y',y,x)               *
* ----------------------------------------------------------------------- *
*  INPUTS:                                                                *
*    f         Equation to integrate                                      *
*    xi, xf    Begin, end values of variable x                            *
*    Yi        Begin values at xi (of f and derivatives)                  *
*    m         Number of points to calculate                              *
*    n         Order of differential equation                             *
*    fi        finesse (number of intermediate points)                    *
*                                                                         *
*  OUTPUTS:                                                               *
*    t         real vector storing m results for function y               *
* ----------------------------------------------------------------------- *      
*  EXAMPLE:    y" = (2 y'y' + yy) / y with y(4) = 2, y'(4) = -2*tan(1)    *
*              Exact solution:  y = (2 cos(1))/cos(x-5)                   *
**************************************************************************}
 PROCEDURE Equadiffn;
 Var h,x,a,b,c,d : real_ar;
     ta,tb,tc,td,y,z : Table;
     i,j,k,ni,n1,n2 : INTEGER;
 Begin
   if (fi<1) OR (n>MAXSIZE) OR (m>SIZE) then exit;
   h := (xf - xi) / fi / (m-1);
   n1:=n-1; n2:=n-2;
   t^[1]:=Yi[0];
   for k:=0 to n1 do
   begin
     y[k]:=Yi[k]; z[k]:=Yi[k]
   end;
   for i:=1 to m do
   begin
     ni:=(i-1)*fi-1;
     for j:=1 to fi do
     begin
       x:=xi+h*(ni+j);
       for k:=0 to n1 do y[k]:=z[k];
       a:=h*fp(x,y);
       for k:=0 to n2 do
       begin
         ta[k]:=h*y[k+1]; y[k]:=z[k]+ta[k]/2
       end;
       y[n1]:=z[n1]+a/2;
       x:=x+h/2;
       b:=h*fp(x,y);
       for k:=0 to n2 do
       begin
         tb[k]:=h*y[k+1]; y[k]:=z[k]+tb[k]/2
       end;
       y[n1]:=z[n1]+b/2;
       c:=h*fp(x,y);
       for k:=0 to n2 do
       begin
         tc[k]:=h*y[k+1]; y[k]:=z[k]+tc[k]
       end;
       y[n1]:=z[n1]+c;
       x:=x+h/2;
       d:=h*fp(x,y);
       for k:=0 to n2 do
         z[k]:=z[k]+(ta[k]+2*tb[k]+2*tc[k]+h*y[k+1])/6;
       z[n1]:=z[n1]+(a+b+b+c+c+d)/6
     end;
     t^[i+1]:=z[0]
   end;
   Affiche(t,m,xi,xf)
 End;

{***************************************************************************
*         SOLVING DIFFERENTIAL SYSTEMS WITH P VARIABLES OF ORDER 1         *
*                 of type yi' = f(y1,y2,...,yn), i=1..n                    *
* ------------------------------------------------------------------------ *
*  INPUTS:                                                                 *
*    f         table of equations to integrate                             *
*    xi, xf    begin, end values of variable x                             *
*    yi        table of begin values of functions at xi                    *
*    p         number of independant variables                             *
*    fi        finesse (number of intermediate points)                     *
*                                                                          *
*  OUTPUTS:                                                                *
*    t1,t2     real vectors storing the results for first two functions,   *
*              y1 and y2.                                                  *
* ------------------------------------------------------------------------ *      
*  EXAMPLE:    y1'=y2+y3-3y1, y2'=y1+y3-3y2, y3'=y1+y2-3y3                 *
*              Exact solution :  y1 = 1/3 (exp(-4x)  + 2 exp(-x))          *
*                                y2 = 1/3 (4exp(-4x) + 2 exp(-x))          *
*                                y3 = 1/3 (-5exp(-4x)+ 2 exp(-x))          *
***************************************************************************}
 PROCEDURE Equadiffp;
 Var h,x : real_ar;
     ta,tb,tc,td,y,z : Table;
     i,j,k,ni : INTEGER;
 Begin
   if (fi<1) OR (p>MAXSIZE) OR (m>SIZE) then exit;
   h := (xf - xi) / fi / (m-1);
   p:=p-1;
   t1^[1]:=Yi[0];
   t2^[1]:=Yi[1];
   t3^[1]:=Yi[2];
   t4^[1]:=Yi[3];
   for k:=0 to p do
   begin
     y[k]:=Yi[k]; z[k]:=Yi[k]
   end;
   for i:=1 to m do
   begin
     ni:=(i-1)*fi-1;
     for j:=1 to fi do
     begin
       x:=xi+h*(ni+j);
       for k:=0 to p do y[k]:=z[k];
       for k:=0 to p do ta[k]:=h*fp(k,x,y);
       for k:=0 to p do y[k]:=z[k]+ta[k]/2;
       x:=x+h/2;
       for k:=0 to p do tb[k]:=h*fp(k,x,y);
       for k:=0 to p do y[k]:=z[k]+tb[k]/2;
       for k:=0 to p do tc[k]:=h*fp(k,x,y);
       for k:=0 to p do y[k]:=z[k]+tc[k];
       x:=x+h/2;
       for k:=0 to p do td[k]:=h*fp(k,x,y);
       for k:=0 to p do
         z[k]:=z[k]+(ta[k]+2*tb[k]+2*tc[k]+td[k])/6;
     end;
     t1^[i+1]:=z[0];
     t2^[i+1]:=z[1];
     t3^[i+1]:=z[2];
     t4^[i+1]:=z[3]
   end;
   Display(t1,t2,t3,t4,m,p+1,xi,xf)
 End;


END.

{End of file eqdif1.pas}