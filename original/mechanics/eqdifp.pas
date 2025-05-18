{*******************************************************************
*          Differential equations of order 1 to N                  *
*             by Runge-Kutta method of order 4                     *
* ---------------------------------------------------------------- *
* Reference: "Analyse en Turbo Pascal versions 5.5 et 6.0" By Marc *
*             DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991     *
*             [BIBLI 03].                                          *
*                                                                  *
*             See demonstration programs teqdif1.pas and           *
*             pendulum.pas.                                        * 
*******************************************************************}
UNIT eqdifp;

INTERFACE

USES WinCrt, Type_def;

CONST MAXSIZE = 9;

TYPE
     Table  = Array[0..MAXSIZE] of real_ar;
     FPrime = FUNCTION(X,Y:real_ar): real_ar;
     Fn     = FUNCTION(X:real_ar;Y:Table): real_ar;
     Fp1    = FUNCTION(K:INTEGER;X:real_ar;Y:Table): real_ar;

PROCEDURE Display(V1,V2:RV; ndata:INTEGER; xi,xf:real_ar);
PROCEDURE Equadiffp( fp:Fp1; VAR t1, t2:RV; xi,xf:real_ar; Yi:Table;m,p,fi:INTEGER );

IMPLEMENTATION

PROCEDURE Aff_reel(x:real_ar);
Begin
  write(x:10:6)
End;


PROCEDURE Display;
Var i:INTEGER;
    h,x:real_ar;
    fp1,fp2: TEXT;  {Optional output files for theta and theta'}
Begin
  Assign(fp1,'theta.asc'); Rewrite(fp1);
  Assign(fp2,'theta1.asc'); Rewrite(fp2);
  writeln(fp1,'Theta');
  writeln(fp1,ndata);
  writeln(fp2,'Theta''');
  writeln(fp2,ndata);
  h:=(xf-xi)/(ndata-1); 
  x:=xi-h;
  Clrscr;
  writeln;
  writeln('    Time      Theta     Theta''   ');
  writeln('----------------------------------');
  for i:=1 to ndata do
  begin
    x:=x+h;
    aff_reel(x); aff_reel(V1^[i]); aff_reel(V2^[i]); writeln;
    writeln(fp1,x:10:6,V1^[i]:10:6);
    writeln(fp2,x:10:6,V2^[i]:10:6)
  end;
  close(fp1); close(fp2);
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
     t2^[i+1]:=z[1]
   end;
   Display(t1,t2,m,xi,xf)
 End;


END.

{End of file eqdifp.pas}