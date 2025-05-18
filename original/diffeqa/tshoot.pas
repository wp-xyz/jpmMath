{*******************************************************************
*             Differential equations of order N                    *
*             by Runge-Kutta method of order 4                     *
*                with boundary conditions                          *
* ---------------------------------------------------------------- *
* Reference: "Analyse en Turbo Pascal versions 5.5 et 6.0" by Marc *
*             DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991.    *
*                                                                  *
*                                Modified By J-P Moreau, Paris     *
*                                 (added boundary condition).      *
*                                      (www.jpmoreau.fr)           *
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
*                                                                  *
* Example: integrate y"=(2y'y'+yy)/y from x=4 to x=6               *
*          with initial conditions: y(4)=2 and y'(4)=-2            *
*          and final condition y(6)=2                              *
*                                                                  *
* Exact solution is:   y = 2cos(1)/cos(x-5)                        *
*                                                                  *
*        DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER N          *
*              of type y(n) = f(y(n-1),...,y'',y,x)                *
*                   with boundary conditions                       *
*                                                                  *
*   order of equation: 2                                           *
*                                                                  *
*   begin value x    : 4                                           *
*   end value x      : 6                                           *
*   y0 value at x0   : 2                                           *
*   y1 value at x0   : -2                                          *
*   Desired value at xf: 2                                         *
*   number of points   : 11                                        *
*   finesse            : 20                                        *
*   Max. iterations    : 200                                       *
*                                                                  *
*   Last iteration:                                                *
*   u  =-3.11481552124023E+0000  du  =-7.62939453125000E-0007      *
*   dyb= 8.19545675767586E-0008  67 steps                          *
*                                                                  *
*       X            Y                                             *
* ---------------------------                                      *
*   4.000000     2.000000                                          *
*   4.200000     1.551018                                          *
*   4.400000     1.309291                                          *
*   4.600000     1.173217                                          *
*   4.800000     1.102583                                          *
*   5.000000     1.080605                                          *
*   5.200000     1.102583                                          *
*   5.400000     1.173217                                          *
*   5.600000     1.309291                                          *
*   5.800000     1.551018                                          *
*   6.000000     2.000000                                          *
*                                                                  *
*******************************************************************}
Program TShoot;
Uses WinCrtMy, Type_def;

TYPE

  Table  = Array[0..9] of real_ar;
  Fn     = FUNCTION(X:real_ar;Y:Table): real_ar;

VAR

  fi,i,ndata,ordre: INTEGER;
  xi,xf : REAL_AR;
  yi    : Table;
  vect  : RV;      {see type_def.pas}

  du,dyb,olddyb,u,yb: Real_ar;
  maxiter,niter:Integer;

  {Example: y" = 4y + 0.1*x^2}
  FUNCTION fp(x:REAL_AR;y:Table):REAL_AR; FAR;
  Begin
    fp:=(2.0*y[1]*y[1]+y[0]*y[0])/y[0]
  End;

  {Example: y" = -10*y^3
  FUNCTION fp(x:REAL_AR;y:Table):REAL_AR; FAR;
  Begin
    fp:=-10.0*y[0]*y[0]*y[0]
  End;  }

PROCEDURE Disp_real(x:real_ar);
Begin
  write(x:10:6)
End;

PROCEDURE Display(V:RV; ndata:INTEGER; xi,xf:real_ar);
Var i:INTEGER;
    h,x:real_ar;
Begin
  h:=(xf-xi)/(ndata-1); 
  x:=xi-h;
  writeln;
  writeln('      X         Y     ');
  writeln('----------------------');
  for i:=1 to ndata do
  begin
    x:=x+h;
    Disp_real(x); Disp_real(V^[i]); writeln
  end
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
*    flag      False: no printing  True: print results                    *
*                                                                         *
*  OUTPUTS:                                                               *
*    t         real vector storing m results for function y               *
* ----------------------------------------------------------------------- *      
*  EXAMPLE:    y" = (2 y'y' + yy) / y with y(4) = 2, y'(4) = -2*tan(1)    *
*              Exact solution:  y = (2 cos(1))/cos(x-5)                   *
**************************************************************************}
 PROCEDURE Equadiffn(fp:Fn; var t:RV; xi,xf:real_ar; Yi:Table;m,n,fi:INTEGER;flag:Boolean );
 Const MAXSIZE=9;
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
   if flag then Display(t,m,xi,xf)
 End;


BEGIN

  New(vect);

  Clrscr;
  writeln;
  writeln('    DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER N');
  writeln('          of type y(n) = f(y(n-1),...,y'',y,x)   ');
  writeln('               with boundary conditions');
  writeln;
  write('    order of equation: '); readln(ordre);
  writeln;
  write('    begin value x    : '); readln(xi);
  write('    end value x      : '); readln(xf);
  for i:=0 to ordre-1 do
  begin
    write('    y',i,' value at x0   : '); readln(yi[i])
  end;
  write('    Desired value at xf: '); readln(yb);
  write('    number of points   : '); readln(ndata);
  write('    finesse            : '); readln(fi);
  write('    Max. iterations    : '); readln(maxiter);

  niter:=1;
  u:=yi[1];         {First derivative initial value}
  du:=Abs(u)/10.0;
  {estimate initial error on desired final value}
  equadiffn(fp,vect,xi,xf,yi,ndata,ordre,fi,False);
  olddyb:=abs(vect^[ndata]-yb);

  {Iterate until desired final value yb is obtained}
  Repeat
    u:=u+du; yi[1]:=u;
    equadiffn(fp,vect,xi,xf,yi,ndata,ordre,fi,False);
    dyb:=abs(vect^[ndata]-yb);
    if abs(dyb)>abs(olddyb) then  {overshoot?}
      if niter=1 then
        du:=-du
      else
        du:=-du/2.0;
    olddyb:=dyb; Inc(niter);
  Until (dyb<0.0000001) or (niter>=maxiter);

  writeln;
  writeln('    Last iteration:');
  writeln('    u  =',u,'  du  =',du);
  writeln('    dyb=',dyb,'  ',niter,' steps');
  readkey;

  {print final results}
  equadiffn(fp,vect,xi,xf,yi,ndata,ordre,fi,True);
   
  Readkey;
  Dispose(vect);
  DoneWinCrt

END.

{End of file tshoot.pas}