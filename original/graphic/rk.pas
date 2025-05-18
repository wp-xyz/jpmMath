{***************************************************************
* View the solutions of a differential equations system having *
* the form:   y'=f(x,y,z),  z'=g(x,y,z)                        *
*                                                              *
*                                  TPW version By J-P Moreau   *
*                                      (www.jpmoreau.fr)       *
* ------------------------------------------------------------ *
* REFERENCE:                                                   *
* "Graphisme dans le plan et dans l'espace avec Turbo Pascal   *
*  4.0 By R. Dony - MASSON, Paris 1990, pages 254-260"         *
*  [BIBLI 12].                                                 *
* ------------------------------------------------------------ *
* TUTORIAL:                                                    *
* It can be proved that any differential equation of order > 1 *
* may be replaced by a 1st order differential equations system.*
* For example, the 2nd order equation y" = f(x,y,y') is equi-  *
* valent to the system:                                        *
*                         y' = z                               *
*                         z' = f(x,y,z)                        *
*                                                              *
* This program draws the solutions (integral curves) of a sys- *
* tem with the form:                                           *
*                        y' = f(x,y,z)                         *
*                        z' = g(x,y,z)                         *
*                                                              *
* The proposed example is:    y' = x*x*y*y - 1                 *
*                             z' = x*x+y*y - 4                 *
*                                                              *
* The steps to follow are:                                     *
*                                                              *
*     1. choose a physical coordinates window x1,x2,y1,y2      *
*     2. choose an increment step, dx                          *
*     3. for each point of plane (x,y) inside previous window, *
*        calculate by Runge-Kutta method and draw the solution *
*        curve y=f(x) passing through this point.              *
*                                                              *
* Obviously, these curves will eventually go out of screen. It *
* can be managed by the clipping capability of unit CrtGr2D.   *
* However to limit computation time, a curve is stopped when:  *
*                                                              *
*     - the following point will exit graph window,            *
*     - the number of segments already drawn > 300,            *
*     - the curve is approaching a critical point where        *
*     - f(x,y) and g(x,y) ==> 0.                               *
*                                                              *
* Other possible examples:                                     *
*                                                              *
*     1) y' = x + y            z' = x*y                        *
*     2) y' = y(y*y-1)         z' = x(x*x-1)                   *
*     3) y' = -x               z' = x*x + y*y - 1              *
*     4) y' = cos(x)-x*sin(y)  z' = x*sin(y) + y*sin(x)        *
*                                                              *
* ------------------------------------------------------------ *
* NOTE: the variable switch (in Runge-Kutta procedure) is used *
*       to "switch" direction. For each start point in plane   *
*       (x,y), the integral curve can be drawwn in two opposite*
*       directions. When switch=1, the curve is drawn in one   *
*       direction until a stop condition is met; at this moment*
*       switch is forced to -1 and the other part of the curve *
*       is drawn until a new stop condition is met. Then switch*
*       is again put to 1 and the process will continue with a *
*       new starting point.                                    *
***************************************************************}
 Program Runge_Kutta;
 Uses WinCrtMy,WinTypes,Type_def,CrtGr2D;

 CONST  nbreseg = 300;

 VAR

        f1,f2,f3,f4     : REAL_AR;
        i,j,h,xpos,ypos : REAL_AR;

       {Example:  y' = x*x*y*y-1,  z' = x*x+y*y-4 }

       FUNCTION F(x,y : REAL_AR): REAL_AR;
       begin
	 F:=x*x*y*y-1.0;
       end;

       FUNCTION G(x,y : REAL_AR): REAL_AR;
       begin
	 G:=x*x+y*y-4.0;
       end;


 PROCEDURE Data;
 Begin
   clrscr;
   writeln;
   writeln('     Integration by Runge-Kutta method');
   writeln('     _________________________________');
   gotoxy(3,10);
   write('  Input physical window x1,x2,y1,y2: ');
   readln(f1,f2,f3,f4);
   gotoxy(3,11);
   write('  Input drawing step dx: '); readln(h)
 End;

 {Integration by Runge-Kutta method}
 PROCEDURE RungeKutta(P:HDC; y0, z0 : REAL_AR);
   
   var y,z,y1,z1 : REAL_AR;
       seg,switch : INTEGER;
       pointcrit : REAL_AR;
       fin : INTEGER;
       k1,k2,k3,k4,l1,l2,l3,l4:REAL_AR;

 Begin
   switch:=1;
   REPEAT
     y:=y0; z:=z0;
     MoveXY(P,y,z);
     y1:=y; z1:=z; seg:=0;
     Repeat
       LineXY(P,y1,z1);
       Inc(seg);
       y:=y1; z:=z1;
       k1:=switch*h*F(y,z);
       l1:=switch*h*G(y,z);
       y:=y1+k1/2.0; z:=z1+l1/2.0;
       k2:=switch*h*F(y,z);
       l2:=switch*h*G(y,z);
       y:=y1+k2/2.0; z:=z1+l2/2.0;
       k3:=switch*h*F(y,z);
       l3:=switch*h*G(y,z);
       y:=y1+k3/2.0; z:=z1+l3/2.0;
       k4:=switch*h*F(y,z);
       l4:=switch*h*G(y,z);
       y1:=y+switch*h/6.0*(k1+2*k2+2*k3+k4);
       z1:=z+switch*h/6.0*(l1+2*l2+2*l3+l4);
       pointcrit:=sqrt(sqr(F(y,z))+sqr(G(y,z)));
       y1:=y+switch*h*F(y,z);
       z1:=z+switch*h*G(y,z);
     Until (y1<f1) OR (y1>f2) OR (z1<f3) OR (z1>f4) OR (seg>nbreseg) OR (pointcrit<1E-4);
     switch:=-switch;
   UNTIL switch=1;
 End; { RungeKutta }


 PROCEDURE Initgraph(P:HDC);
 Begin
   clrscr;
   Fenetre(f1,f2,f3,f4);
   Cloture(40,MaxX-25,95,MaxY-5);
   Gradue(P,(f2-f1)/10,(f4-f3)/10);
   Grille(P,(f2-f1)/10,(f4-f3)/10);
   Axes(P);
   Bordure(P)
 End;


 {main program}
 BEGIN
   WinCrtInit('RUNGE-KUTTA');
   REPEAT
     Data;
     Initgraph(CrtDC);
     i:=0.0;
     While i < 1 do
     begin
       i:=i+h;
       xpos:=f1+i*(f2-f1);
       j:=0;
       While j < 1 do
       begin
         j:=j+h;
         ypos:=f3+j*(f4-f3);
         RungeKutta(CrtDC,xpos,ypos)
       end
     end;
     SortieGraphique
   UNTIL rep='n';
   DoneWinCrt
 END.

{End of file rk.pas}