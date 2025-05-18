 {************************************************************************
 *                          DRAW  2D CURVES                              *
 * --------------------------------------------------------------------- *
 * This program can draw 2d curves that can be defined as:               *
 *              1)  y = f(x) or y = f(t)                                 *
 *              2)  Rho = f(x) or f(t)                                   *
 *              3)  x = f(t) and y = g(t)                                *
 * --------------------------------------------------------------------- *
 * Reference:                                                            *
 *  From "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0  *
 *        de R. Dony - MASSON 1990 page 113" [BIBLI 12].                 *
 *                                                                       *
 *                                    TPW Release By J-P Moreau, Paris.  *
 *                                           (www.jpmoreau.fr)           *
 ************************************************************************}    
 Program MULTIFUNC;
 Uses WinCrtMy,WinTypes,WinProcs,WObjects,Strings,Type_def,CrtGr2D,UFunct;

 Var  a,b,pas,ux,uy,x,y: real_ar;
      fmin,fmax : real_ar;
      c1,c2,c3,c4 : integer;
      dx,dy,valdep,valarr,xdep,ydep,xcour,ycour: real_ar;
      incr,periode,rmax,theta : real_ar; kind,choix: char;
      f,fx,fy,rho : FONC;


 Procedure SeekPeriod;
 var r,diff1,diff2,max: real_ar;
 begin
   max:=200*pi;
   periode:=0;
   case kind of
     '2' : begin
             xdep:=Evaluate(rho,periode);
             ydep:=0
           end;
     '3' : begin
             xdep:=Evaluate(FX,periode);
             ydep:=Evaluate(FY,periode)
           end
   end;
   Repeat
     periode:=periode+2*pi;
     case kind of
       '2' : begin
               r:=Evaluate(rho,periode);
               xcour:=r*cos(periode);
               ycour:=r*sin(periode);
             end;
       '3' : begin
               xcour:=Evaluate(FX,periode);
               ycour:=Evaluate(FY,periode);
             end;
     end;
     diff1:=abs(xcour-xdep);
     diff2:=abs(ycour-ydep);
   Until ((diff1 < 1E-6) and (diff2 < 1E-6)) or (periode > max)
 end;

 Procedure MaxRadius;
 var t,r: real_ar;
 begin
   rmax:=0;
   t:=valdep;
   while t <= valarr do
   begin
     case kind of
       '2' : begin
               r:=Evaluate(rho,periode);
               xcour:=r*cos(periode);
               ycour:=r*sin(periode);
             end;
       '3' : begin
               xcour:=Evaluate(FX,t);
               ycour:=Evaluate(FY,t);
             end;
     end;
     if abs(xcour) > rmax then rmax:=abs(xcour);
     if abs(ycour) > rmax then rmax:=abs(ycour);
     t:=t+incr
   end;
   rmax:=1.1*rmax
 end;

 Procedure Data;
 var coeff: real_ar;
     ch,ch1: CHAINE;
 begin
   clrscr;
   gotoxy(22,2);
   writeln('DRAW 2D CURVES');
   gotoxy(22,3);
   writeln('--------------');
   gotoxy(5,5);
   writeln('Input Angles in radians.');
   writeln;
   writeln('    1 : Y = F(X)');
   writeln('    2 : Rho = F(t)');
   writeln('    3 : X = FX(t) et Y = FY(t)');
   writeln;
   write('    What is the curve kind (1, 2 or 3): '); readln(kind);
   gotoxy(1,13);
   case kind of
     '1' : begin
             write('    Input function f(x): '); readln(ch);
             if Not EnterFunction(ch,f)
               then exit;
             writeln('    Interval [Xmin,Xmax] and step: ');
             write('      Xmin = '); readln(a);
             write('      Xmax = '); readln(b);
             write('      Step = '); readln(pas);
           end;
     '2','3' : begin
                 if kind='3' then
                 begin
                   writeln('    Input functions x(t) and y(t): ');
                   write('        FX(t)= '); readln(ch);
                   write('        FY(t)= '); readln(ch1);
                   if Not EnterFunction(ch,fx) then exit;
                   if Not EnterFunction(ch1,fy) then exit;
                 end
                 else if kind='2' then
                 begin
                   write('    Input function Rho(t): ');
                   readln(ch);
                   if Not EnterFunction(ch,rho) then exit;
                 end;

                 write('    Is it a periodic function (y/n): ');
                 readln(choix);
                 if upcase(choix)='Y' then
                 begin
                   write('    Period (0 if unknown): ');
                   readln(periode);
                   if periode=0 then SeekPeriod;
                   valdep:=0; valarr:=periode;
                   incr:=periode/511
                 end
                 else
                 begin
                   writeln('    Interval [tmin,tmax] and step: ');
                   write('      t1 = '); readln(valdep);
                   write('      t2 = '); readln(valarr);
                   write('      step = '); readln(incr);
                 end;
                 write('    Maximum Radius (0 if unknown): ');
                 readln(rmax);
                 if rmax=0 then MaxRadius
                           else rmax:=1.1*rmax
               end
   end
 end;

 Procedure DrawCurveXY;
 var x,y: real_ar;
 begin
   x:=a; y:=Evaluate(f,x);
   MoveXY(CrtDC,x,y);
   repeat
     y:=Evaluate(f,x);
     LineXY(CrtDC,x,y);
     x:=x+pas
   until x > b
 end;

 Procedure DrawParamCurve;
 var x1,y1,x2,y2,t,r: real_ar;
 begin
   t:=valdep;
   case kind of
     '2' : begin
             r:=Evaluate(rho,t);
             x1:=r*cos(t);
             y1:=1.25*r*sin(t)
           end;
     '3' : begin
             x1:=Evaluate(FX,t);
             y1:=1.25*Evaluate(FY,t)
           end
   end;
   MoveXY(CrtDC,x1,y1);
   while t <= valarr do
   begin
     t:=t+incr;
     case kind of
     '2' : begin
             r:=Evaluate(rho,t);
             x2:=r*cos(t);
             y2:=1.25*r*sin(t)
           end;
     '3' : begin
             x2:=Evaluate(FX,t);
             y2:=1.25*Evaluate(FY,t)
           end
     end;
     LineXY(CrtDC,x2,y2)
   end
 end;


 BEGIN
   Repeat
     WinCrtInit('DRAW A 2D CURVE');
     Data;
     choix:='O';
     Case kind of
        '1' : begin
                clrscr;
                fmin:=1E20; fmax:=-1E20; dx:=(b-a)/10.0;
                x:=a; y:=Evaluate(f,x);
                repeat
                  y:=Evaluate(f,x);
                  if y<fmin then fmin:=y;
                  if y>fmax then fmax:=y;
                  x:=x+pas
                until x > b;
                fmin:=Round(fmin-0.5); fmax:=Round(fmax);
                dy:=(fmax-fmin)/10.0;
                Fenetre(a,b,fmin,fmax);
                Cloture(50,MaxX-50,100,MaxY-10);
                Grille(CrtDC,dx/2,dy/2);
                Bordure(CrtDC);
                Axes(CrtDC);
                Gradue(CrtDC,dx,dy);
                DrawCurveXY;
                {Legendes(10,'Y = F (x)','','X','Y','');  }
              end;
        '2','3' : begin
                    clrscr;
                    Fenetre(-rmax,rmax,-rmax,rmax);
                    Cloture(50,MaxX-50,75,MaxY-10);
                    Bordure(CrtDC);
                    DrawParamCurve
                  end;
     End;
     SortieGraphique
   Until rep='n';
   DoneWinCrt
 END.

 {end of file multifon.pas}