{**************************************************************************
*                   FRACTALS:   HENON'S  ATTRACTORS                       *
* ----------------------------------------------------------------------- *
* TUTORIAL:                                                               *
*                                                                         *
* The french politician and mathematician Poincaré highlighted the com-   *
* plexity of planet trajectories that simple laws of celestial mechanics  *
* do not show at first. A typical example is the movement of the Earth    *
* around the sun, everybody knows that the period is one year long. This  *
* is only true in first approximation, but not if you take into conside-  *
* ration all the small perturbations caused by the other planets !        *
*                                                                         *
* Periodic trajectories are useful for at least two reasons:              *
*                                                                         *
*  1) they can be predicted with accuracy by Kepler's laws,               *
*                                                                         *
*  2) it is possible to describe conveniently the small trajectory        *
*     variations by using the following:                                  *
*                                                                         *
* Let be T the basic periodic trajectory that we intercept by a fixed,    *
* perpendicular plane, and O the point, taken as origin, of the plane     *
* intercepted by the T trajectory. At next periods, the trajectory with   *
* its small variations will meet the same plane at nearby points A0, A1,  *
* A2 etc.                                                                 *
*                                                                         *
* This allows us to replace a complex problem of 3D trajectories by a     *
* simpler problem of 2D representation of points.                         *
*                                                                         *
* The french mathematician Hénon described in 1969 a dynamic system which *
* illustrates the proposed method:                                        *
*                                                                         *
* 		x n+1 = x n cos a - ( y n - x n² ) sin a                  *
*		y n+1 = x n sin a + ( y n - x n² ) cos a                  *
*                                                                         *
* a, being the angle of basic trajectory with respect to the observation  *
* plane (in radians).                                                     *
*                                                                         *
* The program randomly chooses an angle a (usually between 0.5 and 2.0    *
* radians). The number of iterations equals arbitrarily 200.              * 
*                                                                         *
* As one can see, caracreristic figures are obtained, the points have a   *
* noticeable tendancy to agglutinate themselves in particular "zones of   *
* attraction".                                                            *
* ----------------------------------------------------------------------- *
* REFERENCE:                                                              *
*     "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0      *
*      By R. Dony - MASSON, Paris 1990" [BIBLI 12].                       *
*                                          TPW version by J-P Moreau      *
*                                              (www.jpmoreau.fr)          *
**************************************************************************}
 PROGRAM HENON;     {WinCrtMy version with printing capability}
 USES WinCrtMy,WinTypes,WinProcs,Strings,Type_def,CrtGr2D,WinPrint;

 LABEL 100;

 VAR  f1,f2,f3,f4 : real_ar;
      angle,x,y,xpos,ypos,cosinus,sinus,dx : real_ar;
      compte,nbreiter,xcur,ycur : integer;
      s:array[0..20] of char; s5:string[5];
      Pinfo : PPrinterInfo;

 PROCEDURE Data;
 BEGIN
   clrscr;
   writeln;
   writeln('        HENON''S ATTRACTORS');
   writeln('        ------------------');
   writeln;
   write('  Input angle in rad. ( 1.5732 ): '); readln(angle);
   cosinus:=cos(angle); sinus:=sin(angle);
   nbreiter:=250;
   dx:=0.01;
 END;

 PROCEDURE Henon1(P:HDC;x,y: real_ar);
 VAR aux: real_ar;
     n  : integer;
 BEGIN
   FOR n:=0 TO nbreiter DO
   BEGIN
     MoveXY(P,x,y);
     LineXY(P,x+dx,y);
     aux:=x;
     IF x < 1E19 THEN
     BEGIN
       x:=x*cosinus-(y-x*x)*sinus;
       y:=aux*sinus+(y-aux*aux)*cosinus;
     END;
   END;
 END;

 PROCEDURE InitGraph;
 BEGIN
   Fenetre(-1.4,1.4,-1.4,1.4);
   Cloture(42,MaxX-25,95,MaxY-10)
 END;


 {main program}
 BEGIN
   Randomize;
   REPEAT
     compte:=0;
     WinCrtInit('HENON');
     New(Pinfo,Init);
     Data;
     InitGraph;
     Clrscr;
     TextOut(CrtDC,60,15,'HENON''S ATTRACTORS',18);
     Grille(CrtDC,0.2,0.2);
     Bordure(CrtDC);
     Axes(CrtDC);
     Gradue(CrtDC,0.4,0.4);
     Str(angle:5:2,s5);
     StrPCopy(s,'Angle ='+s5);
     TextOut(CrtDC,MaxX-125,MaxY-120,s,strlen(s));
     REPEAT
       xcur:=round(random*MaxX);
       ycur:=round(random*MaxY);
       xpos:= xcur/xrapport+xgfen;
       ypos:= (ycur/yrapport)+ybfen;
       Henon1(CrtDc,xpos,ypos);
       inc(compte);
       IF Keypressed THEN GOTO 100;
     UNTIL compte=300;
100: SortieGraphique;
     if rep='i' then  {send to printer}
     begin
       with Pinfo^ do
       begin
	 MaxX:=3200; MaxY:=2200;
         StartDoc('HENON');
         Cloture(200,MaxX-15,45,MaxY-10);
         TextOut(PrintDC,250,100,'HENON''S ATTRACTORS',18);
         Grille(PrintDC,0.2,0.2);
         Bordure(PrintDC);
         Axes(PrintDC);
         Gradue(PrintDC,0.4,0.4);
         Str(angle:5:2,s5);
         StrPCopy(s,'Angle ='+s5);
	 TextOut(PrintDC,MaxX-500,MaxY-200,s,strlen(s));
	 compte:=0;
         REPEAT
           xcur:=round(random*MaxX);
           ycur:=round(random*MaxY);
           xpos:= xcur/xrapport+xgfen;
           ypos:= (ycur/yrapport)+ybfen;
           Henon1(PrintDc,xpos,ypos);
           inc(compte);
	 UNTIL compte=300;
         NewFrame;
	 EndDoc;
         MaxX:=600; MaxY:=425;
       end;
       SortieGraphique
     end;
   UNTIL rep='n';
   Dispose(Pinfo,Done);
   DoneWinCrt
 END.

{End of file henon.pas}