{*********************************************************
* Program to demonstrate Apollonius circles and clipping *
* capability of unit CrtGr2D.pas                         *
*                                                        *
*                      TPW version by J-P Moreau, Paris  *
*                              (www.jpmoreau.fr)         *
* ------------------------------------------------------ *
* DESCRIPTION:  At each step, the program:               *
*   1. randomly defines 3 circles so that they are at    *
*      least partially visible and each one is tangent   *
*      to the other two.                                 *
*   2. recursively draws all the internal (apollonius)   *
*      circles until radius is too small, using a fast   *
*      procedure to draw circles.                        *
*   The number of steps is limited to 20 and the average *
*   drawing time of a step is constantly displayed.      *
* ------------------------------------------------------ *
* INSTRUCTIONS:  Press space bar to stop program after a *
*                step is over.                           *
*                                                        *
* NOTE:    The print option is not implemented here      *
*********************************************************}
Program Apollo;
Uses WinCrtMy,WinProcs,Maths_2D,Graph_2,Type_def,CrtGr2D,Time;

Type
      {object circle defined by center and 1/radius}
      cercle = Record
                 centre: point2D;
                 cb    : real_ar;
               End;

Var
      X_min,X_max,Y_min,Y_max: real_ar;
      C : Array [1..3] of cercle;
      i : 1..3;
      moyenne : real_ar;
      compte,err,j,k : integer;

  {returns true if circle is not entirely visible}
  Function horsLimites(C: cercle): boolean;
  Var r: real_ar;
  Begin
    With C,centre do
    begin
      r:=1/cb;
      horsLimites:=true;
      if (x<X_min) and (r<X_min-x) then exit;
      if (y<Y_min) and (r<Y_min-y) then exit;
      if (x>X_max) and (r<x-X_max) then exit;
      if (y>Y_max) and (r<y-Y_max) then exit;
      horsLimites:=false
    end
  End;

    {returns the internal circle to 3 given circles}
    Procedure CercleInterne(c1,c2,c3: cercle; var Co: cercle);
    Var u,v,d: real_ar;
        V2,V3: Vecteur2D;
    Begin
      With Co do
      begin
        if c1.cb > c2.cb then begin
                                Co:=c1; c1:=c2; c2:=Co
                              end;
        if c2.cb > c3.cb then begin
                                Co:=c2; c2:=c3; c3:=Co
                              end;
        GetVecteur2D(c1.centre,c2.centre,V2);
        GetVecteur2D(c1.centre,c3.centre,V3);
        d:=V2.x*V3.y-V2.y*V3.x;
        Erreur2D:=abs(d)<eps; if erreur2D then exit;
        u:=sqrt(c3.cb)*sqrt(c1.cb+c2.cb+c1.cb*(c2.cb/c3.cb));
        cb:=c1.cb+c2.cb+c3.cb+2*u;
        u:=(1+c2.cb/c1.cb+(c2.cb-c1.cb)/cb)/c1.cb/c2.cb;
        v:=(1+c3.cb/c1.cb+(c3.cb-c1.cb)/cb)/c1.cb/c3.cb;
        centre.x:=c1.centre.x+(u*V3.y-V2.y*v)/d;
        centre.y:=c1.centre.y+(v*V2.x-V3.x*u)/d
      end
    End;

    {draw a given circle using the fast Circle1 procedure
     defined in unit CrtGr2D} 
    Procedure traceCercle(C: cercle);
    VAR  cx,cy,r : real_ar;
    Begin
      with C do
      begin
        cx:=centre.x; cy:=centre.y; r:=1.0/cb;
        Circle1(CrtDC,cx,cy,r,TRUE);
      end
    End;

    {Recursive procedure to draw all the visible internal
     circles until size is too small}
    Procedure TraceInterne (c1,c2,c3: cercle);
    Var h11,h12,h13: boolean;
        Co         : cercle;
    Begin
      if KeyPressed then exit;
      h11:=horsLimites(c1); h12:=horsLimites(c2);
      if h11 and h12 then exit;
      h13:=horsLimites(c3);
      if h13 and (h11 or h12) then exit;
      if c1.cb+c2.cb+c3.cb < 200 then
      begin
        CercleInterne(c1,c2,c3,Co);
        TraceCercle(Co);
	TraceInterne(C1,C2,Co);
        TraceInterne(C1,C3,Co);
        TraceInterne(C2,C3,Co)
      end
    End;

    {**********************************************************
    * INPUTS:                                                 *
    *            Circle c1 defined by center and 1/radius     *
    *            Circle c2 defined by center, radius t.b.d.   *
    *            Circle c3 defined by 1/radius, center t.b.d. *
    * OUTPUTS:                                                *
    *            Circle c2 with 1/radius defined              *
    *            Circle c3 with center defined                *
    *            so that each circle is tangent to the other  *
    *            two.                                         *
    * ------------------------------------------------------- *
    * NOTE: t.b.d. means To be defined.                       *
    **********************************************************}    
    Procedure troisiemeCercle(c1 : cercle; var c2,c3 : cercle);
    Var  a,b,c: real_ar;
         V    : vecteur2D;
    Begin
      c:=1/Dist2D(c1.centre,c2.centre);
      Erreur2D:=(c >= c1.cb);
      If erreur2D then exit;
      a:=c/c1.cb; b:=1-a; c2.cb:=c/b;
      GetVecteur2D(c1.centre,c2.centre,V);
      c:=c/c3.cb; a:=a+c; b:=b+c;
      c:=(sqr(a)-sqr(b)+1)/2.0;
      a:=sqrt(sqr(a)-sqr(c));
      c3.centre:=c1.centre;
      With c3.centre do
      begin
        x:=x+V.x*c+a*V.y;
        y:=y+V.y*c-a*V.x
      end
    End;


{main program}
BEGIN

  Randomize;
  x_min:=-0.60; y_min:=-0.45; x_max:=0.60; y_max:=0.45;
  WinCrtInit('APOLLONIUS CIRCLES');
  Fenetre(x_min,x_max,y_min,y_max);
  Cloture(70,MaxX-30,95,MaxY-5);
  StartTiming;
  rep:='o';

  Repeat
    compte:=0;
    Repeat
      ClrScr;
      Axes(CrtDC);
      Grille(CrtDC,0.1,0.05);
      Gradue(CrtDC,0.2,0.10);
      Bordure(CrtDC);
      Repeat
        Inc(compte);
        gotoxy(62,2); write(compte:2);
        With C[1] do
        begin
	  centre.x := 1.0*Random-1.0; centre.y := 0.5*Random-1.25;
	  cb := 1.0*Random+0.3;
        end;
        With C[2].centre do
        begin
	  x:=1.0*Random-1.0; y:=0.5*Random+0.5;
        end;
        C[3].cb:=1.0*Random+0.3;
        TroisiemeCercle(C[1],C[2],C[3]);
      Until Not Erreur2D;
      For i:=1 to 3 do traceCercle(C[i]);
      TraceInterne(C[1],C[2],C[3]);
      StopTiming;
      Val(Elapsed,moyenne,err);
      moyenne:=moyenne/compte;
      gotoxy(12,2);
      write('Mean time: ',moyenne:6:1,' s.');
      for j:=1 to 30000 do   {waiting loop}
        for k:=1 to 3000 do;
    Until (compte>19) or Keypressed;
    SortieGraphique
  Until rep='n';
  DoneWinCrt

END.

{end of file apollo.pas}