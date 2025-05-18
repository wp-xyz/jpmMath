{**************************************************************
*  This program draws the Euler's circle of a triangle ABC    *
*  with automatic scales.                                     *
*                                                             *
*                                      By J-P Moreau, Paris.  *
**************************************************************}
Program Circle_Euler;

Uses WinCrt1, Type_def, Graph_2d;

Var
    ax,ay,bx,by,cx,cy: real_ar;
    a1x,a1y,b1x,b1y,c1x,c1y: real_ar;
    a,b,x0,y0,x1,y1,x2,y2: real_ar;


{Solve linear system ax+b=c and a1x+b1y=c1}
Procedure System2D(a,b,c,a1,b1,c1:real_ar; Var x,y:real_ar);
Begin
  if Abs(b*a1-a*b1)<1E-12 then
    writeln(' * Error in System2D - System is singular!')
  else
  begin
    y := (c*a1-c1*a)/(b*a1-a*b1);
    if a<>0.0 then x := (c-b*y)/a
              else y := (c1-b1*y)/a1
  end             
End;

{determine coefficients a, b of y=ax+b straight line passing
 through 2 points (x1,y1) and (x2,y2) }
Procedure Line2D(x1,y1,x2,y2:real_ar; Var a,b:real_ar);
Begin
  System2D(x1,1.0,y1,x2,1.0,y2,a,b)
End;

{Draw line y=ax+b from (x1,y1) to (x2,y2) }
Procedure Draw_Line(x1,x2,a,b:real_ar);
Begin
  MoveXY(CrtDc,x1,a*x1+b); LineXY(CrtDc,x2,a*x2+b);
End;

{Draw a circle passing through 3 points (x1,y1),(x2,y2),(x3,y3) }
Procedure Draw_Circle3P(x1,y1,x2,y2,x3,y3:Real_ar);
Var c,r: real_ar; {coefficients of circle equation
                   x²+y²-2ax-2by+c=0 and radius r}
Begin
  System2D(2.0*(x2-x1),2.0*(y2-y1),(x2*x2+y2*y2-x1*x1-y1*y1),
           2.0*(x3-x1),2.0*(y3-y1),(x3*x3+y3*y3-x1*x1-y1*y1),a,b);
  c:=2.0*a*x1+2.0*b*y1-x1*x1-y1*y1;
  r:=sqrt(a*a+b*b-c);

  Circle(CrtDc,a,b,r,True);
  CroixXY(CrtDc,a,b);

End;


{main}
BEGIN
  
  {define 3 points in plane (Ox,Oy) }
  ax:=0.5; ay:=1.6;
  bx:=0.0; by:=0.0;
  cx:=1.5; cy:=0.0;

  WinCrtInit(' Circle of Euler');
  MaxX:=MaxY;
  InitFenetre(CrtDc,10,-0.4,2.0,-0.4,2.0);

  TextXY(CrtDc,ax-0.05,ay+0.15,'A');
  TextXY(CrtDc,bx-0.05,by-0.05,'B');
  TextXY(CrtDc,cx+0.025,cy-0.025,'C');
  TextXY(CrtDc,bx+0.55,by+0.25,'H');
  TextXY(CrtDc,bx+0.65,by+0.4,'O1');
  TextXY(CrtDc,bx+0.75,by+0.6,'O');
  TextXY(CrtDc,bx+0.5,by-0.05,'HA');
  TextXY(CrtDc,bx+1.1,by+0.8,'HB');
  TextXY(CrtDc,bx+0.2,by+0.55,'HC');

  MoveXY(CrtDc,ax,ay); LineXY(CrtDc,bx,by);
  LineXY(CrtDc,cx,cy); LineXY(CrtDc,ax,ay);

  Draw_Circle3P(ax,ay,bx,by,cx,cy);
  TextXY(CrtDc,ax+0.05,ay+0.3,'Circumscribed Circle');
  x1:=a; y1:=b;

  {define 3 points of Euler's circle}
  a1x:=(ax+bx)/2.0; a1y:=(ay+by)/2.0;
  b1x:=(cx+bx)/2.0; b1y:=(cy+by)/2.0;
  c1x:=(ax+cx)/2.0; c1y:=(ay+cy)/2.0;

  {draw the 3 bisectors of triangle ABC}
  MoveXY(CrtDc,ax,ay); LineXY(CrtDc,b1x,b1y);
  MoveXY(CrtDc,bx,by); LineXY(CrtDc,c1x,c1y);
  MoveXY(CrtDc,cx,cy); LineXY(CrtDc,a1x,a1y);

  {drax euler circle of triangle ABC}
  Draw_Circle3P(a1x,a1y,b1x,b1y,c1x,c1y);
  TextXY(CrtDc,0.3,0.82,'Euler''s Circle');
  x2:=a; y2:=b;

  Line2D(x1,y1,x2,y2,a,b);
  Draw_Line(ax,0.9*cx,a,b);
  TextXY(CrtDc,ax+0.39,ay-0.3,'Euler''s Line');
  
  Legendes(CrtDc,'Euler''s Circle of triangle ABC','X','Y');

  {Draw line C-Hc perpendicular to side AB}
  Line2D(ax,ay,bx,by,a,b);
  System2D(a,-1,-b,(-1.0/a),-1.0,-(cy+cx/a),x0,y0); {Hc is in (x0,y0) }
  MoveXY(CrtDc,cx,cy); LineXY(CrtDc,x0,y0);

  {Draw line B-Hb perpendicular to side AC}
  Line2D(ax,ay,cx,cy,a,b);
  System2D(a,-1,-b,(-1.0/a),-1.0,-(by+bx/a),x0,y0); {Hb is in (x0,y0) }
  MoveXY(CrtDc,bx,by); LineXY(CrtDc,x0,y0);

  {Draw line A-Ha perpendicular to side BC
  Line2D(bx,by,cx,cy,a,b);
  System2D(a,-1,-b,(-1.0/a),-1.0,-(ay+ax/a),x0,y0); {Hb is in (x0,y0)
  MoveXY(CrtDc,ax,ay); LineXY(CrtDc,x0,y0); }
  MoveXY(CrtDc,ax,ay); LineXY(CrtDc,ax,0.0);  {special case a=0}

  SortieGraphique;
  DoneWinCrt

END.
 