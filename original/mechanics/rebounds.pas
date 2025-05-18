{**************************************************************************
*                      THE   BOUNCING  BALL                               *
* ----------------------------------------------------------------------- *
* Explanations:                                                           *
*                                                                         *
* Let us consider an elastic ball dropped with a horizontal speed Vx from *
* a height H.  Here are the main equations to calculate the successive    *
* rebounds on the ground:                                                 *
*                                                                         *
* Initializing:                                                           *
* The first part of the trajectory is a parabolic section, upwards orien- *
* ted. From the physical laws for a falling body, we know that:           *
*      Vy = sqrt(2GH)  and V0 = sqrt(Vx² + Vy²) with theta=arctan(Vx/Vy), *
* where V0 is the speed at ground level and theta = pi - angle of V0 with *
* the Ox axis.                                                            *
* The ball moves with the horizontal speed Vx and the vertical speed -Vy  *
* during time t=Vy/G until it reaches the ground (G=gravity acceleration).*
* Hence the horizontal distance run is x2=Vx*t=Vx*Vy/G.                   *
* The equation of the parabolic section is: y = -(H/x2²)*x + H from x=0   *
* to x=x2.                                                                *
*                                                                         *
* Rebounds:                                                               *
* For each rebound, the initial speed V0 is multiplied by a damping       *
* coefficient, called amort (<1), but the theta initial angle is kept the *
* same. The previous x2 value becomes the next x1 value. The parabolic    *
* equation of the next trajectory section has the complete form:          *
*                         y = a*x² + b*x + c  (1)                         *
* P=2*Vx*Vy/G is the correct length for all the rebounds. Hence the next  *
* x2 = x1 + P. The parabola summit coordinates are the following:         *
*              xs = (x1+x2)/2  and  ys = (V0²*sin²(theta))/(2G)           *
* So to draw the trajectory section of a rebound, we have to know the     *
* equation of a parabola passing through the 3 points: (x1,0), (x2,0) and *
* (xs,ys). This leads to the following linear system to solve of order 3: *
*              a*x1² + b*x1 + c = 0                                       *
*              a*xs² + b*xs + c = 0                                       *
*              a*x2² + b*x2 + c = 0                                       *
* where the unknowns are a, b, c, the coefficients of the parabola (1).   *
* Using the Cramer's method, we calculate the main determinant Dp and the *
* three sub-determinants, Da, Db, Dc. We find:                            *
*          Dp = x1²*(xs-x2) + x2²*(x1-xs) + xs²*(x2-x1)                   *
*          Da = ys*(x2-x1)  Db = ys*(x1²-x2²)  Dc = ys*x1*x2*(x2-x1)      *
* Hence:   a = Da/Dp   b = Db/Dp   c = Dc/Dp                              *
* Now we can draw the parabolic section of the next rebound in [x1,x2].   *
* We continue the same process until the required number of rebounds is   *
* reached.                                                                *
* ----------------------------------------------------------------------- *
* Reference:                                                              *
* From "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0 de  *
* R. Dony - MASSON 1990 page 113" [BIBLI 12].                             *
*                                                                         *
*                                     TPW Release By J-P Moreau, Paris.   *
*                                            (www.jpmoreau.fr)            *
**************************************************************************}
Program Rebounds;

Uses WinCrtMy, WinProcs, Strings, Type_def, CrtGr2D;

Const Height = 12;
      G = 9.81;
      step = 0.005;

Var   amort,vx,vy,v0,a,b,c: real;
      x1,x2,xs,ys,angle,Length: real;
      bounds,nrebounds: byte;

  Procedure Data;
  begin
    clrscr;
    writeln;
    writeln('  The bouncing ball');
    writeln('  -----------------');
    writeln;
    write('  Initial horizontal speed: '); readln(vx);
    writeln;
    write('  Damping coefficient.....: '); readln(amort);
    writeln;
    write('  Number of rebounds......: '); readln(nrebounds)
  end;             

  Procedure Init;
  begin
    bounds:=0;
    x1:=0.0;
    vy:=sqrt(2*G*Height);
    angle:=ArcTan(vy/vx);
    v0:=sqrt(sqr(vx)+sqr(vy));
    x2:=vx*vy/G;
    a:=-Height/sqr(x2);
    b:=0.0;
    c:=Height
  end;

  Procedure DrawParabol(pt1,pt2:real);
  var x,y: real;
  begin
    x:=pt1;
    MoveXY(CrtDc,x,0.0);
    while x<pt2 do
    begin
      y:=a*sqr(x)+b*x+c;
      LineXY(CrtDc,x,y);
      x:=x+step
    end;
    Inc(bounds)
  end;

  Procedure CoeffParabol;
  var detprinc,deta,detb,detc: real;
  begin
    ys:=sqr(v0*sin(angle))/(2*G);
    vy:=sqrt(2*G*ys);
    Length:=2*vx*vy/G;
    x2:=x1+Length;
    xs:=(x1+x2)/2.0;
    detprinc:=sqr(x1)*(xs-x2)+sqr(x2)*(x1-xs)+sqr(xs)*(x2-x1);
    deta:=ys*(x2-x1);
    a:=deta/detprinc; {1st coefficient of parabola}
    detb:=ys*(sqr(x1)-sqr(x2));
    b:=detb/detprinc; {2nd coefficient}
    detc:=ys*x1*x2*(x2-x1);
    c:=detc/detprinc {3rd coefficient}
  end;

  Procedure PrintData;
  var s1,s2,s3,s4,ch: String;
      s:array[0..79] of char;
  begin
    Str(Height:3,s1);
    Str(vx:4:1,s2);
    Str(x2:5:2,s3);
    Str(amort:4:2,s4);
    ch:=Concat('Height=',s1,'  Speed=',s2,'  Total length=',s3,'  Damping=',s4);
    StrPCopy(s,ch);
    TextOut(CrtDc,75,15,s,strlen(s))
  end;   

{main program}
BEGIN

  WinCrtInit(' The bouncing Ball');
  Repeat
    Data;
    Init;
    ClrScr;
    Fenetre(0.0,20,0.0,Height);
    Cloture(50,MaxX-10,100,MaxY-40);
    Axes(CrtDc);
    Grille(CrtDc,1.0,1.0);
    Gradue(CrtDc,2.0,2.0);
    Bordure(CrtDc);
    DrawParabol(x1,x2);
    x1:=x2;
    Repeat
      v0:=v0*amort;
      CoeffParabol;
      DrawParabol(x1,x2);
      x1:=x2
    Until bounds=nrebounds;
    PrintData;
    SortieGraphique
  Until rep='n';
  DoneWinCrt

END.

{end of file rebounds.pas}