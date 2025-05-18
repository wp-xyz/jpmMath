{****************************************************
*            REDUCTION  OF  CONICALS                *
*     of equation: ax^2+2bxy+cy^2+2dx+2ey+f=0       *
* ------------------------------------------------- *
* The program, knowing the coefficients a, b, c, d, *
* e and f of the cartesian equation, finds the      *
* caracteristics of the reduction elements:         *
* type of conical  (hyperbola,parabola,ellipse or   *
* circle) and center position,  focus position(s),  *
* axes and excentricity.                            *
* ------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"  *
* [BIBLI 05].                                       *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
*                                                   *
*  REDUCTION OF CONICALS                            *
*                                                   *
*    a  = 1                                         *
*    2b = 2                                         *
*    c  = 1                                         *
*    2d = -13                                       *
*    2e = -11                                       *
*    f  = 32                                        *
*                                                   *
*  Type: Parabola                                   *
*                                                   *
*  Center:  x=0.999999  y=5.000000                  *
*  Symmetry direction:  x=1  y=-1                   *
*  Focus:   x=1.124999  y=4.875000                  *
*  Parameter:           0.35355339                  *
*                                                   *       
****************************************************}
Program Conicals;
Uses WinCrt;

LABEL continue;

VAR
     {input variables}
     a,b,c,d,e,f: DOUBLE;   {coefficients of cartesian equation}

     {output variables:
        typ: type of conical (integer)
             1: ellipse
             2: hyperbola
             3: parabola
             4: circle
             5: line
             6: two lines
             7: one point
             8: no conical at all.
        ex: excentricity (or parameter for a parabola),
        xc,yc: coordinates of center,
        la,lb: axis half-length for an ellipse (la, radius for a circle),
        xf1,yf1,xf2,yf2: focus coordinates (only one focus for a parabola),
        xs1,ys1,xs2,ys2: summit coordinates for an hyperbola,
        xv1,yv1,xv2,yv2: vector directions of symmetry axes
                         (only one direction for a parabola). 
     }
     ex,xc,yc,la,lb,xf1,yf1,xf2,yf2,xs1,ys1,xs2,ys2,xv1,yv1,xv2,yv2: DOUBLE;

     {other internal variables}
     delta,u,l,m,x2,y2: DOUBLE;
     i,j,typ: INTEGER;


  Procedure Parabola;
  Var x,y: DOUBLE;
  Begin
    typ:=3;
    if (c=0) and (a=0) then {the parabola does not exist or is degenerated into a line}
    begin
      if (d=0) and (e=0) then typ:=8 {no conical}else typ:=5{line};
      exit
    end;
    if a=0 then
    begin
      x2:=1; y2:=0; l:=0; m:=1;
    end
    else
    begin
      if a<0 then
      begin
        f:=-f; e:=-e; d:=-d;
        c:=-c; b:=-b; a:=-a;
      end;
      l:=sqrt(a); m:=sqrt(c);
      if b<0 then m:=-m;
      u:=sqrt(a+c);
      x2:=m/u; y2:=-l/u;
      f:=f*u;c:=(a+c)*u;
      u:=d*m-e*l; e:=d*l+e*m; d:=u;
    end;
    if d=0 then
    begin
      if e*e<c*f then typ:=8 else typ:=6;{two lines}
      exit
    end
    else
    begin
      x:=(e*e-c*f)/2/c/d; y:=-e/c;
      xc:=x*x2-y*y2; yc:=y*x2+x*y2;
      ex:=-d/c;
      xf1:=xc+ex*x2/2; yf1:=yc+ex*y2/2;
      xv1:=m; yv1:=-l;
    end
  End;  {Parabola}
                   
  Procedure Hyperbola;
  Begin
    typ:=2;
    la:=sqrt(-abs(f)/l); lb:=sqrt(abs(f)/m);
    if f<0 then
    begin
      u:=la; la:=lb; lb:=u;
      u:=x2; x2:=-y2; y2:=u;
      u:=xv1; xv1:=-yv1; yv1:=u;
    end;
    xv2:=-yv1; yv2:=xv1;
    u:=sqrt(la*la+lb*lb);
    xs1:=xc+x2*la; ys1:=yc+y2*la;
    xs2:=xc-x2*la; ys2:=yc-y2*la;
    xf1:=xc+x2*u; yf1:=yc+y2*u;
    xf2:=xc-x2*u; yf2:=yc-y2*u;
    ex:=u/la;
  End;

  Procedure Ellipse;
  Begin
    typ:=1;
    if f*l>0 then begin typ:=8;{no conical}exit end;
    la:=sqrt(-f/l); lb:=sqrt(-f/m);
    if l<0 then
    begin
      u:=la; la:=lb; lb:=u;
      u:=x2; x2:=-y2; y2:=u;
      u:=xv1; xv1:=-yv1; yv1:=u;
    end;
    xv2:=-yv1; yv2:=xv1;
    u:=sqrt(la*la-lb*lb);
    xf1:=xc+x2*u; yf1:=yc+y2*u;
    xf2:=xc-x2*u; yf2:=yc-y2*u;
    ex:=u/la;
  End;


{main program}
BEGIN

  writeln;
  writeln('  REDUCTION OF CONICALS');
  writeln;
  write('    a  = '); readln(a);
  write('    2b = '); readln(b);
  write('    c  = '); readln(c);
  write('    2d = '); readln(d);
  write('    2e = '); readln(e);
  write('    f  = '); readln(f);
  writeln;
  b:=b/2; d:=d/2; e:=e/2;
  delta:=a*c-b*b;
  if delta=0 then Parabola
  else
  begin
    xc:=(b*e-d*c)/delta;
    yc:=(b*d-a*e)/delta;
    f := f + a*xc*xc+2*b*xc*yc+c*yc*yc+2*d*xc+2*e*yc;
    if f=0 then
    begin
      if delta>0 then typ:=7{one point}else typ:=6;{two lines}
      goto continue
    end;
    u:=sqrt(sqr(a-c)+4*b*b);
    l:=(a+c-u)/2; m:=(a+c+u)/2;
    if (a=c) and (b=0) then
      if f*a>=0 then
      begin
        typ:=8;{no conical}
        goto continue
      end
      else
      begin
        typ:=4; {circle}
        la:=sqrt(-f/a);
        ex:=1;
        goto continue;
      end;
    if (a<c) and (b=0) then
    begin
      x2:=1; y2:=0; xv1:=1; yv1:=0
    end
    else
    begin
      xv1:=b; yv1:=1-a;
      u:=sqrt(xv1*xv1+yv1*yv1);
      x2:=xv1/u; y2:=yv1/u;
    end;
    if delta<0 then Hyperbola else Ellipse
  end; {else if delta=0}

  {print results}
continue: case typ of
    1:begin
        writeln('  Type: Ellipse');
        writeln;
        writeln('  Center                :  x=',xc:9:6,'  y=',yc:9:6);
        writeln('  Direction big axis    :  x=',xv1:5:2,'  y=',yv1:5:2);
        writeln('  Direction small axis  :  x=',xv2:5:2,'  y=',yv2:5:2);
        writeln('  Half length big axis  : ',la:9:6);
        writeln('  Half length small axis: ',lb:9:6);
        writeln('  First focus           :  x=',xf1:9:6,'  y=',yf1:9:6);
        writeln('  Second focus          :  x=',xf2:9:6,'  y=',yf2:9:6);
        writeln('  Excentricity          :  ',ex:10:8);
      end;
    2:begin
        if a+c=0 then writeln('  Type: Equilateral Hyperbola')
                 else writeln('  Type: Hyperbola');
        writeln;
        writeln('  Center                :  x=',xc:9:6,'  y=',yc:9:6);
        writeln('  Direction first axis  :  x=',xv1:5:2,'  y=',yv1:5:2);
        writeln('  Direction second axis :  x=',xv2:5:2,'  y=',yv2:5:2);
        writeln('  First summit          :  x=',xs1:9:6,'  y=',ys1:9:6);
        writeln('  Second summit         :  x=',xs2:9:6,'  y=',ys2:9:6);
        writeln('  First focus           :  x=',xf1:9:6,'  y=',yf1:9:6);
        writeln('  Second focus          :  x=',xf2:9:6,'  y=',yf2:9:6);
        writeln('  Excentricity          :  ',ex:10:8);
      end;
    3:begin
        writeln('  Type: Parabola');
        writeln;
        writeln('  Center:  x=',xc:9:6,'  y=',yc:9:6);
        writeln('  Symmetry direction:   x=',xv1:5:2,'  y=',yv1:5:2);
        writeln('  Focus:   x=',xf1:9:6,'  y=',yf1:9:6);
        writeln('  Parameter:             ',ex:10:8);
      end;
    4:begin
        writeln('  Type: Circle');
        writeln;
        writeln('  Center         :  x=',xc:9:6,'  y=',yc:9:6);
        writeln('  Radius         : ',la:9:6);
        writeln('  Excentricity   :  ',ex:1:0);
      end;
    5:writeln('  Conical degenerated into a line.');
    6:writeln('  Conical degenerated into two lines.');
    7:writeln('  Conical degenerated into a single point.');
    8:writeln('  Not a conical !');

  end;

  writeln;
  ReadKey; DoneWinCrt

END.

{end of file conical.pas}