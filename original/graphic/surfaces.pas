{******************************************************
*                  3D  SURFACES                       *
* --------------------------------------------------- *
* This program draws 3D surfaces defined by equations *
* Z = F(X,Y) by using unit graph_3d.pas.              *
* --------------------------------------------------- *
* SAMPLE RUN:                                         *
* ( Draw the surface defined by:                      *
*   Z = 8 sin(sqrt(X*X+Y*Y))/sqrt(X*X+Y*Y) )          *
*                                                     *
*  Input intervals [x1,x2] and [y1,y2]:               *
*    X1 X2 = -10 10                                   *
*    Y1 Y2 = -10 10                                   *
*  Number of lines ..........: 100                    *
*  Number of points per line : 100                    *
*  Real or uniform view (r/u): u                      *
*                                                     *
*  Choice of projection type:                         *
*  =========================                          *
*    1. real perspective                              *
*    2. ordinary parallel                             *
*    3. dimetric parallel                             *
*    4. isometric parallel                            *
*                                                     *
*  Your choice (1 to 4): 4                            *
* --------------------------------------------------- *
* Ref.: "Graphisme dans le plan et dans l'espace en   *
*        Turbo Pascal 4.0 By R. Dony, MASSON Paris,   *
*        1990" [BIBLI 12].                            *
*                                                     *
*                   TPW Release By J-P Moreau, Paris  *
*                          (www.jpmoreau.fr)          *
******************************************************}
Program SurfacesZ;
Uses WinCrtMy,WinProcs,Strings,Type_def,Graph_3D;

Const  LimX = 625;             {dimensions in pixels of graphic}
       LimY = 400;             {window open with WinCrtInit.   }

Type   Table = Array[0..LimX] of integer;

Var    Hmax,Hmin     : Table;
       x1,x2,y1,y2   : real_ar;
       incx,incy     : real_ar;
       f1,f2,f3,f4   : real_ar;
       c1,c2,c3,c4   : integer;
       echx,echy,ech : real_ar;
       xg,yg,xd,yd   : integer;
       nbrelignes,
       nbrepoints    : integer;
       vue           : char;

 {Function F(x,y:real_ar) : real_ar;
  begin
    F:=x*x+y*y
  end;

  Function F(x,y:real_ar) : real_ar;
  begin
    F:=5.0*sin(x)*sin(y)
  end;

  Function F(x,y:real_ar) : real_ar;
  var k:real_ar;
  begin
    k:=sqrt(x*x+y*y);
    F:=-k+abs(sin(k))
  end; }

  Function F(x,y:real_ar) : real_ar;
  var k:real_ar;
  begin
    k:=sqrt(x*x+y*y);
    if abs(k)>1e-10 then
      F:=8.0*sin(k)/k
    else
      F:=8.0;
  end;

  Function Signe(x: real_ar): integer;
  begin
    if x > 0 then signe:=1
             else if x < 0 then signe:=-1
             else signe:=0
  end;

  Procedure Data;
  var choix : byte;  ch : String;
  begin
    clrscr;
    writeln;
    writeln('  Input intervals [x1,x2] and [y1,y2]:');
    write('    X1 X2 = '); readln(x1,x2); 
    write('    Y1 Y2 = '); readln(y1,y2); 
    write('  Number of lines ..........: '); readln(nbrelignes);
    write('  Number of points per line : '); readln(nbrepoints);
    write('  Real or uniform view (r/u): '); readln(vue);
   {x1:=-10; x2:=10; y1:=-10; y2:=10;
    nbrelignes:=100; nbrepoints:=100; vue:='u'; }
    writeln;
    writeln('  Choice of projection type:');
    writeln('  =========================  ');
    writeln;
    writeln('  1. real perspective');
    writeln('  2. ordinary parallel');
    writeln('  3. dimetric parallel');
    writeln('  4. isometric parallel');
    writeln;
    write('   Your choice (1 to 4): '); readln(choix);
    writeln;
    DE:=1;
    Projection:=Parallele;
    Case choix of
      1: begin
           write('  Distance Rho : '); readln(rho);
           write('  Angle Theta  : '); readln(theta);
           write('  Angle Phi    : '); readln(phi);
           Projection:=Perspective
         end;
      2: begin
           write('  Angle Theta  : '); readln(theta);
           write('  Angle Phi    : '); readln(phi)
         end;
      3: begin
           theta:=22.20765;
           phi:=20.704811
         end;
      4: begin
           theta:=45;
           phi:=35.26439
         end
    end;
    writeln;
    write('  Calculating surface points...')
  end;   { Data }

  Procedure Init;
  var  aux: real_ar; i: integer;
  begin
    incx:=(x2-x1)/nbrepoints;
    incy:=(y2-y1)/nbrelignes;
    c1:=0; c2:=LimX-1; c3:=1; c4:=LimY;
    f1:=1E10; f2:=-f1; f3:=f1; f4:=-f1;
    xg:=-1; yg:=-1; xd:=-1; yd:=-1;
    FillChar(Hmax,SizeOf(Hmax),0);
    for i:=0 to LimX do Hmin[i]:=LimY;
    if (theta < 0.0) OR (theta > 180.0) then
    begin
      aux:=x1; x1:=x2; x2:=aux; incx:=-incx;
      aux:=y1; y1:=y2; y2:=aux; incy:=-incy
    end
  end;  { Init }

  Procedure Fenetre;
  var x,y,z : real_ar;
      ligne,point : integer;
  begin
    for ligne:=0 to nbrelignes do
    begin
      y:=y2-ligne*incy;
      for point:=0 to nbrepoints do
      begin
        x:=x1+point*incx;
        z:=F(x,y);
        Project(x,y,z);
        if xproj < f1 then f1:=xproj;
        if xproj > f2 then f2:=xproj;
        if yproj < f3 then f3:=yproj;
        if yproj > f4 then f4:=yproj
      end
    end
  end;   { Fenetre }

  Procedure Echelles;
  begin
    echx:=(c2-c1)/(f2-f1);
    echy:=(c4-c3)/(f4-f3);
    if Upcase(vue)='R' then
      if echx < echy then echy:=echx else echx:=echy
  end;

  Procedure Horizon(x1,y1,x2,y2: integer);
  var  x,y,dx : integer;
       pente  : real_ar;

       function Max(x1,x2:integer): integer;
       begin
         if x1>x2 then Max:=x1 else Max:=x2
       end;

       function Min(x1,x2:integer): integer;
       begin
         if x1<x2 then Min:=x1 else Min:=x2
       end;

  begin
    dx:=Signe(x2-x1);
    if dx=0 then
    begin
      Hmax[x2+1]:=Max(Hmax[x2],y2);
      Hmin[x2+1]:=Min(Hmin[x2],y2)
    end
    else
    begin
      pente:=(y2-y1)/(x2-x1);
      for x:=x2+1 to x1 do
      begin
        y:=Round(pente*(x-x1)+y1);
        Hmax[x]:=Max(Hmax[x],y);
        Hmin[x]:=Min(Hmin[x],y)
      end
    end
  end;  { Horizon }

  Procedure Visibilite (x,y:integer; var visi: integer);
  begin
    if ((y < Hmax[x]) AND (y > Hmin[x])) then
      visi:=0
    else if y >= Hmax[x] then
      visi:=1
    else visi:=-1
  end;

  Procedure Inter1(x1,y1,x2,y2:integer; var Tabaux:Table; var xi,yi:integer);
  var  ct1,ct2,p1,p2,xii,yii : real_ar;
  begin
    if x2-x1 = 0 then
    begin
      xii:=x2;
      yii:=Tabaux[x2]
    end
    else
    begin
      p1:=(y2-y1)/(x2-x1); p2:=(Tabaux[x2]-Tabaux[x1])/(x2-x1);
      if (abs(p1) > 1e-10) and (abs(p1-p2) > 1e-10) then
      begin
        ct1:=y1-p1*x1; ct2:=Tabaux[x1]-p2*x1;
        yii:=(p1*ct2-p2*ct1)/(p1-p2);
        xii:=(yii-ct1)/p1
      end
      else
      begin
        xii:=x2;
        yii:=y2
      end
    end;
    xi:=Round(xii);
    yi:=Round(yii)
  end;  {Inter1}

  Procedure AreteFermeture(x,y:integer;var xlateral,ylateral:integer);
  begin
    if xlateral <> -1 then Horizon(xlateral,ylateral,x,y);
    xlateral:=x;
    ylateral:=y
  end;

  {Draw a surface Z=F(X,Y) removing hidden sections}
  Procedure DessinFonction;
  var  xe,ye,ligne,point,xi,yi : integer;
       xprec,yprec,xcour,ycour : integer;
       visicour,visiprec       : integer;
       x,y,z                   : real_ar;
  begin
    for ligne:=0 to nbrelignes do
    begin
      y:=y2-ligne*incy;
      x:=x1;
      z:=F(x,y);
      Project(x,y,z);
      xprec:=Round((xproj-f1)*echx)+c1;
      yprec:=Round((yproj-f3)*echy)+c3;
      AreteFermeture(xprec,yprec,xd,yd);
      MoveTo(CrtDc,xprec,LimY-yprec);
      Visibilite(xprec,yprec,visiprec);
      for point:=0 to nbrepoints do
      begin
        x:=x1+point*incx;
        z:=F(x,y);
        Project(x,y,z);
        xcour:=Round((xproj-f1)*echx)+c1;
        ycour:=Round((yproj-f3)*echy)+c3;
        Visibilite(xcour,ycour,visicour);
        if (Hmax[xcour]=0) OR (Hmin[xcour]=LimY) then
          visicour:=visiprec;
        if visicour=visiprec then
        begin
          if (visicour=1) OR (visicour=-1) then
          begin
            MoveTo(CrtDc,xprec,LimY-yprec);
            LineTo(CrtDc,xcour,LimY-ycour);
            Horizon(xprec,yprec,xcour,ycour)
          end
        end
        else
        begin
          if visicour=0 then
          begin
            if  visiprec=1 then
              Inter1(xprec,yprec,xcour,ycour,Hmax,xi,yi)
            else
              Inter1(xprec,yprec,xcour,ycour,Hmin,xi,yi);
            MoveTo(CrtDc,xprec,LimY-yprec);
            LineTo(CrtDc,xi,LimY-yi);
            Horizon(xprec,yprec,xi,yi)
          end
          else
          begin
            if visicour=1 then
            begin
              if visiprec=0 then
              begin
                Inter1(xprec,yprec,xcour,ycour,Hmax,xi,yi);
                MoveTo(CrtDc,xi,LimY-yi);
                LineTo(CrtDc,xcour,LimY-ycour);
                Horizon(xi,yi,xcour,ycour)
              end
              else
              begin
                Inter1(xprec,yprec,xcour,ycour,Hmin,xi,yi);
                MoveTo(CrtDc,xprec,LimY-yprec);
                LineTo(CrtDc,xi,LimY-yi);
                Horizon(xprec,yprec,xi,yi);
                Inter1(xprec,yprec,xcour,ycour,Hmax,xi,yi);
                MoveTo(CrtDc,xi,LimY-yi);
                LineTo(CrtDc,xcour,LimY-ycour);
                Horizon(xi,yi,xcour,ycour)
              end
            end
            else
            begin
              if visiprec=0 then
              begin
                Inter1(xprec,yprec,xcour,ycour,Hmin,xi,yi);
                MoveTo(CrtDc,xi,LimY-yi);
                LineTo(CrtDc,xcour,LimY-ycour);
                Horizon(xi,yi,xcour,ycour)
              end
              else
              begin
                Inter1(xprec,yprec,xcour,ycour,Hmax,xi,yi);
                MoveTo(CrtDc,xprec,LimY-yprec);
                LineTo(CrtDc,xi,LimY-yi);
                Horizon(xprec,yprec,xi,yi);
                Inter1(xprec,yprec,xcour,ycour,Hmin,xi,yi);
                MoveTo(CrtDc,xi,LimY-yi);
                LineTo(CrtDc,xcour,LimY-ycour);
                Horizon(xi,yi,xcour,ycour)
              end
            end
          end  { if visicour=0 }
        end;
        visiprec:=visicour;
        xprec:=xcour;
        yprec:=ycour
      end;  { for point }
      AreteFermeture(xcour,ycour,xg,yg)
    end   { for ligne }
  end;  { DessinFonction }

  Procedure Affiche;
  var  aux,ch : string;
       s : array[0..60] of char;
  begin
    str(x1:4:1,aux); ch:=Concat('X=[',aux,',');
    str(x2:4:1,aux); ch:=Concat(ch,aux,']   ');
    str(y1:4:1,aux); ch:=Concat(ch,'Y=[',aux,',');
    str(y2:4:1,aux); ch:=Concat(ch,aux,']   ');
    if projection=Perspective then
    begin
      str(rho:4:1,aux); ch:=Concat(ch,'Rho=',aux,'  ')
    end;
    str(theta:5:3,aux); ch:=Concat(ch,'Theta=',aux,'  ');
    str(phi:5:3,aux); ch:=Concat(ch,'Phi=',aux,'  ');
    StrPCopy(s,ch);
    TextOut(CrtDc,20,5,s,strlen(s));
    str(nbrelignes:2,aux); ch:=Concat('( ',aux,' lines','  ');
    str(nbrepoints:3,aux); ch:=Concat(ch,aux,' pts per line )');
    StrPCopy(s,ch);
    TextOut(CrtDc,20,20,s,strlen(s))
  end;

  {main program}
  BEGIN
    WinCrtInit('Surfaces Z=F(X,Y)');
    Data;
    Init;
    InitProj;
    Fenetre;
    Echelles;
    clrscr;
    DessinFonction;
    Affiche;
    Pause;
    DoneWinCrt
  END.

{end of file surfaces.pas}