{******************************************************
*          SURFACES OF PARAMETRIC EQUATIONS           *
* --------------------------------------------------- *
* This program draws 3D surfaces defined by equations *
* x = f(u,v), y = g(u,v), z = h(u,v)  by using unit   *
* graph_3d.pas.                                       *
* --------------------------------------------------- *
* SAMPLE RUNS:                                        *
*                                                     *
*       Surface Data:                                 *
*       =============                                 *
*                                                     *
*   Surface Type: 1. Ellipsoid                        *
*                 2. Sphere                           *
*                 3. Torus                            *
*                 4. Hyperboloid  3                   *
*   Projection kind (1:perspective 2:parallel): 2     *
*   Value of distance to screen: 30                   *
*   Value of angle Theta (deg.): 60                   * 
*   Value of angle Phi   (deg.): 30                   *
*                                                     *
*   Input U begin, U end and step: -1.57 1.57 0.2     *
*   Input V begin, V end and step: -3.14 3.14 0.2     *
*   (An ellipsoid is displayed in parralel mode).     *
*                                                     *
*   Your choice(1 to 4): 2                            *
*   Projection kind (1:perspective 2:parallel): 2     *
*   Value of distance to screen: 25                   *
*   Value of angle Theta (deg.): 45                   * 
*   Value of angle Phi   (deg.): 30                   *
*                                                     *
*   Input U begin, U end and step: -3.14 3.14 0.2     *
*   Input V begin, V end and step: -3.14 3.14 0.2     *
*   (A shere is displayed in parralel mode).          *
*                                                     *
*   Your choice(1 to 4): 3                            *
*   Projection kind (1:perspective 2:parallel): 2     *
*   Value of distance to screen: 20                   *
*   Value of angle Theta (deg.): 30                   * 
*   Value of angle Phi   (deg.): 30                   *
*                                                     *
*   Input U begin, U end and step: -3.14 3.14 0.4     *
*   Input V begin, V end and step: -3.14 3.14 0.2     *
*   (A torus is displayed in parralel mode).          *
*                                                     *
*   Your choice(1 to 4): 4                            *
*   Projection type (1:perspective 2:parallel): 2     *
*   Value of distance to screen: 100                  *
*   Value of angle Theta (deg.): 60                   * 
*   Value of angle Phi   (deg.): 30                   *
*                                                     *
*   Input U begin, U end and step: -1  1  0.1         *
*   Input V begin, V end and step: -1  1  0.1         *
*   (An hyperpoloid is displayed in parralel mode).   *
*                                                     *
*   Your choice(1 to 4): 3                            *
*   Projection type (1:perspective 2:parallel): 1     *
*   Value of Rho               : 150                  *
*   Value of distance to screen: 2500                 *
*   Value of angle Theta (deg.): 30                   * 
*   Value of angle Phi   (deg.): 30                   *
*                                                     *
*   Input U begin, U end and step: -3.14 3.14 0.4     *
*   Input V begin, V end and step: -3.14 3.14 0.2     *
*   (A torus is displayed in perspective mode).       *
*                                                     *
* --------------------------------------------------- *
* Ref.: "Graphisme dans le plan et dans l'espace en   *
*        Turbo Pascal 4.0 By R. Dony, MASSON Paris,   *
*        1990" [BIBLI 12].                            *
*                                                     *
*                   TPW Release By J-P Moreau, Paris. *
*                          (www.jpmoreau.fr)          *
*******************************************************
NOTE: In perspective mode, for a given rho value, you
can increase the dimension of object by increasing the
value of distance to screen (see last sample run).
------------------------------------------------------}
Program Surfaces;
Uses WinCrtMy,Type_def,Graph_3D;

VAR  U,Udebut,Ufin,dU : real_ar;
     V,Vdebut,Vfin,dV : real_ar;
     A,B,C : real_ar;
     typ   : word;
     rep   : char;


  Function FX(u,v : real_ar) : real_ar;
  begin
    Case typ of
    1: FX:= A*cos(u)*cos(v);        { Ellipsoid   }
    2: FX:= A*cos(u)*cos(v);        { sphere      }
    3: FX:=(A+B*cos(u))*cos(v);     { torus       }
    4: FX:=u                        { hyperboloid }
    end
  end;

  Function FY(u,v : real_ar) : real_ar;
  begin
    Case typ of
    1: FY:= B*cos(u)*sin(v);        { Ellipsoid   }
    2: FY:= B*cos(u)*sin(v);        { sphere      }
    3: FY:=(A+B*cos(u))*sin(v);     { torus       }
    4: FY:=v                        { hyperboloid }
    end
  end;

  Function FZ(u,v : real_ar) : real_ar;
  begin
    Case typ of
    1: FZ:=C*sin(u);                { Ellipsoid   }
    2: FZ:=C*sin(u);                { sphere      }
    3: FZ:=C*sin(u);                { torus       }
    4: FZ:=u*u-v*v                  { hyperboloid }
    end
  end;

  Procedure Data;
  var proj : byte;
  begin
    clrscr;
    writeln;
    writeln('    Projection Data:');
    writeln('    =============== ');
    gotoxy(1,5);
    writeln('  Kind of surface: 1. Ellipsoïd     ');
    writeln('                   2. Sphere        ');
    writeln('                   3. Torus         ');
    write  ('                   4. Hyperboloid   '); readln(typ);
    write  ('  Projection kind (1:perspective 2:parallel): '); readln(proj);
    gotoxy(1,10);
    if proj=1 then
    begin
      projection:=Perspective;
      write('  Value of Rho ...........: '); readln(rho);
      write('  Value of screen distance: '); readln(DE)
    end
    else
    begin
      projection:=Parallele;
      Rho:=1E20;
      write('  Value of screen distance: '); readln(DE)
    end;
    write('  Angle Theta ............: '); readln(theta);
    write('  Angle Phi ..............: '); readln(phi);
    writeln;
    write('  Input U begin, U end and step: ');
    readln(Udebut,Ufin,dU);
    write('  Input V begin, V end and step: ');
    readln(Vdebut,Vfin,dV);
    Case typ of
    1: begin A:=6.0; B:=3.0; C:=2.0 end;
    2: begin A:=4.0; B:=4.0; C:=4.5 end;
    3: begin A:=6.0; B:=3.0; C:=3.0 end;
    4: begin A:=1.0; B:=1.0; C:=1.0 end
    end
  end;  { Entrees }

  Procedure CourbesEnU;
  var x,y,z : real_ar;
  begin
    u:=udebut;
    while u <= ufin do
    begin
      v:=vdebut;
      x:=FX(u,v);
      y:=FY(u,v);
      z:=FZ(u,v);
      MoveXYZ(CrtDc,x,y,z);
      while v <= vfin do
      begin
        x:=FX(u,v);
        y:=FY(u,v);
        z:=FZ(u,v);
        DrawXYZ(CrtDc,x,y,z);
        v:=v+dv
      end;
      u:=u+du
    end
  end;

  Procedure CourbesEnV;
  var x,y,z : real_ar;
  begin
    v:=vdebut;
    while v <= vfin do
    begin
      u:=udebut;
      x:=FX(u,v);
      y:=FY(u,v);
      z:=FZ(u,v);
      MoveXYZ(CrtDc,x,y,z);
      while u <= ufin do
      begin
        x:=FX(u,v);
        y:=FY(u,v);
        z:=FZ(u,v);
        DrawXYZ(CrtDc,x,y,z);
        u:=u+du
      end;
      v:=v+dv
    end
  end;

  Procedure Affiche;
  var s1,s2,s3,s4,s5,s6,s7,s8,s9,s10 : string;
  begin
    if projection = Parallele then s1:='Infinity'
                              else str(rho:4:1,s1);
    str(theta:4:1,s2);  str(phi:4:1,s3);  str(DE:4:1,s4);
    str(udebut:5:2,s5); str(ufin:5:2,s6); str(du:5:2,s7);
    str(vdebut:5:2,s8); str(vfin:5:2,s9); str(dv:5:2,s10);
    gotoxy(3,2);
    writeln('    Rho='+s1+' Theta='+s2+' Phi='+s3+' Screen='+s4);
    gotoxy(3,3);
    writeln('    U=['+s5+','+s6+']  Step='+s7);
    gotoxy(3,4);
    writeln('    V=['+s8+','+s9+']  Step='+s10);
    gotoxy(50,23);
    Case typ of
    1: write('Ellipsoid');
    2: write('Sphere');
    3: write('Torus');
    4: write('Hyperboloid')
    end;
    gotoxy(5,24);
    write('Another view (y/n)? '); rep:=readkey
  end;

  {main program}
  BEGIN
    WinCrtInit(' PARAMETRIC SURFACES');
    Repeat
      Data;
      Clrscr;
      Cloture(0,MaxX,0,MaxY);
      Border(CrtDc,10);
      InitProj;
      Axes(CrtDc,1.5*A,1.5*B,1.5*C);
      CourbesEnU;
      CourbesEnV;
      Affiche
    Until rep='n';
    DoneWinCrt
  END.

{end of file surfpara.pas}