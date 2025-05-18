{************************************************************
*                UNIT TO DRAW 3D CURVES                     *
* --------------------------------------------------------- *
* Reference:                                                *
* From Pascal Unit GRAPH3D.PAS By Robert DONY - MASSON 1990 *
* "Graphisme dans le plan et dans l'espace en Turbo Pascal  *
*  4.0" - Adapted to Windows By J-P Moreau, july 1993.      *
*  [BIBLI 12].                                              *
*                                                           *
************************************************************}
Unit Graph_3D;

Interface
USES  WinCrtMy,WinTypes,WinProcs,Strings,Type_Def;
{real_ar type is defined in Type_Def.pas}

const blanc=15;      {white}
      noir=0;        {black}
      bleu=9;        {blue}
      jaune=14;      {yellow}
      echelle=1.35;

      {projection kinds}
type  GenreProjection = (perspective,parallele);

var   rho,theta,phi,de     : real_ar;
      aux1,aux2,aux3,aux4  : real_ar;
      aux5,aux6,aux7,aux8  : real_ar;
      projection           : GenreProjection;
      xobs,yobs,zobs       : real_ar;
      xproj,yproj          : real_ar;
      xecran,yecran        : integer;
      xgclot,xdclot,ybclot,yhclot: integer;
      ferme,ouvert         : string;
      nbrecoul             : byte;

  Procedure Cloture(c1,c2,c3,c4:integer);
  Procedure InitProj;
  Procedure Project(x,y,z: real_ar);
  Procedure DrawXYZ(P:HDC;x,y,z: real_ar);
  Procedure MoveXYZ(P:HDC;x,y,z: real_ar);
  Procedure Border(P:HDC;bord:integer);
  Procedure Axes(P:HDC;x,y,z:real_ar);
  Procedure Pause;
  Procedure WinCrtInit(Nom:PChar);
  Procedure Bip;

Implementation

  Procedure Cloture;
  begin
    xgclot:=c1;
    xdclot:=c2;
    ybclot:=c3;
    yhclot:=c4;
  end;

  Procedure InitProj; { calcul des variables auxiliaires }
  var th,ph: real_ar;
  begin
    th:=pi*theta/180; ph:=pi*phi/180;
    aux1:=sin(th); aux2:=sin(ph); aux3:=cos(th); aux4:=cos(ph);
    aux5:=aux3*aux2; aux6:=aux1*aux2; aux7:=aux3*aux4; aux8:=aux1*aux4;
  end;

  Procedure Project;
  begin
    xobs:=-x*aux1+y*aux3; yobs:=-x*aux5-y*aux6+z*aux4;
    if projection=perspective then
    begin
      zobs:=-x*aux7-y*aux8-z*aux2+rho;
      xproj:=de*xobs/zobs; yproj:=de*yobs/zobs;
    end
    else
    begin
      xproj:=de*xobs; yproj:=de*yobs
    end
  end;

  Procedure DrawXYZ;
  begin
    Project(x,y,z);
    xecran:=round(xproj*echelle+MaxX div 2);
    yecran:=round(MaxY div 2-yproj);
    lineto(P,xecran,yecran)
  end;

  Procedure MoveXYZ;
  begin
    Project(x,y,z);
    xecran:=round(xproj*echelle+MaxX div 2);
    yecran:=round(MaxY div 2-yproj);
    moveto(P,xecran,yecran)
  end;

  Procedure Border;
  begin
    rectangle(P,bord,bord,MaxX-bord,MaxY-bord)
  end;

  Procedure Axes;  { Tracé des axes }
  begin
    MoveXYZ(P,0,0,0); DrawXYZ(P,x,0,0);
    Inc(xecran,5); Inc(yecran,5);
    textout(P,xecran,yecran,'X',1);
    MoveXYZ(P,0,0,0); DrawXYZ(P,0,y,0);
    Inc(xecran,5); Inc(yecran,5);
    textout(P,xecran,yecran,'Y',1);
    MoveXYZ(P,0,0,0); DrawXYZ(P,0,0,z);
    Inc(xecran,10); Dec(yecran,5);
    textout(P,xecran,yecran,'Z',1);
  end;

  Procedure Pause;
  var ch: char;
  begin
    repeat until keypressed;
    ch:=readkey
  end;

  procedure WinCrtInit(Nom:PChar);
  begin
    WindowOrg.X:=100;
    WindowOrg.Y:=100;
    WindowSize.X:=650;
    WindowSize.Y:=450;
    StrCopy(WindowTitle,Nom);
    MaxX:=625; MaxY:=375;
    InitWinCrt;
    CrtDC:=GetDC(CrtWindow);
  end;

  Procedure Bip;
  begin
    MessageBeep(0)
  end;


END.

{end of file graph_3d.pas}