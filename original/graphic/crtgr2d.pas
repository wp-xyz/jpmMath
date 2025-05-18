{********************************************************
*            UNIT  CRTGR2D  (English version)           *
* ----------------------------------------------------- *
* REFERENCE:                                            *
* Adapted from DOS unit GRAPH2D.PAS by Robert DONY -    *
* "MASSON 1990, Graphisme dans le plan et dans l'espace *
* en Turbo Pascal 4.0" to Windows By J-P Moreau, Paris  *
* [BIBLI 12].                                           *
* ----------------------------------------------------- *
* DESCRIPTION:                                          *
*   Procedures to draw a 2D curve with manual scaling.  *
********************************************************}
UNIT CrtGr2D;

INTERFACE
USES  WinCrtMy,WinDos,WinTypes,WinProcs,Strings,Type_Def,Savecrt;
{real_ar is defined in unit Type_Def.pas}

TYPE Table = ARRAY[1..300] OF real_ar;   {used by Polygone}
     Ptab  = ^Table;
     Table1= ARRAY[1..120] OF real_ar;   {used by Circle1}
     S12   = String[12];

VAR xp1,xp2,yp1,yp2,xgfen,xdfen,ybfen,yhfen: real_ar;
    xgclot,xdclot,ybclot,yhclot: integer;
    xrapport,yrapport,XOrig,YOrig : real_ar;
    T1       : Table1;
    S        : ARRAY[0..60] OF char;
    ferme    : STRING[6];
    rep      : char;
    CrtDC    : HDC;    {identification of CRT window}
		       {used by calling program     }


  PROCEDURE Fenetre(f1,f2,f3,f4:real_ar);
  PROCEDURE Cloture(c1,c2,c3,c4:integer);
  PROCEDURE PleinEcran;
  PROCEDURE LineXY(P:HDC;x,y: real_ar);
  PROCEDURE MoveXY(P:HDC;x,y: real_ar);
  PROCEDURE Bordure(P:HDC);
  PROCEDURE Polygone(P:HDC;X,Y: Ptab; Lim: integer; Mode: STRING);
  PROCEDURE Axes(P:HDC);
  PROCEDURE Gradue(P:HDC;UnitX,UnitY: real_ar);
  PROCEDURE Grille(P:HDC;UnitX,UnitY: real_ar);
  PROCEDURE Cercle(P:HDC;xc,yc,r: real_ar; trait: boolean);
  PROCEDURE Circle1(P:HDC;xc,yc,r: real_ar; trait: boolean);
  PROCEDURE WinCrtInit(Nom:PChar);
  PROCEDURE SortieGraphique;

IMPLEMENTATION

  PROCEDURE Fenetre;  {define limits in physical coordinates}
  BEGIN              
    xgfen:=f1;
    xdfen:=f2;
    ybfen:=f3;
    yhfen:=f4
  END;

  PROCEDURE Cloture;  {define limits in window pixels}
  BEGIN
    xgclot:=c1;
    xdclot:=c2;
    ybclot:=c3;
    yhclot:=c4;
    xrapport:=(xdclot-xgclot)/(xdfen-xgfen);
    yrapport:=(yhclot-ybclot)/(yhfen-ybfen);
  END;

  PROCEDURE PleinEcran;        {use all application CRT window}
  BEGIN                   
    Cloture(0,MaxX,0,MaxY)
  END;

  PROCEDURE Decoupage( P:HDC;x1,y1,x2,y2: real_ar);
  {stop drawing at limits defined by Cloture (clipping) }
  TYPE  region = (gauche,droite,basse,haute);
        Code   = SET OF Region;
  VAR   c,c1,c2: Code;
        x,y    : real_ar;
        xx1,yy1,xx2,yy2 : Integer;

    PROCEDURE CodeBin( x,y: real_ar; VAR c: code);
    BEGIN
      c:=[];
      IF x < xgfen THEN c:=[gauche]
                   ELSE IF x > xdfen THEN c:=[droite];
      IF y < ybfen THEN c:=c+[basse]
                   ELSE IF y > yhfen THEN c:=c+[haute]
    END;

  BEGIN
    CodeBin(x1,y1,c1);
    CodeBin(x2,y2,c2);
    WHILE (c1 <> []) OR (c2 <> []) DO
    BEGIN
      IF (c1*c2) <> [] THEN exit;
      IF c1 = [] THEN c:=c2 ELSE c:=c1;
      IF gauche IN C THEN
      BEGIN
        x:=xgfen;
        y:=y1+(y2-y1)*(xgfen-x1)/(x2-x1)
      END
      ELSE IF droite IN C THEN
      BEGIN
        x:=xdfen;
        y:=y1+(y2-y1)*(xdfen-x1)/(x2-x1)
      END
      ELSE IF basse IN C THEN
      BEGIN
        y:=ybfen;
        x:=x1+(x2-x1)*(ybfen-y1)/(y2-y1)
      END
      ELSE IF haute IN C THEN
      BEGIN
        y:=yhfen;
        x:=x1+(x2-x1)*(yhfen-y1)/(y2-y1)
      END;
      IF c=c1 THEN
      BEGIN
        x1:=x; y1:=y;
        CodeBin(x,y,c1)
      END
      ELSE
      BEGIN
        x2:=x; y2:=y;
        CodeBin(x,y,c2)
      END
    END;
    xx1:=round((x1-xgfen)*xrapport);
    yy1:=round((yhfen-y1)*yrapport);
    xx2:=round((x2-xgfen)*xrapport);
    yy2:=round((yhfen-y2)*yrapport);
    MoveTo(P,xgclot+xx1,MaxY-yhclot+yy1);
    LineTo(P,xgclot+xx2,MaxY-yhclot+yy2)
  END;

  PROCEDURE LineXY;  {in physical coordinates}
  BEGIN
    xp2:=x; yp2:=y;
    Decoupage(P,xp1,yp1,xp2,yp2);
    xp1:=xp2; yp1:=yp2
  END;

  PROCEDURE MoveXY;  {in physical coordinates}
  BEGIN
    xp1:=x; yp1:=y;
  END;

  PROCEDURE Bordure;  {draw a frame around drawing zone}
  BEGIN
    MoveXY(P,xgfen,ybfen);
    LineXY(P,xdfen,ybfen);
    LineXY(P,xdfen,yhfen);
    LineXY(P,xgfen,yhfen);
    LineXY(P,xgfen,ybfen);
  END;

  PROCEDURE Polygone;  {draw a polygon defined byr X(i), Y(i) }
  VAR  i: integer;     {in physical coordinates closed or open}
  BEGIN
    MoveXY(P,X^[1],Y^[1]);
    FOR i:=2 TO Lim DO
       LineXY(P,X^[i],Y^[i]);
    IF Mode = ferme THEN LineXY(P,X^[1],Y^[1])
  END;

  PROCEDURE Axes;  {draw axis Ox and Oy}
  BEGIN
    IF (XgFen < 0) AND (XdFen > 0) THEN XOrig:=0
                                   ELSE XOrig:=XgFen;
    IF (YbFen < 0) AND (YhFen > 0) THEN YOrig:=0
                                   ELSE YOrig:=YbFen;
    MoveXY (P,XgFen,YOrig);
    LineXY (P,XdFen,YOrig);
    MoveXY (P,XOrig,YbFen);
    LineXY (P,XOrig,YhFen);
  END;

    PROCEDURE CorrectX(xorig,xgfen,unitx:real_ar; VAR corX:real_ar);
    {used by procedure Gradue}
    VAR ntir: real_ar;
    BEGIN
      IF xorig=0 THEN
      BEGIN
        ntir:=(xorig-xgfen)/unitx;
	corx:=(ntir-trunc(ntir))*unitx
      END
      ELSE IF xorig > 0 THEN corx:=trunc(xorig/unitx+1)*unitx-xorig
			ELSE corx:=abs(xorig)+trunc(xorig/unitx)*unitx
    END;

    PROCEDURE CorrectY(yorig,ybfen,unity:real_ar; VAR corY:real_ar);
    {used by procedure Gradue}
    VAR ntir: real_ar;
    BEGIN
      IF yorig=0 THEN
      BEGIN
        ntir:=(yorig-ybfen)/unity;
	cory:=(ntir-trunc(ntir))*unity
      END
      ELSE IF yorig > 0 THEN cory:=trunc(yorig/unity+1)*unity-yorig
			ELSE cory:=abs(yorig)+trunc(yorig/unity)*unity
    END;

  PROCEDURE Gradue;    {graduate axis using steps unitx and unity}
  VAR  CorX,CorY,X,Y : real_ar;
       tiretx,tirety : real_ar;
       esp,xx,yy: integer;
       mot: ARRAY[0..10] OF char;
  BEGIN
    tiretx:=(xdfen-xgfen)/100; tirety:=(yhfen-ybfen)/100;
    IF UnitX > 0 THEN
    BEGIN
      CorrectX(xorig,xgfen,unitx,corX);
      x:=xgfen+corX;
      REPEAT
        MoveXY(P,X,YOrig+tirety);
        LineXY(P,X,YOrig-2*tirety);
        Str(X:5:2,mot);
        xx:=xgclot+round((X-xgfen)*xrapport);
        if MaxX > 1000 then esp := 50 else esp:=5;
        TextOut(P,xx-15,MaxY-ybclot+esp,mot,strlen(mot));
        X:=X + unitX
      UNTIL X > 1.05*XdFen
    END;

    IF UnitY > 0 THEN
    BEGIN
      CorrectY(yorig,ybfen,unity,corY);
      y:=ybfen+corY;
      REPEAT
        MoveXY(P,XOrig-tiretx,Y);
        LineXY(P,XOrig+2*tiretx,Y);
        Str(Y:5:2,mot);
        yy:=MaxY-yhclot+round((yhfen-Y)*yrapport);
	TextOut(P,5,yy-10,mot,strlen(mot));
        Y:=Y + unitY
      UNTIL Y > YhFen
    END;
  END;

  PROCEDURE Grille;     {draw optional grid in light red}
  VAR i: real_ar;
      HPPen,RedP,OldP: HPen;  {Pen handles}
  BEGIN
    if MaxX < 1000 then
    begin
      RedP:=CreatePen(ps_Dot,1,RGB(255,0,127));
      OldP:=SelectObject(P,RedP)  {dotted line in red, thickness=1}
    end
    else
    begin
      HPPen:=CreatePen(ps_Dot,1,RGB(0,0,0));
      OldP:=SelectObject(P,HPPen)  {dotted line in black, thickness=1}
    end;
    i:=xgfen;                     {vertical grid}
    WHILE i <= xdfen DO
    BEGIN
      MoveXY (P,i,ybfen);
      LineXY (P,i,yhfen);
      i:=i+unitX
    END;
    i:=yhfen;                     {horizontal grid}
    WHILE i >= ybfen DO
    BEGIN
      MoveXY (P,xgfen,i);
      LineXY (P,xdfen,i);
      i:=i-unitY
    END;
    SelectObject(P,OldP);
    if MaxX < 1000 then DeleteObject(RedP)         {free pen memory}
                   else DeleteObject(HPPen)
  END;

  PROCEDURE Cercle;
  {algorithm to draw a circle in physical
   coordinates (dotted line or normal line} 
  VAR  dx,s,c,x,y,aux : real_ar;
       n : integer;
  BEGIN
    s:=sin(pi/36); c:=cos(pi/36); dx:=r/50;
    x:=xc+r; y:=yc;
    MoveXY(P,x,y);
    FOR n:=2 TO 74 DO
    BEGIN
      aux:=xc+(x-xc)*c-(y-yc)*s;
      y  :=yc+(y-yc)*c+(x-xc)*s;
      x  :=aux;
      IF NOT trait THEN      {dotted line}
      BEGIN
         MoveXY(P,x,y);
         LineXY(P,x+dx,y)
      END
      ELSE                   {normal line}
        LineXY(P,x,y)
    END
  END;  

  Procedure Init_Table;
  { initialize table T1 used by procedure Circle1}
  var i: integer;
      dt,t: real_ar;
  begin
    dt:=pi/30; t:=0.0;
    for i:=1 to 60 do
    begin
      t:=t+dt;
      T1[i]:=cos(t);
      T1[i+60]:=sin(t)
    end
  end;

  PROCEDURE Circle1;
  { Faster algorithm to draw a circle in 
    physical coordinates in dotted line
    or normal line; read the values of
    cos(kpi/30) and sin(kpi/30) in T1 table} 
  VAR  c,dx,x,y : real_ar;
       i : integer;
  BEGIN
    dx:=r/30.0;
    c:=1.0;
    x:=xc+r; y:=yc; 
    MoveXY(P,x,y);
    for i:=1 to 60 do
    begin
      x:=xc+r*T1[i];
      y:=yc+c*r*T1[i+60];
      IF NOT trait THEN      {dotted line}
      BEGIN
        MoveXY(P,x,y);
        LineXY(P,x+dx,y)
      END
      ELSE                   {normal line}
        LineXY (P,x,y)
    end
  END;

  PROCEDURE WinCrtInit(Nom:PChar);
  { Open a CRT window with graphic capability (use unit
    WinCrtMy instead of standard WinCrt to have parameter
    CrtWindow visible). }
  BEGIN
    WindowOrg.X:=200;    {position of upper left corner}
    WindowOrg.Y:=200;    {and dimensions ofCRT window  }
    WindowSize.X:=600;
    WindowSize.Y:=500;
    StrCopy(WindowTitle,Nom);  {title of CRT window}
    MaxX:=590; MaxY:=475;      {useful drawing zone}
    InitWinCrt;                {call Borland procedure}
    CrtDC:=GetDC(CrtWindow)    {used by calling program}
  END;

  { options to exit graph}
  PROCEDURE SortieGraphique;
  VAR s:string; c: char;
  BEGIN
    REPEAT
      gotoxy(2,28); clreol;
      gotoxy(10,29);
      write('S:Save - R:Read - P:Print - C:Continue - E:Exit : ');
      c:=readkey; s:=''; rep:=#0;
      CASE Upcase(c) OF
	'S' : BEGIN    {Save picture to disk in B & W}
                Repeat
         	  gotoxy(2,28); clreol;
		  write('Input file name: '); Readln(s)
                Until length(s)>0;
		IF Copy(s,length(s)-3,1)<>'.' THEN
		  s:=s+'.IMG';
                gotoxy(2,28); clreol;
                write('Printing to disk...');
		if length(s)>0 then WCrtToFile(CrtDC,s);
	      END;
        'P' : rep:='i';  {print screen to printer}
	'R' : BEGIN      {read a B & W picture from disk}
		gotoxy(2,28); clreol;
		write('Input file name (without .img): '); Readln(s);
		IF Copy(s,length(s)-3,1)<>'.' THEN
		  s:=s+'.IMG';
		Clrscr;
                if length(s)>0 then WLoadCrt(CrtDC,s);
	      END;
	'C' : rep:='o'; {continue program}
	'E' : rep:='n'  {exit program}
      END
    UNTIL rep IN ['i','o','n']
  END;


BEGIN
  Init_Table;  { initialize the table of cos(kpi/30) et sin(kpi/30) }
               { values used by procedure Circle1.                  }
END.           { english version last modified 06/21/2000           }

{End of file crtgr2d.pas}