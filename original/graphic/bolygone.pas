{**************************************************************************
*                          THE  BOLYGONES                                 *
* A "bolygone" is a mathematical equivalent of a game for kids consisting *
* in creating nice figures with a piece of string around  nails on a      *
* wooden board. Let be a circle with a radius=1, centered at point 0,0.   *
* Let us consider an angle at centre i varying from  0 to 360 degrees     *
* with a given step p. The bolygone of order n is the envelope of all the *
* cords joining each angle i to the angle n.i, where n is an integer > 1. *
* ----------------------------------------------------------------------- *
* From "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0     *
* By R. Dony - MASSON, Paris 1990 page 49" [BIBLI 12].                    *
* ----------------------------------------------------------------------- *
* The program asks for order n (usually from 2 to 30) and step of         * 
* increase in degrees (usually from 1 to 4).                              *
**************************************************************************}
program Test_Bolygone;
{wincrt version with printing capability}
uses WinCrtMy,WinTypes,WinProcs,WObjects,WinPrint,Strings,Type_def,CrtGr2D;

 var ordre,pas: integer;
     adeg,arad: real;
     f1,f2,f3,f4: real;
     Pinfo : PPrinterInfo;
     oldmaxx,oldmaxy : integer;

    {read a string at screen location x,y}
    procedure ReadXY(x,y:integer; var s1:string);
    var c: char;
        S: array[0..40] of char;
        i,x0: integer;
    begin
      s1:=''; i:=1; x0:=x;
      repeat
        CursorTo(x,y);
        c:=readkey;
        if c<>#13 then
        begin
          s1:=s1+c;
          StrPCopy(S,s1);
          TextOut(CrtDC,8*x0,20*y-10,S,strlen(S));
          inc(x); inc(i)
        end
      until c=#13
    end;

    {input order and step for bolygone}
    procedure Data;
    var s1:string;
        S :array[0..60] of char;
        err:integer;
    begin
      clrscr;
      StrPCopy(S,'Order (integer) of bolygone (0 = end): ');
      TextOut(CrtDC,10,10,S,strlen(S));
      ReadXY(32,1,s1); val(s1,ordre,err);
      if ordre<>0 then
      begin
        StrPCopy(S,'Step value (integer) in degrees: ');
        TextOut(CrtDC,10,30,S,strlen(S));
        ReadXY(28,2,s1); val(s1,pas,err);
        f1:=-1.0; f2:=-f1; f3:=-1.0; f4:=-f3;
      end;
    end;

    {flag: FALSE=screen TRUE=printer}
    procedure Bolygone(P:HDC; flag:boolean);
    var c1,c2,pi:real;
        ch:array[0..2] of char;
        FenPen,FenPen1:HPen;
    begin
      if not flag then  { écran }
      begin
        FenPen:=CreatePen(ps_Solid,1,RGB(0,0,255));      { bleu  }
        FenPen1:=CreatePen(ps_Solid,1,RGB(180,50,50));   { rouge }
        SelectObject(P,FenPen1)
      end;
      Fenetre(f1,f2,f3,f4);
      Cloture(15,MaxX-30,90,MaxY-10);
      Bordure(P);
      if not flag then SelectObject(P,FenPen);
      adeg:=0; pi:=3.1416; c1:=0.7; c2:=0.9;
      repeat
        arad:=pi*adeg/180.0;
        MoveXY(P,c1*cos(arad),c2*sin(arad));
        LineXY(P,c1*cos(ordre*arad),c2*sin(ordre*arad));
        adeg:=adeg+pas
      until adeg > 360;
      TextOut(P,Round(0.825*MaxX),Round(0.05*MaxY),'N =',3);
      Str(ordre:2,ch);
      TextOut(P,Round(0.88*MaxX),Round(0.05*MaxY),ch,2);
      TextOut(P,Round(0.075*MaxX),Round(0.05*MaxY),'BOLYGONE',8);
    end;

  {main program}
  begin
    WinCrtInit('Bolygone');
    New(Pinfo, Init);
    repeat
      Data;
      clrscr;
      if ordre<>0 then
      begin
        Bolygone(CrtDC,FALSE);
        SortieGraphique;
        if rep='i' then           {send to printer}
          with Pinfo^ do
          begin
            oldmaxx:=MaxX; oldmaxy:=MaxY;
            MaxX:=3300; MaxY:=2350;  {HP laser landscape 300 DPI }
            StartDoc('Bolygone');
            Bolygone(PrintDC,TRUE);
            NewFrame;
            EndDoc;
            { restitution valeurs écran SVGA }
            MaxX:=oldmaxx; MaxY:=oldmaxy;
            SortieGraphique
          end
      end
    until (rep='n') OR (ordre=0);
    Dispose(Pinfo, Done);
    DoneWinCrt
  end.

  {Jean-Pierre Moreau - last modified august 1998

End of file bolygone.pas}