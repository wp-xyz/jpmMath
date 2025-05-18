{***************************************************************
*               THE PROBLEM OF CHESS KNIGHT                    *
* A question, in chess game, is to know wether a knight can    *
* go through all the board squares passing only once on each   *
* square. The answer is yes for a 8 x 8 board, as this small   *
* program shows. We use here the law givan by Warnsdorff in    *
* 1823 that almost always allows finding a solution:           *
* "In each strike, you must play the knight to the square from *
* which there are the less exits towards the squares not yet   *
* full."                                                       *
* ------------------------------------------------------------ *
*  After "Graphisme dans le plan et dans l'espace avec Turbo   *
*  Pascal 4.0 de R. Dony - MASSON 1990 page 227".              *
*                                                              *
*                            TPW Release By J-P Moreau, Paris. *
***************************************************************}
program Knights;

uses WinCrt,WinTypes,WinProcs,Strings,CrtGr2D;

const M = 335;      {set size of chessboard}
      N = 38;

var  echiquier      : array[-1..10,-1..10] of byte;
     Xcasefuite     : array[1..8] of byte;
     Ycasefuite     : array[1..8] of byte;
     Nbrecasefuite  : array[1..8] of byte;
      casefuite     : boolean;
     lig,col,ct     : integer;
     duree          : longint;
     compteur,cpt2,ordre,ligne,colonne: byte;
     FenPen,FenPen1:HPen;

const pasy=N;

  procedure ReadOrder;
  begin
    Clrscr;
    writeln;
    repeat
      write('  Size of chessboard (4 to 8): '); readln(ordre)
    until (ordre > 3) and (ordre < 9);
    duree:=100000;  {for speed regulation}
  end;

  procedure ReadStart;
  begin
    writeln;
    writeln('  Starting coordinates: ');
    repeat
      write('    column n° '); readln(colonne)
    until (colonne > 0) and (colonne <= ordre);
    writeln;
    repeat
      write('    line n° '); readln(ligne)
    until (ligne > 0) and (ligne <= ordre);
  end;

  procedure DisplaySolution;
  var  ch,s1,s2: string;
       t: array[0..20] of char;
  begin
    str(ligne,s1);
    str(colonne,s2);
    ch:=concat('(',s2,',',s1,')');
    StrPCopy(t,ch);
    inc(ct);
    col:=ct mod 3;
    Case col of
      0:  textout(CrtDc,col*55+420,lig,t,strlen(t));
      1:  textout(CrtDc,col*55+420,lig,t,strlen(t));
      2:  begin
	    textout(CrtDc,col*55+420,lig,t,strlen(t));
            lig:=lig+15
          end
    end
  end;

  procedure Init;
  var l,c: integer;
  begin
    fillchar(echiquier,sizeof(echiquier),1);
    for l:=1 to ordre do
      for c:=1 to ordre do Echiquier[l,c]:=0;
    compteur:=0;
    casefuite:=true; ct:=-1;
    lig:=10+N*(8-ordre); col:=0;
    FenPen:=CreatePen(ps_Solid,1,RGB(0,0,255));
    FenPen1:=CreatePen(ps_Solid,1,RGB(255,255,255));
  end;

  procedure DrawKnight(P:HDC; x,y:integer);
  var xc,yc: integer;
  begin
    moveto(P,x+9,y);
    xc:=x+24; yc:=y+33; moveto(P,xc,yc);
    xc:=xc+25; lineto(P,xc,yc);
    xc:=xc+5;  yc:=yc-13; lineto(P,xc,yc);
    xc:=xc-5;  yc:=yc-10; lineto(P,xc,yc);
    xc:=xc-10; yc:=yc-2;  lineto(P,xc,yc);
    yc:=yc-3;  lineto(P,xc,yc);
    xc:=xc-5;  yc:=yc+3;  lineto(P,xc,yc);
    xc:=xc-2;  yc:=yc-3;  lineto(P,xc,yc);
    xc:=xc-3;  yc:=yc+3;  lineto(P,xc,yc);
    xc:=xc-5;  yc:=yc+3;  lineto(P,xc,yc);
    xc:=xc-10; yc:=yc+7;  lineto(P,xc,yc);
    yc:=yc+3;  lineto(P,xc,yc);
    xc:=xc+7;  yc:=yc-2;  lineto(P,xc,yc);
    xc:=xc-5; yc:=yc+4;  lineto(P,xc,yc);
    xc:=xc+8; yc:=yc-2;  lineto(P,xc,yc);
    xc:=xc+10; lineto(P,xc,yc);
    xc:=xc-10; yc:=yc+9;  lineto(P,xc,yc);
    yc:=yc+3;  lineto(P,xc,yc);
    Ellipse(P,x+30,y+10,x+36,y+16)  {eye circle}
  end;

  procedure PutKnights(P:HDC);
  var ch:array[0..2] of char;
      x,y:integer;
      i:longint;
      z: real;
  begin
    z:=1;
    x:=50*(colonne-1); y:=M-N*ligne;
    SelectObject(CrtDC,FenPen);     {blue pen}
    DrawKnight(P,x,y);
    MessageBeep(1);
    for i:=1 to duree do z:=z*z*z;  {delay}
    SelectObject(CrtDC,FenPen1);    {white pen to erase}
    DrawKnight(P,x,y);
    inc(compteur);
    str(compteur:2,ch);
    textout(P,x+28,y+13,ch,2);
    Echiquier[ordre-ligne+1,colonne]:=compteur;
    DisplaySolution
  end;

  procedure DrawBoard(P: HDC);
  var i: byte;
  begin
    SelectObject(CrtDC,FenPen);    { plume bleue }
    Rectangle(P,5,5,570,350);
    for i:=0 to ordre do
    begin
      moveto(P,10,M-i*N);
      lineto(P,10+50*ordre,M-i*N)
    end;
    for i:=0 to ordre do
    begin
      moveto(P,10+i*50,M);
      lineto(P,10+50*i,M-N*ordre)
    end
  end;

  procedure FreeSquares;
  var  cpt,i,ll,cc,ind: byte;

    procedure GoAround(lg,cl,ind:byte);
    var  l,c: integer;

      procedure Test;
      begin
        if Echiquier[ordre-l+1,c]=0 then
        begin
          inc(cpt);
          if ind=1 then
          begin
            Xcasefuite[cpt]:=c;
            Ycasefuite[cpt]:=l
          end
        end
      end;

    begin
      cpt:=0;
      casefuite:=true;
      l:=lg+2; c:=cl-1; Test;
      l:=lg+2; c:=cl+1; Test;
      l:=lg+1; c:=cl+2; Test;
      l:=lg-1; c:=cl+2; Test;
      l:=lg-2; c:=cl+1; Test;
      l:=lg-2; c:=cl-1; Test;
      l:=lg-1; c:=cl-2; Test;
      l:=lg+1; c:=cl-2; Test;
    end;

  begin  {FreeSquares}
    FillChar(Xcasefuite,sizeof(Xcasefuite),0);
    FillChar(Ycasefuite,sizeof(Ycasefuite),0);
    GoAround(ligne,colonne,1);
    if cpt <> 0 then
    begin
      cpt2:=cpt;
      for i:=1 to cpt2 do
      begin
        cc:=Xcasefuite[i];
        ll:=Ycasefuite[i];
        GoAround(ll,cc,2);
        nbrecasefuite[i]:=cpt;
      end
    end
    else casefuite:=false
  end;

  procedure NextChoice;
  var  min,pos,i: byte;
  begin
    min:=255; pos:=0;
    for i:=1 to cpt2 do
      if nbrecasefuite[i] < min then
      begin
        min:=nbrecasefuite[i];
        pos:=i
      end;
    ligne:=Ycasefuite[pos];
    colonne:=Xcasefuite[pos]
  end;

  procedure SeekSolution(P:HDC);
  begin
    while (compteur < ordre*ordre) and (casefuite) do
    begin
      if KeyPressed then exit;
      PutKnights(P);
      FreeSquares;
      if casefuite then NextChoice
    end
  end;

  BEGIN    {main program}
    WinCrtInit(' PROBLEM OF CHESS KNIGHT');
    repeat
      ReadOrder;
      ReadStart;
      Init;
      Clrscr;
      DrawBoard(CrtDc);
      gotoxy(16,2); write('PROBLEM of CHESS KNIGHT');
      SeekSolution(CrtDc);
      Sortiegraphique  {exit menu}
    until rep='n';
    DeleteObject(FenPen);
    DeleteObject(FenPen1);
    DoneWinCrt
  END.

{last update - August, 1998}