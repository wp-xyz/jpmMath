{*************************************************************
*             Chess: Problem of the Eight Queens             *
* ---------------------------------------------------------- *
* Description:                                               *
* This famous problem consists in finding locations of eight *
* queens on a chess board,  such as no one can be taken by   *
* another one. After Arthur Engel Cedic (1979), this problem *
* was posed for the first time in 1848 by Max Bezzel in a    *
* Chess Journal. This publication gave way to an extraordi-  *
* nary passion and on sept. 21th, 1850 a certain Dr Nauck    *
* gave all the solutions when the famous mathematician Gauss *
* had only found 72 out of 92. But it is true that he had    *
* better to attend to... This small program gives all the    *
* solutions from n=4 to n=8 (size of board) by using the so  *
* called "backtracking method".                              *
* ---------------------------------------------------------- *
* REFERENCE:                                                 *
* "Graphisme dans le plan et dans l'espace avec Turbo Pascal *
*  4.0 By R. Dony - MASSON, Paris 1990".                     *
*                                                            *
*                          TPW Release By J-P Moreau, Paris. *
*************************************************************}
Program Qweens;

Uses WinCrt1,WinTypes,WinProcs,Strings;

const M = 365;
      N = 38;
      Blue  =  9;
      White = 15;

type Rt = array[0..8] of byte; 

var  R        : Rt;                {table storing one solution}
     S        : array[1..92] of Rt;      {maximum 92 solutions}
     Xr,Yr    : array[1..8] of integer;
     i,j,order: byte;
     ok,fin   : boolean;
     x,y,k,l  : integer;
     num, old : word;
     rep      : char;
     Numsol   : array[1..8] of word;
     CrtPen,CrtPen1: HPen;                  {blue & white pens}
     CrtDC: HDC;                               {Device context}

const pasy=N;

  procedure ReadOrder;
  begin
    clrscr;
    writeln;
    repeat
      write('  Size of Board (4 to 8): '); readln(order)
    until (order > 3) and (order < 9);
    fillchar(R,sizeof(R),0)
  end;

  procedure DrawQween(P:HDC; x,y:integer);
  var xc,yc:integer;
  begin
    xc:=x+15; yc:=y+35; moveto(P,xc,yc);
    xc:=xc+30; lineto(P,xc,yc);
    yc:=yc-3;  lineto(P,xc,yc);
    xc:=xc-30; lineto(P,xc,yc);
    yc:=yc+3;  lineto(P,xc,yc);
    xc:=xc+5; yc:=yc-5; moveto(P,xc,yc);
    yc:=yc+2;  lineto(P,xc,yc);
    xc:=xc+20; lineto(P,xc,yc);
    yc:=yc-2;  lineto(P,xc,yc);
    xc:=xc-20; lineto(P,xc,yc);
    xc:=xc-5; yc:=yc-5; lineto(P,xc,yc);
    xc:=xc-5; yc:=yc-8; lineto(P,xc,yc);
    xc:=xc+2; yc:=yc-1; moveto(P,xc,yc);
    xc:=xc+3; yc:=yc+4; lineto(P,xc,yc);
    xc:=xc+5; yc:=yc-2; lineto(P,xc,yc);
    yc:=yc-4;   lineto(P,xc,yc);
    xc:=xc+1;   moveto(P,xc,yc);
    xc:=xc+1; yc:=yc+3; lineto(P,xc,yc);
    xc:=xc+6; yc:=yc-2; lineto(P,xc,yc);
    xc:=xc+1; yc:=yc-3; lineto(P,xc,yc);
    xc:=xc+2;   moveto(P,xc,yc);
    xc:=xc+1; yc:=yc+3; lineto(P,xc,yc);
    xc:=xc+6; yc:=yc+2; lineto(P,xc,yc);
    xc:=xc+1; yc:=yc-3; lineto(P,xc,yc);
    xc:=xc+1;   moveto(P,xc,yc);
    yc:=yc+4;   lineto(P,xc,yc);
    xc:=xc+5; yc:=yc+2; lineto(P,xc,yc);
    xc:=xc+3; yc:=yc-4; lineto(P,xc,yc);
    xc:=xc+2; yc:=yc+1; moveto(P,xc,yc);
    xc:=xc-5; yc:=yc+8; lineto(P,xc,yc);
    xc:=xc-5; yc:=yc+5; lineto(P,xc,yc);
    Ellipse(P,x+ 8,y+13,x+12,y+17);
    Ellipse(P,x+18,y+10,x+22,y+14);
    Ellipse(P,x+28,y+ 8,x+32,y+12);
    Ellipse(P,x+38,y+10,x+42,y+14);
    Ellipse(P,x+48,y+13,x+52,y+17);
  end;

  procedure Init;
  var c: integer;
  begin
    clrscr;
    CrtPen:=CreatePen(ps_Solid,1,RGB(0,0,255));
    CrtPen1:=CreatePen(ps_Solid,1,RGB(255,255,255));
    Xr[1]:=0; Yr[1]:=M-42;
    for i:=2 to order do Xr[i]:=Xr[i-1]+50;
    for i:=2 to order do Yr[i]:=Yr[i-1]-N;
    fin:=FALSE;
    i:=0; num:=0; k:=-1;
    l:=(9-order)*N-25;
    SelectObject(CrtDc,CrtPen);
    DrawQween(CrtDc,6,0);
    gotoxy(5,2); write(order);
    Numsol[4]:= 2;
    Numsol[5]:=10;
    Numsol[6]:= 4;
    Numsol[7]:=40;
    Numsol[8]:=92;
  end;

  procedure DrawBoard(P:HDC);
  var i: byte;
  begin
    for i:=0 to order do
    begin
      moveto(P,5,M-i*N);
      lineto(P,5+50*order,M-i*N)
    end;
    for i:=0 to order do
    begin
      moveto(P,5+i*50,M);
      lineto(P,5+i*50,M-N*order)
    end
  end;

  procedure TestPreviousQweens;
  begin
    ok:=TRUE;
    j:=1;
    While (j <= i-1) AND (ok) do
    begin
      if (R[i]=R[j]) OR (abs(R[i]-R[j])=i-j) then ok:=FALSE;
      Inc(j)
    end
  end;

  procedure ColumnDown;
  begin
    R[i]:=0;
    Dec(i);
    if i=0 then fin:=TRUE
           else Dec(i)
  end;

  procedure DrawSolution(n:word;coul:word);
  var i : integer;
  begin
    if coul=Blue then
      SelectObject(CrtDc,CrtPen)
    else
      SelectObject(CrtDc,CrtPen1);
    for i:=1 to order do
      DrawQween(CrtDc,Xr[i],Yr[S[n][i]])
  end;

  procedure Display;
  var ch,chaine : string[8];
      s1: array[0..8] of char;
      j:byte;
  begin
    chaine:='';
    inc(num); S[num]:=R;
    for j:=1 to order do
    begin
      Str(R[j]:1,ch);
      chaine:=Concat(chaine,ch)
    end;
    Inc(k);
    StrPCopy(s1,chaine);
    Case (k MOD 3) of
      0: textout(CrtDc,410,l,s1,strlen(s1));
      1: textout(CrtDc,480,l,s1,strlen(s1));
      2: begin
	   textout(CrtDc,550,l,s1,strlen(s1));
           Inc(l,15)
         end
    End;
  end;

  procedure SeekSolution;
  begin
    While Not(fin) do
    begin
      Inc(i); 
      While (i<=order) AND (Not(fin)) do
      begin
        ok:=FALSE;
        While (R[i]<order) AND (Not(ok)) do
        begin
          R[i]:=R[i]+1;
          TestPreviousQweens
        end;
        if Not(ok) then ColumnDown;
	Inc(i);
	If KeyPressed then exit
      end;
      if Not(fin) then
      begin
        Display;
        R[order]:=0;
        i:=order-2
      end
    end;
    gotoxy(21,2);
    MessageBeep(0);
    writeln(num:2,' solutions.')
  end;


  PROCEDURE WinInit(Nom:PChar);
  { open a graphic and/or CRT window under the condition
    to use WinCrtMy or WinCrt1 instead of WinCrt to make 
    CrtWindow visible.                 }
  BEGIN
    WindowOrg.X:=0;       {upper left corner position}
    WindowOrg.Y:=20;    {and dimensions of Crt window}
    WindowSize.X:=640;
    WindowSize.Y:=440;
    StrCopy(WindowTitle,Nom);         {window caption}
    InitWinCrt;
    CrtDC:=GetDC(CrtWindow);
  END;

  Function Valid: Boolean;
  begin
    if (num>=0) AND (num < Numsol[order]+1) then Valid:=TRUE
					    else Valid:=FALSE
  end;

  {main program}
  BEGIN
    WinInit(' PROBLEM OF QUEENS (4 TO 8)');
    rep:='o';
    Repeat
      ReadOrder;
      Clrscr;
      Init;
      DrawBoard(CrtDc);
      SeekSolution;
      Repeat
        Repeat
	  gotoxy(13,3);
	  write('Number of solution to draw:    ');
          gotoxy(41,3); if Valid then old:=num; read(num)
        Until Valid;
        if num<>0 then
        begin
          DrawSolution(old,White);
	  DrawSolution(num,Blue)
        end
      Until num=0;
      gotoxy(20,4); write(' Continue (o/n) ? ');
      rep:=readkey
    Until rep='n';
    DeleteObject(CrtPen);
    DeleteObject(CrtPen1);
    DoneWinCrt
  END.

{end of file queens.pas}