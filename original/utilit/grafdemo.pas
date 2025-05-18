{*******************************************************************
*                   DEMO 1 OF UNIT GRAPH_2D.PAS                    *
* ---------------------------------------------------------------- *
* This program demonstrates the use of the unit Graph_2d to draw a *
* 2D curve with automatic scaling.                                 *
*                                                                  *
*                               TPW version by J-P Moreau, Paris   *
*                                      (www.jpmoreau.fr)           *
*******************************************************************}
PROGRAM GRAPHDEMO;
Uses WinCrtMy, Winprocs, Type_def, Graph_2d;

{Draw curve y=sin(x) from x=0 with a dx step of 0.1 (250 pts) }
Procedure Graph_demo;
Var
  dx,x   : REAL_AR;
  i,nwin : INTEGER;
Begin
  dx:=0.1; x:=0; nwin:=10;
  {open graphic zone n° 10 (physical coordinates) }
  InitFenetre(CrtDC,nwin,0,10,-1.0,1.0);
  {main drawing loop}
  MoveXY(CrtDC,x,sin(x));
  for i:=1 to 250 do
  begin
    x:=x+dx;
    LineXY(CrtDC,x,sin(x));
  end;
  {write captions}
  Legendes(CrtDC,' FUNCTION  Y = SIN(X) ','X','Y');
  TextOut(CrtDC,MaxX-100,40,'250 points',10);
End;

{main program}
BEGIN
  WinCrtInit('DEMO 1 OF UNIT GRAPH_2D');
  Repeat
    Clrscr;
    Graph_demo;
    Sortiegraphique
  Until rep='n';
  DoneWinCrt
END.

{end of file grafdemo.pas}