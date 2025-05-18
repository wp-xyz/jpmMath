{*******************************************************************
*                  DEMO 2 OF UNIT GRAPH_2D.PAS                     *
* ---------------------------------------------------------------- *
* This program demonstrates the use of the unit Graph_2d to draw a *
* 2D curve with automatic scaling with use of procedure CourbeXY.  *
*                                                                  *
*                               TPW version by J-P Moreau, Paris   *
*                                      (www.jpmoreau.fr)           *
*******************************************************************}
PROGRAM Grafdemo1;
Uses WinCrtMy, Winprocs, Type_def, Graph_2d;

Var  Y: RV;

{Draw curve y=sin(x)+2cos(2x)-3sin(4x) from x1=0 to x2=5*pi
 (256 pts) using automatic function CourbeXY()  }
Procedure Graph_demo;
Var
  dx,x,x1,x2: REAL_AR;
  i,ndata,nwin:INTEGER;
Begin
  ndata:=1024;
  nwin:=10;
  {Choose linear scale for axes Ox and Oy}
  Log_X:=FALSE; Log_Y:=FALSE;
  {store curve in table Y}
  x1:=0.0; x2:=21.25; dx:=(x2-x1)/(ndata-1);
  x:=x1-dx;
  for i:=1 to ndata do
  begin
    x:=x+dx;
    Y^[i] := sin(x)+2*cos(2*x)-3*sin(4*x);
  end;
  {draw curve (scaling adapted to window size) }
  CourbeXY(CrtDC,ndata,nwin,Y,x1,x2);
  {write captions}
  Legendes(CrtDC,' Y = SIN(X) + 2 COS(2X) - 3 SIN(4X) ','X','Y');
  TextOut(CrtDC,MaxX-150,40,'1024 points',11);
  TextOut(CrtDC,MaxX-185,MaxY-130,'Linear scaling in X and Y',25);
End;

{main program}
BEGIN
  WinCrtInit('DEMO 2 OF UNIT GRAPH_2D');
  New(Y);
  Repeat
    Clrscr;
    Graph_demo;
    Sortiegraphique
  Until rep='n';
  Dispose(Y);
  DoneWinCrt
END.

{end of file grafdem1.pas}