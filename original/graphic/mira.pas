{*************************************************************
*                   BEAUTY  OF  CHAOS                        *
* ---------------------------------------------------------- *
*   The two scientists Gumowski and Mira proposed in their   *
* publications the following iterative model:                *
*                                                            *
*           x n+1 = b.y n + F(x n)                           *
*           y n+1 =  -x n + F(x n+1)                         *
*                                                            *
*     with  F(x)  = a.x + (2(1-a)x²)/(1+x²)                  *
*                                                            *
*   According to the value of parameters a and b, and the    *
* coordinates of the starting poit, the achieved drawing can *
* have a surprising beauty.                                  *
* ---------------------------------------------------------- *
* Some interesting combinations:                             *
* a=0.1 b=0.99 x0=3 y0=0 iter=30000 f1=-8 f2=8 f3=-6 f4=6    *
* a=0.01 b=0.96 x0=5 y0=0 iter=25000 f1=-12 f2=12 f3=-8 f4=8 *
* a=0.32 b=1 x0=-5 y0=0 iter=30000 f1=-16 f2=16 f3=-12 f4=12 *
* a=0.31 b=1 x0=12 y0=0 iter=30000 f1=-40 f2=40 f3=-25 f4=25 *
* Can you find other ones ?                                  *
* ---------------------------------------------------------- *
* From "Graphisme dans le plan et dans l'espace avec Turbo   *
*       Pascal 4.0 de R. Dony - MASSON 1990, page 198"       *
*       [BIBLI 12].                                          *
*                                                            *
*                        TPW version by J-P Moreau, Paris    *
*                                (www.jpmoreau.fr)           *
*************************************************************}
Program mira;   {WINCRT version - print option inactive}
Uses WinCrtMy,WinTypes,WinProcs,WObjects,Strings,Type_Def,CrtGr2D;

VAR
        a,b,dx,x0,y0,x1,y1 : REAL_AR;
        f1,f2,f3,f4 : REAL_AR;
        n, niter    : LONGINT;
        CrtPen      : HPen;


    Function F(x:REAL_AR): REAL_AR;
    begin
      F:=a*x+(1.0-a)*(2.0*x*x/(1.0+x*x))
    end;


{main program}
BEGIN

  WinCrtInit('MIRA');

  Repeat
    {data section}
    clrscr;
    writeln;
    writeln('           MIRA''S AESTHETIC CHAOS');
    writeln('           ----------------------');
    writeln; writeln;
    write('  Input A: '); readln(a);
    write('  Input B: '); readln(b);
    writeln;
    write('  Input X0, Y0: '); readln(x0,y0);
    writeln;
    write('  Input number of iterations: '); readln(niter);
    writeln;
    write('  Input f1,f2,f3,f4: '); readln(f1,f2,f3,f4);
    dx:=(f2-f1)/MaxX;

    {drawing section}
    clrscr;
    CrtPen:=CreatePen(ps_Solid,1,RGB(0,0,255));
    SelectObject(CrtDC,CrtPen);
    Fenetre(f1,f2,f3,f4);
    Cloture(15,MaxX-25,80,MaxY-15);
    Bordure(CrtDc);
    TextOut(CrtDC,20,20,'MIRA''S CHAOS',12);

    {main drawing loop}
    for n:=1 to niter do
    begin
      MoveXY(CrtDc,x0,y0);
      LineXY(CrtDc,x0+dx,y0);
      x1:=b*y0 + F(x0);
      y1:=F(x1)- x0;
      x0:=x1; y0:=y1
    end;

    SortieGraphique {menu to exit graph}
  Until rep='n';
  DoneWinCrt

END.

{end of file mira.pas}