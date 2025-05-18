{**********************************************************************
*               SIMULATION OF AN ELLIPTICAL BILLARD                   *
* ------------------------------------------------------------------- *
* In one of his books (*), the polish mathematician Hugo Steinhaus    *
* outlines the idea of an elliptical billard table! He distinguishes  *
* three cases of trajectory, depending on throwing the ball:          *
* 1) between the two focuses of the ellipse 2) between one focus and  *
* the edge of the billard  3) passing through a focus.  This program  *
* allows visualizing the envelope of trajectories in the three cases  *
* with a  limitation of 100 rebounds.                                 *
* ------------------------------------------------------------------- *
* From "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0 *
* By R. Dony - MASSON 1990, page 113" [BIBLI 12].                     *
* ------------------------------------------------------------------- *
* (*) Mathématiques en instantanés By Hugo Steinhaus - Flammarion 1964*
**********************************************************************}
Program billard;    {WinCrtMy version with print option}
Uses WinCrtMy,WinTypes,WinProcs,Strings,Type_Def,CrtGr2D,WinPrint;

 Const bord = 10;         {margin in pixels           }
       A = 5;  B = 3;     {parameters of ellipse      }
       limite = 100;      {maximum number of rebounds }
       increment = 0.1;   {drawing accuracy of ellipse}

 Var f1,f2,f3,f4,m,n,x1,y1,x2,y2: real_ar;
     nbrerebonds: integer;
     intersection,aborted:boolean;
     oldmaxx,oldmaxy : integer;
     FenPen:HPen;
     Pinfo:PPrinterInfo;
     ch : char;

    Procedure Data;
    begin
      clrscr;
      writeln;
      writeln('         SIMULATION OF AN ELLIPTICAL BILLARD');
      writeln('         -----------------------------------');
      gotoxy(3,10);
      write('La ball is thrown along a line Y = MX + N');
      gotoxy(3,13);
      write('Input M et N: '); readln(m,n);
    end;

    {intersection is TRUE if the line meets the ellipse}
    procedure Test;
    var delta: real_ar;
    begin
      delta:=A*A*B*B*(A*A*m*m+B*B-n*n);
      intersection:=(delta>0)
    end;

    {draw billard elliptic contour in blue}
    procedure DessineEllipse(P:HDC);
    var angle,x1,y1,x2,y2: real_ar;
    begin
      FenPen:=CreatePen(ps_Solid,1,RGB(0,0,255));
      SelectObject(P,FenPen);
      Cloture(bord,MaxX-3*bord,100,MaxY-bord);
      Bordure(P);
      angle:=0; x1:=A; Y1:=0;
      MoveXY(P,x1,y1);
      while angle < 2*pi do
      begin
        angle:=angle+increment;
        x2:=A*cos(angle);
        y2:=B*sin(angle);
        LineXY(P,x2,y2)
      end
    end;

    {draw a cross at each focus of the ellipse}
    procedure TraceFoyers(P:HDC);
    var xf,l: real_ar;

      procedure crossatfocus(P:HDC; x,y: real_ar);
      begin
        MoveXY(P,x-l,y);
        LineXY(P,x+l,y);
        MoveXY(P,x,y-l);
        LineXY(P,x,y+l)
      end;

    begin
      xf:=sqrt(a*a-b*b); l:=0.03*a;
      crossatfocus(P,-xf,0);
      crossatfocus(P,xf,0)
    end;

    procedure swap;
    var aux:real_ar;
    begin
      aux:=x1; x1:=x2; x2:=aux;
      aux:=y1; y1:=y2; y2:=aux
    end;

    {seek intersection of line with ellipse}
    procedure ChercheIntersection;
    var delta,denom: real_ar;
    begin
      delta:=a*a*b*b*(a*a*m*m+b*b-n*n);
      denom:=a*a*m*m+b*b;
      if delta>0 then x1:=(-a*a*m*n+sqrt(delta))/denom;
      y1:=m*x1+n;
      if delta>0 then x2:=(-a*a*m*n-sqrt(delta))/denom;
      y2:=m*x2+n;
      if m>0 then if y1>y2 then swap
      else if m<0 then if y1<y2 then swap
      else if x1>x2 then swap
    end;

    {draw the 100 rebounds}
    procedure Trajectory(P:HDC);
    var xprec,yprec,tgphi,tgteta,angleinc,anglephi,oldm,oldn: real_ar;
    begin
      nbrerebonds:=0;
      oldm:=m;  oldn:=n;
      MoveXY(P,f1,m*f1+n);
      LineXY(P,x2,y2);
      Inc(nbrerebonds);
      Repeat
        if x2<>0 then tgphi:=a*a*y2/(b*b*x2)
                 else tgphi:=1E20;
        tgteta:=m;
        angleinc:=arctan((tgphi-tgteta)/(1+tgphi*tgteta));
        anglephi:=arctan(tgphi);
        m:=sin(angleinc+anglephi)/cos(angleinc+anglephi);
        if abs(n) < 1E-3 then begin m:=0; n:=0 end
                         else n:=y2-m*x2;
        xprec:=x2; yprec:=y2;
        ChercheIntersection;
        if abs(xprec-x1) > 1E-6 then swap;
        MoveXY(P,xprec,yprec);
        LineXY(P,x2,y2);
        Inc(nbrerebonds);
      Until nbrerebonds > limite;
      m:=oldm; n:=oldn
    end;


    BEGIN  {main program}
      WinCrtInit('ELLIPTICAL BILLARD');
      f1:=-1.1*A; f2:=-f1; f3:=-1.1*B; f4:=-f3;
      Fenetre(f1,f2,f3,f4);
      New(Pinfo,Init);
      Repeat
	rep:='n';
        Data;
        Test;
        if intersection then  {la droite coupe bien l'ellipse}
        begin
          Clrscr;
          DessineEllipse(CrtDc);
          TraceFoyers(CrtDc);
          ChercheIntersection;
	  Trajectory(CrtDc);
	  TextOut(CrtDc,25,MaxY-90,'ELLIPTICAL',10);
	  TextOut(CrtDc,475,MaxY-90,'BILLARD',7);
        end
	else  {the line does not meet the ellipse}
	begin
	  MessageBeep(0);
	  MessageBox(CrtWindow,'The line does not meet the ellipse !',
		     'Warning',mb_Ok OR mb_IconInformation);
	  rep:='o'
        end;
	if rep<>'o' then Sortiegraphique;
        {section impression}
	if rep='i' then
	begin
	  With Pinfo^ do
	  begin
	    oldmaxx:=MaxX; oldmaxy:=MaxY;
	    MaxX:=3300; MaxY:=2350;    {HP laser landscape}
            StartDoc('BILLARD');
	    DessineEllipse(PrintDc);
            TraceFoyers(PrintDc);
            ChercheIntersection;
	    Trajectory(PrintDc);
	    TextOut(PrintDc,50,MaxY-150,'ELLIPTICAL',10);
	    TextOut(PrintDc,MaxX-350,MaxY-150,'BILLARD',7);
            NewFrame;
	    EndDoc;
	    MaxX:=oldmaxx; MaxY:=oldmaxy
	  end;
	  SortieGraphique
        end
      Until rep='n';
      DeleteObject(FenPen);
      Dispose(Pinfo, Done);
      DoneWinCrt
    END.

    {Author: Jean-Pierre Moreau - last modified 06/21/2000

End of file billard.pas}