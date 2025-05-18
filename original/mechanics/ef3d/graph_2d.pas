{*******************************************************************
*        UNIT Graph_2D.Pas  English Release 1.0 (04/18/2000)       *
*                                   By J-P Moreau, Paris (1)       *
* ---------------------------------------------------------------- *
* Set of procedures to draw a 2D curve y=f(x) in a window defined  *
* by procedure Fenetre() and to manually or automatically adjust   *
* scales of axes in linear or logarithmic mode.                    *
*                                                                  *
* If input of printer DPI by calling program is correct, the graph *
* will have centimetric graduations in axes.                       *      
*                                                                  *
* Main global variables:                                           *
*                                                                  *
* xgclot,xdclot,ybcot,yhclot : limits in pixels of the window defi-*
*                              ned by Fenetre().                   *
* MaxX, MaxY   : screen resolution in X et Y                       *
*                (Ex.: MaxX=639  MaxY=479 for VGA).                *
* Log_X, Log_Y : log. scales if TRUE.                              *
* Ech_auto     : automatic scales if TRUE.                         *
*                                                                  *
*         (1) from a pascal dos program by J.-P. Dumont, France.   *
*******************************************************************}
UNIT GRAPH_2D;

INTERFACE
USES  WinCrt1,WinTypes,WinProcs,Strings,Type_Def,SaveCrt;

 Const
    XrIBM=31.10;   {Pixels by screen centimeter in x and y}
    YrIBM=31.10;   
    XrNEC=70.87;   {Pixels per cm for a NEC printer}
    YrNEC=70.87;
    XrEPS=142.9;   {Pixels per cm for EPSON Stylus color 760}
    YrEPS=142.9;
    XrHP=118.1;    {Pixels per cm for a HP laser 300 dpi printer}
    YrHP=118.1;

    XMini=0;   YMini=0;     Bord=5;

 Var cm_par_grad_x,cm_par_grad_y,win : integer;
     xmin,xmax,ymin,ymax,Cxmx,Cymx,Xratio,Yratio : real_ar;
     dx,dy,Echx,Echy,x0,y0,Xc,Yc,wl,xcm,ycm      : real_ar;
     Ixmn,Ixmx,Iymn,Iymx,xmin1,xmax1,ymin1,ymax1 : integer;
     MAXLIG   : integer;
     MAXCOL   : integer;
     ROND     : real_ar;
     CrtDC    : HDC;        { Handle of WinCrtMy window }
     fen,fen10,fen11,id_imprim : boolean;

    Function Log10 (x : real_ar) : real_ar;

    Function Power(x:real_ar; n:integer): real_ar;

    Procedure Fenetre(P: HDC; num: word);

    Procedure PleinEcran;

    Procedure TracerLesAxes(P:HDC);

    Procedure GraduerlesAxes(P:HDC);

    Procedure InitFenetre(P,fntr:integer;xmn,xmx,ymn,ymx:real_ar);

    Procedure MoveXY(P:HDC; x,y: real_ar);

    Procedure LineXY(P:HDC; x,y: real_ar);

    Procedure TextXY(P:HDC;x,y:real_ar;text:PChar);

    Procedure CroixXY(P:HDC;x,y:real_ar);

    Procedure Legendes (P:HDC; titre,titrex,titrey:Pchar);

    Procedure MinMax(n:integer;Y:RV; VAR ymin,ymax:real_ar);

    Procedure CourbeXY(P:HDC;n,fntr: integer; Y: RV; xn,xm: real_ar);

    Procedure TracerXY(P:HDC;n: integer; Y: RV; xn,xm: real_ar);

    Procedure Circle(P:HDC;xc,yc,r: real_ar; trait: boolean);

    Procedure WinCrtInit(Nom:PChar);

    Procedure SortieGraphique;


 IMPLEMENTATION


    Function Log10 (x : real_ar) : real_ar;
    begin
      if x>0 then
        log10:= ln(x)/ln(10)
      else
        log10:= -1E12
    end;

    {x power n}
    Function Power(x:real_ar; n:integer): real_ar;
    var i,m : integer;
        result :real_ar;
    begin
      result := 1.0;
      if n=0 then
      begin
        Power:=result;
        exit;
      end;
      m:=  n;
      if n<0 then m:=-n;
      for i:=1 to m do result :=x*result;
      Power :=result;
      if n<0 then Power:=1.0/result;
    end;

    {open one of 11 possible partial windows inside current calling program window}
    Procedure Fenetre(P: HDC; num: word);
    begin
      if (num < 1) or (num >11) then exit;
      Case num of
      1: begin          {upper left quarter}
           Ixmn:=Round(0.01*(MaxX+1));
           Ixmx:=Round(0.49*(MaxX+1));
           Iymn:=Round(0.01*(MaxY+1));
           Iymx:=Round(0.44*(MaxY+1));
         end;
      2: begin          {upper right quarter}
           Ixmn:=Round(0.50*(MaxX+1));
           Ixmx:=Round(0.97*(MaxX+1));
           Iymn:=Round(0.01*(MaxY+1));
           Iymx:=Round(0.44*(MaxY+1));
         end;
      3: begin          {lower half}
           Ixmn:=Round(0.01*(MaxX+1));
           Ixmx:=Round(0.97*(MaxX+1));
           Iymn:=Round(0.45*(MaxY+1));
           Iymx:=Round(0.85*(MaxY+1));
         end;
      4: begin          {all current window with internal graduations and title}
           Ixmn:=Round(0.01*(MaxX+1));
           Ixmx:=Round(0.97*(MaxX+1));
           Iymn:=Round(0.01*(MaxY+1));
           Iymx:=Round(0.89*(MaxY+1));
         end;
      5: begin          {upper half}
           Ixmn:=Round(0.01*(MaxX+1));
           Ixmx:=Round(0.97*(MaxX+1));
           Iymn:=Round(0.02*(MaxY+1));
           Iymx:=Round(0.43*(MaxY+1));
         end;
      6: begin          {lower left quarter}
           Ixmn:=Round(0.01*(MaxX+1));
           Ixmx:=Round(0.49*(MaxX+1));
           Iymn:=Round(0.45*(MaxY+1));
           Iymx:=Round(0.89*(MaxY+1));
         end;
      7: begin          {lower right quarter}
           Ixmn:=Round(0.50*(MaxX+1));
           Ixmx:=Round(0.97*(MaxX+1));
           Iymn:=Round(0.45*(MaxY+1));
           Iymx:=Round(0.89*(MaxY+1));
         end;
      8: begin          {left half}
           Ixmn:=Round(0.01*(MaxX+1));
           Ixmx:=Round(0.49*(MaxX+1));
           Iymn:=Round(0.01*(MaxY+1));
           Iymx:=Round(0.84*(MaxY+1));
         end;
      9: begin          {right half}
           Ixmn:=Round(0.50*(MaxX+1));
           Ixmx:=Round(0.97*(MaxX+1));
           Iymn:=Round(0.01*(MaxY+1));
           Iymx:=Round(0.84*(MaxY+1));
         end;
     10: begin          {screen with external graduations and upper title}
           Ixmn:=Round(0.08*(MaxX+1));
           Ixmx:=Round(0.98*(MaxX+1));
	   Iymn:=Round(0.07*(MaxY+1));
           Iymx:=Round(0.80*(MaxY+1));
           fen10:=true;
         end;
     11: begin          {same as number 10 adapted to HP laser printer}
           Ixmn:=Round(0.08*(MaxX+1));
           Ixmx:=MaxX;
           Iymn:=Round(0.10*(MaxX+1));
           Iymx:=Round(0.92*(MaxY+1));
           fen11:=true;
         end;

      end;
      {Draw a window frame}
      Rectangle(P,Ixmn,Iymn,Ixmx,Iymx+3);
      fen:=true;
      if num<>10 then fen10:=false;
    end;

    {Return to most used number 10 window}
    Procedure PleinEcran;
    begin
      Ixmn:=0;
      Ixmx:=MaxX;
      Iymn:=0;
      Iymx:=MaxY;
      fen:=False;
      fen10:=True
    end;

    {*********************************************************
    * You can choose a scale factor among to sets of values: *
    * [10, 5, 2.5, 2, 1.25, 1] or [10, 5, 2.5, 1].           *
    * ------------------------------------------------------ *
    * NOTE:                                                  *
    * The price to pay for allways centimetric graduations   *
    * is that the curve may not fill all available surface.  * 
    *********************************************************}
    PROCEDURE Echelle (VAR ech :real_ar);
    LABEL fin;
    CONST nval = 6;
          wval :  ARRAY[1..nval] OF real_ar
		     = (10.0,5.0,2.5,2.0,1.25,1.0);
          nval1 = 4;
          wval1:  ARRAY[1..nval1] OF real_ar
                     = (10.0,5.0,2.5,1.0);
    VAR
          i,ixpo              : integer;
          wlp,wlm,LogEch      : real_ar;
    BEGIN
      wlp:=wval[1];
      LogEch:=Log10(ech);
      ixpo := Trunc(LogEch);
      IF(LogEch<0) THEN ixpo:=ixpo-1;
      wl   := Ech/Power(10,ixpo);
      if (Not Log_X) and (Not Log_Y) then
        FOR i:= 1 TO nval-1 DO
        BEGIN
          wlm  := wlp;
          wlp  := wval[i+1];
          IF (wl-wlm)*(wl-wlp) <= 0 THEN GOTO fin;
        END
      else
        FOR i:= 1 TO nval1-1 DO
        BEGIN
          wlm  := wlp;
          wlp  := wval1[i+1];
          IF (wl-wlm)*(wl-wlp) <= 0 THEN GOTO fin;
        END;
      wlm  := 1;
 fin: wl   := wlm;
      ech  := wl * Power(10,ixpo);
    END;           {Echelle}

   {**********************************************
    Find scale for X axis                        }
    Procedure EchelleX(Var Echx:real_ar);
    Var ech : real_ar;
    Begin
      dx:=xmax-xmin;
      if dx <= 0 then Exit
                 else
                 begin
                   Ech:=dx/Cxmx;
                   Echelle(Ech);
                 end;
      Echx:=Ech;
    End;            

   {**********************************************
    Find scale for Y axis                        }
    Procedure EchelleY(Var Echy:real_ar);
    Var ech :real_ar;
    Begin
      dy := ymax- ymin;
      if dy <= 0 then exit
                 else
                 begin
                   Ech:=dy/Cymx;
                   Echelle(Ech);
                 end;
      Echy:=Ech;
    End;            

 {*************************************************************
    Conversion of physical coordinates into integer pixels in
    window defined by procedure InitFenetre. Used by MoveXY and
    LineXY.                                                   }  
    Procedure Conversion(xx,yy:real_ar; var Ix,Iy:integer);
    var xcm,ycm: real_ar;
    begin
      xcm:=x0; ycm:=y0;
      if Not Log_X then xcm:=x0+(xx-xmin)/Echx
                   else if abs(xx) > macheps then xcm:=x0+(log10(abs(xx))-xmin)/Echx;
      if Not Log_Y then ycm:=y0+(yy-ymin)/Echy
                   else if abs(yy) > macheps then ycm:=y0+(log10(abs(yy))-ymin)/Echy;
      Ix:=Ixmn+Round(XRatio*xcm);
      Iy:=Iymx-Round(YRatio*ycm)
    end;

  {************************************************************
    Draw X and Y axis                                         }
    PROCEDURE TracerLesAxes(P:HDC);
    VAR ix1,iy1  : integer;
        xcm,ycm  : real_ar;

        procedure Line(x1,y1,x2,y2:integer);
        begin
          MoveTo(P,x1,y1);
          LineTo(P,x2,y2);
        end;

    BEGIN
      xcm := x0-Xmin/Echx;
      ix1 := Ixmn +Round(Xratio*xcm);
      ycm := y0-Ymin/Echy;
      iy1 := Iymx -Round(Yratio*ycm);
      if (ix1>Ixmn) and (ix1<Ixmx) then Line(ix1,Iymn,ix1,Iymx);
      if (iy1>Iymn) and (iy1<Iymx) then Line(Ixmn,iy1,Ixmx,iy1);
    END; (* Tracer_les_Axes *)


    {***********************************************************
    * Adjustment of the number of digits displayed for gradua- *
    * tion of axis.                                            *
    * -------------------------------------------------------- *
    *                  Version dated 02/15/1987 by J-P Dumont  *
    ***********************************************************}
    Procedure Ajuster_Format(xs,dx_par_graduation,xmax:real_ar;Var n,m:Integer);
    Var test_zero,test_precision  :real_ar;
        m0,ng                     :Integer;
    begin
      if wl=7.5 then m0:=2 else m0:=1;    {2 digits for scale in 7.5}
      if(Abs(xs)<0.000001*Abs(xmax)) then test_zero:=0
                                          else test_zero:=Log10(Abs(xs));
      if(Test_zero=0) then
      begin
        n:=1;
        m:=0;
        exit;
      end;
      if (test_zero<0) then ng:=1
                       else ng:=Trunc(test_zero)+1;
      test_precision:= Log10(Abs(dx_par_graduation));
      if (test_precision>=0) then m:=0
                             else m:=m0+Trunc(Abs(test_precision)-0.01);
      n:=ng+m+2;       {2 for sign and dot}
    end;

    {*******************************************************************
    Graduate axis                                                      }
    Procedure GraduerLesAxes(P:HDC);
    Const
          grad_x    =4;
          grad_y    =4;
          longueur =16;

    Var   nx,comptx, ny,compty,sous_grad,i,ii,imin,imax,jj,
          ixxg,ixxg0,ixxg1,ixxgmx,iyyg,iyyg0,iyyg1,iyygmx,
          npar,mpar,trait,colonne,ligne,depart,arrivee,raccord: integer;
          xcm,ycm,xg,yg,test,Xlim,Ylim,dx_par_graduation,
          dy_par_graduation :real_ar;
          mot : string[16];
          mot1,mot2,flag:string[1];
          ChNum : array[0..longueur] of char;
          { Xlim,Ylim = physical coord. of left bottom corner of window}
          RedPen,OldPen : HPen;

          Procedure vidage;
          {Empty ChNum string}
          var i:word;
          begin
            for i:=0 to longueur do ChNum[i]:=#0
          end;

    begin
      dx_par_graduation:=Echx*cm_par_grad_x;
      dy_par_graduation:=Echy*cm_par_grad_y;

      RedPen:=CreatePen(ps_Solid,1,RGB(255,0,0));       {red pen thickness=1}

    if Not Log_X then   {linear graduations for X axis}
    begin
      Xlim:=xmin-x0*Echx;
      Ylim:=ymin-y0*Echy;
      nx:= Round(Xlim/Echx)-1;
      comptx := nx mod cm_par_grad_x;
      if comptx<0 then comptx:= cm_par_grad_x + comptx;
      Repeat
        comptx := comptx +1;
        if comptx = cm_par_grad_x then comptx :=0;
        Inc(nx);
        xg := nx *Echx;
        xcm := x0 +(xg-xmin)/Echx;
        ixxg := Ixmn+Round(xratio * xcm);
        trait :=grad_x;
        if comptx=0 then trait := 2*grad_x;
        If ixxg > Ixmn-1 then
        begin
          raccord:=1;
          MoveTo(P,ixxg,iymx+raccord);
          LineTo(P,ixxg,iymx+raccord-trait);
          If comptx =0 then
          begin
            Ajuster_Format(xg,dx_par_graduation,xmax,npar,mpar);
            Str(xg:npar:mpar,mot);
            if (abs(xg-1)<0.0001) then mot:='  1';
            if (abs(xg+1)<0.0001) then mot:=' -1';
            if (abs(xg-1000)<0.01) then mot:='  1000';
            if (abs(xg+1000)<0.01) then mot:=' -1000';
            if (abs(xg)<0.0001) then mot:='   0';
            vidage; { effacement de ChNum }
            for jj:=1 to length(mot) do ChNum[jj-1]:=mot[jj];
            if fen10 then
              TextOut(P,ixxg-10,iymx+5,ChNum,strlen(ChNum))
            else if fen11 then
              TextOut(P,ixxg-50,iymx+15,ChNum,strlen(ChNum))
            else
              TextOut(P,ixxg-10,iymx-30,ChNum,strlen(ChNum))
          end;
          MoveTo(P,ixxg,iymn);
          LineTo(P,ixxg,iymn+trait);

          OldPen:=SelectObject(P,RedPen);  {select blue pen}
          for i:=1 to (iymx+2-iymn) div 2 do   {draw vertical grid}
          begin
          jj:=2*i;
          MoveTo(P,ixxg,Iymn+jj);
          LineTo(P,ixxg,Iymn+jj+1);
          end;
          SelectObject(P,OldPen)
        end
      until xcm > Cxmx-0.95;
    end  { Fin Log_X = faux }
    else
    begin          {Log. graduations for X axis}
      Xlim:=xmin-x0*Echx;
      Ylim:=ymin-y0*Echy;
      nx:= Round(Xlim/Echx)-1;
      comptx := nx mod cm_par_grad_x;
      if comptx<0 then comptx:= cm_par_grad_x + comptx;
      flag:='*';
      Repeat
        comptx := comptx +1;
        sous_grad:=sous_grad+1;
        if comptx = cm_par_grad_x then begin
        comptx :=0; sous_grad:=0; mot1:=' '; flag:=' '; end;
        nx := nx +1;
        xg := nx *Echx;
        xcm := x0 +(xg-xmin)/Echx;
        ixxg := ixmn +Round(xratio * xcm);
        if comptx=0 then ixxg0:=ixxg;
        if sous_grad=1 then
        begin
          ixxg:=ixxg0+Round((ixxg-ixxg0)*1.204);
          mot1:='2'; mot2:='3';
          ixxg1:=ixxg0+Round((ixxg-ixxg0)*1.585);
        end;
        if sous_grad=2 then
        begin
          ixxg:=ixxg0+Round((ixxg-ixxg0)*1.364);
          mot1:='5'
        end;
        if sous_grad=3 then
        begin
          ixxg:=ixxg0+Round((ixxg-ixxg0)*1.14);
          mot1:='7'
        end;
        trait :=grad_x;
        if comptx=0 then trait := 2*grad_x;
        ixxgmx:=Ixmx-25;
        If (flag<>'*') and (ixxg<ixxgmx) then
        begin
          raccord:=1;
          MoveTo(P,ixxg,Iymx+raccord);
          LineTo(P,ixxg,Iymx+raccord-trait);
          vidage; { effacement de ChNum }
          for jj:=1 to length(mot1) do ChNum[jj-1]:=mot1[jj];
          TextOut(P,ixxg-10,Iymx-15,ChNum,strlen(ChNum));

          If sous_grad=1 then
          begin
            MoveTo(P,ixxg1,Iymx+raccord);
            LineTo(P,ixxg1,Iymx+raccord-trait);
            vidage; { effacement de ChNum }
            for jj:=1 to length(mot2) do ChNum[jj-1]:=mot2[jj];
            TextOut(P,ixxg1-10,Iymx-15,ChNum,strlen(ChNum))
          end;
          if comptx =0 then
          begin
            If Abs(Round(xg)-(xg))<0.001 then
            begin
              xg:=Power(10,Round(xg));
              Ajuster_Format(xg,dx_par_graduation,xmax,npar,mpar);
              If xg>1 then mpar:=0;
              If xg<1 then mpar:=Round(log10(1.0/xg));
              If xg>0.5 then npar:=npar-1;
              Str(xg:npar:mpar,mot);
              if (abs(xg-1)<0.0001) then mot:=' 1';
              vidage; { effacement de ChNum }
              for jj:=1 to length(mot) do ChNum[jj-1]:=mot[jj];
              if not fen10 then
                TextOut(P,ixxg-10,iymx-30,ChNum,strlen(ChNum))
              else
                TextOut(P,ixxg-10,iymx+5,ChNum,strlen(ChNum));
            end;
          end;
          raccord:=1;
          MoveTo(P,ixxg,Iymn+raccord);
          LineTo(P,ixxg,Iymn+raccord+trait);
          for i:=1 to (iymx+2-iymn) div 2 do   {draw vertical grid}
          begin
            jj:=2*i;
            MoveTo(P,ixxg,Iymn+jj);
            LineTo(P,ixxg,Iymn+jj+1);
          end;

          If sous_grad=1 then
          begin
            MoveTo(P,ixxg1,Iymn+raccord);
            LineTo(P,ixxg1,Iymn+raccord+trait);
            for i:=1 to (iymx+2-iymn) div 2 do   { additional vertical grid}
            begin                                
              jj:=2*i;
              MoveTo(P,ixxg1,Iymn+jj);
              LineTo(P,ixxg1,Iymn+jj+1);
            end;
         end;
        end;
      until xcm > Cxmx
    end; {else}

    if Not Log_Y then   {Linear graduations for Y axis}
    begin
      ny:= Round(Ylim/Echy)-1;
      compty := ny mod cm_par_grad_y;
      if compty<0 then compty:=cm_par_grad_y+compty;
      Repeat
        compty := compty +1;
        if compty = cm_par_grad_y then compty :=0;
        ny := ny +1;
        yg := ny *Echy;
        ycm := y0 +(yg-ymin)/Echy;
        iyyg := Iymx-Round(Yratio * ycm);
        trait :=grad_y;
        if compty = 0 then trait:=2*grad_y;
        if iyyg < Iymx then
        begin
          raccord:=1;
          MoveTo(P,ixmn+raccord,iyyg);
          LineTo(P,ixmn+raccord+trait,iyyg);
          if compty=0 then
          begin
            Ajuster_Format(yg,dy_par_graduation,ymax,npar,mpar);
            Str(yg:npar:mpar,mot);
            if (abs(yg-1)<0.0001) then mot:='  1';
            if (abs(yg+1)<0.0001) then mot:=' -1';
            if (abs(yg-1000)<0.01) then mot:='  1000';
            if (abs(yg+1000)<0.01) then mot:=' -1000';
            if (abs(yg)<0.0001) then mot:='  0';
            vidage; { effacement de ChNum }
            for jj:=1 to length(mot) do ChNum[jj-1]:=mot[jj];
            if fen10 then
              TextOut(P,5,iyyg-10,ChNum,strlen(ChNum))
            else if fen11 then
              TextOut(P,0,iyyg-10,ChNum,strlen(ChNum))
            else
              TextOut(P,Ixmn+10,iyyg-10,ChNum,strlen(ChNum))
          end;
          MoveTo(P,ixmx,iyyg);
          LineTo(P,ixmx-trait,iyyg);

          OldPen:=SelectObject(P,RedPen);  {select blue pen}
          for i:=1 to (ixmx-ixmn) div 3 do  {draw horizontal grid}
          begin
            jj:=3*i;
            MoveTo(P,Ixmn+jj,iyyg);
            LineTo(P,Ixmn+jj+1,iyyg);
          end;
          SelectObject(P,OldPen);
        end;
      until ycm > Cymx-1.0
      {Until iyyg < Iymn - 10  }
    end
    else   {Log. graduations for Y axis}
    begin
      dy_par_graduation:=Echy*cm_par_grad_x;
      ny:= Round(Ylim/Echy)-1;
      compty := ny mod cm_par_grad_x;
      if compty<0 then inc(compty,cm_par_grad_x);
      flag:='*';
      Repeat
        inc(compty); inc(sous_grad); inc(ny);
        yg := ny *Echy;
        ycm := y0 +(yg-ymin)/Echy;
        iyyg := iymx -Round(Yratio * ycm);
        if compty = cm_par_grad_x then
        begin
          compty :=0; sous_grad:=0; mot1:=' '; flag:=' '
        end;
        if compty=0 then iyyg0:=iyyg;
        if sous_grad=1 then
        begin
          iyyg:=iyyg0+Round((iyyg-iyyg0)*1.204);
          mot1:='2'; mot2:='3';
          iyyg1:=iyyg0+Round((iyyg-iyyg0)*1.585);
        end;
        if sous_grad=2 then begin
          iyyg:=iyyg0+Round((iyyg-iyyg0)*1.364);
          mot1:='5'
        end;
        if sous_grad=3 then begin
          iyyg:=iyyg0+Round((iyyg-iyyg0)*1.14);
          mot1:='7'
        end;
        trait :=grad_y;
        if compty=0 then trait := 2*grad_y;
        iyygmx:=Iymx-10;
        if (iyyg > 15) and (flag<>'*') and (iyyg < iyygmx) then
        begin
          raccord:=1;
          MoveTo(P,Ixmn+raccord,iyyg);
          LineTo(P,Ixmn+raccord+trait,iyyg);
          vidage; { effacement de ChNum }
          for jj:=1 to length(mot1) do ChNum[jj-1]:=mot1[jj];
          TextOut(P,Ixmn+10,iyyg-10,ChNum,strlen(ChNum));
          If sous_grad=1 then
          begin
            MoveTo(P,Ixmn+raccord,iyyg1);
            LineTo(P,Ixmn+raccord+trait,iyyg1);
            vidage; { effacement de ChNum }
            for jj:=1 to length(mot2) do ChNum[jj-1]:=mot2[jj];
            TextOut(P,Ixmn+10,iyyg1-10,ChNum,strlen(ChNum))
          end;
          if compty=0 then
          begin
            If Abs(Round(yg)-(yg))<0.001 then
            begin
              yg:=Power(10,Round(yg));
              Ajuster_Format(yg,dy_par_graduation,ymax,npar,mpar);
              If yg>1 then mpar:=0;
              If yg<1 then mpar:=Round(log10(1.0/yg));
              If yg>0.5 then npar:=npar-1
                        else npar:=npar+1;
              Str(yg:npar:mpar,mot);
              if (abs(yg-1)<0.0001) then mot:='  1';
              vidage; { effacement de ChNum }
              for jj:=1 to length(mot) do ChNum[jj-1]:=mot[jj];
              if not fen10 then
                TextOut(P,Ixmn+10,iyyg-10,ChNum,strlen(ChNum))
              else
                TextOut(P,5,iyyg-10,ChNum,strlen(ChNum));
            end;
          end;
          MoveTo(P,ixmx,iyyg);
          LineTo(P,ixmx-trait,iyyg);
          for i:=1 to (ixmx-ixmn) div 3 do  {draw horizontal grid}
          begin
            jj:=3*i;
            MoveTo(P,Ixmn+jj,iyyg);
            LineTo(P,Ixmn+jj+1,iyyg);
          end;
          If sous_grad=1 then
          begin
            MoveTo(P,ixmx,iyyg1);
            LineTo(P,ixmx-trait,iyyg1);

            for i:=1 to (ixmx-ixmn) div 3 do  {additional horizontal grid}
            begin
              jj:=3*i;
              MoveTo(P,Ixmn+jj,iyyg1);
              LineTo(P,Ixmn+jj+1,iyyg1);
            end;
          end;
        end;
      until ycm>Cymx-2
    end; {else}

    DeleteObject(RedPen)

    End;     {GraduerLesAxes}

   {******************************************************************************
    Open virtual window in physical coordinates                                  }
    Procedure InitFenetre(P,fntr:integer;xmn,xmx,ymn,ymx:real_ar);
    begin

      Fenetre(P,fntr);     {call one of the 11 predefined windows}

      {centimeters in X and Y available for drawing}
      Cxmx:=(Ixmx-Ixmn)/Xratio;
      Cymx:=(Iymx-Iymn)/Yratio; 

      if Ech_auto then  {automatic scaling}
      begin
        xmin:=xmn; xmax:=xmx; ymin:=ymn; ymax:=ymx
      end
      else     {X_mini,X_maxi,Y_mini,Y_maxi:     }
      begin    {to be defined by calling program!}
        xmin:=X_mini; xmax:=X_maxi; ymin:=Y_mini; ymax:=Y_maxi
      end;    

      if Log_X then  {Log. scaling for X axis}
      begin
        if (Ixmx-Ixmn) > MaxX div 2 then xmin:=xmax/10000.0
                                    else xmin:=xmax/100.0;
        xmax:=Round(Log10(xmax));
        xmin:=Round(Log10(xmin))
      end;

      if Log_Y then  {Log. scaling for Y axis}
      begin
        if (Iymx-Iymn) > MaxY div 2 then ymin:=ymax/1000.0
                                    else ymin:=ymax/10.0;
        ymax:=Round(Log10(ymax));
        ymin:=Round(Log10(ymin))
      end;

      {Define scales in X and Y}
      EchelleX(Echx);
      EchelleY(Echy);

      {Determine upper left point (Xc,Yc) of drawing}
      x0 := 0.5 * (Cxmx - dx/Echx);
      if xmin=0 then
      begin
        Xc:=xmin-x0*Echx;
        if ( Xc<0 ) then x0:=xmin/Echx;
      end;
      y0 := 0.5 * (Cymx - dy/Echy);
      if ymin=0 then
      begin
        Yc:=ymin-y0*Echy;
        if (Yc<0 ) then y0:=ymin/Echy;
      end;

      {Graduate axis}
      GraduerLesAxes(P);
      {Draw axis}
      TracerLesAxes(P);
    end;                  {InitFenetre}

 {*************************************************************
    Move pen to current physical point (x,y)                  }
    Procedure MoveXY(P:HDC; x,y: real_ar);
    var ix,iy: integer;
    begin
      if Not Log_X then xcm:=x0+(x-xmin)/Echx
                   else if abs(x) > macheps then xcm:=x0+(log10(abs(x))-xmin)/Echx;
      if Not Log_Y then ycm:=y0+(y-ymin)/Echy
                   else if abs(y) > macheps then ycm:=y0+(log10(abs(y))-ymin)/Echy;
      Ix:=Ixmn+Round(XRatio*xcm);
      Iy:=Iymx-Round(YRatio*ycm);
      if (ix>Ixmn-1) and (ix<Ixmx+1) and (iy>Iymn-1) and (iy<Iymx+1) then
        MoveTo(P,ix,iy)
    end;

 {*******************************************************************
    Draw a line from current physical point to physical point (x,y) }
    Procedure LineXY(P:HDC; x,y: real_ar);
    var ix,iy: integer;                        
    begin
      if Not Log_X then xcm:=x0+(x-xmin)/Echx
                   else if abs(x) > macheps then xcm:=x0+(log10(abs(x))-xmin)/Echx;
      if Not Log_Y then ycm:=y0+(y-ymin)/Echy
                   else if abs(y) > macheps then ycm:=y0+(log10(abs(y))-ymin)/Echy;
      ix:=Ixmn+Round(XRatio*xcm);
      iy:=Iymx-Round(YRatio*ycm);
      if (ix>Ixmn-1) and (ix<Ixmx+1) and (iy>Iymn-1) and (iy<Iymx+1) then
        LineTo(P,ix,iy)
    end;

 {********************************************************************************
    Write a text at physical point (x,y)                                         }
    Procedure TextXY;
    var ix,iy: integer;
    begin
      if Not Log_X then xcm:=x0+(x-xmin)/Echx
                   else if abs(x) > macheps then xcm:=x0+(log10(abs(x))-xmin)/Echx;
      if Not Log_Y then ycm:=y0+(y-ymin)/Echy
                   else if abs(y) > macheps then ycm:=y0+(log10(abs(y))-ymin)/Echy;
      ix:=Ixmn+Round(XRatio*xcm);
      iy:=Iymx-Round(YRatio*ycm);
      if (ix>Ixmn-1) and (ix<Ixmx+1) and (iy>Iymn-1) and (iy<Iymx+1) then
        TextOut(P,ix,iy,text,strlen(text))
    end;

    {Draw a cross at physical point (x,y)                                         }
    Procedure CroixXY;
    var ix,iy: integer;
    begin
      Conversion(x,y,ix,iy);
      MoveTo(P,ix-5,iy); LineTo(P,ix+5,iy);
      MoveTo(P,ix,iy-5); LineTo(P,ix,iy+5)
    end;
 {********************************************************************************
  Write graph caption, names of X and Y axis                                       }
  Procedure Legendes;
  VAR ic,ix,lh: word;  rx: real;
      OldFont,TitleFont : HFont;
      TitleLogFont: TLogFont;
  BEGIN
    {Define title font}
    IF Fen11 THEN BEGIN lh:=100; ic:=50 END
    ELSE IF Id_imprim THEN BEGIN lh:=50; ic:=25 END
    ELSE BEGIN lh:=25; ic:=12 END;  {case screen}
    WITH TitleLogFont DO  
    BEGIN
      lfHeight:=lh;
      lfwidth:=0;
      lfEscapement:=0;
      lfOrientation:=0;
      lfWeight:=fw_UltraBold;
      lfItalic:=1;
      lfUnderline:=0;
      lfStrikeOut:=0;
      lfcharSet:=ANSI_CharSet;
      lfOutPrecision:=Out_Default_Precis;
      lfClipPrecision:=Clip_Default_Precis;
      lfQuality:=Default_Quality;
      lfPitchAndFamily:=Default_Pitch OR (ff_DontCare SHL 4);
      StrCopy(@lfFaceName,'Police Titre');
    END;
    {Draw centered graph title}
    TitleFont:=CreateFontIndirect(TitleLogFont);
    OldFont:=SelectObject(P,TitleFont);
    ix:=(Ixmx-Ixmn-ic*strlen(titre)) DIV 2;
    if Ixmx-Ixmn < MaxX DIV 2 then ix:=ix+40;
    {small windows: title shifted 40 pixels to the right}
    IF Fen10 THEN       {screen window n° 10}
      TextOut(P,Ixmn+ix,1,titre,strlen(titre))
    ELSE IF Fen11 THEN  {printer window n° 11 }
      TextOut(P,Ixmn+(ix DIV 2),Round(0.5*Iymn),titre,strlen(titre))
    ELSE                {other cases}
      TextOut(P,Ixmn+ix,Iymn+10,titre,strlen(titre));
    {back to normal font}
    SelectObject(P,OldFont);
    DeleteObject(TitleFont);

    {Write X axis name}
    if Ixmx-Ixmn > MaxX DIV 2 then rx:=0.85 else rx:=0.70;
    TextOut(P,Ixmn+Round(rx*(Ixmx-Ixmn)),Iymn+Round(0.90*(Iymx-Iymn)),
            titreX,strlen(titreX));
     {Write Y axis name}
    TextOut(P,Ixmn+Round(0.05*(Ixmx-Ixmn)),Iymn+Round(0.03*(Iymx-Iymn)),
            titreY,strlen(titreY));
  END;   {Legendes}

 {**********************************************************
  Searching minimum and maximum of a dynamic Table Y       } 
  Procedure MinMax;
  Var i: integer;
  Begin
    ymin:=Y^[1]; ymax:=ymin;
    for i:=2 to n do
    begin
      if Y^[i] < ymin then ymin:=Y^[i];
      if Y^[i] > ymax then ymax:=Y^[i]
    end
  End;

  {*********************************************************
  * Draw a signal made of n tabulated values (n<=2048) to  *
  * device context P, with a call to InitFenetre.          *
  * ------------------------------------------------------ *
  * INPUTS:  n = number of signal points.                  *
  *          fntr = predefined window number (1 to 11)     *
  *          Y = table  of n values to draw (REAL type)    * 
  *          xn, xm = values of begin and end abscissas    *
  *          X sampling is constant.                       *
  *                        --------                        *
  * OUTPUTS: Display signal in desired window.             *
  * ------------------------------------------------------ *
  * Nota :   To superpose several curves, use TracerXY to  *
  *          draw the other curves.                        *
  *********************************************************}
  Procedure CourbeXY;
  Var i: integer;
      x,dx,ymin,ymax: real_ar;
  Begin
    ymin:=Y^[1]; ymax:=ymin;
    for i:=2 to n do
    begin
      if Y^[i] < ymin then ymin:=Y^[i];
      if Y^[i] > ymax then ymax:=Y^[i]
    end;
    dx:=(xm-xn)/(n-1);
    InitFenetre(P,fntr,xn,xm,ymin,ymax);
    x:=xn;
    MoveXY(P,x,Y^[1]);
    for i:=2 to n do
    begin
      x:=x+dx;
      LineXY(P,x,Y^[i])
    end
  End;

  {*********************************************************
  * Draw a signal made of n tabulated values (n<=2048) to  *
  * device context P, without a call to InitFenetre.       *
  * ------------------------------------------------------ *
  * INPUTS:  n = number of signal points.                  *
  *          Y = table  of n values to draw (REAL type)    * 
  *          xn, xm = values of begin and end abscissas    *
  *          X sampling is constant.                       *
  *                        --------                        *
  * OUTPUTS: Display signal in desired window.             *
  * ------------------------------------------------------ *
  * Nota :   To draw a single (or the first) curve, USE    *
  *          CourbeXY with a call to InitFenetre.          *
  *********************************************************}
  Procedure TracerXY;  
  Var i   : integer;
      x,dx: real_ar;
  Begin
    dx:=(xm-xn)/(n-1);
    x:=xn;
    MoveXY(P,x,Y^[1]);
    for i:=2 to n do
    begin
      x:=x+dx;
      LineXY(P,x,Y^[i])
    end
  End;

  {*******************************************
  * algorithm to draw a circle in physical   *
  *  coordinates (dotted line or normal line *
  * ---------------------------------------- *
  * INPUTS:                                  *
  *          xc,yc : center of circle        *
  *              r : radius of circle        *
  *                                          *
  *          trait : TRUE for normal line    *
  *                  FALSE for dotted line   * 
  *******************************************}
  PROCEDURE Circle(P:HDC;xc,yc,r: real_ar; trait: boolean);
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


  {**********************************************************
  * Open a CRT window with caption with possibility of text *
  * and/or graph with the condition to uses WinCrt1 unit    *
  * instead of WinCrt (else CrtWindow is not visible).      *
  **********************************************************}
  PROCEDURE WinCrtInit;
  BEGIN
    WindowOrg.X:=200;    {upper left corner position}
    WindowOrg.Y:=200;    {and sizes of CRT window   }
    WindowSize.X:=600;
    WindowSize.Y:=500;
    StrCopy(WindowTitle,Nom);  {window title}
    MaxX:=590; MaxY:=475;      {client zone for drawing}
    InitWinCrt;                {call standard Borland procedure}
    CrtDC:=GetDC(CrtWindow);   {define device context of window}
  END;

  {options to exit graph}
  PROCEDURE SortieGraphique;
  VAR s:string; c,rep: char;
  BEGIN
    REPEAT
      gotoxy(2,28); clreol;
      gotoxy(15,29);
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

{initialize unit} 
BEGIN
  MaxX:=800; MaxY:=600;          {SVGA by default}
  XRatio:=XrIBM; YRatio:=YrIBM;  
  cm_par_grad_x:=4;
  cm_par_grad_y:=2;
  PleinEcran                     {window n° 10 by default}
END. {of unit wgraph_2d.pas}
