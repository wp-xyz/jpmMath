{*******************************************************************
*                          FRACTAL CURVES                          *
* ---------------------------------------------------------------- *
* This program generates a fractal curve using the formula:        *
*              x' = x(x2-x1) - y(y2-y1) + x1                       *
*              y' = x(y2-y1) + y(x2-x1) + y1                       *
* using a basis composed of segments (such as a single segment, a  *
* square, a triangle etc.) forming a closed figure, and a genera-  *
* trix simple open figure, such as two bits of line at 60°.        *
* The basis is memorized by its summit coordinates stored in two   *
* tables Xbasis and Ybasis. The generatrix is stored in the same   *
* way in tables Xgen and Ygen. Also generates a H Fractal.         *
* ---------------------------------------------------------------- *
* REFERENCE:                                                       *
* "Graphisme dans le plan et dans l'espace avec Turbo-Pascal 4.0   *
*  By R. Dony - MASSON, Paris 1990" [BIBLI 12].                    *
*                                                                  *
*                                TPW Release By J-P Moreau, Paris. *
*******************************************************************}
Program Fractal;

Uses WinCrtMy, WinTypes, WinProcs, Strings, Type_def, CrtGr2D;

Const limite = 4096;

Type Ptb = ^Tab1;
     Tab1 = Array[0..limite] of real_ar;

Var  Xbase,Ybase  : array[0..10] of real_ar;
     Xgen,Ygen    : array[0..20] of real_ar;
     Xfract,Yfract: Ptb;
     Example, Order, PointsBase, PointsGen: Byte;
     f1,f2,f3,f4: real_ar;
     NomCourbe: string[80];

     BluePen, RedPen: HPen;

    {specific HFractal}
    X1,Y1,X2,Y2,X3,Y3,X4,Y4: array[1..20] of Real_ar;
    A: Real_ar;
    M,N:Integer;

    {integer power base^exponant}
    function IPower(base,exponant:integer):word;
    var i,prod:word;
    begin
      prod:=1;
      for i:=1 to exponant do prod:=prod*base;
      IPower:=prod
    end;

    function Power(base:real_ar;exponant:integer):real_ar;
    var i:word; prod:real_ar;
    begin
      prod:=1.0;
      for i:=1 to exponant do prod:=prod*base;
      Power:=prod
    end;

    Procedure Data;
    begin
      Writeln;
      Write(' Input # Example (1 to 4): '); Readln(Example);
      if Example<1 then Example:=1;
      if Example>4 then Example:=4;
      Case Example of
      1: begin
         NomCourbe:='Von Koch''s Triangular Snowflakes';
         PointsBase:=3;
         XBase[0]:=1;       YBase[0]:=0;
         XBase[1]:=-0.5;    YBase[1]:=0.866025;
         XBase[2]:=-0.5;    YBase[2]:=-0.866025;
         XBase[3]:=1;       YBase[3]:=0;
         PointsGen:=5;
         XGen[0]:=0;        YGen[0]:=0;
         XGen[1]:=0.3333;   YGen[1]:=0;
         XGen[2]:=0.5;      YGen[2]:=0.2887;
         XGen[3]:=0.6666;   YGen[3]:=0;
         XGen[4]:=1;        YGen[4]:=0;
         f1:=-1.1; f2:=1.1; f3:=-1.1; f4:=1.1;
         Writeln;
         Write(' Input Fractal Order (1 to 6): '); Readln(order)
         end;
      2: begin
         NomCourbe:='Example #2';
         PointsBase:=5;
         XBase[0]:= 0.0;    YBase[0]:= 0.0;
         XBase[1]:= 0.0;    YBase[1]:= 1.0;
         XBase[2]:= 1.0;    YBase[2]:= 1.0;
         XBase[3]:= 1.0;    YBase[3]:= 0.0;
         XBase[4]:= 0.0;    YBase[4]:= 0.0;
         PointsGen:=3;
         XGen[0]:= 0.0;     YGen[0]:=0.0;
         XGen[1]:= 0.5;     YGen[1]:=0.5;
         XGen[2]:= 1.0;     YGen[2]:=0.0;
         f1:=-1.1; f2:=2.1; f3:=-1.1; f4:=2.25;
         Writeln;
         Write(' Input Fractal Order (1 to 12): '); Readln(order)
         end;
      3: begin
         NomCourbe:='Example #3';
         PointsBase:=5;
         XBase[0]:= 0.0;    YBase[0]:= 0.0;
         XBase[1]:= 0.0;    YBase[1]:= 1.0;
         XBase[2]:= 1.0;    YBase[2]:= 1.0;
         XBase[3]:= 1.0;    YBase[3]:= 0.0;
         XBase[4]:= 0.0;    YBase[4]:= 0.0;
         PointsGen:=5;
         XGen[0]:=0;        YGen[0]:=0;
         XGen[1]:=0.4;      YGen[1]:=0;
         XGen[2]:=0.5;      YGen[2]:=0.4;
         XGen[3]:=0.6;      YGen[3]:=0;
         XGen[4]:=1;        YGen[4]:=0;
         f1:=-0.5; f2:=1.5; f3:=-0.5; f4:=1.6;
         Writeln;
         Write(' Input Fractal Order (1 to 6): '); Readln(order)
         end;
      4: begin
         NomCourbe:='H Fractal';
         Writeln;
         Write(' Input fractal order (1 to 8): '); Readln(order);
         Fenetre(-1,1,-1,1);
         Cloture(10,MaxX-10,15,MaxY);
         end
      end
    end;

    procedure CalculateCurve(O,B,G: byte);
    var i,aux1,aux2,j,m1,m2,k: word;
	x1, y1: real_ar;
    begin
      Xfract^[0]:=0;
      Yfract^[0]:=0;
      Xfract^[IPower(G,O)]:=1;
      Yfract^[IPower(G,O)]:=0;
      for i:=0 to O-1 do
      begin
        aux1:=IPower(G,O-i);
        aux2:=IPower(G,O-i-1);
        j:=0;
        Repeat
          m1:=j+aux1;
          x1:=Xfract^[m1]-Xfract^[j];
          y1:=Yfract^[m1]-Yfract^[j];
          for k:=1 to G-1 do
          begin
            m2:=j+k*aux2;
            Xfract^[m2]:=x1*Xgen[k]-y1*Ygen[k]+Xfract^[j];
            Yfract^[m2]:=y1*Xgen[k]+x1*Ygen[k]+Yfract^[j];
          end;
          j:=j+aux1;
        Until j > IPower(G,O)-1
      end
    end;

    procedure DrawBase;
    var i: integer;
    begin
      BluePen:=CreatePen(ps_Solid,1,RGB(0,0,255));     {blue pen}
      RedPen :=CreatePen(ps_Solid,1,RGB(180,50,50));   {light red pen}
      SelectObject(CrtDC,RedPen);
      Fenetre(-2,2,-2,2);
      Cloture(10,Round(0.33*(MaxX+1)),Round(0.55*(MaxY+1)),MaxY-10);
      Bordure(CrtDc);
      MoveXY(CrtDc,Xbase[0],Ybase[0]);
      for i:=1 to PointsBase do LineXY(CrtDc,Xbase[i],Ybase[i]);
      TextOut(CrtDc,25,18,'Base',4)
    end;

    procedure DrawGeneratrix;
    const r=-0.1;
    var i: integer;
    begin
      Fenetre(-1,2,-1.5,1.5);
      Cloture(10,Round(0.33*(MaxX+1)),80,Round(0.53*(MaxY+1)));
      Bordure(CrtDc);
      MoveXY(CrtDc,r+Xgen[0],r+Ygen[0]);
      for i:=1 to PointsGen-1 do LineXY(CrtDc,r+Xgen[i],r+Ygen[i]);
      TextOut(CrtDc,25,Round(0.49*(MaxY+1)),'Generatrix',11)
    end;

    procedure DrawCurve(O,B,G: byte);
    var diff1,diff2,xx,yy : real_ar;
        ix,m,n : word;
        chaine: string;
        nom : array[0..60] of char;

    begin
      SelectObject(CrtDC,BluePen);
      Fenetre(f1,f2,f3,f4);
      Str(order,chaine);
      chaine:=Concat(NomCourbe,'         Order = ',chaine);
      StrPCopy(nom,chaine);
      ix:=130 + ((MaxX-8*strlen(nom)) DIV 2);
      TextOut(CrtDc,ix,18,nom,Strlen(nom));
      Cloture(Round(0.34*(MaxX+1)),MaxX-13,80,MaxY-10);
      Bordure(CrtDc);
      MoveXY(CrtDc,Xbase[0],Ybase[0]);
      for m:=0 to B-1 do
      begin
        diff1:=Xbase[m+1] - Xbase[m];
        diff2:=Ybase[m+1] - Ybase[m];
        for n:=0 to IPower(G,O) do
        begin
          xx:=diff1*Xfract^[n] - diff2*Yfract^[n] + Xbase[m];
          yy:=diff2*Xfract^[n] + diff1*Yfract^[n] + Ybase[m];
          LineXY(CrtDc,xx,yy)
        end
      end
    end;

    {called by DrawHFractal}
    Procedure DrawSegments(S:byte);
    Var j:integer;
        X,Y,B,C: Real_ar;
    Begin
      For j:=S to Order do
      begin
        X:=X1[j-1];
        Y:=Y1[j-1];
        B:=Power(A,j);
        C:=A*B*1.5;
        X1[j]:=X+B; Y1[j]:=Y+C;
        X2[j]:=X+B; Y2[j]:=Y-C;
        X3[j]:=X-B; Y3[j]:=Y+C;
        X4[j]:=X-B; Y4[j]:=Y-C;
        MoveXY(CrtDc,X-B,Y); LineXY(CrtDc,X+B,Y);
        MoveXY(CrtDc,X1[j],Y1[j]); LineXY(CrtDc,X2[j],Y2[j]);
        MoveXY(CrtDc,X3[j],Y3[j]); LineXY(CrtDc,X4[j],Y4[j])
      end
    End;

    Procedure DrawHFractal;  {option #4}
    Var S: byte;
        chaine: string;
        nom : array[0..60] of char;
    Begin
      A:=0.5;
      X1[1]:=0;
      Y1[1]:=0;
      S:=1;
      DrawSegments(S);
      For M:=1 to Round(Power(4,Order-1))-1 do
      begin
        N:=M;
        S:=Order;
        While (N Mod 4) = 0 do
        begin
          N:=N Div 4;
          S:=S-1
        end;
        X1[S-1]:=X2[S-1]; X2[S-1]:=X3[S-1];
        X3[S-1]:=X4[S-1]; Y1[S-1]:=Y2[S-1];
        Y2[S-1]:=Y3[S-1]; Y3[S-1]:=Y4[S-1];
        DrawSegments(S)
      end;
      Str(order,chaine);
      chaine:=Concat(NomCourbe,'           Order = ',chaine);
      StrPCopy(nom,chaine);
      TextOut(CrtDc,190,18,nom,Strlen(nom));
    End;

{main program}
BEGIN
  WinCrtInit(' Fractals');
  New(Xfract); New(Yfract);
  Repeat
    ClrScr;
    Data;
    ClrScr;
    if Example < 4 then
    begin
      CalculateCurve(Order, PointsBase, PointsGen-1);
      DrawBase;
      DrawGeneratrix;
      DrawCurve(Order, PointsBase, PointsGen-1)
    end
    else
      DrawHFractal;
    SortieGraphique
  Until rep='n';
  Dispose(XFract); Dispose(YFract);
  DoneWinCrt
END.

{end of file fractal.pas}