{*************************************************************************
*                       THE  HANOI  TOWERS                               *
* ---------------------------------------------------------------------- *
* Explanations:                                                          *
*   This famous game was invented by a French matematician, Lucas. who   *
* was then a teacher in St Louis College of Paris. Three vertical needles*
* are regularly spaced on a board. At the beginning, N disks lie on the  *
* left needle, the largest at the bottom and the other ones, smaller and *
* smaller up to the top of the heap. The game consists in moving all the *
* N disks, respecting the two following rules:                           *
*                                                                        *
*          1) move only one disk at a time,                              *
*          2) never put a disk on a smaller one.                         *
*                                                                        *
* We can prove that the number of necessary disk shifts is 2^N-1.        *
* This program using recursivity, shows the disk shifts for a number N   *
* of disks given by the user.                                            *
* ---------------------------------------------------------------------- *
* Reference:                                                             *
* From "Graphisme dans le plan et dans l'espace avec Turbo Pascal 4.0    *
*       By R. Dony - MASSON 1990 page 113" [BIBLI 12].                   *
*                                                                        *
*                                    TPW Release By J-P Moreau, Paris.   *
*                                            (www.jpmoreau.fr)           *
*************************************************************************}
Program Hanoi;

Uses WinCrtMy, WinTypes, WinProcs, Strings, Type_def, CrtGr2D;


Const
      IncHor=  10;
      IncVer=   1;
      White =   0;
      Black =   1;

Type  Tpile = Array[1..10] of byte;
      TcadrePile = Array[1..10,1..4] of integer;

Var   Pile1,Pile2,Pile3: Tpile;
      CadrePile1,CadrePile2,CadrePile3: TcadrePile;
      Pointeur1,Pointeur2,Pointeur3: byte;
      NbDisk, NoDisk: byte;
      HPile,Hmax: byte;
      TrajHor, Notrans: integer;
      C1,C2,C3,C4: integer;
      BlackPen,WhitePen: HPen;

      Procedure Delay;
      var k,duree:longint;
          zz:double;
      Begin
        duree:=125000;
        for k:=1 to duree do  zz:=sin(zz)+cos(zz);
      End;

      Procedure DrawFrame(C:byte);
      Begin
        if C=White then
          SelectObject(CrtDc,WhitePen)
        else
          SelectObject(CrtDc,BlackPen);
        MoveXY(CrtDc,C1,C3);
        LineXY(CrtDc,C2,C3);
        LineXY(CrtDc,C2,C4);
        LineXY(CrtDc,C1,C4);
        LineXY(CrtDc,C1,C3)
      End;

      Procedure Transfer(N,Origin,Dest,Interm:integer);

        procedure MoveFrameUpwards;
        begin
          DrawFrame(White);
          SelectObject(CrtDc,BlackPen);
          MoveXY(CrtDc,C1,C3);
          LineXY(CrtDc,C2,C3);  {redraw lower side}
          repeat
            C3:=C3+IncVer;
            C4:=C4+IncVer;
            DrawFrame(Black);
            Delay;
            DrawFrame(White)
          until C4 >= Hmax
        end;

        procedure MoveFrameHoriz;
        var inc: integer;
        begin
          if Dest-Origin > 0 then
            inc:=IncHor
          else
            inc:=-IncHor;
          repeat
            DrawFrame(White);
            C1:=C1+inc;
            C2:=C2+inc;
            DrawFrame(Black);
            Delay
          until C1=TrajHor
        end;

        procedure MoveFrameDownwards;
        begin
          repeat
            DrawFrame(White);
            C3:=C3-IncVer;
            C4:=C4-IncVer;
            DrawFrame(Black);
            Delay
          until C4 <= HPile
        end;

        procedure MoveDisk(Origin,Dest:integer);
        var s1,s2,s3: String;
            ch:array[0..79] of char;
        begin
          Inc(Notrans);
          Str(Notrans:4,s1);
          Str(Origin:1,s2);
          Str(Dest:1,s3);
          StrPCopy(ch,s1+')  Move a disk from tower '+s2+' towards tower '+s3);
          TextOut(CrtDc,5,25,ch,strlen(ch));
          case Origin of
            1: begin
                 NoDisk:=Pile1[Pointeur1];
                 Dec(Pointeur1);
                 C1:=CadrePile1[NoDisk,1];
                 C2:=CadrePile1[NoDisk,2];
                 C3:=CadrePile1[NoDisk,3];
                 C4:=CadrePile1[NoDisk,4]
               end;
            2: begin
                 NoDisk:=Pile2[Pointeur2];
                 Dec(Pointeur2);
                 C1:=CadrePile2[NoDisk,1];
                 C2:=CadrePile2[NoDisk,2];
                 C3:=CadrePile2[NoDisk,3];
                 C4:=CadrePile2[NoDisk,4]
               end;
            3: begin
                 NoDisk:=Pile3[Pointeur3];
                 Dec(Pointeur3);
                 C1:=CadrePile3[NoDisk,1];
                 C2:=CadrePile3[NoDisk,2];
                 C3:=CadrePile3[NoDisk,3];
                 C4:=CadrePile3[NoDisk,4]
               end
          end;
          TrajHor:=(Dest-Origin)*30 + C1;
          MoveFrameUpwards;
          MoveFrameHoriz;

          case Dest of
            1: begin
                 Inc(Pointeur1);
                 HPile:=Pointeur1
               end;
            2: begin
                 Inc(Pointeur2);
                 HPile:=Pointeur2
               end;
            3: begin
                 Inc(Pointeur3);
                 HPile:=Pointeur3
               end
          end;

          MoveFrameDownwards;

          case Dest of
            1: begin
                 Pile1[Pointeur1]:=NoDisk;
                 CadrePile1[NoDisk,1]:=C1;
                 CadrePile1[NoDisk,2]:=C2;
                 CadrePile1[NoDisk,3]:=C3;
                 CadrePile1[NoDisk,4]:=C4
               end;
            2: begin
                 Pile2[Pointeur2]:=NoDisk;
                 CadrePile2[NoDisk,1]:=C1;
                 CadrePile2[NoDisk,2]:=C2;
                 CadrePile2[NoDisk,3]:=C3;
                 CadrePile2[NoDisk,4]:=C4
               end;
            3: begin
                 Pile3[Pointeur3]:=NoDisk;
                 CadrePile3[NoDisk,1]:=C1;
                 CadrePile3[NoDisk,2]:=C2;
                 CadrePile3[NoDisk,3]:=C3;
                 CadrePile3[NoDisk,4]:=C4
               end
          end
        end;   {MoveDisk}

      Begin    {Transfert}
        if N>0 then
        begin
          Transfer(N-1,Origin,Interm,Dest);
          MoveDisk(Origin,Dest);
          Transfer(N-1,Interm,Dest,Origin)
        end
      End;

      Procedure InputData;
      Begin
        ClrScr;
        writeln;
        write('  Number of disks: ');
        readln(NbDisk)
      End;

      Procedure Init;
      Var L:byte;
      Begin
        Hmax := NbDisk+2;
        Notrans:=0;
        {Init heap 1}
        CadrePile1[1,1]:=0;
        CadrePile1[1,2]:=20;
        CadrePile1[1,3]:=0;
        CadrePile1[1,4]:=1;
        For L:=2 to NbDisk do
        begin
          CadrePile1[L,1]:=CadrePile1[L-1,1]+1;
          CadrePile1[L,2]:=CadrePile1[L-1,2]-1;
          CadrePile1[L,3]:=CadrePile1[L-1,3]+1;
          CadrePile1[L,4]:=CadrePile1[L-1,4]+1
        end;
        {Init heap pointers}
        Pointeur1:=NbDisk;
        Pointeur2:=0;
        Pointeur3:=0;
        {Init Disk numbers}
        For L:=1 to NbDisk do
        begin
          Pile1[L]:=L;
          Pile2[L]:=0;
          Pile3[L]:=0
        end
      End;

      Procedure InitDrawing;
      Var i:byte;
      Begin
        BlackPen:=CreatePen(ps_Solid,1,RGB(0,0,0));
        WhitePen:=CreatePen(ps_Solid,1,RGB(255,255,255));
        Fenetre(0,80,0,20);
        Cloture(20,MaxX-40,80,MaxY);
        SelectObject(CrtDc,BlackPen);
        MoveXY(CrtDc,0,0);
        LineXY(CrtDc,80,0);
        For i:=1 to NbDisk do
        begin
          C1:=CadrePile1[I,1];
          C2:=CadrePile1[I,2];
          C3:=CadrePile1[I,3];
          C4:=CadrePile1[I,4];
          DrawFrame(Black)
        end
      End;


BEGIN
  WinCrtInit(' HANOI TOWERS');
    InputData;
    ClrScr;
    Init;
    InitDrawing;
    Transfer(NbDisk,1,3,2);
  SortieGraphique;
  DoneWinCrt
END.

{end of file hanoi.pas}