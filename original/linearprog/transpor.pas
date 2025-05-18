{************************************************************
*                 THE TRANSPORT  MODEL                      *
*                 --------------------                      *
*                                                           *
* LIST OF MAIN VARIABLES:                                   *
*                                                           *
* ISOU          : NUMBER OF SOURCES                         *
* IDES          : NUMBER OF DESTINATIONS                    *
* XM(ISOU,IDES) : RESOURCES DISTRIBUTION                    *
* C(ISOU,IDES)  : UNITARY TRANSPORT COSTS                   *
* R(ISOU,IDES)  : PATH SEARCH                               *
* P(4,ISOU+IDES): TABLE OF FOUND PATHS                      *
* DAO(ISOU)     : AVAILABLE SOURCE QUANTITIES               *
* RAD(IDES)     : REQUIRED DESTINATION QUANTITIES           *
* IX, IY        : AUXILIARY VARIABLES AT TRANSPORT PATH     *
* CC            : PATH COST                                 *
* XINFCC        : MINIMUM OF PATH COSTS                     *
* NR            : NUMBER OF KNOWN MARKS                     *
* IF1           : FLAG FOR MARKS STILL TO ELIMINATE         *
* ID            : DISPLACEMENT LENGTH                       *
* TT            : TRANSPORT TOTAL QUANTITY                  *
* LT            : TOTAL PATH LENGTH                         *
* CT            : TOTAL TRANSPORT COST                      *
* IOPTIMAL      : FLAG=1, IF OPTIMALITY IS FOUND            *
* QT            : TRANSPORT QUANTUTY FOR ONE PATH           *
* --------------------------------------------------------- *
* PROBLEM DESCRIPTION                                       *
* A papermaker owns three factories respectively producing  *
* monthly 5000, 10000 and 6000 tons. He has to deliver the  *
* following quantities respectively to four customers in    *
* different cities: 2000 T (Paris), 11000 T (Lyon), 4000 T  *
* (Marseille) and 4000 T (Grenoble) at the lowest total     *
* transportation cost. The transport costs grid is the      *
* following:                                                *
*               PARIS  LYON  MARSEILLE  GRENOBLE            *
*   FACTORY 1    200    700     800     very high           *
*   FACTORY 2    400    400     500       500               *
*   FACTORY 3    400    500     600       600               *
* --------------------------------------------------------- *
* SAMPLE RUN:                                               *
*                                                           *
* TRANSPORT MODEL                                           *
*                                                           *
* NUMBER OF SOURCES ? 3                                     *
*                                                           *
* NUMBER OF DESTINATIONS ? 4                                *
*                                                           *
* INPUT THE AVAILABLE SOURCE QUANTITIES:                    *
*  SOURCE #1 ? 5000                                         *
*  SOURCE #2 ? 10000                                        *
*  SOURCE #3 ? 6000                                         *
*                                                           *
* INPUT THE REQUIRED DESTINATION QUANTITIES:                *
*  DESTINATION #1 ? 2000                                    *
*  DESTINATION #2 ? 11000                                   *
*  DESTINATION #3 ? 4000                                    *
*  DESTINATION #4 ? 4000                                    *
*                                                           *
* INPUT TRANSPORT COSTS MATRIX:                             *
*  FROM SOURCE #1 TO DESTINATION #1 ? 200                   *
*  FROM SOURCE #1 TO DESTINATION #2 ? 700                   *
*  FROM SOURCE #1 TO DESTINATION #3 ? 800                   *
*  FROM SOURCE #1 TO DESTINATION #4 ? 9999                  *
*  FROM SOURCE #2 TO DESTINATION #1 ? 400                   *
*  FROM SOURCE #2 TO DESTINATION #2 ? 400                   *
*  FROM SOURCE #2 TO DESTINATION #3 ? 500                   *
*  FROM SOURCE #2 TO DESTINATION #4 ? 500                   *
*  FROM SOURCE #3 TO DESTINATION #1 ? 400                   *
*  FROM SOURCE #3 TO DESTINATION #2 ? 500                   *
*  FROM SOURCE #3 TO DESTINATION #3 ? 600                   *
*  FROM SOURCE #3 TO DESTINATION #4 ? 600                   *
*                                                           *
* TRANSPORT MODEL                                           *
*                                                           *
* TRANSPORTS:                                               *
*    FROM SOURCE #1 TO DESTINATION #1:  2000.00             *
*    FROM SOURCE #1 TO DESTINATION #2:  3000.00             *
*    FROM SOURCE #2 TO DESTINATION #2:  8000.00             *
*    FROM SOURCE #2 TO DESTINATION #3:  2000.00             *
*    FROM SOURCE #3 TO DESTINATION #3:  2000.00             *
*    FROM SOURCE #3 TO DESTINATION #4:  4000.00             *
*                                                           *
* TOTAL TRANSPORT COST: 10300000.0                          *
*                                                           *
* --------------------------------------------------------- *
* REFERENCE:                                                *
* Modèles pratiques de décision Tome 2, By Jean-Pierre      *
* Blanger, PSI Editions, France, 1982.                      *
*                                                           *
*                 Pascal Release 1.0 By J-P Moreau, Paris.  *
*                          (www.jpmoreau.fr)                *
************************************************************}
PROGRAM TRANSPORT;

Uses Wincrt;

Const
     NMAX=10;

Type
     TAB  = Array[1..NMAX,1..NMAX] of REAL;
     TAB0 = Array[0..NMAX,0..NMAX] of REAL;
     VEC  = Array[1..NMAX] of REAL;

Var
     C, XM: TAB;
     R: TAB0;
     DAO, RAD: VEC;

     P: Array[1..4,1..25] of INTEGER;

     ID,IDES,IF1,IOPTIMAL,ISOU,IX,IY,LT,NR: Integer;
     CC,CT,QT,TT,XINFCC: REAL;


Procedure Init;
Var I,J: Integer;
Begin
  writeln;
  writeln(' TRANSPORT MODEL');
  writeln;
  write(' NUMBER OF SOURCES ? '); readln(ISOU);
  writeln;
  write(' NUMBER OF DESTINATIONS ? '); readln(IDES);
  writeln;
  writeln(' INPUT THE AVAILABLE SOURCE QUANTITIES:');
  FOR I := 1 TO ISOU DO
  begin
    write('  SOURCE #',I,' ? '); readln(DAO[I])
  end;
  writeln;
  writeln(' INPUT THE REQUIRED DESTINATION QUANTITIES:');
  FOR I := 1 TO IDES DO
  begin
    write('  DESTINATION #',I,' ? '); readln(RAD[I])
  end;
  writeln;
  writeln(' INPUT TRANSPORT COSTS MATRIX:');
  FOR I:=1 TO ISOU DO
    FOR J:=1 TO IDES DO
    begin
      write('  FROM SOURCE #',I,' TO DESTINATION #',J,' ? ');
      readln(C[I,J])
    end;
  writeln;
  writeln(' TRANSPORT MODEL')
End;

  Procedure Corner; Forward;
  Procedure Optimal; Forward;
  Procedure TCost; Forward;
  Procedure Seek_Path(I,J:Integer); Forward;
  Procedure Increase(I,J:Integer); Forward;


Procedure Main;  {STEPPING STONE ALGORITHM}
Begin
  Corner;
  Optimal;
  Tcost
End;

Procedure Corner;  {N-W CORNER}
Label 20,50,60;
Var I,J: Integer;
Begin
     I:=1; J:=1;
20:  IF DAO[I]<=RAD[J] THEN GOTO 50;
     XM[I,J]:=XM[I,J]+RAD[J];DAO[I]:=DAO[I]-RAD[J];
     RAD[J]:=0; Inc(J);
     GOTO 60;  {ELSE}
50:  XM[I,J]:=XM[I,J]+DAO[I];RAD[J]:=RAD[J]-DAO[I];
     DAO[I]:=0; Inc(I);
60:  IF (I<=ISOU) AND (J<=IDES) THEN GOTO 20
End;

Procedure Optimal;
Label 10,70,140,150;
Var I,J: Integer;
Begin
10:  XINFCC:=0.0;
     FOR I:=1 TO ISOU DO
       FOR J:=1 TO IDES DO
       begin
         IF XM[I,J]<>0.0 THEN GOTO 70;
         Seek_Path(I,J);
         Increase(I,J);
70:    end;
     IF XINFCC>=0.0 THEN
     begin
       IOPTIMAL:=1;
       GOTO 150;
     end;
     FOR I:=1 TO LT DO
     begin
       IX:=P[3,I]; IY:=P[4,I];
       IF (I MOD 2)=0 THEN
       begin
         XM[IX,IY]:=XM[IX,IY]-TT;
         GOTO 140
       end;
       XM[IX,IY]:=XM[IX,IY]+TT;
140: end;
150: IF IOPTIMAL=0 THEN GOTO 10
End;

Procedure Seek_Path(I,J:Integer);
Label 70,160,260;
Var I1,I2: Integer;
Begin
     FOR I1:=1 TO ISOU DO
       FOR I2:=1 TO IDES DO
         R[I1,I2]:=XM[I1,I2];
     FOR I1:=1 TO ISOU DO R[I1,0]:=0.0;
     FOR I2:=1 TO IDES DO R[0,I2]:=0.0;
     R[I,J]:=1.0;
70:  FOR I2:=1 TO IDES DO
     begin
       IF R[0,I2]=1 THEN GOTO 160;
       NR:=0;
       FOR I1:=1 TO ISOU DO
         IF R[I1,I2]<>0.0 THEN Inc(NR);
       IF NR<>1 THEN GOTO 160;
       FOR I1:=1 TO ISOU DO R[I1,I2]:=0.0;
       R[0,I2]:=1.0; IF1:=1;
160: end;
     FOR I1:=1 TO ISOU DO
     begin
       IF R[I1,0]=1 THEN GOTO 260;
       NR:=0;
       FOR I2:=1 TO IDES DO
         IF R[I1,I2]<>0.0 THEN Inc(NR);
       IF NR<>1 THEN GOTO 260;
       FOR I2:=1 TO IDES DO R[I1,I2]:=0.0;
       R[I1,0]:=1.0; IF1:=1;
260: end;
     IF IF1=1 THEN
     begin
       IF1:=0; GOTO 70
     end
End;

Procedure Increase(I,J:Integer);
Label 20,70,130,170,180,230;
Var I1,I2: Integer;
Begin
     P[1,1]:=I; P[2,1]:=J; IX:=I; IY:=J; ID:=1; CC:=0.0; QT:=999999.0;
20:  Inc(ID); IF1:=0;
     FOR I1:=1 TO ISOU DO
     begin
       IF (R[I1,IY]=0) OR (I1=IX) THEN GOTO 70;
       P[1,ID]:=I1; P[2,ID]:=IY; IX:=I1; CC:=CC-C[IX,IY];
       IF1:=1; I1:=ISOU;
       IF (XM[IX,IY] < QT) AND (XM[IX,IY] > 0.0) THEN QT:=XM[IX,IY];
70:  end;
     IF IF1=0 THEN GOTO 170;
     Inc(ID); IF1:=0;
     FOR I2:=1 TO IDES DO
     begin
       IF (R[IX,I2]=0) OR (I2=IY) THEN GOTO 130;
       P[1,ID]:=IX; P[2,ID]:=I2; IY:=I2; CC:=CC+C[IX,IY];
       IF1:=1; I2:=IDES;
130: end;
     IF IF1=0 THEN GOTO 170;
     IF (IX<>I) OR (IY<>J) THEN GOTO 20;
     GOTO 180;
170: writeln(' DEGENERATE SOLUTION !');
     Halt(0);
180: IF (CC>0) OR (CC>XINFCC) THEN GOTO 230;
     TT:=QT; XINFCC:=CC; ID:=ID-1; LT:=ID;
     FOR I1:=1 TO ID DO
     begin
       P[3,I1]:=P[1,I1]; P[4,I1]:=P[2,I1]
     end;
230:End;

Procedure TCost;
Label 10;
Var I,J:Integer;
Begin
     CT:=0.0;
     writeln;
     writeln(' TRANSPORTS:');
     FOR I:=1 TO ISOU DO
       FOR J:=1 TO IDES DO
       begin
         CT:=CT+XM[I,J]*C[I,J];
         IF XM[I,J]=0 THEN GOTO 10;
         writeln('    FROM SOURCE #',I,' TO DESTINATION #',J,': ', XM[I,J]:8:2);
10:    end;
     writeln;
     writeln(' TOTAL TRANSPORT COST: ', CT:10:1);
     writeln
End;

{main program}
BEGIN
  ClrScr;
  Init;
  Main;
  ReadKey;
  DoneWinCrt
END.

{end of file transpor.pas}