{*********************************************************
*                  DANTZIG'S  MODEL                      *
*                                                        *
* LIST OF MAIN VARIABLES:                                *
*                                                        *
*   NR       : TOTAL NUMBER OF MARKS                     *
*   V(NR,NR) : PATH VALUES                               *
*   XL(NR)   : USED FOR WEIGHTING OF EACH PATH           *
*   IE(NR*2) : SET OF ADOPTED MARKS                      *
*   V        : PATH VALUE                                *
*   IR       : NUMBER OF MARKS COMING FROM A SAME MARK   *
*   RA       : INDEX OF AN ARRIVAL MARK                  *
*   R        : = "YES" ON OPTIMAL ROUTE                  *
*   XM1,XM2  : USED TO SEEK MINIMUM, MAXIMUM PATH        *
*   IASSO    : ARRIVAL MARK ASSOCIATED WITH SET IE       *
* ------------------------------------------------------ *
* PROBLEM DESCRIPTION:                                   *
* A mountaineer wants to climb up a stiff face. He has   *
* preliminarily estimated that he can follow different   *
* paths from 0 (bottom) to 8 (summit) with different     *
* lengths:                                               *
*              Path 1:  0 - 1 - 4 - 7 - 8                *
*              Path 2:  0 - 3 - 5 - 8                    *
*              Path 3:  0 - 2 - 6 - 5 - 8                *
*              Path 3a: 0 - 2 - 6 - 7 - 8                *
*              Path 4:  0 - 1 - 2 - 4 - 7 - 8            *
*              Path 4a: 0 - 1 - 2 - 4 - 6 - 7 - 8        *
* The mountaineer wants to know the maximum length rope  *
* to take,  knowing the following distance matrix        *
* (meters):                                              *
*     START/ARR/DIST  START/ARR/DIST  START/ARR/DIST     *
*       0    1   5      0    2   11     0    3   10      *
*       1    2   4      1    4   6                       *
*       2    4   5      2    6   9                       *
*       3    5   2      4    7   5                       *
*       5    6   4      5    8   15                      *                                                
*       6    7   6      7    8   8                       *
* ------------------------------------------------------ *
* SAMPLE RUN:                                            *
*                                                        *
* DO YOU WANT THE MINIMUM PATH (Y/N) ? N                 *
*                                                        *
* NUMBER OF MARKS ? 9                                    *
*                                                        *
* NUMBER OF ARCS FROM MARK #0 ? 3                        *
* ARRIVAL MARK, VALUE ? 1 5                              *
* ARRIVAL MARK, VALUE ? 2 11                             *
* ARRIVAL MARK, VALUE ? 3 10                             *
*                                                        *
* NUMBER OF ARCS FROM MARK #1 ? 2                        *
* ARRIVAL MARK, VALUE ? 2 4                              *
* ARRIVAL MARK, VALUE ? 4 6                              *
*                                                        *
* NUMBER OF ARCS FROM MARK #2 ? 2                        *
* ARRIVAL MARK, VALUE ? 4 5                              *
* ARRIVAL MARK, VALUE ? 6 9                              *
*                                                        *
* NUMBER OF ARCS FROM MARK #3 ? 1                        *
* ARRIVAL MARK, VALUE ? 5 2                              *
*                                                        *
* NUMBER OF ARCS FROM MARK #4 ? 1                        *
* ARRIVAL MARK, VALUE ? 7 5                              *
*                                                        *
* NUMBER OF ARCS FROM MARK #5 ? 2                        *
* ARRIVAL MARK, VALUE ? 6 4                              *
* ARRIVAL MARK, VALUE ? 8 15                             *
*                                                        *
* NUMBER OF ARCS FROM MARK #6 ? 1                        *
* ARRIVAL MARK, VALUE ? 7 6                              *
*                                                        *
* NUMBER OF ARCS FROM MARK #7 ? 1                        *
* ARRIVAL MARK, VALUE ? 8 8                              *
*                                                        *
* RESULTS:                                               *
*                                                        *
* THE MAXIMUM VALUE PATH IS:                             *
*                                                        *
*     <- 0 - 2 - 6 - 7 - 8 ->                            *
*                                                        *
* PATH VALUE =  34.0000                                  *
*                                                        *
* ------------------------------------------------------ *
* REFERENCE:                                             *
*  Modèles pratiques de décision Tome 2, By Jean-Pierre  *
*  Blanger, PSI Editions, France, 1982.                  *
*                                                        *
*               Pascal Release 1.0 By J-P Moreau, Paris. *
*                           (www.jpmoreau.fr)            *
*********************************************************}
PROGRAM DANTZIG;

Uses WinCrt;

Const NMAX = 25;

Type
     pTAB = ^TAB;
     TAB = Array[0..NMAX,0..NMAX] of REAL;
     pVEC = ^VEC;
     VEC = Array[0..NMAX] of REAL;
     pIVEC = ^IVEC;
     IVEC = Array[0..2*NMAX] of INTEGER;

Var
     R:  Char;
     K,NR: Integer;
     V:  pTAB;
     XL: pVEC;
     IE: pIVEC;
     XM1,XM2,XV: REAL;

Procedure Data;  {INPUT DATA}
Var I,IR,J,RA: Integer;
    XV:REAL;
Begin
  For i:=0 to NMAX do
    For j:=0 to NMAX do
      V^[I,j]:=0;
  Writeln;
  Write(' DO YOU WANT THE MINIMUM PATH (Y/N) ? '); Readln(R);
  Writeln;
  Write(' NUMBER OF MARKS ? '); Readln(NR);
  Writeln;
  FOR I := 0 TO NR - 2 DO
  begin
    Write(' NUMBER OF ARCS FROM MARK #', I,' ? ');
    Readln(IR);
    FOR J := 1 TO IR DO
    begin
      Write(' ARRIVAL MARK, VALUE ? '); Readln(RA, XV);
      V^[I, RA] := XV
    end;
    Writeln
  end;
  Writeln;
  Writeln(' RESULTS');
  Writeln;
  IF UpCase(R) = 'N' THEN
    Writeln(' THE MAXIMUM VALUE PATH IS:')
  ELSE
    Writeln(' THE MINIMUM VALUE PATH IS:');
  Writeln
End;

Procedure DANTZIG1;
Label 30, 50, 70;
Var F,IASSO,IO,ID,J: Integer;
Begin
     XL^[0] := 0;
     K := 1; IE^[K] := 0;
30:  IF UpCase(R) = 'N' THEN
     begin
       XM1 := -999999; XM2 := -999999
     end
     ELSE
     begin
       XM1 := 999999; XM2 := 999999
     end;
     FOR IO := 1 TO K DO
     begin
       FOR ID := 1 TO NR DO
       begin
         IF UpCase(R) = 'N' THEN
         begin
           IF (V^[IE^[IO], ID] < XM1) OR (V^[IE^[IO], ID] = 0) THEN GOTO 70
         end
         ELSE
         begin
           IF (V^[IE^[IO], ID] > XM1) OR (V^[IE^[IO], ID] = 0) THEN GOTO 70
         end;
         FOR J := 1 TO K DO
         begin
           IF IE^[J] = ID THEN
           begin
             F := 1;
             goto 50
           end
         end;
50:      IF F = 1 THEN begin F := 0; GOTO 70 end;
         XM1 := V^[IE^[IO], ID]; IASSO := ID;
70:    end;
       XV := XL^[IE^[IO]] + XM1;
       IF UpCase(R) = 'N' THEN
         IF (XV > XM2) THEN begin XM2 := XV; IE^[K + 1] := IASSO end
       ELSE
         IF (XV < XM2) THEN begin XM2 := XV; IE^[K + 1] := IASSO end
     end;
     XL^[IE^[K + 1]] := XM2;
     K := K + 1; IF IE^[K] <> NR - 1 THEN GOTO 30
End;

Procedure Results;  {EDIT OPTIMAL PATH}
Label 30, 50;
Var IO,ID,K: Integer;
    L:REAL;
Begin
     IO := IE^[1]; ID := IE^[2]; K := 1; L:=0;
     Write(' <- 0 - ');
30:  IF V^[IO, ID] <> 0 THEN
     begin
       L:=L+V^[IO,ID];
       Write(ID,' - '); IO := ID; K := K + 1; ID := IE^[K];
       GOTO 50
     end;
     K := K + 1; ID := IE^[K];
50:  IF ID <> NR - 1 THEN GOTO 30;
     L:=L+V^[IO,ID];
     Writeln(ID,' ->');
     Writeln;
     writeln(' Path Value = ', L:8:4);
End;


{main Program}
BEGIN
  ClrScr;
  New(V); New(XL); New(IE);
  Data;      {INPUT DATA}
  Dantzig1;  {DANTZIG MODEL}
  Results;   {EDIT OPTIMAL PATH};
  ReadKey;
  Dispose(V); Dispose(XL); Dispose(IE);
  DoneWinCrt
END.

{end of file dantzig.pas}