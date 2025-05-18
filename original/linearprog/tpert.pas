{********************************************************
'*                 TIME  PERT  MODEL                    *
'*                                                      *
'* LIST OF MAIN VARIABLES:                              *
'*                                                      *
'* N        : NUMBER OF ACTIVITIES                      *
'* NA       : NUMBER OF PRECEDING ACTIVITIES            *
'* D(N)     : PROBABLE DURATION FOR EACH ACTIVITY       *
'* PRED(N)  : INDEX OF PRECEDING ACTIVITIES             *
'* ANTE(NA) : PRECEDING ACTIVITIES OF EACH ACTIVITY     *
'* MARQ(2,N): ACTIVITY LEVEL                            *
'* TOT(N)   : SOONEST DATE FOR EACH ACTIVITY            *
'* TARD(N)  : LATEST DATE FOR EACH ACTIVITY             *
'* DE       : ESTIMATED DURATION                        *
'* DOP      : OPTIMISTIC DURATION (DO IS NOT VALID)     *
'* DP       : PESSIMISTIC DURATION                      *
'* IANTE    : NUMBER OF PREVIOUS ACTIVITIES             *
'* K        : NUMBER OF ACTIVITY LEVEL                  *
'* R        : USED TO SEEK A LEVEL                      *
'* TOT      : USED TO SEEK SOONEST DATE                 *
'* TARD     : USED TO SEEK LATEST DATE                  *
'* FTA      : LATEST ENDING DATE                        *
'* FOT      : SOONEST ENDING DATE                       *
'* LO       : MAXIMUM WIDTH OG GANT'S DIAGRAM           *
'* MAXT     : PROBABLE MAXIMUM TIME                     *
'* X        : NORMALIZATION COEFFICIENT FOR GANT'S      *
'*            DIAGRAM                                   *
'* Y        : NUMBER OF UNITS FOR A SIGN * OR + IN      *
'*            GANT'S DIAGRAM                            *
'* ---------------------------------------------------- *
'* PROBLEM DESCRIPTION:                                 *
'* A Company, owning one mainframe computer, wants to   *
'* buy and install a second one. Before a false floor   *
'* must be laid to hide connection cables, and an air-  *
'* conditioning unit must be put in. The activities to  *
'* perform are given in the following table (durations  *
'* in hours):                                           *
'*   -------------------------------------------------  *
'*   Number  Task Description   Preceding  DE  DO  DP   *
'*                              Activities              *
'*   -------------------------------------------------  *
'*     1     Set the Dates of                           *
'*           All the Tasks          -       3   2   4   *
'*     2     Install Electric                           *
'*           Connections            1       7   4  10   *
'*     3     Install Air Condit.    2      16   8  24   *
'*     4     Disconnect Old                             *
'*           Computer               2       4   2   6   *
'*     5     Start New Air Cond.    3       2   1   3   *
'*     6     Lay False Floor       3,4     24  16  32   *
'*     7     Clean Up the Room      6       3   2   4   *
'*     8     Install New Computer  5,7     16   8  24   *
'*     9             ---            8       0   0   0   *
'*   -------------------------------------------------  *
'* SAMPLE RUN:                                          *
'*                  TIME PERT MODEL                     *
'*               ACTIVITY DESCRIPTION                   *
'*                                                      *
'* NUMBER OF ACTIVITIES ? 9                             *
'*                                                      *
'* NUMBER OF PRECEDING ACTIVITIES ? 11                  *
'*                                                      *
'* ACTIVITY #1:                                         *
'*   ESTIMATED DURATION   ? 3                           *
'*   OPTIMISTIC DURATION  ? 2                           *
'*   PESSIMISTIC DURATION ? 4                           *
'*   NUMBER OF PRECEDING ACTIVITIES ? 0                 *
'*                                                      *
'* ACTIVITY #2:                                         *
'*   ESTIMATED DURATION   ? 7                           *
'*   OPTIMISTIC DURATION  ? 4                           *
'*   PESSIMISTIC DURATION ? 10                          *
'*   NUMBER OF PRECEDING ACTIVITIES ? 1                 *
'*   PREC. ACTIVITY #1 ? 1                              *
'*                                                      *
'* ACTIVITY #3:                                         *
'*   ESTIMATED DURATION   ? 16                          *
'*   OPTIMISTIC DURATION  ? 8                           *
'*   PESSIMISTIC DURATION ? 24                          *
'*   NUMBER OF PRECEDING ACTIVITIES ? 1                 *
'*   PREC. ACTIVITY #1 ? 2                              *
'*                                                      *
'* ACTIVITY #4:                                         *
'*   ESTIMATED DURATION   ? 4                           *
'*   OPTIMISTIC DURATION  ? 2                           *
'*   PESSIMISTIC DURATION ? 6                           *
'*   NUMBER OF PRECEDING ACTIVITIES ? 1                 *
'*   PREC. ACTIVITY #1 ? 2                              *
'*                                                      *
'* ACTIVITY #5:                                         *
'*   ESTIMATED DURATION   ? 2                           *
'*   OPTIMISTIC DURATION  ? 1                           *
'*   PESSIMISTIC DURATION ? 3                           *
'*   NUMBER OF PRECEDING ACTIVITIES ? 1                 *
'*   PREC. ACTIVITY #1 ? 3                              *
'*                                                      *
'* ACTIVITY #6:                                         *
'*   ESTIMATED DURATION   ? 24                          *
'*   OPTIMISTIC DURATION  ? 16                          *
'*   PESSIMISTIC DURATION ? 32                          *
'*   NUMBER OF PRECEDING ACTIVITIES ? 2                 *
'*   PREC. ACTIVITY #1 ? 3                              *
'*   PREC. ACTIVITY #2 ? 4                              *
'*                                                      *
'* ACTIVITY #7:                                         *
'*   ESTIMATED DURATION   ? 3                           *
'*   OPTIMISTIC DURATION  ? 2                           *
'*   PESSIMISTIC DURATION ? 4                           *
'*   NUMBER OF PRECEDING ACTIVITIES ? 1                 *
'*   PREC. ACTIVITY #1 ? 6                              *
'*                                                      *
'* ACTIVITY #8:                                         *
'*   ESTIMATED DURATION   ? 16                          *
'*   OPTIMISTIC DURATION  ? 8                           *
'*   PESSIMISTIC DURATION ? 24                          *
'*   NUMBER OF PRECEDING ACTIVITIES ? 2                 *
'*   PREC. ACTIVITY #1 ? 5                              *
'*   PREC. ACTIVITY #2 ? 7                              *
'*                                                      *
'* ACTIVITY #9:                                         *
'*   ESTIMATED DURATION   ? 0                           *
'*   OPTIMISTIC DURATION  ? 0                           *
'*   PESSIMISTIC DURATION ? 0                           *
'*   NUMBER OF PRECEDING ACTIVITIES ? 1                 *
'*   PREC. ACTIVITY #1 ? 8                              *
'*                                                      *
'*                  TIME PERT MODEL                     *
'*                      RESULTS                         *
'*                                                      *
'*   CA  = CRITICAL ACTIVITY                            *
'*   SBD = SOONEST BEGIN DATE                           *
'*   LBD = LATEST BEGIN DATE                            *
'*   SED = SOONEST END DATE                             *
'*   LED = LATEST END DATE                              *
'*                                                      *
'*   ACTIVITY #1  PROBABLE DURATION: 3  * CA *          *
'*     SBD = 0   LBD = 0                                *
'*     SED = 3   LED = 3                                *
'*     TOTAL MARGIN = 0                                 *
'*                                                      *
'*   ACTIVITY #2  PROBABLE DURATION: 7  * CA *          *
'*     SBD = 3   LBD = 3                                *
'*     SED = 10  LED = 10                               *
'*     TOTAL MARGIN = 0                                 *
'*                                                      *
'*   ACTIVITY #3  PROBABLE DURATION: 16  * CA *         *
'*     SBD = 10  LBD = 10                               *
'*     SED = 26  LED = 26                               *
'*     TOTAL MARGIN = 0                                 *
'*                                                      *
'*   ANY KEY TO CONTINUE... ?                           *
'*                                                      *                                                  
'*   ACTIVITY #4  PROBABLE DURATION: 4                  *
'*     SBD = 10  LBD = 22                               *
'*     SED = 14  LED = 26                               *
'*     TOTAL MARGIN = 12                                *
'*                                                      *
'*   ACTIVITY #5  PROBABLE DURATION: 2                  *
'*     SBD = 26  LBD = 51                               *
'*     SED = 28  LED = 53                               *
'*     TOTAL MARGIN = 25                                *
'*                                                      *
'*   ACTIVITY #6  PROBABLE DURATION: 24  * CA *         *
'*     SBD = 26  LBD = 26                               *
'*     SED = 50  LED = 50                               *
'*     TOTAL MARGIN = 0                                 *
'*                                                      *
'*   ANY KEY TO CONTINUE... ?                           *
'*                                                      *
'*   ACTIVITY #7  PROBABLE DURATION: 3  * CA *          *
'*     SBD = 50  LBD = 50                               *
'*     SED = 53  LED = 53                               *
'*     TOTAL MARGIN = 0                                 *
'*                                                      *
'*   ACTIVITY #8  PROBABLE DURATION: 16  * CA *         *
'*     SBD = 53  LBD = 53                               *
'*     SED = 69  LED = 69                               *
'*     TOTAL MARGIN = 0                                 *
'*                                                      *
'*   ACTIVITY #9  PROBABLE DURATION: 0  * CA *          *
'*     SBD = 69  LBD = 69                               *
'*     SED = 69  LED = 69                               *
'*     TOTAL MARGIN = 0                                 *
'*                                                      *
'*   ANY KEY TO CONTINUE... ?                           *
'*                                                      *
'*   GANT'S DIAGRAM                                     *
'*                                                      *
'*   THE SIGNS * OR + ARE WORTH 3 TIME UNITS            *
'*                                                      *
'*   -----------------------------------------          *
'*    * 3                                               *
'*    *** 10                                            *
'*        ******** 26                                   *
'*        **++++++ 26                                   *
'*                 *++++++++++++ 53                     *
'*                 ************ 50                      *
'*                             * 53                     *
'*                              ******** 69             *
'*   -----------------------------------------          *
'*   * = PROBABLE DURATION    + = TOTAL MARGIN          *
'*                                                      *
'* ---------------------------------------------------- *
'* REFERENCE:                                           *
'* Modèles pratiques de décision Tome 2, By Jean-Pierre *
'* Blanger, PSI Editions, France, 1982.                 *
'*                                                      *
'*             Pascal Release 1.0 By J-P Moreau, Paris. *
*                        (www.jpmoreau.fr)              *
'********************************************************
'NOTE: PERT is an acronym for Program Evaluation & Review Technique.}

PROGRAM TPERT;
Uses WinCrt;

Const NMAX = 12;

Type TAB = Array[0..NMAX] of REAL;
     pITAB = ^ITAB;
     ITAB = Array[0..NMAX] of INTEGER;

Var  K, N, NA: Integer;
     D, TOT, TARD: TAB;
     ANTE, PRED: pITAB;
     MARQ: Array[1..2,0..NMAX] of Integer;
     XTOT,XTARD: REAL;

Procedure Data;  {DESCRIBE ACTIVITIES}
Label 20;
Var I,IANTE,J: Integer;
    DE,DOP,DP: REAL;
Begin
     writeln;
     writeln('                  TIME PERT MODEL');
     writeln('               ACTIVITY DESCRIPTION');
     writeln;
     write(' NUMBER OF ACTIVITIES ? '); readln(N);
     writeln;
     write(' NUMBER OF PRECEDING ACTIVITIES ? '); readln(NA);
     writeln;
     PRED^[0]:=0; ANTE^[0]:=0; ANTE^[1]:=0;
     FOR I := 1 TO N DO
     begin
       writeln(' ACTIVITY #', I,':');
       write('   ESTIMATED DURATION   ? '); readln(DE);
       write('   OPTIMISTIC DURATION  ? '); readln(DOP);
       write('   PESSIMISTIC DURATION ? '); readln(DP);
       D[I] := INT(10 * (DOP + 4 * DE + DP) / 6) / 10;
       write('   NUMBER OF PRECEDING ACTIVITIES ? '); readln(IANTE);
       IF IANTE = 0 THEN
       begin
         IANTE := 1; GOTO 20
       end;
       FOR J := 1 TO IANTE DO
       begin
         write('   PREC. ACTIVITY #', J,' ? ');
         readln(ANTE^[PRED^[I - 1] + J])
       end;
20:    writeln; PRED^[I] := PRED^[I - 1] + IANTE
     end;
     writeln;
     writeln('                  TIME PERT MODEL');
     writeln('                      RESULTS');
     writeln
End;

Procedure Level;   {Determine Level}
Label 30;
Var F,I,J,R: Integer; 
Begin
     MARQ[2, 0] := 1;
     FOR I := 1 TO N DO
     begin
       K := K + 1;
       FOR J := 1 TO N DO
       begin
         IF MARQ[2, J] <> 0 THEN GOTO 30;
         FOR R := PRED^[J - 1] + 1 TO PRED^[J] DO
           IF MARQ[2, ANTE^[R]] = 0 THEN F := 1;
         IF F = 1 THEN begin F := 0; GOTO 30 end;
         FOR R := PRED^[J - 1] + 1 TO PRED^[I] DO
           MARQ[1, J] := K;
30:    end;
       FOR J := 1 TO N DO
         IF MARQ[1, J] = K THEN MARQ[2, J] := K
     end;
     K := MARQ[2, 1];
     FOR I := 1 TO N DO
       IF MARQ[2, I] > K THEN K := MARQ[2, I]
End;

Procedure SDates;   {SOONEST DATES}
Label 10;
Var I,I1,J:Integer;
Begin
     FOR I := 1 TO K DO
       FOR J := 1 TO N DO
       begin
         IF MARQ[2, J] <> I THEN GOTO 10;
         FOR I1 := PRED^[J - 1] + 1 TO PRED^[J] DO
         begin
           XTOT := INT(10 * (TOT[ANTE^[I1]] + D[ANTE^[I1]])) / 10;
           IF TOT[J] = 0 THEN TOT[J] := XTOT;
           IF XTOT > TOT[J] THEN TOT[J] := XTOT
         end;
10:    end
End;

Procedure LDates;   {LATEST DATES}
Label 10,20;
Var I,I1,J: Integer;
Begin
     FOR I := 1 TO N DO
       IF MARQ[2, I] = K THEN
       begin
         TARD[I] := TOT[I]; I := N
       end;
     FOR I := K DOWNTO 1 DO
     begin
       FOR J := 1 TO N DO
       begin
         IF MARQ[2, J] <> I THEN GOTO 20;
         FOR I1 := PRED^[J - 1] + 1 TO PRED^[J] DO
         begin
           XTARD := INT(10 * (TARD[J] - D[ANTE^[I1]])) / 10;
           IF TARD[ANTE^[I1]] = 0 THEN
           begin
             TARD[ANTE^[I1]] := XTARD; GOTO 10;
           end;
           IF TARD[ANTE^[I1]] > XTARD THEN TARD[ANTE^[I1]] := XTARD;
10:      end;
20:    end
     end
End;

Procedure EDates;   {EDIT DATES}
Label 20,30;
Var I: Integer;
    FOT,FTA,XMT: REAL;
Begin
     writeln('   CA  := CRITICAL ACTIVITY');
     writeln('   SBD := SOONEST BEGIN DATE');
     writeln('   LBD := LATEST BEGIN DATE');
     writeln('   SED := SOONEST END DATE');
     writeln('   LED := LATEST END DATE');
     writeln;
     FOR I := 1 TO N DO
     begin
       write(' ACTIVITY #', I, '  PROBABLE DURATION: ', D[I]:4:0);
       FTA := TARD[I] + D[I]; FOT := TOT[I] + D[I];
       XMT := TARD[I] - TOT[I];
       IF XMT = 0 THEN writeln('  * CA *'); GOTO 20;
       writeln;
20:    writeln('  SBD = ', TOT[I]:4:0,'  LBD = ', TARD[I]:4:0);
       writeln('  SED = ', FOT:4:0,'  LED = ', FTA:4:0);
       writeln('  TOTAL MARGIN = ', XMT:4:0);
       writeln;
       IF ((I MOD 3) <> 0) OR (I = N) THEN GOTO 30;
       ReadKey;
       writeln;
30:  end;
     ReadKey;
     writeln
End;

Procedure Gant;   {GANT'S DIAGRAM}
Label 100,170,200;
Var XLO, XMAXT, Y: Integer;
    I, J, L: Integer;
    X: REAL;
Begin
     XLO := 35; XMAXT := 0;
     writeln(' GANT''S DIAGRAM');
     writeln;
     FOR I := 1 TO N DO
       IF TOT[I] + D[I] > XMAXT THEN XMAXT := Trunc(TOT[I] + D[I]);
     IF XMAXT < XLO THEN
     begin
       X := 1; GOTO 100
     end;
     X := Trunc(10*XLO / XMAXT)/10;
     Y := Trunc((XMAXT / XLO) / X);
100: writeln(' THE SIGNS * OR + ARE WORTH ', Y, ' TIME UNITS');
     FOR I := 1 TO Round(XLO + 8) DO write('-');
     writeln;
     FOR I := 1 TO N - 1 DO
     begin
       L := Trunc(TOT[I] * X);
       FOR J := 1 TO L DO  write(' ');
       L := Trunc(D[I] * X); IF L = 0 THEN GOTO 170;
       FOR J := 1 TO L DO  write('*');
170:   L := Trunc((TARD[I] - TOT[I]) * X);
       IF L = 0 THEN GOTO 200;
       FOR J := 1 TO L DO  write('+');
200:   writeln((TARD[I] + D[I]):4:0)
     end;
     FOR I := 1 TO Round(XLO + 8) DO write('-');
     writeln;
     writeln(' * = PROBABLE DURATION    + = TOTAL MARGIN');
End;


{main program}
BEGIN
  ClrScr;
  New(ANTE); New(PRED);
  Data;    {DESCRIBE ACTIVITIES}
  Level;   {DETERMINE LEVEL}
  SDates;  {SOONEST DATES}
  LDates;  {LATEST DATES}
  EDates;  {EDIT DATES}
  Dispose(ANTE); Dispose(PRED);
  Gant;    {GANT'S DIAGRAM}
  ReadKey;
  DoneWinCrt
END.

{end of file tpert.pas}