{*********************************************************
'*                   SIMPLEX METHOD                      *
'*                   --------------                      *
'*                                                       *
'* LIST OF MAIN VARIABLES:                               *
'*                                                       *
'*  R:          MAXIMIZE = Y, MINIMIZE = N               *
'*  NV:         NUMBER OF VARIABLES OF ECONOMIC FUNCTION *
'*              (TO MAXIMIZE OR MINIMIZE).               *
'*  NC:         NUMBER OF CONSTRAINTS                    *
'*  TS:         SIMPLEX TABLE OF SIZE NC+1 x NV+1        *
'*  R1:         =1 TO MAXIMIZE, =-1 TO MINIMIZE          *
'*  R2:         AUXILIARY VARIABLE FOR INPUTS            *
'*  NOPTIMAL    BOOLEAN IF FALSE, CONTINUE ITERATIONS    *
'*  XMAX:       STORES GREATER COEFFICIENT OF ECONOMIC   *
'*              FUNCTION.                                *
'*  RAP         STORES SMALLEST RATIO > 0                *
'*  V:          AUXILIARY VARIABLE                       *
'*  P1,P2:      LINE, COLUMN INDEX OF PIVOT              *
'*  XERR:       BOOLEAN IF TRUE, NO SOLUTION.            *
'* ----------------------------------------------------- *
'* PROBLEM DESCRIPTION:                                  *
'* A builder of houses can make 3 kinds of them with     *
'* various profits: 15000$, 17000$ and 20000$.           *
'* Each kind must respect following conditions:          *
'* 1) for supply reasons, the number of houses of kind 2 *
'*    built each month must not exceed the number of     *
'*    kind 3 by more than two.                           *
'* 2) for staff reasons, the buider can make each month  *
'*    up to 5, 5 and 3, respectively of kinds 1, 2 and 3.*
'* 3) for organisation reasons, the builder can make up  *
'*    to 8 houses monthly of kinds 1,2 and 3, respecti-  *
'*    vely in the proportions of 3, 2 and 1.             *
'* The builder, having these data, wants to maximize his *
'* monthly profit by knowing the number oh houses of     *
'* each kind to build monthly.                           *
'* ----------------------------------------------------- *
'* SAMPLE RUN:                                           *
'* (Maximize 15 X1 + 17 X2 + 20 X3 with conditions:      *
'*                    X2 -   X3 <= 2                     *
'*           3 X1 + 3 X2 + 5 X3 <= 15                    *
'*           3 X1 + 2 X2 +   X3 <= 8     )               *
'*                                                       *
'* LINEAR PROGRAMMING                                    *
'*                                                       *
'* MAXIMIZE ? Y                                          *
'*                                                       *
'* NUMBER OF VARIABLES OF ECONOMIC FUNCTION ? 3          *
'*                                                       *
'* NUMBER OF CONSTRAINTS ? 3                             *
'*                                                       *
'* INPUT COEFFICIENTS OF ECONOMIC FUNCTION:              *
'*       #1 ? 15                                         *
'*       #2 ? 17                                         *
'*       #3 ? 20                                         *
'*       Right hand side ? 0                             *
'*                                                       *
'* CONSTRAINT #1:                                        *
'*       #1 ? 0                                          *
'*       #2 ? 1                                          *
'*       #3 ? -1                                         *
'*       Right hand side ? 2                             *
'*                                                       *
'* CONSTRAINT #2:                                        *
'*       #1 ? 3                                          *
'*       #2 ? 3                                          *
'*       #3 ? 5                                          *
'*       Right hand side ? 15                            *
'*                                                       *
'* CONSTRAINT #3:                                        *
'*       #1 ? 3                                          *
'*       #2 ? 2                                          *
'*       #3 ? 1                                          *
'*       Right hand side ? 8                             *
'*                                                       *
'* RESULTS:                                              *
'*       VARIABLE #1:   0.333333                         *
'*       VARIABLE #2:   3.000000                         *
'*       VARIABLE #3:   1.000000                         *
'*                                                       *
'*       ECONOMIC FUNCTION:   76.000000                  *
'*                                                       *
'* (Building monthly 1/3, 3 and 1 house(s) of kinds 1,2, *
'*  3, the builder can make a monthly profit of 76000$). *
'* ----------------------------------------------------- *
'* REFERENCE:                                            *
'* Modèles pratiques de décision Tome 2, By Jean-Pierre  *
'* Blanger, PSI Editions, France, 1982.                  *
'*                                                       *
'*            Pascal Release 1.0 By J-P Moreau, Paris.   *
'*                         (www.jpmoreau.fr)             *
'********************************************************}
PROGRAM SIMPLEX;
Uses WinCrt;

Const
      CMAX = 10;  {max. number of variables in economic function}
      VMAX = 10;  {max. number of constraints}

Var
      NC, NV, NOPTIMAL,P1,P2,XERR: Integer;
      TS: Array[0..CMAX,0..VMAX] of Double;


Procedure Data;
Var R1,R2: Double;
    R: Char;
    I,J: Integer;
Begin
  writeln;
  writeln(' LINEAR PROGRAMMING');
  writeln;
  write(' MAXIMIZE (Y/N) ? '); readln(R);
  writeln;
  write(' NUMBER OF VARIABLES OF ECONOMIC FUNCTION ? '); readln(NV);
  writeln;
  write(' NUMBER OF CONSTRAINTS ? '); readln(NC);
  writeln;
  IF Upcase(R) = 'Y' THEN
    R1 := 1
  ELSE
    R1 := -1;
  writeln(' INPUT COEFFICIENTS OF ECONOMIC FUNCTION:');
  FOR J := 1 TO NV DO
  begin
    write('       #',J,' ? '); readln(R2);
    TS[1, J + 1] := R2 * R1
  end;
  write('       Right hand side ? '); readln(R2);
  TS[1, 1] := R2 * R1;
  FOR I := 1 TO NC DO
  begin
    writeln;
    writeln(' CONSTRAINT #', I);
    FOR J := 1 TO NV DO
    begin
      write('       #',J,' ? '); readln(R2);
      TS[I + 1, J + 1] := -R2
    end;
    write('       Right hand side ? '); readln(TS[I + 1, 1])
  end;
  writeln;
  writeln(' RESULTS:');
  writeln;
  FOR J := 1 TO NV DO TS[0, J + 1] := J;
  FOR I := NV + 1 TO NV + NC DO TS[I - NV + 1, 0] := I
End;

  Procedure Pivot; Forward;
  Procedure Formula; Forward;
  Procedure Optimize; Forward;

Procedure SIMPLEX1;
Label 10;
Begin
10: PIVOT;
    FORMULA;
    OPTIMIZE;
    IF NOPTIMAL = 1 THEN GOTO 10
End;

Procedure PIVOT;
Label 100;
Var RAP,V,XMAX: Double;
    I,J: Integer;
Begin
  XMAX := 0.0;
  FOR J := 2 TO NV + 1 DO
  begin
    IF (TS[1, J] > 0) AND (TS[1, J] > XMAX) THEN
    begin
      XMAX := TS[1, J];
      P2 := J
    end
  end;
  RAP := 999999.0;
  FOR I := 2 TO NC + 1 DO
  begin
    IF TS[I, P2] >= 0 THEN GOTO 100;
    V := ABS(TS[I, 1] / TS[I, P2]);
    IF V < RAP THEN
    begin
      RAP := V;
      P1 := I
    end;
100: end;
  V := TS[0, P2]; TS[0, P2] := TS[P1, 0]; TS[P1, 0] := V
End;

Procedure FORMULA;
Label 60,70,100,110;
Var I,J: Integer;
Begin
     FOR I := 1 TO NC + 1 DO
     begin
       IF I = P1 THEN GOTO 70;
       FOR J := 1 TO NV + 1 DO
       begin
         IF J = P2 THEN GOTO 60;
         TS[I, J] := TS[I, J] - TS[P1, J] * TS[I, P2] / TS[P1, P2];
60:    end;
70:  end;
     TS[P1, P2] := 1.0 / TS[P1, P2];
     FOR J := 1 TO NV + 1 DO
     begin
       IF J = P2 THEN GOTO 100;
       TS[P1, J] := TS[P1, J] * ABS(TS[P1, P2]);
100: end;
     FOR I := 1 TO NC + 1 DO
     begin
       IF I = P1 THEN GOTO 110;
       TS[I, P2] := TS[I, P2] * TS[P1, P2];
110: end
End;

Procedure OPTIMIZE;
Label 10;
Var I,J: Integer;
Begin
  FOR I := 2 TO NC + 1 DO
    IF TS[I, 1] < 0 THEN XERR := 1;
  NOPTIMAL := 0;
  IF XERR = 1 THEN GOTO 10;
  FOR J := 2 TO NV + 1 DO
    IF TS[1, J] > 0 THEN NOPTIMAL := 1;
10: End;

Procedure RESULTS;
Label 30,70,100;
Var I,J: Integer;
Begin
    IF XERR = 0 THEN GOTO 30;
    writeln(' NO SOLUTION.'); GOTO 100;
30: FOR I := 1 TO NV DO
      FOR J := 2 TO NC + 1 DO
      begin
        IF TS[J, 0] <> I THEN GOTO 70;
        writeln('       VARIABLE #', I,': ', TS[J, 1]:10:6);
70:   end;
    writeln;
    writeln('       ECONOMIC FUNCTION: ', TS[1, 1]:10:6);
100:writeln; writeln
End;

{main program}
BEGIN
  ClrScr;
  Data;
  Simplex1;
  Results;
  ReadKey;
  DoneWinCrt
END.

{end of file simplex.pas}