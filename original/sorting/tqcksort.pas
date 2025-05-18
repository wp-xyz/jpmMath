{**************************************************************
*         Sorting an array with the Quicksort method          *
* ----------------------------------------------------------- *
* REFERENCE:                                                  *
*                                                             *
*   "NUMERICAL RECIPES by W.H. Press, B.P. Flannery,          *
*    S.A. Teukolsky and W.T. Vetterling, Cambridge            *
*    University Press, 1986".                                 *
*                                                             *
*                              Pascal version by J-P Moreau   *
* ----------------------------------------------------------- *
* SAMPLE RUN:                                                 *
*                                                             *
* Table to be sorted:                                         *
*                                                             *
* 102.2 837.1 270.3 484.0 164.8 236.6 805.3 142.6 988.2  23.2 *
* 352.8   7.4 745.3 691.9 907.8 674.6 918.7 854.5 894.6 845.9 *
* 207.4 116.6 992.9 291.4 974.7 494.5 115.0 262.6 831.1 554.4 *
* 135.3 937.0 383.4 567.8 640.7 894.0 147.0 754.9 266.1 673.7 *
* 579.2 271.0 345.8 927.9 229.9  31.0 663.1 295.8 823.6  30.0 *
*  82.1 987.1   1.7 661.4 580.2 415.0 553.6 577.1 593.7 334.0 *
* 853.8 183.1 255.2 793.3 692.3 137.7 820.9 871.2 224.9 814.9 *
* 335.0 509.0  53.5 141.2 725.4 462.4 805.6 652.7 641.2 609.3 *
*                                                             *
* Sorted table (Quicksort method):                            *
*                                                             *
*   1.7   7.4  23.2  30.0  31.0  53.5  82.1 102.2 115.0 116.6 *
* 135.3 137.7 141.2 142.6 147.0 164.8 183.1 207.4 224.9 229.9 *
* 236.6 255.2 262.6 266.1 270.3 271.0 291.4 295.8 334.0 335.0 *
* 345.8 352.8 383.4 415.0 462.4 484.0 494.5 509.0 553.6 554.4 *
* 567.8 577.1 579.2 580.2 593.7 609.3 640.7 641.2 652.7 661.4 *
* 663.1 673.7 674.6 691.9 692.3 725.4 745.3 754.9 793.3 805.3 *
* 805.6 814.9 820.9 823.6 831.1 837.1 845.9 853.8 854.5 871.2 *
* 894.0 894.6 907.8 918.7 927.9 937.0 974.7 987.1 988.2 992.9 *
*                                                             *
**************************************************************}
PROGRAM QUICKSORT;
Uses WinCrt;

CONST   SIZE = 100;                   {maximum size of table}

TYPE
        Table = Array[1..SIZE] of real;

VAR
        A         : Table;               {Table to be sorted}
        MAX_VALUE : real;            {Maximum value of table}
        i,N       : integer;

{**********************************************************
* Sorts an array ARR of length N into ascending numerical *
* order using the Quicksort algorithm. N is input, ARR is *
* replaced on output by its sorted rearrangement.         *
**********************************************************}
Procedure QCKSRT (N:Integer; Var ARR:Table);

Const M = 7; NSTACK = 50;
{ Here M is the size of subarrays sorted by straight insertion,
  NSTACK is the required auxuliary storage. }
Label 10, 12, 20, 21, 22, 30, 40;
Var ISTACK:Array[1..NSTACK] of integer;
    I,IQ,IR,J,JSTACK,L: Integer;
    A:real;
Begin
    JSTACK := 0;
    L := 1; IR := N;
10: IF (IR - L < M) THEN
    begin
      FOR J := L + 1 TO IR DO
      begin
        A := ARR[J];
        FOR i := J - 1 DOWNTO 1 DO
        begin
          IF ARR[i] <= A THEN GOTO 12;
          ARR[i + 1] := ARR[i]
        end;
        i := 0;
12:     ARR[i + 1] := A
      end;
      IF (JSTACK = 0) THEN GOTO 40;
      IR := ISTACK[JSTACK];
      L := ISTACK[JSTACK - 1];
      JSTACK := JSTACK - 2
    end
    ELSE
    begin
      i := L; J := IR;
    { Generate a random integer IQ between L and IR inclusive }
      IQ := Round(L + (IR - L) * Random);
      A := ARR[IQ];
      ARR[IQ] := ARR[L];
20:   IF (J > 0) THEN
21:     IF A < ARR[J] THEN
        begin
          J := J - 1;
          GOTO 21
        end;
      IF J <= i THEN
      begin
        ARR[i] := A;
        GOTO 30
      end;
      ARR[i] := ARR[J];
      i := i + 1;
22:   IF i <= N THEN
        IF A > ARR[i] THEN
        begin
          i := i + 1;
          GOTO 22
        end;
      IF J <= i THEN
      begin
        ARR[J] := A;
        i := J;
        GOTO 30
      end;
      ARR[J] := ARR[i];
      J := J - 1;
      GOTO 20;
30:   JSTACK := JSTACK + 2;
      IF JSTACK > NSTACK THEN
      begin
        Writeln(' NSTACK too small!');
        goto 40
      end;
      IF (IR - i >= i - L) THEN
      begin
        ISTACK[JSTACK] := IR;
        ISTACK[JSTACK - 1] := i + 1;
        IR := i - 1
      end
      ELSE
      begin
        ISTACK[JSTACK] := i - 1;
        ISTACK[JSTACK - 1] := L;
        L := i + 1
      end
    end;
    GOTO 10;
40: END;


{write table of size N to standard output}
Procedure TWRIT(N:INTEGER;ARR:Table);
Var i:integer;
begin
  writeln;
  for i:=1 to N do
  begin
    write(ARR[i]:6:0);
    if (i MOD 10)=0 then writeln
  end
end;


{main program}
BEGIN
  {initialize random number generator}
  Randomize;

  N:=80;    {initialize size of table}
  MAX_VALUE := 1000.0;

  {generate random table of numbers (between 0 and MAX_VALUE) }
  for i:=1 to N do A[i]:=MAX_VALUE*Random;

  writeln;
  writeln(' Table to be sorted:');
  TWRIT(N,A);
 
  {call sorting subroutine}
  QCKSRT (N, A);

  writeln;
  writeln(' Sorted table (Quicksort method):');
  TWRIT(N,A);

  ReadKey; DoneWinCrt

END.

{end of file tqcksrt.pas}
