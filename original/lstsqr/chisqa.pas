{********************************************************
*      Program to demonstrate Chi-square Statistic      *
* ----------------------------------------------------- *
*    Reference: BASIC Scientific Subroutines, Vol. II   *
*    By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                       *
*                Pascal version by J-P Moreau, Paris    *
*                         (www.jpmoreau.fr)             *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
*                                                       *
*  CHI-SQUARE CUMULATIVE DISTRIBUTION APPROXIMATION     *
*                                                       *
*   P(X)        X2       REAL X2                        *
*  -----------------------------                        *
*  0.05     124.5347     124.3                          *
*  0.10     118.6326     118.5                          *
*  0.15     114.5999                                    *
*  0.20     111.4342                                    *
*  0.25     108.7913     109.1                          *
*  0.30     106.5055                                    *
*  0.35     104.4821                                    *
*  0.40     102.6610                                    *
*  0.45     101.0018                                    *
*  0.50      99.1944      99.3                          *
*  0.55      97.6863                                    *
*  0.60      96.0812                                    *
*  0.65      94.3595                                    *
*  0.70      92.4934                                    *
*  0.75      90.4428      90.1                          *
*  0.80      88.1446                                    *
*  0.85      85.4893                                    *
*  0.90      82.2522      82.4                          *
*  0.95      77.7867      77.9                          *
*  1.00       0.0000                                    *
*  -----------------------------                        *
*                                                       *
********************************************************}
PROGRAM CHISQA;
Uses WinCrt;

Var  x,y : DOUBLE;
     i,m : integer;
     X1  : array[1..20] of string[5];



{***************************************************
* Chi-square cumulative distribution approximation *
* Good for m > 100.                                *
* Reference: Statistics Manual, Crow, Maxfield and *
* Davis (Dover, 1960).                             *
* The input value is y, the probability, the output* 
* value is the corresponding Chi-square statistic. *
***************************************************}
PROCEDURE Chi_Square;
Label 100, 200, 210;
Var z : DOUBLE;
Begin

  x := y;
  {Guard against zero discontinuity}
  IF x <= 0 THEN x := EXP(-100);
  IF x > 0.5 THEN GOTO 100;
  x := -Ln(x);
  {Regressed table correction}
  z := -0.803 + 1.312 * x - 0.2118 * x * x + 0.016 * x * x * x;
  GOTO 200;
100: x := 1.0 - x;
  {Guard against zero discontinuity}
  IF x <= 0 THEN GOTO 210;
  x := -Ln(x);
  z := 0.803 - 1.312 * x + 0.2118 * x * x - 0.016 * x * x * x ;
200: x := 2.0 / (9.0 * m);
  x := 1.0 - x + z * SQRT(x);
210: x := m * x * x * x

End;


{main program}
BEGIN

  x1[1]  := '124.3'; x1[2]  := '118.5';
  x1[5]  := '109.1'; x1[10] := ' 99.3';
  x1[15] := ' 90.1'; x1[18] := ' 82.4'; x1[19] := ' 77.9';
  Writeln(' CHI-SQUARE CUMULATIVE DISTRIBUTION APPROXIMATION');
  Writeln('   P(X)         X²      REAL X² ');
  Writeln('   ---------------------------- ');
  m := 100; i := 0;
  y := 0.05;
  While y < 1.01 do
  begin
    i := i + 1;
    Chi_Square;
    Writeln('  ',y:5:2,'     ',x:8:4,'    ',x1[i]:5);
    y := y + 0.05;
  end;
  Write('   ---------------------------- ');
  Readln; DoneWinCrt

END.

{End of file chisqa.pas}
