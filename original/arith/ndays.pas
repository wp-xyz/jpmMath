{*******************************************
* Compute number of days between two dates *
* ---------------------------------------- *
* Ref.: "Mathématiques par l'informatique  *
*        individuelle (1) par H. Lehning   *
*        et D. Jakubowicz, Masson, Paris,  *
*        1982" [BIBLI 06].                 *
*                                          *
*     Pascal version By J-P Moreau, Paris. *
*              (www.jpmoreau.fr)           *
*******************************************}                     

{*****************************************************************************
 Note: The two dates must be between 01/01/1901
       and 12/31/2099 to avoid treatment of special
       "secular" years, such as 1800, 1900, 2100...
       (2000, a multiple of 400, is a normal year!)

 The first date is M, D, Y
 N, number of days between first date and fictitious
 starting date 01/00/1901, is given by formula:

       N = INT(365,25*(Y-1901)+A(M)+D)

 where A(M) is given by table below:

 M    1  2  3     4      5      6      7      8      9     10     11     12
 A(M) 0 31 59,25 90,25 120,25 151,25 181,25 212,25 243,25 273,25 304,25 334,25

       First N is saved in R.

 The second date is again M, D, Y
 N is again computed with same formula.

 Final result is N - R.
******************************************************************************}
PROGRAM NDAYS;
Uses WinCrt;

VAR
        A : Array[1..12] of REAL;
        M : WORD;                     {month number}
        D : WORD;                       {day number}
        Y : INTEGER;                   {year number}
        N,R : REAL;                 {number of days}

Procedure Calcul;
Begin
  Readln(M, D, Y);
  N := 365.25 * (Y - 1901) + A[M] + D;
  N := Round(N)
End;


BEGIN
  A[1] := 0; A[2] := 31; A[3] := 59.25; A[4] := 90.25; A[5] := 120.25; A[6] := 151.25;
  A[7] := 181.25; A[8] := 212.25; A[9] := 243.25; A[10] := 273.25; A[11] := 304.25;
  A[12] := 334.25;
  writeln;
  writeln(' Compute number of days between two dates');
  writeln('     (from 01/01/1901 to 12/31/2099)');
  writeln;
  write(' Enter first date (M D Y): ');
  Calcul;
  R := N;
  write(' Enter second date (M D Y): ');
  Calcul;
  R := N - R;
  writeln;
  writeln(' Number of days: ',R:6:0);
  writeln;
  ReadKey;
  DoneWinCrt
END.

{end of file ndays.pas}