{*****************************************************
*        Computing the statistical functions         *
*             for one or two variables               *
*                                                    *
* -------------------------------------------------- *
* REFERENCE:  "Mathematiques et statistiques by H.   *
*              Haut, PSI Editions, France, 1981"     *
*              [BIBLI 13].                           *
*                                                    *
*              Pascal version by J-P Moreau, Paris   *
*                       (www.jpmoreau.fr)            *
* -------------------------------------------------- *
* SAMPLE RUN:                                        *
*                                                    *
* TUTORIAL                                           *
*                                                    *
* 1. Define type of calculus:"                       *
*    1: Statistical functions for a set X(i)         *
*    2: Statistical functions for a set X(i), Y(i)   *
*                                                    *
* 2. Input number n of data                          *
*                                                    *
* 3. Input successively the n values X(i) [and Y(i)] *
*                                                    *
* Type of calculus (1 or 2): 2                       *
*                                                    *
* Number of data: 5                                  *
*                                                    *
*   1  1 12                                          *
*   2  2 9                                           *
*   3  3 7                                           *
*   4  4 15                                          *
*   5  5 6                                           *
*                                                    *
*                                                    *
* Mean of X(i)....................: 3.00000000       *
*                                                    *
* (n-1) standard deviation of X(i): 1.58113883       *
*   (n) standard deviation of X(i): 1.41421356       *
*                                                    *
* (n-1) standard dev. of X mean...: 0.70710678       *
*   (n) standard dev. of X mean...: 0.63245553       *
*                                                    *
* Mean of Y(i)....................: 9.80000000       *
*                                                    *
* (n-1) standard deviation of Y(i): 3.70135111       *
*   (n) standard deviation of Y(i): 3.31058908       *
*                                                    *
* (n-1) standard dev. of Y mean...: 1.65529454       *
*   (n) standard dev. of Y mean...: 1.48054045       *
*                                                    *
*                                                    *
* (n-1) covariance of X,Y.........: -1.50000000      *
*   (n) covariance of X,Y.........: -1.20000000      *
*                                                    *
* Correlation coefficient.........: -0.25630730      *
*                                                    *
******************************************************
 Description: This program computes the usual statistical functions
              for a set of n variables X(i) or for a set of n pairs
              X(i), Y(i).   }
PROGRAM FStat;
Uses WinCrt;

Label   fin;

Const   SIZE = 100;

Var
        X, Y : Array[1..SIZE] of DOUBLE;
        a3,a4,a5,a6,a7 : DOUBLE;
        v1,v2,v3,v4,v5,v6,v7,v8: DOUBLE;
        i, n, nt : INTEGER;


{****************************************************
*       Subroutine of statistical functions         *
*              for 1 or 2 variables                 *
* ------------------------------------------------- *
* INPUTS:                                           *
*         nt: type of calculus =1 for one variable  *
*                              =2 for two variables *
*          n: number of data X(i) [and Y(i) ]       *
* OUTPUTS:                                          *
*         v1: mean of X(i)                          *
*         v2: mean of Y(i)                          *
*         v3: (n-1) standard deviation of v1        *
*         v4: (n-1) standard deviation of v2        *
*         v5:   (n) standard deviation of v1        *
*         v6:   (n) standard deviation of v2        *
*         v7: (n-1) covariance of X,Y               *
*         v8: correlation coefficient               *
*         a3: (n-1) standard deviation of X         *
*         a4: (n-1) standard deviation of Y         *
*         a5:   (n) standard deviation of X         *
*         a6:   (n) standard deviation of Y         *
*         a7:   (n) covariance of X,Y               *
****************************************************}
Procedure Stat_functions;
Label 100;
Var i,n1:INTEGER;
    xnr:DOUBLE;
Begin
  v1:=0.0; v2:=0.0; v3:=0.0; v4:=0.0; v5:=0.0;
  n1:=n-1; xnr:=SQRT(n);
  {choose type of calculus}
  if nt=2 then goto 100;
  {case of one set X(i)}
  for i:=1 to n do
  begin
    v1:=v1+X[i];
    v3:=v3+SQR(X[i])
  end;
  v6:=v3-v1*v1/n;
  v1:=v1/n;
  a3:=SQRT(v6/n1); a5:=SQRT(v6/n);
  v5:=a5/xnr; v3:=a3/xnr;
  exit;
100: {case of a set X(i), Y(i) }
  for i:=1 to n do
  begin
    v1:=v1+X[i];
    v2:=v2+Y[i];
    v3:=v3+SQR(X[i]);
    v4:=v4+SQR(Y[i]);
    v5:=v5+X[i]*Y[i]
  end;
  v6:=v3-v1*v1/n;
  v7:=v4-v2*v2/n;
  v8:=v5-v1*v2/n;
  v1:=v1/n; v2:=v2/n;
  a3:=SQRT(v6/n1); a4:=SQRT(v7/n1);
  a5:=SQRT(v6/n); a6:=SQRT(v7/n);
  v7:=v8/n1; a7:=v8/n;
  v3:=a3/xnr; v4:=a4/xnr;
  v5:=a5/xnr; v6:=a6/xnr;

  if (a3=0) or (a4=0) then exit; {correlation coefficient not defined}
  v8:=v8/(n*a5*a6)
End;


{main program}
BEGIN

  clrscr;
  writeln;
  writeln(' TUTORIAL');
  writeln;
  writeln(' 1. Define type of calculus:');
  writeln('    1: Statistical functions for a set X(i)');
  writeln('    2: Statistical functions for a set X(i), Y(i)');
  writeln;
  writeln(' 2. Input number n of data');
  writeln;
  writeln(' 3. Input successively the n values X(i) [and Y(i)]');
  writeln;
  write(' Type of calculus (1 or 2): '); read(nt);
  writeln;
  write(' Number of data: '); read(n);
  writeln;
  {read n data}
  for i:=1 to n do
  begin
    write('  ',i,'  ');
    if nt=1 then readln(X[i]);
    if nt=2 then readln(X[i], Y[i])
  end;

  {call stat subroutine}
  Stat_functions;

  {writeln; results}
  writeln;
  writeln(' Mean of X(i)....................: ', v1:12 :8);
  writeln;
  writeln(' (n-1) standard deviation of X(i): ', a3:12:8);
  writeln('   (n) standard deviation of X(i): ', a5:12:8);
  writeln;
  writeln(' (n-1) standard dev. of X mean...: ', v3:12:8);
  writeln('   (n) standard dev. of X mean...: ', v5:12:8);
  writeln;
  if nt=1 then goto fin;
  writeln(' Mean of Y(i)....................: ', v2:12:8);
  writeln;
  writeln(' (n-1) standard deviation of Y(i): ', a4:12:8);
  writeln('   (n) standard deviation of Y(i): ', a6:12:8);
  writeln;
  writeln(' (n-1) standard dev. of Y mean...: ', v4:12:8);
  writeln('   (n) standard dev. of Y mean...: ', v6:12:8);
  writeln;
  writeln;
  writeln(' (n-1) covariance of X,Y.........: ', v7:12:8);
  writeln('   (n) covariance of X,Y.........: ', a7:12:8);
  writeln;
  writeln(' Correlation coefficient.........: ', v8:12:8);
  writeln;
fin: Readkey; DoneWinCrt

END.

{end of file fstat.pas}