{***************************************************
*      Program to demonstrate least squares        *
*         polynomial fitting subroutine            *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal version by J-P Moreau, Paris.   *
*                    (www.jpmoreau.fr)             *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* LEAST SQUARES POLYNOMIAL FITTING                 *
*                                                  *
* What is the order of the fit      : 4            *
* What is the error reduction factor: 0            *
* How many data pooints are there   : 11           *
*                                                  *
* Input the data points as prompted:               *
*                                                  *
*  1   X, Y = 1,1                                  *
*  2   X, Y = 2,2                                  *
*  3   X, Y = 3,3                                  *
*  4   X, Y = 4,4                                  *
*  5   X, Y = 5,5                                  *
*  6   X, Y = 6,6                                  *
*  7   X, Y = 7,7                                  *
*  8   X, Y = 8,8                                  *
*  9   X, Y = 9,9                                  *
*  10   X, Y = 0,0                                 *
*  11   X, Y = 1,1                                 *
*                                                  *
* Coefficients are:                                *
*  0    -0.000000                                  *
*  1    1                                          *
*  2    -0.000000                                  *
*  3    0                                          *
*  4    -0.000000                                  *
*                                                  *
* Standard deviation =  0.00000000                 *
*                                                  *
***************************************************}
PROGRAM LSQRPOLY;
Uses WinCrt;

TYPE    TAB = Array[0..25] of DOUBLE;

VAR
        i,l,m,n : integer;
        dd,e1,vv : DOUBLE;

        x,y,v,a,b,c,d,c2,e,f : TAB;



{****************************************************************
*        LEAST SQUARES POLYNOMIAL FITTING PROCEDURE             *
* ------------------------------------------------------------- *
* This program least squares fits a polynomial to input data.   *
* Forsythe orthogonal polynomials are used in the fitting.      *
* The number of data points is n.                               *
* The data is input to the subroutine in x[i], y[i] pairs.      *
* The coefficients are returned in c[i],                        *
* the smoothed data is returned in v[i],                        *
* the order of the fit is specified by m.                       *
* The standard deviation of the fit is returned in d.           *
* There are two options available by use of the parameter e:    *
*  1. if e = 0, the fit is to order m,                          *
*  2. if e > 0, the order of fit increases towards m, but will  *
*     stop if the relative standard deviation does not decrease *
*     by more than e between successive fits.                   *
* The order of the fit then obtained is l.                      *
****************************************************************}
PROCEDURE LS_POLY;
Label 10,15,20,30,50,fin;
Var i,l2,n1 : integer;
    a1,a2,b1,b2,c1,d1,f1,f2,v1,v2,w : DOUBLE;
Begin
  n1 := m + 1;
  v1 := 1e7;
  {Initialize the arrays}
  FOR i := 1 TO n1 DO
  begin
    a[i] := 0; b[i] := 0; f[i] := 0
  end;
  FOR i := 1 TO n DO
  begin
    v[i] := 0; d[i] := 0
  end;
  d1 := SQRT(n); w := d1;
  FOR i := 1 TO n DO
  begin
    e[i] := 1 / w
  end;
  f1 := d1; a1 := 0;
  FOR i := 1 TO n DO
  begin
    a1 := a1 + x[i] * e[i] * e[i]
  end;
  c1 := 0;
  FOR i := 1 TO n DO
  begin
    c1 := c1 + y[i] * e[i]
  end;
  b[1] := 1 / f1; f[1] := b[1] * c1;
  FOR i := 1 TO n DO
  begin
    v[i] := v[i] + e[i] * c1
  end;
  m := 1;
10: {Save latest results}
  FOR i := 1 TO l DO c2[i] := c[i];
  l2 := l; v2 := v1; f2 := f1; a2 := a1; f1 := 0;
  FOR i := 1 TO n DO
  begin
    b1 := e[i];
    e[i] := (x[i] - a2) * e[i] - f2 * d[i];
    d[i] := b1;
    f1 := f1 + e[i] * e[i]
  end;
  f1 := SQRT(f1);
  FOR i := 1 TO n DO e[i] := e[i] / f1;
  a1 := 0;
  FOR i := 1 TO n DO a1 := a1 + x[i] * e[i] * e[i];
  c1 := 0;
  FOR i := 1 TO n DO c1 := c1 + e[i] * y[i]; 
  m := m + 1; i := 0;
15: l := m - i; b2 := b[l]; d1 := 0;
  IF l > 1 THEN d1 := b[l - 1];
  d1 := d1 - a2 * b[l] - f2 * a[l];
  b[l] := d1 / f1; a[l] := b2; i := i + 1;
  IF i <> m THEN GOTO 15;
  FOR i := 1 TO n DO v[i] := v[i] + e[i] * c1; 
  FOR i := 1 TO n1 DO
  begin
    f[i] := f[i] + b[i] * c1;
    c[i] := f[i]
  end;
  vv := 0;
  FOR i := 1 TO n DO vv := vv + (v[i] - y[i]) * (v[i] - y[i]); 
  {Note the division is by the number of degrees of freedom}
  vv := SQRT(vv / (n - l - 1)); l := m;
  IF e1 = 0 THEN GOTO 20;
  {Test for minimal improvement}
  IF ABS(v1 - vv) / vv < e1 THEN GOTO 50;
  {If error is larger, quit}
  IF e1 * vv > e1 * v1 THEN GOTO 50;
  v1 := vv;
20: IF m = n1 THEN GOTO 30;
  GOTO 10;
30: {Shift the c[i] down, so c(0) is the constant term}
  FOR i := 1 TO l DO c[i - 1] := c[i];
  c[l] := 0;
  {l is the order of the polynomial fitted}
  l := l - 1; dd := vv;
  goto fin;
50: {Aborted sequence, recover last values}
  l := l2; vv := v2;
  FOR i := 1 TO l DO c[i] := c2[i];
  GOTO 30;
fin: End;


BEGIN
  Clrscr;
  Writeln(' LEAST SQUARES POLYNOMIAL FITTING');
  Writeln;
  write(' What is the order of the fit      : '); Read(m);
  write(' What is the error reduction factor: '); Read(e1);
  write(' How many data pooints are there   : '); Read(n);
  Writeln;
  Writeln(' Input the data points as prompted:');
  Writeln;
  FOR i := 1 TO n DO
  begin
    Write(' ',i,'  X  Y = '); Read(x[i], y[i])
  end;
  Writeln;

  LS_POLY;

  Writeln(' Coefficients are:');
  FOR i := 0 TO l DO Writeln('  ',i,'  ',c[i]:9:6);
  Writeln;
  Writeln(' Standard deviation: ',dd:9:8);
  writeln;
  ReadKey; DoneWinCrt
END.

{End lsqrpoly.pas}




