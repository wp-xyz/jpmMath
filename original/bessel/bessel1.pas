{*****************************************************************
*  Program to demonstrate the Bessel function asymptotic series  *
* -------------------------------------------------------------- *
*       Reference: BASIC Scientific Subroutines, Vol. II         *
*    By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [BIBLI 01]     *
*                                                                *
*                        Pascal Version by J.-P. Moreau, Paris.  *
*                                  (www.jpmoreau.fr)             *
* -------------------------------------------------------------- *
* SAMPLE RUN:                                                    *
*                                                                *
*  X            J0(X)            J1(X)        N            E     *
* -------------------------------------------------------------- *
*  1           0.733562         0.402234      1       0.1121521  *
*  2           0.221488         0.578634      1       0.0070095  *
*  3          -0.259956         0.340699      2       0.0007853  *
*  4          -0.396826        -0.065886      2       0.0001398  *
*  5          -0.177486        -0.327696      3       0.0000155  *
*  6           0.150635        -0.276760      3       0.0000036  *
*  7           0.300051        -0.004696      4       0.0000004  *
*  8           0.171638         0.234650      5       0.0000000  *
*  9          -0.090332         0.245324      5       0.0000000  *
* 10          -0.245930         0.043476      6       0.0000000  *
* 11          -0.171187        -0.176788      5       0.0000000  *
* 12           0.047689        -0.223450      5       0.0000000  *
* 13           0.206925        -0.070319      4       0.0000000  *
* 14           0.171072         0.133376      4       0.0000000  *
* 15          -0.014225         0.205105      4       0.0000000  *
* -------------------------------------------------------------- *
*                                                                *
*****************************************************************}
PROGRAM BESSEL1;
Uses WinCrtMy, Typedef, Utilit;

Var
        e, e3, J0, J1, x  : ar_reel;
        m1, m2, n1, n2 : ar_reel;
        i,n : integer;


{*********************************************************
* Bessel function asymptotic series subroutine. This     *
* program calculates the zeroth and first order Bessel   *
* functions using an asymptotic series expansion.        *
* The required input are X and a convergence factor e.   *
* Returned are the two Bessel functions, J0(X) and J1(X) *
* ------------------------------------------------------ *
* Reference; Algorithms for RPN calculators, by Ball,    *
* L.A. Wiley and sons.                                   *
*********************************************************}
PROCEDURE Calcul_Bessel;
Label 100, 200;
Var   a,a1,a2,b,c,e1,e2,e4,x1 : ar_reel;
      m : integer;
Begin
  {Calculate P and Q polynomials; P0(x):=m1, P1(x):=m2,
   Q0(x):=n1, Q1(x):=n2  }
  a := 1; a1 := 1; a2 := 1; b := 1; c := 1;
  e1 := 1e6;
  m := -1;
  x1 := 1.0 / (8.0 * x);
  x1 := x1 * x1;
  m1 := 1.0; m2 := 1.0;
  n1 := -1.0 / (8.0 * x);
  n2 := -3.0 * n1;
  n := 0;
100:  m := m + 2; a := a * m * m;
  m := m + 2; a := a * m * m;
  c := c * x1;
  a1 := a1 * a2; a2 := a2 + 1.0; a1 := a1 * a2;
  a2 := a2 + 1.0; e2 := a * c / a1;
  e4 := 1.0 + (m + 2) / m + (m + 2) * (m + 2) / (a2 * 8 * x) + (m + 2) * (m + 4) / (a2 * 8 * x);
  e4 := e4 * e2;
  {Test for divergence}
  IF ABS(e4) > e1 THEN GOTO 200;
  e1 := ABS(e2);
  m1 := m1 - e2;
  m2 := m2 + e2 * (m + 2) / m;
  n1 := n1 + e2 * (m + 2) * (m + 2) / (a2 * 8 * x);
  n2 := n2 - e2 * (m + 2) * (m + 4) / (a2 * 8 * x);
  n := n + 1;
  {Test for convergence criterion}
  IF e1 < e3 THEN GOTO 200;
  GOTO 100;
200: a := PI;
  e := e2;
  b := SQRT(2.0 / (a * x));
  J0 := b * (m1 * COS(x - a / 4) - n1 * SIN(x - a / 4.0));
  J1 := b * (m2 * COS(x - 3 * a / 4.0) - n2 * SIN(x - 3 * a / 4.0))
End;


{main program}
BEGIN
  CLRSCR;
  e3:=1e-10;
  Writeln;
  Writeln('   X            J0(X)            J1(X)        N            E      ');
  Writeln('------------------------------------------------------------------');
  FOR i := 1 TO 15 DO
  begin
    x := i;
    Calcul_Bessel;
    Write('  ',x:2:0,'        '); Aff_reel(J0); Write('       '); Aff_reel(J1);
    Write('       ',n); Write('       '); Aff_reel(e); Writeln
  end;
  Writeln('------------------------------------------------------------------');
  Readkey; DoneWinCrt
END.

{End of file bessel1.pas}




