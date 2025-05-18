{*******************************************************
*             Equations of degree 2 to 4               *
* ---------------------------------------------------- *
* This program calculates the real or complex roots of *
* algebraic equations of degree 2, 3 and 4.            *
* ---------------------------------------------------- *
* Reference: Mathématiques et statitiques - Programmes *
*            en BASIC. Editions du P.S.I., 1981        *
*            [BIBLI 13].                               *
*                                                      *
*                 Pascal Version By J-P Moreau, Paris  *
*                          (www.jpmoreau.fr)           *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
*                                                      *
* ROOTS OF EQUATIONS OF DEGREE 2, 3 OU 4               *
* Example: X^4+4X^2-3X+3 = 0                           *
*                                                      *
* INSTRUCTIONS:                                        *
*                                                      *
* 1. Input degree of equation                          *
* 2. Input coefficients beginning with highest exponent*
* 3. Solutions are given under the form:               *
*     (real part) + I (imaginary part)                 *
*                                                      *
* Degree of equation (2 to 4): 4                       *
*                                                      *
* Coefficient of X^4 = 1                               *
* Coefficient of X^3 = 0                               *
* Coefficient of X^2 = 4                               *
* Coefficient of X^1 = -3                              *
* Coefficient of X^0 = 3                               *
*                                                      *
* Root 1 = 0.449702657 + I 0.731070322                 *
* Root 2 = 0.449702657 - I 0.731070322                 *
* Root 1 = -0.449702657 + I 1.967231848                *
* Root 2 = -0.449702657 - I 1.967231848                *
*                                                      *
*******************************************************}
Program Root4;
Uses WinCrt;

Type    TAB = Array[0..4] of DOUBLE;
        MAT = Array[0..4,0..4] of DOUBLE;

Var
        a : MAT;
        r : TAB;
        i,n : INTEGER;
        aa,b,c,d,ii,im,k,l,m,q,rr,s,sw : DOUBLE;


Function f(x:double): double;
begin
  f := x * x * x + a[3, 2] * x * x + a[3, 1] * x + a[3, 0]
end;


{******************************************
* This Subroutine calculates the roots of *
*    X^2 + A(2,1)*X + A(2,0) = 0          *
*    W = Determinant of equation          *
******************************************}
Procedure Root_2;
Var q1, q2, w : double;
begin
  w := a[2, 1] * a[2, 1] - 4 * a[2, 0];
  q1 := -a[2, 1] / 2.0;
  if w < 0 then
  begin
    q2 := SQRT(-w) / 2.0;
    im := q2;
    r[1] := q1; r[2] := q1
  end
  else if w=0 then
  begin
    r[1] := q1; r[2] := q1; im := 0.0
  end
  else if w > 0 then
  begin
    q2 := SQRT(w) / 2.0; im := 0.0;
    r[1] := q1 + q2; r[2] := q1 - q2
  end
end;


{*******************************************
* This subroutine calculates the roots of  *
* X^3 + A(3,2)*X^2 + A(3,1)*X + A(3,0) = 0 *
*******************************************}
Procedure Root_3;
Label 100,200,500;
Var am, er, te, tt, xa, xb, xc, y1, y2 : double;
begin
  {one root equals zero}
  IF a[3, 0] = 0 THEN
  BEGIN
    xc := 0; GOTO 500
  END;
  {looking for maximum coefficient in absolute magnitude}
  am := ABS(a[3, 0]);
  FOR i := 1 TO 3 DO
  begin
    tt := ABS(a[3, i]);
    IF am < tt THEN am := tt
  end;
  {Define interval where a real root exists
   according to sign of A(3,0)  }
  IF a[3, 0] > 0 THEN
  begin
    aa := -am - 1;
    b := 0;
    GOTO 100
  end;
  aa := 0; b := am + 1;

100:
  {Searching for xc = real root in interval (a,b)
   (by Bisection method)  }
  
  {Define segment (xa,xb) including the root}
  xa := aa; xb := b; y1 := f(aa); y2 := f(b);

  {Desired précision}
  er := 0.000001;
  IF y1 = 0 THEN
  begin
    xc := aa; GOTO 500
  end;
  IF y2 = 0 THEN
  begin
    xc := b; GOTO 500
  end;
  {xc : middle of segment}
  200: xc := (xa + xb) / 2.0; te := f(xc);
  IF te = 0 THEN GOTO 500;
  IF (xb - xa) < er THEN GOTO 500;
  IF f(xa) * te > 0 THEN
  begin
    xa := xc; GOTO 200
  end;
  xb := xc; GOTO 200;
  {r[3] is a real root}
  500: r[3] := xc;
  IF sw = -1 THEN exit;
  {Calculates the roots of remaining quadratic equation
   Define the equation coefficients }
  a[2, 1] := a[3, 2] + xc;
  a[2, 0] := a[3, 1] + a[3, 2] * xc + xc * xc;
  Root_2
end;

{Root search main subroutine}
Procedure Root_4;
LABEL 100,200,300;
begin
  {Normalization of coefficients }
  FOR i := 0 TO n - 1 DO
  begin
    a[n,i] := a[n,i] / a[n,n]
  end;
  {Branching according to degree n}
  if n=2 then
  begin
    Root_2;
    exit
  end
  else if n=3 then
  begin
    Root_3;
    exit
  end
  else if n=4 then
  begin
    aa := a[4, 3];
    b := a[4, 2];
    c := a[4, 1];
    d := a[4, 0];
    q := b - (3.0 * aa * aa / 8.0);
    rr := c - (aa * b / 2) + (aa * aa * aa / 8.0);
    s := d - (aa * c / 4) + (aa * aa * b / 16.0) - (3 * aa * aa * aa * aa / 256.0);
    {Défine coefficients of cubic equation}
    a[3, 2] := q / 2.0;
    a[3, 1] := (q * q - 4 * s) / 16.0;
    a[3, 0] := -(rr * rr / 64.0);
    {Calculate a real root of this equation}
    IF ((rr <> 0) OR (a[3, 1] >= 0)) THEN GOTO 100;
    {Particular case when this equation is of 2nd order}
    a[2, 1] := a[3, 2]; a[2, 0] := a[3, 1];
    Root_2;
    {One takes the positive root}
    r[3] := r[1];
    GOTO 200;
    {Calling Root_3 with sw=-1 to calculate one root only}
    100: sw := -1;
    Root_3;
    {real root of above equation}
    200: k := SQRT(r[3]);
    {Calculate L and M if k=0}
    IF k = 0 THEN
    begin
      rr := SQRT(q * q - 4 * s); GOTO 300
    end;
    q := q + (4 * r[3]);
    rr := rr / (2 * k);
    300: l := (q - rr) / 2.0; m := (q + rr) / 2.0;
    {Solving two equations of degree 2}
    a[2, 1] := 2 * k; a[2, 0] := l;
    {1st equation}
    Root_2;
    {Transferring solutions in r[3], r[4], ii}
    r[3] := r[1] - (a[4, 3] / 4); r[4] := r[2] - (a[4, 3] / 4.0); ii := im;
    a[2, 1] := -2 * k; a[2, 0] := m;
    {2nd equation}
    Root_2;
    r[2] := r[2] - (a[4, 3] / 4.0); r[1] := r[1] - (a[4, 3] / 4.0)
  end
end;


{main program}
BEGIN
  Writeln( ' ROOTS OF EQUATIONS OF DEGREE 2, 3 OR 4');
  Writeln( ' EXAMPLE : X^4+4X^2-3X+3 = 0');
  Writeln;
  Writeln( ' INSTRUCTIONS:');
  Writeln;
  Writeln( ' 1. Input degree of equation');
  Writeln( ' 2. Input coefficients beginning by highest exponent');
  Writeln( ' 3. Solutions are given under the form:');
  Writeln( '    (real part) + I (imaginary part)');
  Writeln;

  Write(' DEGREE OF EQUATION [2 to 4]: '); Readln(n);
  IF n < 2 THEN n := 2;
  IF n > 4 THEN n := 4;
  Writeln;
  FOR i := n DOWNTO 0 DO
  begin
    Write( ' COEFFICIENT OF X^',i,' = '); Readln(a[n,i])
  end;
  Writeln;
  {calling root search main subroutine}
  Root_4;
  {writing results}
  IF n=2 THEN
  BEGIN
    {cas n=2}
    IF im = 0 THEN
    begin
      Writeln( ' Root 1 = ',r[1]:12:9);
      Writeln( ' Root 2 = ',r[2]:12:9)
    end
    else
    begin
      Writeln( ' Root 1 = ',r[1]:12:9,'  +I  ',im:12:9);
      Writeln( ' Root 2 = ',r[2]:12:9,'  -I  ',im:12:9)
    end;
    Readkey; DoneWinCrt
  END
  ELSE IF n=3 THEN
  BEGIN
    {cas n:=3}
    Writeln( ' Root 1 = ',r[3]:12:9);
    IF im = 0 THEN
    begin
      Writeln( ' Root 2 = ',r[1]:12:9);
      Writeln( ' Root 3 = ',r[2]:12:9)
    end
    else
    begin
      Writeln( ' Root 2 = ',r[1]:12:9,'  +I  ',im:12:9);
      Writeln( ' Root 3 = ',r[2]:12:9,'  -I  ',im:12:9)
    end;
    Readkey; DoneWinCrt
  END
  ELSE IF n=4 THEN
  BEGIN
    {cas n:=4}
    IF im=0 THEN
    begin
      Writeln( ' Root 1 = ',r[1]:12:9);
      Writeln( ' Root 2 = ',r[2]:12:9)
    end
    else
    begin
      Writeln( ' Root 1 = ',r[1]:12:9,'  +I  ',im:12:9);
      Writeln( ' Root 2 = ',r[2]:12:9,'  -I  ',im:12:9)
    end;
    IF ii=0 THEN
    begin
      Writeln( ' Root 3 = ',r[3]:12:9);
      Writeln( ' Root 4 = ',r[4]:12:9)
    end
    else
    begin
      Writeln( ' Root 3 = ',r[3]:12:9,'  +I  ',ii:12:9);
      Writeln( ' Root 4 = ',r[4]:12:9,'  -I  ',ii:12:9)
    end;
    Readkey; DoneWinCrt
  END
END.

{End of file root4.pas}