{*****************************************************
*   ROOTS OF A SECOND DEGREE EQUATION WITH COMPLEX   *
*   COEFFICIENTS  AZ2 + BZ + C = 0                   *
* -------------------------------------------------- *
* SAMPLE RUN:                                        *
* (Find complex roots of equation:                   *
* (-4+i)z2 + (2-3i)z +(5-i) = 0).                    *
*                                                    *
* SOLVING A COMPLEX 2ND DEGREE EQUATION              *
*     az2 + bz + c = 0                               *
*                                                    *
*  a real part      = -4                             *
*  a imaginary part = 1                              *
*  b real part      = 2                              *
*  b imaginary part = -3                             *
*  c real part      = 5                              *
*  c imaginary part = -1                             *
*                                                    *
*  Root 1 =   1.444644 -  0.352759 i                 *
*  Root 2 =  -0.797586 -  0.235476 i                 *
*                                                    *
* -------------------------------------------------- *
* Ref.: "Mathématiques en Turbo-Pascal By M. Ducamp  *
* and A. Reverchon (vol 2), Eyrolles, Paris, 1988"   *
* [BIBLI 03].                                        *
*****************************************************}
Program TEQUA2C;
Uses WinCrt;

Type
     COMPLEX = Record
       x,y: Double;
     End;

Var
     a,b,c,s1,s2: COMPLEX;

Procedure equa2c(a,b,c: COMPLEX; VAR s1, s2: COMPLEX);
Var u,v,w:COMPLEX; r: Double;
Begin
  u.x := b.x * b.x - b.y * b.y - 4 * a.x * c.x + 4 * a.y * c.y;
  u.y := 2 * b.x * b.y - 4 * a.x * c.y - 4 * a.y * c.x;
  r := sqrt(u.x*u.x + u.y*u.y);
  v.x := sqrt((r+u.x)/2);
  v.y := sqrt((r-u.x)/2);
  if u.y<0 then v.y:=-v.y;
  w.x := (-b.x - v.x)/2;
  w.y := (-b.y - v.y)/2;
  u.x := (-b.x + v.x)/2;
  u.y := (-b.y + v.y)/2;
  r := a.x * a.x + a.y * a.y;
  s1.x := (a.x * w.x + a.y * w.y)/r;
  s1.y := (a.x * w.y - a.y * w.x)/r;
  s2.x := (a.x * u.x + a.y * u.y)/r;
  s2.y := (a.x * u.y - a.y * u.x)/r;
End;

BEGIN
  writeln;
  writeln(' SOLVING A COMPLEX 2ND DEGREE EQUATION');
  writeln('      az2 + bz + c := 0');
  writeln;
  write('    a      real part = '); readln(a.x);
  write('    a imaginary part = '); readln(a.y);
  write('    b      real part = '); readln(b.x);
  write('    b imaginary part = '); readln(b.y);
  write('    c      real part = '); readln(c.x);
  write('    c imaginary part = '); readln(c.y);

  equa2c(a,b,c,s1,s2);

  writeln;
  write(' Root 1 = ', s1.x:10:6);
  if s1.y > 0 then writeln(' + ', abs(s1.y):10:6,' i')
  else writeln(' - ', abs(s1.y):10:6,' i');
  write(' Root 2 = ', s2.x:10:6);
  if s2.y > 0 then writeln(' + ', abs(s2.y):10:6,' i')
  else writeln(' - ', abs(s2.y):10:6,' i');
  writeln;
  ReadKey;
  DoneWinCrt
END.

{end of file tequa2.pas}