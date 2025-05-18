{**************************************************************************
*               Test program for Bauhuber's method                        *
* ----------------------------------------------------------------------- *
*  This program uses Bauhuber's Method to find all real or complex        *
*  roots of a polynomial of degree n:                                     *
*                                               n-1           n           *
*      P(x) = a[0] + a[1] * x + ... + a[n-1] * x    + a[n] * x ,          *
*                                                                         *
*  where a[i], i=0, ..., n, can be complex.                               *
* ----------------------------------------------------------------------- *
* SAMPLE RUN:                                                             *
* (Find all roots of polynomial of degree 4:                              *
*  P(x) = 0.000089248 - 0.04368 X + 2.9948 X2 - 6.798 X3 + X4)            *
*                                                                         *
*  ----------------------------------------------------------             *
*  Complex and Real Roots of a Polynomial (Bauhuber's method)             *
*  ----------------------------------------------------------             *
*  Polynomial coefficients:                                               *
*  a[0] =     0.000089    0.000000                                        *
*  a[1] =    -0.043680    0.000000                                        *
*  a[2] =     2.994800    0.000000                                        *
*  a[3] =    -6.798000    0.000000                                        *
*  a[4] =     1.000000    0.000000                                        *
*  Roots (without scaling)                                                *
*    No    Real part   Imaginary part    Function value                   *
*     0    0.002454    0.000000     2.05562294698295E-0015                *
*     1    0.012573    0.000000     2.05556701633865E-0015                *
*     2    0.457319    0.000000     2.07678617766666E-0015                *
*     3    6.325654   -0.000000     5.01993719840853E-0014                *
*                                                                         *
*  Roots (with scaling)                                                   *
*    No    Real part   Imaginary part    Function value                   *
*     0    0.002454    0.000000     2.05563872297159E-0015                *
*     1    0.012573    0.000000     2.05571547209755E-0015                *
*     2    0.457319    0.000000     2.20598139995895E-0015                *
*     3    6.325654   -0.000000     5.01993719840853E-0014                *
*                                                                         *
*  ----------------------------------------------------------             *
*                                                                         *
* ----------------------------------------------------------------------- *
* Ref.: "Numerical algorithms with C, By Gisela Engeln-Muellges and       *
*        Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].                  *
*                                                                         *
*                                       TPW Release By J-P Moreau, Paris. *
*                                               (www.jpmoreau.fr)         *
**************************************************************************}
Program TBauhube;

Uses WinCrt, FBauhube;

Var
    ar, ai, rootr, rooti, val: VEC;
    i, j, n, rc, skala: Integer;
  

BEGIN


  Writeln;
  Writeln('  ----------------------------------------------------------');
  Writeln('  Complex and Real Roots of a Polynomial (Bauhuber''s method)');
  Writeln('  ----------------------------------------------------------');

  n := 4;  {order of polynomial}

  {define ar vector}
  ar[0] :=  0.000089248;
  ar[1] := -0.04368;
  ar[2] :=  2.9948;
  ar[3] := -6.798;
  ar[4] :=  1.0;

  {ai vector is null}
  For i:=0 to n do ai[i]:=0.0;

  Writeln('  Polynomial coefficients:');
  for i := 0 to n do
  begin
    write('  a[',i,'] = ');
    write(ar[i]:12:6);
    writeln(ai[i]:12:6)
  end;

  for j := 0 to 1 do
  begin
    skala := j;
    rc := bauhub (0, skala, n, ar, ai, rootr, rooti, val);
    if rc = 0 then
    begin
      if skala = 0 then
        writeln('  Roots (without scaling)')
      else
        writeln('  Roots (with scaling)');

      writeln('    No    Real part   Imaginary part    Function value');

      for i := 0 to n-1 do
      begin
        write('  ', i:4);
        write(rootr[i]:12:6);
        write(rooti[i]:12:6);
        writeln('    ',val[i])
      end;
      writeln
    end
    else
      writeln('  *** Error in bauhube, rc=', rc)
  end;

  Writeln('  ----------------------------------------------------------');
  ReadKey;
  DoneWinCrt

END.

{end of file tbauhube.pas}