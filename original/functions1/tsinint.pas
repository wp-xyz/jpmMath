{*************************************************************
* Program to demonstrate the use of functions sinintegral(x) *
* or cosintegral(x) using the Simpson's method               *
* ---------------------------------------------------------- *
*                                                            *
* SAMPLE RUNS:                                               *
*                                                            *
*         Kind (1 or 2): 1                                   *
*                                                            *
* Test of Function sinintegral(x) by Simpson's method        *
*                                                            *
* Input value for x variable:                                *
*                                                            *
*         X = 2                                              *
*                                                            *
* Result:  1.60541297680275E+0000                            *
*                                                            *
*         X = 5                                              *
*                                                            *
* Result:  1.54993124494464E+0000                            *
*                                                            *
*         X = 80                                             *
*                                                            *
* Result:  1.57233088691248E+0000                            *
*                                                            *
*                                                            *
*         Kind (1 or 2): 2                                   *
*                                                            *
* Test of Function cosintegral(x) by Simpson's method        *
*                                                            *
* Input value for x variable:                                *
*                                                            *
*         X = 2                                              *
*                                                            *
* Result:  4.22980828774809E-0001                            *
*                                                            *
*         X = 5                                              *
*                                                            *
* Result: -1.90029749656719E-0001                            *
*                                                            *
*         X = 80                                             *
*                                                            *
* Result: -1.24025011551291E-0002                            *
*                                                            *
*                                TPW Version By J-P Moreau.  *
*                                    (www.jpmoreau.fr)       *
*************************************************************}
Program Test_Sinintegral;

Uses WinCrt,  {for writeln, readln, etc.}
     Sinint;  {for sinintegral and cosintegral}

Var
    x   : double;        {function argument}
    kind: integer;       {=1 for sinintegral, =2 for cosintegral}

Begin

  writeln;
  write('         Kind (1 or 2): '); Readln(kind);
  writeln;
  if kind = 1 then
    writeln(' Test of Function sinintegral(x) by Simpson''s method')
  else
    writeln(' Test of Function cosintegral(x) by Simpson''s method');
  writeln;
  writeln(' Input value for x variable:');
  writeln;
  write('         X = '); readln(x);
  writeln;
  if (kind=2) and (X<=0.0) then
    writeln(' Warning! X must be positive for kind=2.')
  else
    if kind = 1 then
      writeln(' Result: ', sinintegral(x))
    else
      writeln(' Result: ', cosintegral(x));
  writeln;
  ReadKey;
  DoneWinCrt

End.

{End of file tsinint.pas}