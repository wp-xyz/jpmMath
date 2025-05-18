{************************************************************************
*           EVALUATE A COMPLEX POLYNOMIAL BY HORNER'S METHOD            *
* --------------------------------------------------------------------- *
* SAMPLE RUN:                                                           *
*                                                                       *
*  EVALUATE A CPMPLEX POLYNOMIAL BY HORNER'S METHOD                     *
*                                                                       *
*  Example: P(Z) = (1+i) Z^4 + (2-i) Z^3 +(3+0.5i) Z^2 + (4-2i) Z       *
*                + (5-0.75i)                                            *
*                                                                       *
*  For Z = ( 1.50000000000000E+0000, 2.00000000000000E+0000)            *
*                                                                       *
*  P(Z) = (-2.89375000000000E+0001,-6.93750000000000E+0000)             *
*                                                                       *
*                                                                       *
*                                TPW Release 1.0 By J-P Moreau, Paris   *
*                                         (www.jpmoreau.fr)             *
************************************************************************}  
PROGRAM TEST_CHORNER;

Uses WinCrt;

Const  NMAX = 25;

Type
       {complex number}
       COMPLEX = Record
                   x,y: REAL;  {algebraic form}
                 End;

       pVEC = ^CVEC;
       CVEC  = Array[1..NMAX] of COMPLEX;

Var
       A: pVEC;
       Z, Y0: COMPLEX;
       I,J,K,N:Integer;

    {multiply two complex numbers}
    Procedure ZMult(z1,z2: COMPLEX; VAR z: COMPLEX);
    Begin
      z.x:=(z1.x*z2.x)-(z1.y*z2.y);
      z.y:=(z1.x*z2.y)+(z1.y*z2.x)
    End;

    Procedure CHorner(A:pVec; N: integer; X0: COMPLEX; Var Y0:COMPLEX);
    Var i: integer; tmp: COMPLEX; 
    Begin
      Y0.x := A^[N+1].x; Y0.y := A^[N+1].y;
      For i:=N Downto 1 do
      begin
        ZMult(Y0,X0,tmp);
        Y0.x := tmp.x + A^[i].x;
        Y0.y := tmp.y + A^[i].y;
      end
    End;


{main program}
BEGIN

  N := 4;      {ORDER OF POLYNOMIAL}

  New(A);

  {define complex coefficients}
  A^[5].x := 1.0; A^[5].y :=  1.0;    {power 4}
  A^[4].x := 2.0; A^[4].y := -1.0;    {power 3}
  A^[3].x := 3.0; A^[3].y :=  0.5;    {power 2}
  A^[2].x := 4.0; A^[2].y := -2.0;    {power 1}
  A^[1].x := 5.0; A^[1].y := -0.75;   {power 0}

  Z.x := 1.5;   {complex argument}
  Z.y := 2.0;

  writeln;
  writeln('  EVALUATE A COMPLEX POLYNOMIAL BY HORNER''S METHOD');
  writeln;
  
  CHORNER(A, N, Z, Y0);

  writeln('  For Z = (',Z.x,',',Z.y,')');
  writeln;
  writeln('  P(Z) = (',Y0.x,',',Y0.y,')');
  writeln;

  Dispose(A);
  ReadKey;
  DoneWinCrt

END.

{end of file chorner.pas}