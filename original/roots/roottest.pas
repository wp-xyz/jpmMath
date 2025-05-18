{****************************************************
*       ROOT TESTING PROGRAM OF A POLYNOMIAL        *
* ------------------------------------------------- *
*   Réf.: Basic Scientific Subroutines Vol. II      *
*   By F.R. Ruckdeschel. BYTE/McGRAW-HILL, 1981     *
*   [BIBLI 01].                                     *
*                                                   *
*                   Pascal Version By J-P Moreau.   *
*                        (www.jpmoreau.fr)          *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
* (Example: Where to look for roots of polynomial   *
*           X^5+5*X^4+10*X^3+10*X^2+5*X+1)          *
*                                                   *
*   This program will help you determine            *
*  where to look for roots of a polynomial.         *
*                                                   *
*  What is the degree of the polynomial: 5          *
*                                                   *
*  Input the polynomial coefficients as prompted:   *
*                                                   *
*     A( 0) = 1                                     *
*     A( 1) = 5                                     *
*     A( 2) = 10                                    *
*     A( 3) = 10                                    *
*     A( 4) = 5                                     *
*     A( 5) = 1                                     *
*                                                   *
*  There are 5 roots.                               *
*                                                   *
*  The magnitude of the largest root is <= 2.236068 *
*                                                   *
*  There is at least one negative real root.        *
*                                                   *
*  There are at most 5 negative real roots.         *
*                                                   *
****************************************************}
PROGRAM RootTest;
Uses WinCrt;

LABEL   100,200,300,fin;

VAR
        A  : Array[0..10] of DOUBLE;
        aa : DOUBLE;
        b,c,d,i,n : INTEGER;

BEGIN
  clrscr;
  writeln;
  writeln('  This program will help you determine   ');
  writeln(' where to look for roots of a polynomial.');
  writeln;
  write(' What is the degree of the polynomial: '); read(n);
  writeln;
  writeln(' Input the polynomial coefficients as prompted:');
  writeln;
  for i:=0 to n do
  begin
    write('   A(',i,') = '); read(A[i])
  end;
  writeln;
  writeln;
  writeln(' There are ',n,' roots.');
  writeln;
  {Find the maximum value of root}
  aa:=A[n-1]*A[n-1]-2.0*A[n-2];
  if aa>0 then goto 100;
  writeln(' There are at least two complex roots.');
  goto fin;
100: aa:=sqrt(aa);
  writeln;
  writeln(' The magnitude of the largest root is <= ',aa:9:6);
  writeln;
  aa:=-1.0;
  if (n MOD 2)=0 then aa:=-aa;
  aa:=A[0]/aa;
  {b will flag a negative root}
  b:=0;
  if aa>0 then goto 200;
  writeln (' There is at least one negative real root.');
  b:=1;
200: writeln;
  {Test for Descartes rule nø 1}
  c:=0;
  for i:=1 to n do
    if A[i-1]*A[i]<0 then c:=c+1;
  if c=1 then write(' There is at most one positive ');
  if c>1 then write(' There are at most ',c,' positive ');
  if c=1 then writeln('real root.');
  if c>1 then writeln('real roots.');
  writeln;
  {Test for Descartes rule nø 2}
  d:=0;
  for i:=1 to n do
    if A[i-1]*A[i]>0 then d:=d+1;
  if d=1 then write(' There is at most one negative ');
  if d>1 then write(' There are at most ',d,' negative ');
  if d=1 then writeln('real root.');
  if d>1 then writeln('real roots.');
  writeln;
  if (c MOD 2)=0 then goto 300;
  writeln(' There is at least one positive real root.');
300: if (d MOD 2)=0 then goto fin;
  if b=0 then
    writeln(' There is at least one negative real root.');
fin: writeln;
  ReadKey; DoneWinCrt
END.

{End of file roottest.pas}