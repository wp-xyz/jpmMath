{***************************************************
*    Program to demonstrate Evaluating elliptic    *
*  integrals of first and second kinds (complete)  *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*           Pascal Version by J.-P. Moreau, Paris. *
*                     (www.jpmoreau.fr)            *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
*    K         K(K)          E(K)      STEPS       *
*  -------------------------------------------     *
*   0.00     1.5707963     1.5707963     1         *
*   0.05     1.5717795     1.5698141     3         *
*   0.10     1.5747456     1.5668619     3         *
*   0.15     1.5797457     1.5619230     3         *
*   0.20     1.5868678     1.5549685     3         *
*   0.25     1.5962422     1.5459573     3         *
*   0.30     1.6080486     1.5348335     3         *
*   0.35     1.6225281     1.5215252     3         *
*   0.40     1.6399999     1.5059416     4         *
*   0.45     1.6608862     1.4879683     4         *
*   0.50     1.6857504     1.4674622     4         *
*   0.55     1.7153545     1.4442435     4         *
*   0.60     1.7507538     1.4180834     4         *
*   0.65     1.7934541     1.3886864     4         *
*   0.70     1.8456940     1.3556611     4         *
*   0.75     1.9109898     1.3184721     4         *
*   0.80     1.9953028     1.2763499     4         *
*   0.85     2.1099355     1.2281083     4         *
*   0.90     2.2805490     1.1716970     4         *
*   0.95     2.5900112     1.1027216     5         *
*   1.00     INFINITY      1.0000000     0         *
*                                                  *
***************************************************}
PROGRAM CLIPTIC;
Uses WinCrt;

Var
        A, B : Array[0..100] of DOUBLE;
        e,e1,e2,xk : DOUBLE;
        i,n : integer;

{*****************************************************
* Complete elliptic integral of the first and second *
* kind. The input parameter is xk, which should be   *
* between 0 and 1. Technique uses Gauss' formula for *
* the arithmogeometrical mean. e is a measure of the * 
* convergence accuracy. The returned values are e1,  *
* the elliptic integral of the first kind, and e2,   *
* the elliptic integral of the second kind.          *
* -------------------------------------------------- *
* Reference: Ball, algorithms for RPN calculators.   *
*****************************************************}
PROCEDURE CElliptic;
Label 100,fin;
Var   j, m : integer;
Begin
  A[0]:=1.0+xk ; B[0]:=1-xk;
  n:=0;
  if xk < 0 then goto fin;
  if xk > 1 then goto fin;
  if e <= 0 then goto fin;
100: n:=n+1;
  {Generate improved values}
  A[n]:=(A[n-1]+B[n-1])/2.0;
  B[n]:=SQRT(A[n-1]*B[n-1]);
  if ABS(A[n]-B[n]) > e then goto 100;
  e1:=pi/2.0/A[n];
  e2:=2.0;
  m:=1;
  for j:=1 to n do
  begin
    e2:=e2-m*(A[j]*A[j]-B[j]*B[j]);
    m:=m*2
  end;
  e2:=e2*e1/2.0;
fin: End;


{main program}
BEGIN
  clrscr;
  writeln;
  e:=1e-7;
  writeln('  K         K(K)          E(K)      STEPS ');
  writeln('------------------------------------------');
  xk:=0.0;
  for i:=1 to 20 do
  begin
    CElliptic;
    writeln(' ',xk:4:2,'     ',e1:9:7,'     ',e2:9:7,'     ',n);
    xk := xk + 0.05
  end;
  writeln(' 1.00     INFINITY      1.0000000     0');
  ReadKey; DoneWinCrt
END.

{End of file cliptic.pas}