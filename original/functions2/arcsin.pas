{*******************************************************
*       Program to demonstrate arcsine recursion       *
* ---------------------------------------------------- *
*   Reference: BASIC Scientific Subroutines, Vol. II   *
*   By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                      *
*             Pascal Version by J.-P. Moreau, Paris    *
*                       (www.jpmoreau.fr)              *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
*                                                      *
*   X       ARCSIN(X)       STEPS        ERROR         *
*  -----------------------------------------------     *
*  0.00     0.0000000         0       0.0000000000     *
*  0.05     0.0500209         7      -0.0000000013     *
*  0.10     0.1001674         8      -0.0000000025     *
*  0.15     0.1505683         9      -0.0000000021     *
*  0.20     0.2013579        10      -0.0000000013     *
*  0.25     0.2526803        10      -0.0000000025     *
*  0.30     0.3046927        11      -0.0000000011     *
*  0.35     0.3575711        11      -0.0000000017     *
*  0.40     0.4115168        11      -0.0000000025     *
*  0.45     0.4667653        12      -0.0000000009     *  
*  0.50     0.5235988        12      -0.0000000012     *
*  0.55     0.5823642        12      -0.0000000016     *
*  0.60     0.6435011        12      -0.0000000021     *
*  0.65     0.7075844        13      -0.0000000007     *
*  0.70     0.7753975        13      -0.0000000008     *
*  0.75     0.8480621        13      -0.0000000010     *
*  0.80     0.9272952        13      -0.0000000012     *
*  0.85     1.0159853        13      -0.0000000014     *
*  0.90     1.1197695        14      -0.0000000004     *
*  0.95     1.2532359        14      -0.0000000004     *
*  1.00     1.5707963         0       0.0000000000     *
*                                                      *
*******************************************************}
PROGRAM Arcsinus;
Uses WinCrt;

VAR
        e,x,y : DOUBLE;
        i,m   : integer;


{************************************************
*       Arcsin(x) recursion subroutine          *
* Input is x (-1<x<1), output is y=arcsin(x),   *
* convergence criteria is e.                    *
* --------------------------------------------- *
* Reference: Computational Analysis by Henrici. *
************************************************}
PROCEDURE Arcsin;
Label 100, 200, 300, fin;
Var  u0,u1,u2 : DOUBLE;
Begin
  m:=0;
  {Guard against failure}
  if e<=0 then goto fin;
  if x<>0 then goto 100;
  y:=0.0;
  goto fin;
  {Check range}
100: if abs(x)>1 then goto 300;
  u0:=x*sqrt(1.0-x*x);
  u1:=x;
200: u2:=u1*sqrt(2.0*u1/(u1+u0));
  y:=u2;
  m:=m+1;
  if abs(u2-u1)<e then goto fin;
  u0:=u1; u1:=u2;
  goto 200;
300: if abs(x-1) < 1e-10 then y:=pi/2.0;
fin: End;


{main program}
BEGIN
  clrscr;
  writeln;
  e:=1e-8;
  writeln(' X       ARCSIN(X)       STEPS        ERROR    ');
  writeln('-----------------------------------------------');
  x:=0.0;
  for i:=1 to 21 do
  begin
    Arcsin;
    writeln(x:4:2,'     ',y:9:7,'        ',m:2,'      ',(sin(y)-x):13:10);
    x:=x+0.05
  end;
  ReadKey; DoneWinCrt
END.

{End of file arcsin.pas}
