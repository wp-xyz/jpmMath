{******************************************************
* Least Squares of order 1 or 2 Demonstration Program *
* --------------------------------------------------- *
*   Reference: BASIC Scientific Subroutines, Vol. II  *
*   By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1]. *
*                                                     *
*               Pascal version by J-P Moreau, Paris   *
*                        (www.jpmoreau.fr)            *
* --------------------------------------------------- *
* FIRST SAMPLE RUN:                                   *
*                                                     *
* This program calculates a linear or parabolic least *
* squares fit to a given data set.                    *
*                                                     *
* INSTRUCTIONS                                        *
* ------------                                        *
*                                                     *
* The number of data coordinates provided must be     *
* greater than three. Otherwise, a divide by zero     *
* error may result.                                   *
*                                                     *
* Order of fit (1 or 2): 1                            *
* Number of data points: 4                            *
*                                                     *
* There are two input options:                        *
* 1. input coordinate pairs                           *
* 2. first input the independant variable values,     *
*    then input dependant ones.                       *
*    Your choice (1 or 2): 1                          *
*                                                     *
* ? 1 1                                               *
* ? 2 2                                               *
* ? 3 3                                               *
* ? 5 5.01                                            *
*                                                     *
* Fitted equation is:                                 *
*   Y = -0.004571 + 1.002571 X                        *
*                                                     *
* Standard deviation of fit: 0.002928                 *
*                                                     *
* Do you want a chart (y/n): n                        *
*                                                     *
* SECOND SAMPLE RUN:                                  * 
*                                                     *
* Order of fit (1 or 2): 2                            *
* Number of data points: 4                            *
*                                                     *
* There are two input options:                        *
* 1. input coordinate pairs                           *
* 2. first input the independant variable values,     *
*    then input dependant ones.                       *
*    Your choice (1 or 2): 1                          *
*                                                     *
* ? 1 1                                               *
* ? 2 4                                               *
* ? 3 9                                               *
* ? 5 24.95                                           *
*                                                     *
* Fitted equation is:                                 *
*   Y = -0.017727 + 0.022045 X + 0.994318 X^2         *
*                                                     *
* Standard deviation of fit: 0.004767                 *
*                                                     *
* Do you want a chart (y/n): y                        *
*                                                     *
*    (A chart is displaid)                            *
*                                                     *
******************************************************}
PROGRAM Least_squares;
Uses WinCrtMy, Strings, Type_def, Graph_2D;

VAR
        X,Y : RV;
        a,b,c,d : real_ar;
        dx,xmn,xmx,xx,ymn,ymx : real_ar;
        iorder,iz,m,n : integer;
        s1,s2 : STRING[12];
        buf   : STRING;
        titre : Array[0..40] of char;


{********************************************************
*          Linear least squares subroutine              *
* ----------------------------------------------------- *
* In;  integer n := number of points                    *
*      n values x[i), y(i) shared with main             *
* Out; coefficients a,b of fit (a+b*x) shared with main *
*      standard deviation d shared with main            *
********************************************************}
PROCEDURE lstsqr1 (n : INTEGER);
Var a1,a2,b0,b1,d1 : real_ar;
    m : INTEGER;
Begin
  a1 := 0;a2 := 0;b0 := 0;b1 := 0;
  FOR m := 1 TO n DO
  begin
    a1 := a1 + x^[m];
    a2 := a2 + x^[m] * x^[m];
    b0 := b0 + y^[m];
    b1 := b1 + y^[m] * x^[m]
  end;
  a1 := a1 / n;a2 := a2 / n;b0 := b0 / n;b1 := b1 / n;
  d := a1 * a1 - a2;
  a := a1 * b1 - a2 * b0;a := a / d;
  b := a1 * b0 - b1;b := b / d;
  {Evaluation of standard deviation d (unbiased estim%ate)}
  d := 0;
  FOR m := 1 TO n DO
  begin
    d1 := y^[m] - a - b * x^[m];
    d := d + d1 * d1
  end;
  d := SQRT(d / (n - 2))
End;

{****************************************************************
*           Parabolic least squares subroutine                  *
* ------------------------------------------------------------- *
* In;  integer n := number of points                            *
*      n values x[i), y(i) shared with main                     *
* Out; coefficients a,b,c of fit (a+b*x+c*x^2) shared with main *
*      standard deviation d shared with main                    *
****************************************************************}
PROCEDURE lstsqr2 (n : INTEGER);
Var a0,a1,a2,a3,a4,b0,b1,b2,d1 : real_ar;
    m : INTEGER;
Begin
  a0 := 1;a1 := 0;a2 := 0;a3 := 0;a4 := 0;
  b0 := 0;b1 := 0;b2 := 0;
  FOR m := 1 TO n DO
  begin
    a1 := a1 + x^[m];
    a2 := a2 + x^[m] * x^[m];
    a3 := a3 + x^[m] * x^[m] * x^[m];
    a4 := a4 + x^[m] * x^[m] * x^[m] * x^[m];
    b0 := b0 + y^[m];
    b1 := b1 + y^[m] * x^[m];
    b2 := b2 + y^[m] * x^[m] * x^[m]
  end;
  a1 := a1 / n;a2 := a2 / n;a3 := a3 / n;a4 := a4 / n;
  b0 := b0 / n;b1 := b1 / n;b2 := b2 / n;
  d := a0 * (a2 * a4 - a3 * a3) - a1 * (a1 * a4 - a2 * a3) + a2 * (a1 * a3 - a2 * a2);
  a := b0 * (a2 * a4 - a3 * a3) + b1 * (a2 * a3 - a1 * a4) + b2 * (a1 * a3 - a2 * a2);
  a := a / d;
  b := b0 * (a2 * a3 - a1 * a4) + b1 * (a0 * a4 - a2 * a2) + b2 * (a1 * a2 - a0 * a3);
  b := b / d;
  c := b0 * (a1 * a3 - a2 * a2) + b1 * (a2 * a1 - a0 * a3) + b2 * (a0 * a2 - a1 * a1);
  c := c / d;
  {Evaluation of standard deviation d}
  d := 0;
  FOR m := 1 TO n DO
  begin
    d1 := y^[m] - a - b * x^[m] - c * x^[m] * x^[m];
    d := d + d1 * d1
  end;
  d := SQRT(d / (n - 3))
End;

{main program}
BEGIN
  WinCrtInit('LSTSQR2');
  New(X); New(Y);
  Repeat
    CLRSCR;
    Writeln('LEAST SQUARES CURVE FIT ROUTINE');
    Writeln;
    Writeln(' This program calculates a linear or parabolic least');
    Writeln(' squares fit to a given data set');
    Writeln;
    Writeln('INSTRUCTIONS');
    Writeln('------------');
    Writeln;
    Writeln(' The number of data coordinates provided must be greater');
    Writeln(' than three. Otherwise, a divide by zero error may result.');
    Writeln;
    Write  (' Order of fit (1 or 2): '); read(iorder);
    IF iorder < 1 THEN iorder := 1;
    IF iorder > 2 THEN iorder := 2;
    Write  (' Number of data points: '); read(n)
  Until n > 3;

  Repeat
    CLRSCR;
    Writeln(' There are two input options:');
    Writeln(' 1. input coordinate pairs (Example;1,2.5)');
    Writeln(' 2. first input the independant variable values, then input dependant ones.');
    Write  (' Your choice (1 or 2): '); read(iz);
    Writeln
  Until iz IN [1..2];

  IF iz = 2 THEN
  begin
    {read data from screen (option 2)}
    FOR m := 1 TO n DO
    begin
      Write(' ',m,' '); read(x^[m])
    end;
    Writeln;
    FOR m := 1 TO n DO
    begin
      Write(' ',m,' '); read(y^[m])
    end
  end
  else
  begin
    {read data from screen (option 1) }
    FOR m := 1 TO n DO
    begin
      Write(' ',m,' '); read(x^[m], y^[m])
    end
  end;

  {Call linear or parabolic least squares subroutine }
  IF iorder = 1 THEN lstsqr1(n);
  IF iorder = 2 THEN lstsqr2(n);

  {Write results to screen taking care of zero, one or minus one values}
  Writeln;
  Writeln(' FITTED EQUATION IS:');
  Write('   Y = ');
  IF a <> 0 THEN write(a:9:6);
  IF b <> 0 THEN
  begin
    IF (b > 0) AND (a <> 0) THEN Write(' +');
    IF ABS(b) <> 1 THEN
    begin
      write(b:9:6); write(' X')
    end
    ELSE IF b > 0 THEN
      Write(' X')
    ELSE
      Write(' -X')
  end;
  IF iorder > 1 THEN
  begin
    IF c <> 0 THEN
    begin
      IF (c > 0) AND (a <> 0) AND (b <> 0) THEN Write(' +');
      IF ABS(c) <> 1 THEN
      begin
        write(c:9:6); write(' X^2')
      end
      ELSE IF c > 0 THEN
        Write(' X^2')
      ELSE
        Write(' -X^2')
    end
  end; 
  Writeln; Writeln;
  Writeln(' STANDARD DEVIATION OF FIT: ',d:9:6); 
  Writeln;
  Write(' Do you want a chart (y/n): ');
  if ReadKey='y' then
  begin
    Clrscr;
    MinMax(n,X,xmn,xmx);
    MinMax(n,Y,ymn,ymx);
    InitFenetre(CrtDC,10,xmn,xmx,ymn,ymx);
    {Represent data points by crosses}
    FOR m := 1 TO n DO
      CroixXY(CrtDC,X^[m],Y^[m]);
    if iorder=1 then
    begin
      {Draw least squares line} 
      MoveXY(CrtDC,xmn,b*xmn+a);
      LineXY(CrtDC,xmx,b*xmx+a)
    end
    else  {iorder=2}
    begin
      {Draw least squares parabola a+bx+cx² with 50 points}
      MoveXY(CrtDC,xmn,c*xmn*xmn+b*xmn+a);
      dx:=(xmx-xmn)/49; xx:=xmn;
      FOR m:=1 to 50 DO
      begin
        xx:=xx+dx;
        LineXY(CrtDC,xx,c*xx*xx+b*xx+a)
      end
    end;
    {Prepare title}
    Str(a:12:6,s1); Str(b:12:6,s2);
    if b>0 then
      buf := ' Y = '+s1+' +'+s2+' X'
    else
      buf := ' Y = '+s1+s2+' X';
    if iorder=2 then
    begin
      Str(c:12:6,s2);
      if c>0 then
        buf := buf+' +'+s2+' X^2'
      else
        buf := buf+s2+' X^2'
    end;
    StrPCopy(titre,buf);

    Legendes(CrtDC,titre,'X','Y');
    SortieGraphique
  end;
  Dispose(X); Dispose(Y);
  DoneWinCrt
END.

{End lstsqr2.pas }