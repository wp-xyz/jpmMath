{******************************************************
*  PROGRAM TO DEMONSTRATE ONE DIMENSIONAL OPERATION   *
*    OF THE MULTI-NONLINEAR REGRESSION SUBROUTINE     *
* --------------------------------------------------- *
*   Ref.: BASIC Scientific Subroutines Vol. II,       *
*   By F.R. Ruckdeschel, Byte/McGRAW-HILL, 1981 [1].  *
*                                                     *
*           Pascal version by J-P Moreau, Paris       *
*           with possibility to have a graph.         *
*                   (www.jpmoreau.fr)                 *
* --------------------------------------------------- *
* SAMPLE RUN:                                         *
*                                                     *
* MULTI-NON LINEAR REGRESSION                         *
*                                                     *
* Number of data points (mini 3, maxi 25)........: 11 *
* Degree of the polynomial to be fitted (maxi 10): 4  *
*                                                     *
* Input the data in (x y) pairs as prompted:          *
*                                                     *
*  1   X Y = 1 1                                      *
*  2   X Y = 2 2                                      *
*  3   X Y = 3 3                                      *
*  4   X Y = 4 4                                      *
*  5   X Y = 5 5                                      *
*  6   X Y = 6 6                                      *
*  7   X Y = 7 7                                      *
*  8   X Y = 8 8                                      *
*  9   X Y = 9 9                                      *
*  10   X, Y = 0 0                                    *
*  11   X, Y = 1 1                                    *
*                                                     *
* The calculated coefficients are:                    *
*                                                     *
*  0    -8.09e-13                                     *
*  1     1.000000                                     *
*  2    -7.18e-13                                     *
*  3     9.11e-14                                     *
*  4    -5.13e-15                                     *
*                                                     *
* Standard deviation:  5.84e-12                       *
*                                                     *
* Do you want a graph (y/n): y                        *
*                                                     *
* (a graph is displayed)                              *
*                                                     *
******************************************************}
PROGRAM LEASTSQR;
USES WinCrtMy, WinProcs, Strings, Type_def, Graph_2D;

CONST   MMAX = 25; NMAX = 10;

VAR
        D, X, Y : RV;        {pointers to real vectors size 2048}
        i, m, n : INTEGER;
        dd,dx,xmn,xmx,xx,ymn,ymx,yy : real_ar;

        Z : Array[1..MMAX,1..NMAX+1] of real_ar;   {maxi m=25 data points}
        A : Array[1..MMAX,1..MMAX] of real_ar;     {maxi n=10 order of regression}
        B : Array[1..MMAX,1..2*MMAX] of real_ar;
        C : Array[1..MMAX,1..MMAX] of real_ar;

        s1,s2 : STRING[12];
        buf   : STRING;
        titre : Array[0..40] of char;


{Matrix transpose}
PROCEDURE Transpose;
Var i,j : integer;
Begin
  FOR i := 1 TO n DO
    FOR j := 1 TO m DO
      B[i, j] := A[j, i]
End;

{Matrix save (A in C) }
PROCEDURE A_IN_C(n1,n2:INTEGER);
Var i1,i2 : integer;
Begin
  IF n1 * n2 <> 0 THEN 
    FOR i1 := 1 TO n1 DO
      FOR i2 := 1 TO n2 DO
        C[i1, i2] := A[i1, i2]
End;

{Matrix save (B in A) }
PROCEDURE B_IN_A(n1,n2:INTEGER);
Var i1,i2 : integer;
Begin
  IF n1 * n2 <> 0 THEN 
    FOR i1 := 1 TO n1 DO
      FOR i2 := 1 TO n2 DO
        A[i1, i2] := B[i1, i2]
End;

{Matrix save (C in B) }
PROCEDURE C_IN_B(n1,n2:INTEGER);
Var i1,i2 : integer;
Begin
  IF n1 * n2 <> 0 THEN 
    FOR i1 := 1 TO n1 DO
      FOR i2 := 1 TO n2 DO
        B[i1, i2] := C[i1, i2]
End;

{Matrix save (C in A) }
PROCEDURE C_IN_A(n1,n2:INTEGER);
Var i1,i2 : integer;
Begin
  IF n1 * n2 <> 0 THEN 
    FOR i1 := 1 TO n1 DO
      FOR i2 := 1 TO n2 DO
        A[i1, i2] := C[i1, i2]
End;

{Matrix multiplication}
PROCEDURE Matmult(m1,n1,n2:INTEGER);
Var i,j,k : integer;
Begin
  FOR i := 1 TO m1 DO
    FOR j := 1 TO n2 DO
    begin
      C[i, j] := 0;
      FOR k := 1 TO n1 DO
        C[i, j] := C[i, j] + A[i, k] * B[k, j]
    end
End;

{Matrix inversion}
PROCEDURE Matinv;
Label 10,20,30;
Var i,j,k : integer; bb:real_ar;
Begin
  FOR i := 1 TO n DO
  begin
    FOR j := 1 TO n DO
    begin
      B[i, j + n] := 0;
      B[i, j] := a[i, j]
    end;
    B[i, i + n] := 1
  end;
  FOR k := 1 TO n DO
  begin
    IF k = n THEN GOTO 10;
    m := k;
    FOR i := k + 1 TO n DO
      IF ABS(B[i, k]) > ABS(B[m, k]) THEN m := i;
    IF m = k THEN GOTO 10;
    FOR j := k TO 2 * n DO
    begin
      bb := B[k, j];
      B[k, j] := B[m, j];
      B[m, j] := bb
    end;
10: FOR j := k + 1 TO 2 * n DO
      B[k, j] := B[k, j] / B[k, k];
    IF k = 1 THEN GOTO 20;
    FOR i := 1 TO k - 1 DO
      FOR j := k + 1 TO 2 * n DO
        B[i, j] := B[i, j] - B[i, k] * B[k, j];
    IF k = n THEN GOTO 30;
20: FOR i := k + 1 TO n DO
      FOR j := k + 1 TO 2 * n DO
        B[i, j] := B[i, j] - B[i, k] * B[k, j]
  end; {k loop}
30: FOR i := 1 TO n DO
      FOR j := 1 TO n DO
        B[i, j] := B[i, j + n]
End; {Matinv}

{Coefficient matrix generation}
PROCEDURE Gene_Coeffs;
Var i,j : integer;
    b : real_ar;
Begin
  FOR i := 1 TO m DO
  begin
    b := 1;
    FOR j := 1 TO n + 1 DO
    begin
      Z[i, j] := b;
      b := b * X^[i]
    end
  end
END;

{Least squares multidimensional fitting subroutine}
PROCEDURE LstSqrN;
Var  i,j,m4,n4 : integer;

Begin
  m4 := m;
  n4 := n;
  FOR i := 1 TO m DO
    FOR j := 1 TO n DO
      A[i, j] := Z[i, j];
  Transpose;                             {b:=Transpose(a) }
  A_IN_C(m,n);                           {move A to C}
  B_IN_A(n,m);                           {move B to A}
  C_IN_B(m,n);                           {move C to B}
  Matmult(n,m,n);                        {multiply A and B}
  C_IN_A(n,n);                           {move C to A}
  Matinv;                                {b:=Inverse(a) }
  m := m4;                               {restore m}
  B_IN_A(n,n);                           {move B to A}
  FOR i := 1 TO m DO
    FOR j := 1 TO n DO
      B[j, i] := Z[i, j];
  Matmult(n,n,m);                        {multiply A and B}
  C_IN_A(n,m);                           {move C to A}
  FOR i := 1 TO m DO B[i, 1] := Y^[i];
  Matmult(n,m,1);                        {multiply A and B}
  {Product C is N by 1 - Regression coefficients are in C(I,1) }
  FOR i := 1 TO n DO D^[i] := C[i, 1]
End;

{Standard deviation calculation for a polynomial fit}
PROCEDURE Standard_deviation;
Var i, j : integer; b, yy : real_ar;
Begin
  dd := 0;
  FOR i := 1 TO m DO
  begin
    yy := 0; b := 1;
    FOR j := 1 TO n + 1 DO
    begin
      yy := yy + D^[j] * b;
      b := b * X^[i]
    end;
    dd := dd + (yy - Y^[i]) * (yy - Y^[i])
  end;
  IF m - n - 1 > 0 THEN
    dd := dd / (m - n - 1)
  else
    dd := 0;
  dd := SQRT(dd)
End;


{main program}
BEGIN
  WinCrtInit('LSTSQR');        {open window with title}
  ClrScr;                      {clear screen}
  New(D); New(X); New(Y);      {allocate real vectors}
  for i:=1 to 11 do D^[i]:=0;
  Writeln(' MULTI-NONLINEAR REGRESSION');
  Writeln;
  Write(' Number of data points (mini 3, maxi 25)........: '); read(m);
  if m<3 then m:=3; if m>25 then m:=25;
  Write(' Degree of the polynomial to be fitted (maxi 10): '); read(n);
  if n<1 then n:=1; if n>10 then n:=10;
  if m<n+2 then m:=n+2;
  Writeln;
  Writeln(' Input the data in (x y) pairs as prompted');
  FOR i := 1 TO m  DO
  begin
    Write(' ',i,'  X  Y = '); read(X^[i], Y^[i])
  end;
  Writeln;
  {Call coefficients generation subroutine}
  Gene_Coeffs;
  {Call regression subroutine}
  n := n + 1;
  LstSqrN;
  Writeln(' The calculated coefficients are:');
  FOR i := 1 TO n DO
    Writeln(' ',i-1,'  ',D^[i]:9:6);
  n := n - 1;
  {Call standard deviation subroutine}
  Standard_deviation;
  Writeln;
  Writeln(' Standard deviation: ',dd:9:6);
  Writeln;
  Write(' Do you want a graph (y/n) : ');
  if ReadKey='y' then
  begin     {graph section}
    Clrscr;
    MinMax(m,X,xmn,xmx);
    MinMax(m,Y,ymn,ymx);
    InitFenetre(CrtDC,10,xmn,xmx,ymn,ymx);
    {Data points are represented by crosses}
    FOR m := 1 TO m DO
      CroixXY(CrtDC,X^[m],Y^[m]);
    if n=1 then
    begin
      {Draw regression line} 
      MoveXY(CrtDC,xmn,D^[2]*xmn+D^[1]);
      LineXY(CrtDC,xmx,D^[2]*xmx+D^[1])
    end
    else if n=2 then  
    begin
      {Draw parabola d(1)+D(2)*x+D(3)*x² using 2*MMAX points}
      MoveXY(CrtDC,xmn,D^[3]*xmn*xmn+D^[2]*xmn+D^[1]);
      dx:=(xmx-xmn)/(2*MMAX-1); xx:=xmn;
      FOR m:=1 to 2*MMAX DO
      begin
        xx:=xx+dx;
        LineXY(CrtDC,xx,D^[3]*xx*xx+D^[2]*xx+D^[1])
      end
    end
    else
    begin
      {Draw curve D(1)+D(2)*x+D(3)*x^2+D(4)*x^3+...+D(11)*x^10 using 4*MMAX points}
      yy:=D^[11]*Power(xmn,10)+D^[10]*Power(xmn,9)+D^[9]*Power(xmn,8);
      yy:=yy+D^[8]*Power(xmn,7)+D^[7]*Power(xmn,6)+D^[6]*Power(xmn,5);
      yy:=yy+D^[5]*Power(xmn,4)+D^[4]*Power(xmn,3)+D^[3]*xmn*xmn+D^[2]*xmn+D^[1];
      MoveXY(CrtDC,xmn,yy);
      dx:=(xmx-xmn)/(4*MMAX-1); xx:=xmn;
      FOR m:=1 to 4*MMAX DO
      begin
        xx:=xx+dx;
        yy:=D^[11]*Power(xx,10)+D^[10]*Power(xx,9)+D^[9]*Power(xx,8);
        yy:=yy+D^[8]*Power(xx,7)+D^[7]*Power(xx,6)+D^[6]*Power(xx,5);
        yy:=yy+D^[5]*Power(xx,4)+D^[4]*Power(xx,3)+D^[3]*xx*xx+D^[2]*xx+D^[1];
        LineXY(CrtDC,xx,yy)
      end
    end;
    {Prepare title of graph}
    Str(D^[1]:10:4,s1); Str(D^[2]:10:4,s2);
    if D^[2] > 0 then
      buf := ' Y = '+s1+' +'+s2+' X'
    else
      buf := ' Y = '+s1+s2+' X';
    if n > 1 then
    begin
      Str(D^[3]:10:4,s2);
      if D^[3] > 0 then
        buf := buf+' +'+s2+' X^2'
      else
        buf := buf+s2+' X^2'
    end;
    if n > 2 then buf:=buf+' +...';
    StrPCopy(titre,buf);
    Legendes(CrtDC,titre,'X','Y');
    Str(n:2,s1); buf:='N='+s1; StrPCopy(titre,buf);
    TextOut(CrtDC,550,40,titre,4);
    SortieGraphique
  end;
  {free memory and exit}
  Dispose(D); Dispose(X); Dispose(Y);
  DoneWinCrt
END.

{End of file lstsqr.pas}













