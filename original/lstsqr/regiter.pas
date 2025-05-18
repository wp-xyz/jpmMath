{******************************************************
*  Program to demonstrate multidimensional operation  *
*  of the multi-nonlinear regression subroutine with  *
*  iterative error reduction                          *                                                 *
* --------------------------------------------------- *
*  Ref.:  BASIC Scientific Subroutines Vol. II, by    *
*         F.R. Ruckdeschel, Byte/McGRAW-HILL, 1981    *
*         [BIBLI 01].                                 *
*                                                     *
*              Pascal Version By J-P Moreau, Paris    *
*                        (www.jpmoreau.fr)            *
* --------------------------------------------------- *
* SAMPLE RUN:                                         *
*                                                     *
* MULTI-DIMENSIONAL AND MULTI-NONLINEAR REGRESSION    *
*                                                     *
* How many data points ? 10                           *
*                                                     *
* How many dimensions ? 2                             *
*                                                     *
* What is the fit for dimension 1 ? 2                 *
* What is the fit for dimension 2 ? 1                 *
*                                                     *
* Input the data as prompted:                         *
*                                                     *
* Y( 1) = ? 7                                         *
* X( 1, 1) = ? 1                                      *
* X( 1, 2) = ? 6                                      *
*                                                     *
* Y( 2) = ? 7                                         *
* X( 2, 1) = ? 6                                      *
* X( 2, 2) = ? 1                                      *
*                                                     *
* Y( 3) = ? 6                                         *
* X( 3, 1) = ? 3                                      *
* X( 3, 2) = ? 3                                      *
*                                                     *
* Y( 4) = ? 8                                         *
* X( 4, 1) = ? 2                                      *
* X( 4, 2) = ? 6                                      *
*                                                     *
* Y( 5) = ? 9                                         *
* X( 5, 1) = ? 1                                      *
* X( 5, 2) = ? 8                                      *
*                                                     *
* Y( 6) = ? 9                                         *
* X( 6, 1) = ? 7                                      *
* X( 6, 2) = ? 2                                      *
*                                                     *
* Y( 7) = ? 6                                         *
* X( 7, 1) = ? 3                                      *
* X( 7, 2) = ? 3                                      *
*                                                     *
* Y( 8) = ? 7                                         *
* X( 8, 1) = ? 3                                      *
* X( 8, 2) = ? 4                                      *
*                                                     *
* Y( 9) = ? 7                                         *
* X( 9, 1) = ? 4                                      *
* X( 9, 2) = ? 3                                      *
*                                                     *
* Y( 10) = ? 2                                        *
* X( 10, 1) = ? 0                                     *
* X( 10, 2) = ? 2                                     *
*                                                     *
*                                                     *
* The calculated coefficients are:                    *
*                                                     *
*  1    0.000000                                      *
*  2    1.000000                                      *
*  3    0.000000                                      *
*  4    1.000000                                      *
*  5    0.000000                                      *
*  6    -0.000000                                     *
*                                                     *
* Standard deviation:  0                              *
*                                                     *
* Number of iterations: 1                             *
*                                                     *
******************************************************}
USES WinCrt;

CONST   LMAX = 9; MMAX = 25; NMAX = 25;

TYPE    pRV = ^Real_vector;
        Real_vector = Array[0..MMAX] of DOUBLE;

VAR
        D, D1, Y, Y1 : pRV;     {pointers to real vectors size MMAX}
        i,j,l,l1,m,n : INTEGER;
        bb,cc,dd,dd1 : DOUBLE;

        X : Array[1..MMAX,1..LMAX] of DOUBLE;     {maxi l:=9  number of dimensions}
        Z : Array[1..MMAX,1..NMAX] of DOUBLE;     {maxi m:=25 number of data points}
        A : Array[1..MMAX,1..MMAX] of DOUBLE;     {maxi n:=25 order of regression}
        B : Array[1..MMAX,1..2*MMAX] of DOUBLE;
        C : Array[1..MMAX,1..MMAX] of DOUBLE;

        M1: Array[1..LMAX] of integer;

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
Var i,j,k : integer; bb:DOUBLE;
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

{*********************************************************
* Least squares fitting subroutine, general purpose      *
* subroutine for multidimensional, nonlinear regression. *
* The equation fitted has the form:                      *
*     Y = D(1)X1 + D(2)X2 + ... + D(n)Xn                 *
* The coefficients are returned by the program in D(i).  *  
* The X(i) can be simple powers of x, or functions.      *      
* Note that the X(i) are assumed to be independent.      *   
* The measured responses are Y(i), there are m of them.  * 
* Y is a m row column vector, Z(i,j) is a m by n matrix. *
* m must be >= n+2. The subroutine inputs are m, n, Y(i) *
* and Z(i,j) previously calculated. The subroutine calls *
* several other matrix routines during the calculation.  *
*********************************************************}
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


{************************************************************
*         Coefficient matrix generation subroutine          *
*            for multiple non-linear regression.            *
* --------------------------------------------------------- *
* Also calculates the standard deviation d, even though     *
* there is some redundant computing.                        *
* The maximum number of dimensions is 9.                    *
* The input data set consists of m data sets of the form:   *
*   Y(i),X(i,1),X(i,2) ... X(i,l)                           *
* The number of dimensions is l.                            *
* The order of the fit to each dimension is M(j).           *
* The result is an (m1+1)(m2+1)...(ml+1)+1 column by m row  *
* matrix, Z. This matrix is arranged as follows             *
* (Ex.:l=2,M(1)=2,M(2)=2):                                  *
* 1 X1 X1*X1 X2 X2*X1 X2*X1*X1 X2*X2 X2*X2*X1 X2*X2*X1*X1   *
* This matrix should be dimensioned in the calling program  *
* as should also the X(i,j) matrix of data values.          *
************************************************************}
PROCEDURE Gene_Coeffs_SD;
Label fin;
Var i,j,k : integer; yy : DOUBLE;

  {internal procedures}
  Procedure S40;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[1] do
    begin
      j := j + 1; Z[i, j] := bb; bb := bb * X[i, 1]
    end;
    bb := cc
  end;

  Procedure S50;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[2] do
    begin
      S40; bb := bb * X[i, 2]
    end;
    bb := cc
  end;

  Procedure S60;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[3] do
    begin
      S50; bb := bb * X[i, 3]
    end;
    bb := cc
  end;

  Procedure S70;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[4] do
    begin
      S60; bb := bb * X[i, 4]
    end;
    bb := cc
  end;

  Procedure S80;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[5] do
    begin
      S70; bb := bb * X[i, 5]
    end;
    bb := cc
  end;

  Procedure S90;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[6] do
    begin
      S80; bb := bb * X[i, 6]
    end;
    bb := cc
  end;

  Procedure S100;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[7] do
    begin
      S90; bb := bb * X[i, 7]
    end;
    bb := cc
  end;

  Procedure S110;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[8] do
    begin
      S100; bb := bb * X[i, 8]
    end;
    bb := cc
  end;

  Procedure S120;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[9] do
    begin
      S110; bb := bb * X[i, 9]
    end;
    bb := cc
  end;


Begin
  {Calculate the total number of dimensions}
  n := 1;
  FOR i := 1 TO l DO
    n := n * (M1[i] + 1);
  dd := 0;
  FOR i := 1 TO m DO
  begin
    {Branch according to dimension l (return if l > 9) }
    IF l < 1 THEN
    begin
      l := 0; goto fin
    end;
    IF l > 9 THEN
    begin
      l := 0; goto fin
    end;

    j := 0;
    IF l = 1 THEN begin bb:=1; S40 end;
    IF l = 2 THEN begin bb:=1; S50 end;
    IF l = 3 THEN begin bb:=1; S60 end;
    IF l = 4 THEN begin bb:=1; S70 end;
    IF l = 5 THEN begin bb:=1; S80 end;
    IF l = 6 THEN begin bb:=1; S90 end;
    IF l = 7 THEN begin bb:=1; S100 end;
    IF l = 8 THEN begin bb:=1; S110 end;
    IF l = 9 THEN begin bb:=1; S120 end;

    yy := 0;
    FOR k := 1 TO n DO
      yy := yy + D^[k] * Z[i, k];
    dd := dd + (Y^[i] - yy) * (Y^[i] - yy)
  end; {i loop}
  {Calculate standard deviation (if m > n) }
  IF m - n < 1 THEN
    dd := 0
  else
  begin
    dd := dd / (m - n);
    dd := SQRT(dd)
  end;

fin: End; {Gene_Coeffs_SD}

{*******************************************************************
*   Multi-dimensional polynomial regression iteration subroutine   *
* ---------------------------------------------------------------- *
* This routine supervises the calling of several other subroutines *
* in order to iteratively fit least squares polynomials in more    * 
* one dimension.                                                   *
* The routine repeatedly calculates improved coefficients until    *
* the standard deviation is no longer reduced. The inputs to the   *
* subroutine are the number of dimensions l, the degree of fit for *
* each dimension m(i), and the input data, x(i) and y(i).          *
* The coefficients are returned in d(i), with the standard devia-  *
* tion in d. Also returned is the number of iterations tried, l1.  *
* y1(i), d1(i) and d1 are used respectively to save the original   *
* values of y(i) and the current values of d(i) and d.             *
*******************************************************************}
PROCEDURE Iter_Supervisor;
Label 50,100,fin;

  {internal procedures}
  Procedure S160;
  var i1:integer;
  begin
    cc := bb;
    FOR i1 := 0 TO M1[1] DO
    begin
      j := j + 1;
      z[i, j] := bb; bb := bb * x[i, 1]
    end;
    bb := cc
  end;

  Procedure S170;
  var i2:integer;
  begin
    cc := bb;
    FOR i2 := 0 TO M1[2] DO
    begin
      S160; bb := bb * x[i, 2]
    end;
    bb := cc
  end;

  Procedure S180;
  var i3:integer;
  begin
    cc := bb;
    FOR i3 := 0 TO M1[3] DO
    begin
      S170; bb := bb * x[i, 3]
    end;
    bb := cc
  end;

  Procedure S190;
  var i4:integer;
  begin
    cc := bb;
    FOR i4 := 0 TO M1[4] DO
    begin
      S180; bb := bb * x[i, 4]
    end;
    bb := cc
  end;

  Procedure S200;
  var i5:integer;
  begin
    cc := bb;
    FOR i5 := 0 TO M1[5] DO
    begin
      S190; bb := bb * x[i, 5]
    end;
    bb := cc
  end;

  Procedure S210;
  var i6:integer;
  begin
    cc := bb;
    FOR i6 := 0 TO M1[6] DO
    begin
      S200; bb := bb * x[i, 6]
    end;
    bb := cc
  end;

  Procedure S220;
  var i7:integer;
  begin
    cc := bb;
    FOR i7 := 0 TO M1[7] DO
    begin
      S210; bb := bb * x[i, 7]
    end;
    bb := cc
  end;

  Procedure S230;
  var i8:integer;
  begin
    cc := bb;
    FOR i8 := 0 TO M1[8] DO
    begin
      S220; bb := bb * x[i, 8]
    end;
    bb := cc
  end;

  Procedure S240;
  var i9:integer;
  begin
    FOR i9 := 0 TO M1[9] DO
    begin
      S230; bb := bb * x[i, 9]
    end;
    bb := cc
  end;

  Procedure S250;
  Var i,j,k : integer; yy : DOUBLE;
  Begin
    FOR i := 1 TO m DO
    begin
      j := 0;
      IF l = 1 THEN begin bb:=1; S160 end;
      IF l = 2 THEN begin bb:=1; S170 end;
      IF l = 3 THEN begin bb:=1; S180 end;
      IF l = 4 THEN begin bb:=1; S190 end;
      IF l = 5 THEN begin bb:=1; S200 end;
      IF l = 6 THEN begin bb:=1; S210 end;
      IF l = 7 THEN begin bb:=1; S220 end;
      IF l = 8 THEN begin bb:=1; S230 end;
      IF l = 9 THEN begin bb:=1; S240 end;
      {Array generated for row i}
      yy := 0;
      FOR k := 1 TO n DO yy := yy + d^[k] * z[i, k];
      y^[i] := y^[i] - yy
    end
  End; {S200}

Begin {Iter_Supervisor}
  l1 := 0;
  {Save the y(i)}
  FOR i := 1 TO m DO y1^[i] := y^[i];
  {Zero d1(i) }
  FOR i := 1 TO n DO d1^[i] := 0;
  {Set the initial standard deviation high}
  dd1 := 1e7;
  {Call coefficients subroutine}
50: Gene_Coeffs_SD;
  {Call regression subroutine}
  LstSqrN;
  {Get standard deviation}
  Gene_Coeffs_SD;
  {If standard deviation is decreasing, continue}
  IF dd1 > dd THEN GOTO 100;
  {Terminate iteration}
  FOR i := 1 TO n DO d^[i] := d1^[i];
  {Restore the y(i)}
  FOR i := 1 TO m DO y^[i] := y1^[i];
  {Get the final standard deviation}
  Gene_Coeffs_SD;
  goto fin;
  {Save the standard deviation}
100: dd1 := dd; l1 := l1 + 1;
  {Augment coefficient matrix}
  FOR i := 1 TO n DO
  begin
    d^[i] := d1^[i] + d^[i];
    d1^[i] := d^[i]
  end;
  {Restore y(i)}
  FOR i := 1 TO m DO y^[i] := y1^[i];
  {Reduce y(i) according to the d(i) }
  S250;
  {We now have a set of error values}
  GOTO 50;

fin: End;  {Iter_Supervisor}


{main program}
BEGIN
  New(D); New(D1); New(Y); New(Y1);
  Writeln(' MULTI-DIMENSIONAL AND MULTI-NONLINEAR REGRESSION');
  Writeln;
  Write(' How many data points ? '); read(m);
  Writeln;
  Write(' How many dimensions ? '); read(l);
  Writeln;
  FOR i := 1 TO l DO
  begin
    Write(' What is the fit for dimension ',i,' ? '); read(M1[i])
  end;
  n := 1;
  FOR i := 1 TO l DO n := n * (M1[i] + 1);
  IF m < n + 2 THEN m := n + 2;

  Writeln;
  Writeln(' Input the data as prompted:');
  Writeln;
  FOR i := 1 TO m DO
  begin
    Write(' Y(',i,') = ? '); read(y^[i]);
    FOR j := 1 TO l DO
    begin
      Write(' X(',i,',',j,') =  ? '); read(X[i, j])
    end;
    Writeln
  end;

  {Call iterations supervisor subroutine}
  Iter_Supervisor;

  Writeln;
  Writeln(' The calculated coefficients are:');
  Writeln;
  FOR i := 1 TO n DO Writeln('  ',i,'  ',D^[i]:9:6); 
  Writeln;
  Write(' Standard deviation: ',dd:11:8);
  Writeln;
  Write(' Number of iterations: ',l1);
  Writeln;
  Dispose(D); Dispose(D1); Dispose(Y); Dispose(Y1);
  Readkey; DoneWinCrt
END.

{End of file REGITER.PAS}
