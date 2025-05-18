{******************************************************
* Program to demonstrate multi-dimensional operation  *
* of the multi-nonlinear regression subroutine.       *
* --------------------------------------------------- *
*  Ref.:  BASIC Scientific Subroutines Vol. II,       *     
*         By F.R. Ruckdeschel, Byte/McGRAW-HILL, 1981 *
*         [BIBLI 01].                                 *
*                                                     *
*              Pascal version by J-P Moreau, Paris    *
*                       (www.jpmoreau.fr)             *
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
* The calculated coefficients are:                    *
*                                                     *
*  1    1.95e-14                                      *
*  2    1.000000                                      *
*  3    1.83e-14                                      *
*  4    1.000000                                      *
*  5    1.69e-14                                      *
*  6    -2.80e-15                                     *
*                                                     *
* Standard deviation:  6.49e-13                       *
*                                                     *
******************************************************}
USES WinCrtMy,Type_def,Graph_2D;

CONST   LMAX = 9; MMAX = 25; NMAX = 25;

VAR
        D, Y : RV;     {pointers to real vectors size 2048}
        i,j,l,m,n : INTEGER;
        bb,cc,dd : real_ar;

        X : Array[1..MMAX,1..LMAX] of real_ar;     {maxi l:=9  number of dimensions}
        Z : Array[1..MMAX,1..NMAX] of real_ar;     {maxi m:=25 number of data points}
        A : Array[1..MMAX,1..MMAX] of real_ar;     {maxi n:=25 order of regression}
        B : Array[1..MMAX,1..2*MMAX] of real_ar;
        C : Array[1..MMAX,1..MMAX] of real_ar;
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

{*********************************************************
*    Least squares fitting subroutine, general purpose   *
*  subroutine for multidimensional, nonlinear regression *
* ------------------------------------------------------ *
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

{*************************************************************
* Coefficient matrix generation subroutine for multiple non- *
* linear regression. Also calculates the standard deviation  *
* d, even though there is some redundant computing.          *
* The maximum number of dimensions is 9.                     *
* The input data set consists of m data sets of the form:    *
*   Y(i),X(i,1),X(i,2) ... X(i,l)                            *
* The number of dimensions is l.                             *
* The order of the fit to each dimension is M(j).            *
* The result is an (m1+1)(m2+1)...(ml+1)+1 column by m row   *
* matrix, Z. This matrix is arranged as follows              *
* (Ex.: l=2,M(1)=2,M(2)=2):                                  *
* 1 X1 X1*X1 X2 X2*X1 X2*X1*X1 X2*X2 X2*X2*X1 X2*X2*X1*X1    *
* This matrix should be dimensioned in the calling program   *
* as should also the X(i,j) matrix of data values.           *
*************************************************************}
PROCEDURE Gene_Coeffs_SD;
Label fin;
Var i,j,k : integer; yy : real_ar;

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

{main program}
BEGIN
  WinCrtInit('MLTNLREG');
  New(D); New(Y);
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

  {Call coefficients generation subroutine}
  Gene_Coeffs_SD;

  {Call regression subroutine}
  LstSqrN;

  Writeln(' The calculated coefficients are:');
  Writeln;
  FOR i := 1 TO n DO
    Writeln('  ',i,'  ',D^[i]:9:6); 
  {Get standard deviation}
  n := n - 1;
  Gene_Coeffs_SD;
  Writeln;
  Writeln(' Standard deviation: ',dd:9:6);
  Writeln;
  Dispose(D); Dispose(Y);
  Readkey; DoneWinCrt
END.

{End of file mltnlreg.pas}











