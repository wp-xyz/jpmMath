{****************************************************
* Calculate the determinant of a real square matrix *
* A(n,n) by Gauss method with full pivoting.        *
* ------------------------------------------------- *
* Ref.: "Algèbre - Algorithmes et programmes en     *
*        Pascal By Jean-Louis Jardrin, Dunod -      *
*        Bordas Paris, 1988 p. 76-79" [BIBLI 10].   *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
* (Calculate the determinant of matrix:             *
*            10 18  1 14 22                         *
*             4 12 25  8 16                         *
*            23  6 19  2 15                         *
*            17  5 13 21  9                         *
*            11 24  7 20  3  )                      *
*                                                   *
* Input size of square real matrix: 5               *
*                                                   *
* Line 1                                            *
* Element 1: 10                                     *
* Element 2: 18                                     *
* Element 3: 1                                      *
* Element 4: 14                                     *
* Element 5: 22                                     *
*                                                   *
* Line 2                                            *
* Element 1: 4                                      *
* Element 2: 12                                     *
* Element 3: 25                                     *
* Element 4:  8                                     *
* Element 5: 16                                     *
*                                                   *
* Line 3                                            *
* Element 1: 23                                     *
* Element 2:  6                                     *
* Element 3: 19                                     *
* Element 4:  2                                     *
* Element 5: 15                                     *
*                                                   *
* Line 4                                            *
* Element 1: 17                                     *
* Element 2:  5                                     *
* Element 3: 13                                     *
* Element 4: 21                                     *
* Element 5:  9                                     *
*                                                   *
* Line 5                                            *
* Element 1: 11                                     *
* Element 2: 24                                     *
* Element 3:  7                                     *
* Element 4: 20                                     *
* Element 5:  3                                     *
*                                                   *
* Determinant = -4.68000000000000E+006              *
*                                                   *
*             TPW version with dynamic allocations  *
*                 By Jean-Pierre Moreau, Paris.     *
*                      (www.jpmoreau.fr)            *
****************************************************}
PROGRAM Calculate_deter;
Uses WinCrt;

Const MAXSIZE = 100;   {this value can be changed by user}

Type
      pMat  = ^Mat;    {pointer to a real square matrix of type Mat}
      Mat   = Array[1..MAXSIZE,1..MAXSIZE] of real;
      pVecI = ^VecI;   {pointer to an integer vector of type VecI}
      VecI  = Array[1..MAXSIZE] of integer;

Var
      n  : integer;      {size of matrix A}
      A  : pMat;         {pointer to input matrix}
      det: real;         {determinant of matrix A}
      eps: real;         {desired precision

The procedure TSRGT applies to input real square matrix A(n,n) the upper
triangularization algorithm of Gauss method with full pivoting and keeps
trace of successive transformations done in integer vectors KP and LP.
----------------------------------------------------------------------------
  Input parameters:
  eps        precision (real)
  n          size of A matrix (integer)
  A          pointer to input real square matrix (pMat)
  Output parameters:
  it         flag=1 if A matrix ok, =0 if A matrix is singular (integer)
  C          pointer to table storing main diagonal elements and supra-
             diagonal elements of upper triangular matrix and the multi-
             plying coefficients used during triangularization process (pMat)
  KP         table storing informations concerning the column exchanges
             during process (pVecI)
  LP         table storing informations concerning the line exchanges
             during process (pVecI)
----------------------------------------------------------------------------
The table C is first initialized to A matrix, then receives at each step k
of the triangularization process, usefull elements of A matrix at step k for
k=1,2,...n.
The variables po(real), lo and ko(integer) store respectively pivot at step k,
its line number and its column number.
Note: types pMat and pVecI must be declared and allocated for by main program.
------------------------------------------------------------------------------}
Procedure TSRGT(eps:real;n:integer;A:pMat;VAR it:integer;VAR C:pMat;VAR Kp,Lp:pVecI);
Var i,j,k,ko,lo:integer; po,t0:real;
Begin
  C:=A; it:=1; k:=1;
  While (it=1) and (k<n) do
  begin
    po:=C^[k,k]; lo:=k;ko:=k;
    for i:=k to n do
      for j:=k to n do
        if abs(C^[i,j])>abs(po) then
        begin
          po:=C^[i,j];lo:=i;ko:=j
        end;
    Lp^[k]:=lo; Kp^[k]:=ko;
    if abs(po)<eps then
      it:=0
    else
    begin
      if lo<>k then
        for j:=k to n do
        begin
          t0:=C^[k,j]; C^[k,j]:=C^[lo,j]; C^[lo,j]:=t0
        end;
      if ko<>k then
        for i:=1 to n do
        begin
          t0:=C^[i,k]; C^[i,k]:=C^[i,ko]; C^[i,ko]:=t0
        end;
      for i:=k+1 to n do
      begin
        C^[i,k]:=C^[i,k]/po;
        for j:=k+1 to n do
          C^[i,j]:=C^[i,j]-C^[i,k]*C^[k,j]
      end;
      Inc(k)
    end
  end;
  if (it=1) and (abs(C^[n,n])<eps) then it:=0
End; {TSRGT}

{The function DMGT returns the determinant of a real square matrix
A(n,n) by Gauss method with full pivoting.
----------------------------------------------------------------------------
  Input parameters:
  eps        precision (real)
  n          size of A matrix (integer)
  A          pointer to input real square matrix (pMat)
  Output parameters:
  None
-----------------------------------------------------------------------------
The procedure TSRGT is used to reduce A matrix to an upper triangular matrix.
Output variables are it(integer), C(pMat), Kp and Lp(pVecI).
If it=0, matrix A is singular, if it=1, matrix A is regular. Table C contains
at location i,j (j>=i) the corresponding element of the upper triangular matrix.
Tables Lp and Kp contain informations relative to exchanges of line or column
that occured during the process. For instance, the element number k of Lp is
an integer <> k if an exchange of line has been made at step k (k=1,2,...,n).
The number of exchanges of lines and columns is stored in L(integer). the
determinant of A matrix is stored in d0(real).
Note: types pMat and pVecI must be declared and allocated for by main program,
except local variables C,Kp,Lp allocated (and disposed of) here.
-----------------------------------------------------------------------------}
Function DMGT(eps:real;n:integer;VAR A:pMat):real;
var it,k,l:integer;d0:real;
    C:pMat; Kp,Lp:pVecI;
Begin
  New(C); New(Kp); New(Lp);   {allocate local matrix and vectors}
  TSRGT(eps,n,A,it,C,Kp,Lp);  {call triangularization procedure}
  if it=0 then
    d0:=0.0  {matrix singular, det=0}
  else
  begin  {matrix regular, det<>0}
    d0:=1.0;
    for k:=1 to n do d0:=d0*C^[k,k];
    l:=0;
    for k:=1 to n-1 do
    begin
      if Lp^[k]<>k then Inc(l);
      if Kp^[k]<>k then Inc(l)
    end;
    if Odd(l) then d0:=-d0
  end;
  DMGT:=d0;  {return determinant}
  Dispose(C); Dispose(Kp); Dispose(Lp)    {free memory}
End;

{input size n, precision eps and elements of matrix A}
Procedure Input_Data;
Var i,j: integer;
Begin
  writeln;
  write(' Input size of square real matrix: '); readln(n);
  write(' Input desired precision: '); readln(eps);
  writeln;
  for i:=1 to n do
  begin
    writeln(' Line ',i);
    for j:=1 to n do
    begin
      write(' Element ',j,': '); readln(A^[i,j])
    end;
    writeln
  end;
End;

{main program}
BEGIN
  writeln;
  writeln(' Calculate the determinant of a real square matrix');
  writeln('        by Gauss method with full pivoting.');
  New(A);     {allocate memory space for input matrix}
  Input_Data; {read n, eps and A matrix}
   
  det:=DMGT(eps,n,A);

  writeln;
  writeln(' Determinant = ',det);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file deter.pas}