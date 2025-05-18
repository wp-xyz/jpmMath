{************************************************************
* This program calculates the determinant of a real square  *
* matrix using a recursive function Deter based on Kramer's *
* formula.                                                  *
* --------------------------------------------------------- *
* Method:                                                   *
*                                                           *
* Let us take the example:                                  *
* n=5                                                       *
* Input matrix is:   ( 3 2 1 6 3 )                          *
*                    ( 6 4 4 0 3 )                          *
*                    ( 4 0 1 0 2 )                          *
*                    ( 3 4 2 4 2 )                          *
*                    ( 4 6 1 5 2 )                          *
*                                                           *
* According to Kramer's rule, the determinant det(A) is     *
* calculated by the following:                              *
*                                                           *
*              (4 4 0 3)         (6 4 0 3)       (6 4 0 3)  *
* det = 3 * det(0 1 0 2) -2 * det(4 1 0 2) + 1 * (4 0 0 2)  *
*              (4 2 4 2)         (3 2 4 2)       (3 4 4 2)  *
*              (6 1 5 2)         (4 1 5 2)       (4 6 5 2)  *
*                                                           *
*              (6 4 4 3)         (6 4 4 0)                  *
*      -6 * det(4 0 1 2) +3 * det(4 0 1 0)                  *
*              (3 4 2 2)         (3 4 2 4)                  *
*              (4 6 1 2)         (4 6 1 5)                  *
*                                                           *
* As can be seen, 3, -2, 1, -6 and 3 are the elements of    *
* the 1st row alternatively multiplied by 1 or -1 and the   *
* the submatrices are obtained by successively eliminating  *
* row=1 and column=1, then row=1 and col=2, etc. Then you   *
* can recursively calculate the detrminants of these sub-   *
* matrices by using the same Function Deter with size n-1   *
* instead of n, and so forth...                             *
* When n=2, the result is immediate by the formula:         *
* A(1,1)*A(2,2)-A(2,1)*A(1,2).                              *
*                                                           *
* Program organization:                                     *
* After input of given square matrix A, here of size n=5,   *
* we call the function Deter(n,A) which returns the value   *
* of determinant. This function calls itself recursively    *
* with n decreasing by one at each step, until n=2.         *
* Auxiliary function is pow(n) returning -(-1)^n and auxi-  *
* liary procedure is Submatrix(n,A,B) that extracts sub-    *
* matrix B from A by eliminating row 1 and column col. Note *
* that B is dynamically allocated and disposed of at each   *
* step to spare memory.                                     *
*                                                           *
* SAMPLE RUN:                                               *
*                                                           *
* Size of matrix: 5                                         *
* Input matrix  :                                           *
*   3.0000  2.0000  1.0000  6.0000  3.0000                  *
*   6.0000  4.0000  4.0000  0.0000  3.0000                  *
*   4.0000  0.0000  1.0000  0.0000  2.0000                  *
*   3.0000  4.0000  2.0000  4.0000  2.0000                  *
*   4.0000  6.0000  1.0000  5.0000  2.0000                  *
*                                                           *
* Determinant = -7.00000000000000E+0001                     *
*                                                           *
* --------------------------------------------------------- *
*                                                           *
*                TPW release by Jean-Pierre Moreau, Paris.  *
*                           (www.jpmoreau.fr)               *
************************************************************}
Program Deter2;
Uses WinCrt;

Const SIZE = 25;

Type
      pMat  = ^Mat;
      Mat   = Array[1..SIZE,1..SIZE] of real;

Var
      n  : Integer;   {size of matrix}
      A  : pMat;      {pointer to matrix SIZExSIZE}
      det: real;      {value of determinant}


{return -(-1)^n}
Function pow(n:integer): integer;
Begin
  If Odd(n) then pow:=1 else pow:=-1
End;

{return submatrix without 1st line and colth column of A}
Procedure Submatrix(n,col:integer;A:pMat;VAR B:pMat);
var i,j:integer;
    flag:boolean;
Begin
  if n<3 then exit;
  for i:=1 to n-1 do
  begin
    flag:=FALSE;
    for j:=1 to n-1 do
      if j<>col then
      begin
        if Not flag then B^[i,j]:=A^[i+1,j]
                    else B^[i,j]:=A^[i+1,j+1]
      end
      else
      begin
        Flag:=TRUE;
        B^[i,j]:=A^[i+1,j+1]
      end
  end
End;

{print square matrix A(n,n) }
Procedure Matprint(n:integer;A:pMat);
Var i,j:integer;
Begin
  for i:=1 to n do
  begin
    for j:=1 to n do
      write(' ',A^[i,j]:8:4);
    writeln
  end
End;

{return determinant of matrix A}
Function Deter(n:integer; A:pMat): real;
Var i:integer;
    B:pMat;
    d:real;
Begin
  d:=0.0;
  if n=2 then
    Deter:=A^[1,1]*A^[2,2]-A^[1,2]*A^[2,1]
  else
  begin
    For i:=1 to n do
    begin
      New(B);   {to store submatrix of size n-1 x n-1}
      if n>2 then
      begin
        Submatrix(n,i,A,B);
        {recursive call of Deter}
        d := d + pow(i)*A^[1][i]*Deter(n-1,B)
      end;
      Dispose(B)
    end;
    Deter:=d
  end
End;

{main program}
BEGIN

  New(A);
  {define size n and input matrix A}
  n:=5;
  A^[1,1]:=3.0; A^[1,2]:=2.0; A^[1,3]:=1.0; A^[1,4]:=6.0; A^[1,5]:=3.0;
  A^[2,1]:=6.0; A^[2,2]:=4.0; A^[2,3]:=4.0; A^[2,4]:=0.0; A^[2,5]:=3.0;
  A^[3,1]:=4.0; A^[3,2]:=0.0; A^[3,3]:=1.0; A^[3,4]:=0.0; A^[3,5]:=2.0;
  A^[4,1]:=3.0; A^[4,2]:=4.0; A^[4,3]:=2.0; A^[4,4]:=4.0; A^[4,5]:=2.0;
  A^[5,1]:=4.0; A^[5,2]:=6.0; A^[5,3]:=1.0; A^[5,4]:=5.0; A^[5,5]:=2.0;

  writeln;
  writeln(' Size of matrix: ', n);
  writeln(' Input matrix  :');
  Matprint(n, A);

  det := Deter(n,A);

  writeln;
  writeln(' Determinant = ', det);

  ReadKey; DoneWinCrt

END.

{end of file deter2.pas}