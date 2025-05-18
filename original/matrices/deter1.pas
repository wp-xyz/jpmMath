{****************************************************
* Calculate the determinant of a real square matrix *
* by the LU decomposition method                    *
* ------------------------------------------------- *
* SAMPLE RUN:                                       *
* (Calculate the determinant of real square matrix: *
*        | 1  0.5  0  0    0  0    0  0    0 |      *
*        |-1  0.5  1  0    0  0    0  0    0 |      *
*        | 0 -1    0  2.5  0  0    0  0    0 |      *
*        | 0  0   -1 -1.5  4  0    0  0    0 |      *
*        | 0  0    0 -1   -3  5.5  0  0    0 |      *
*        | 0  0    0  0   -1 -4.5  7  0    0 |      *
*        | 0  0    0  0    0 -1   -6  8.5  0 |      *
*        | 0  0    0  0    0  0   -1 -7.5 10 |      *
*        | 0  0    0  0    0  0    0 -1   -9 |  )   *
*                                                   *
*                                                   *
*  Determinant =  1.00000000000000E+0000            *
*                                                   *
*                                                   *
*                 TPW version by J-P Moreau, Paris. *
*                         (www.jpmoreau.fr)         *
* ------------------------------------------------- *
* (In this version, the matrix A(i,j) is stored in  *
*  a vector A of size (n+1)(n+1). The element (i,j) *
*  of the matrix A equals A(i*(n+1)+j) of vector A. *
*  The access to the ith element of a vector is     *
*  faster than the access to element (i,j) of a     *
*  matrix ).                                        * 
*****************************************************
 Exact value is: 1
----------------------------------------------------}
PROGRAM Determinant;
Uses WinCrt, Basis, Lu;

Var
    A: pVect;     {input matrix A(n,n) }
    INDX:pIVECT;  {dummy vector not used here}
    i,icode,id,j,n,n1,n2: Integer;
    D:double;

BEGIN
  New(A); New(INDX);
  n:=9; n1:=n+1; n2:=n1*n1;
  {define matrix A(i,j) column by column stored in vector A line by line}
  for i:=0 to n2 do A^[i]:=0.0;
  j:=1; i:=1; A^[i*n1+j]:=1.0; i:=2; A^[i*n1+j]:=-1.0;
  j:=2; i:=1; A^[i*n1+j]:=0.5; i:=2; A^[i*n1+j]:=0.5; i:=3; A^[i*n1+j]:=-1.0;
  j:=3; i:=2; A^[i*n1+j]:=1.0; i:=4; A^[i*n1+j]:=-1.0;
  j:=4; i:=3; A^[i*n1+j]:=2.5; i:=4; A^[i*n1+j]:=-1.5; i:=5; A^[i*n1+j]:=-1.0;
  j:=5; i:=4; A^[i*n1+j]:=4.0; i:=5; A^[i*n1+j]:=-3.0; i:=6; A^[i*n1+j]:=-1.0;
  j:=6; i:=5; A^[i*n1+j]:=5.5; i:=6; A^[i*n1+j]:=-4.5; i:=7; A^[i*n1+j]:=-1.0;
  j:=7; i:=6; A^[i*n1+j]:=7.0; i:=7; A^[i*n1+j]:=-6.0; i:=8; A^[i*n1+j]:=-1.0;
  j:=8; i:=7; A^[i*n1+j]:=8.5; i:=8; A^[i*n1+j]:=-7.5; i:=9; A^[i*n1+j]:=-1.0;
  j:=9; i:=8; A^[i*n1+j]:=10.0; i:=9; A^[i*n1+j]:=-9.0;

  {call LU decomposition}
  LUDCMP(A, n, INDX, id, icode);

  {calculate determinant D and display result}
  writeln;
  if icode=0 then  {LU Ok}
  begin
    D:=id;
    for i:=1 to n do D:=D*A^[i*n1+i];
    writeln(' Determinant = ', D)
  end
  else
    writeln(' Error in LU decomposition of matrix A.');
  writeln;
  Readkey;
  Dispose(A); Dispose(INDX);
  DoneWinCrt
END.

{end of file deter1.pas}