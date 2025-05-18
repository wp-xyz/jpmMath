{********************************************************
* Inversion of a real square matrix by LU decomposition *
* with dynamic allocations                              *
*                                                       *
*                   Pascal version by J-P Moreau, Paris *
*                           (www.jpmoreau.fr)           *
* ----------------------------------------------------- *
* Uses:  units Wincrt, Basis.pas, Lu.pas                *
*                                                       *
* SAMPLE RUN:                                           *
*                                                       *  
* Input file (inv_lu.dat):                              *
*                                                       *
*  4                                                    *
*  8  2    3  12                                        *
*  2  4    7   0.25                                     *
*  3  7    3   5                                        *
* 12  0.25 5   2                                        *
*                                                       *
* Output file (inv_lu.lst):                             *
*                                                       *
* ----------------------------------------------------- *
*  INVERSION OF A REAL SQUARE MATRIX:                   *
* ----------------------------------------------------- *
*  N=4                                                  *
*                                                       *
*  8.000000  2.000000  3.000000  12.00000               *
*  2.000000  4.000000  7.000000  0.250000               *
*  3.000000  7.000000  3.000000  5.000000               *
*  12.00000  0.250000  5.000000  2.000000               *
*                                                       *
*  Inverted matrix Y:                                   *
*                                                       *
* -0.040155 -0.085788  0.056557  0.110262               *
* -0.085788 -0.062018  0.202201  0.016978               *
*  0.056557  0.202201 -0.130316 -0.038826               *
*  0.110262  0.016978 -0.038826 -0.066631               *
*                                                       *
*  Verification A*Y = I:                                *
*                                                       *
*  1.000000 -4.55e-13 -4.55e-13 -4.55e-13               *
*  1.71e-13  1.000000 -1.02e-12 -3.69e-13               *
* -6.82e-13  1.14e-13  1.000000  5.68e-13               *
*  0.000000  2.27e-13 -2.27e-13  1.000000               *
* ----------------------------------------------------- * 
*                                                       *
********************************************************}
  Program Inversion_LU;
  Uses WinCrt,Basis,Lu;

  Var
  A : pVECT;      { matrix 0:n x 0:n stored in a vector (see NOTA2) }
  A1: pVect;      { copy of matrix A }
  Y : pVECT;      { matrix 0:n x 0:n stored in a vector }
  temp : pVECT;   { vector 0:n }
  INDX : pIVECT;  { integer vector 0:n }

  {NOTA1: index zero is not used here.
   NOTA2: The element i,j of matrix A(n,n) is the element i*(n+1)+j
          of the vector A, n being the 2nd dimension of A(n,n) }

  d, i, j, n, rc : integer;
  input, output, s : STRING;
  F1, F2 : TEXT;

  {******************************************                                     
  *       MULTIPLY TWO REAL MATRICES        *
  * --------------------------------------- *                                     
  * INPUTS:    A  MATRIX N*N                *                                     
  *            B  MATRIX N*N                *                                     
  *            N  INTEGER                   *                                     
  * --------------------------------------- *                                     
  * OUTPUT:    C  MATRIX N*N, PRODUCT A*B   *                                     
  *                                         *                                     
  ******************************************}
  Procedure MatMult(A,B:pVect; VAR C:pVect; n:Integer); 
  VAR SUM: REAL;
      I,J,K: INTEGER;
  Begin                                           
    for I:=1 to N do
    begin                                                                  
      for J:=1 to N do
      begin                                                                
        SUM:=0.0;                                                                
        for K:=1 to N do
          SUM:=SUM+A^[I*(n+1)+K]*B^[K*(n+1)+J];                                               
        C^[I*(n+1)+J]:=SUM                                                            
      end                                                                   
    end                                                                     
  End;

{main program}
BEGIN

  Writeln;
  Write(' Data file name (without .dat): '); read(s);
  input := s + '.dat';
  output := s + '.lst';

  Assign(F1,input); Reset(F1);
  Assign(F2,output); Rewrite(F2);

  readln(F1,n);     {size of matrix}

  New(A); New(A1); New(Y); New(temp); New(INDX);

  WriteHead(F2,' INVERSION OF A REAL SQUARE MATRIX:');
  writeln(F2,'  N=',n);
  writeln(F2);

  for i:=1 to n do
  begin
    ReadVec1(F1,n,n,temp);  {read a line}
    for j:=1 to n do
    begin
      A^[i*(n+1)+j] := temp^[j];
      A1^[i*(n+1)+j] := temp^[j];
      Y^[i*(n+1)+j] := 0.0
    end;
    Y^[i*(n+1)+i] := 1.0;
    WriteVec1(F2,n,n,temp)   {write a line}
  end;
  close(F1);

{call LU decomposition routine (only once) }
  LUDCMP(A,n,INDX,D,rc);

{call solver if previous return code is ok
 to obtain inverse of A one column at a time}
  if rc=0 then
    for j:=1 to n do
    begin
      for i:=1 to n do temp^[i]:=Y^[i*(n+1)+j];
      LUBKSB(A,n,INDX,temp);
      for i:=1 to n do Y^[i*(n+1)+j]:=temp^[i]
    end;
  {inverse of matrix A is now in matrix Y,
   matrix A is destroyed. }

  if rc=1 then
    writeln(F2,' The matrix is singular !')
  else
  begin
    writeln(F2);
    writeln(F2,'  Inverted matrix Y:');
    writeln(F2);  
    for i:=1 to n do
    begin
      for j:=1 to n do
        f_aff_reel(F2,Y^[i*(n+1)+j]);
      writeln(F2)
    end
  end;
  {verify A1 x Y = I (result put in A) }
  MatMult(A1,Y,A,n);
  {A should now contain identity matrix}
  writeln(F2);
  writeln(F2,'  Verification A*Y = I:');
  writeln(F2);
  for i:=1 to n do
  begin
    for j:=1 to n do
      f_aff_reel(F2,A^[i*(n+1)+j]);
    writeln(F2)
  end;
  WriteEnd(F2);
  Close(F2);
  Dispose(A); Dispose(A1); Dispose(Y); Dispose(temp); Dispose(INDX);
  Writeln;
  Writeln(' Results in file ',output,'.');
  Writeln;
  Readkey;
  DoneWinCrt

END.

{End of file test_lu.pas}