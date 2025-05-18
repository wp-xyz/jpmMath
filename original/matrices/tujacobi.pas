{*********************************************************************
*        Eigenvalues and eigenvectors of a real symmetric            *
*                square matrix by JACOBI'S METHOD                    *
* ------------------------------------------------------------------ *
*   Uses procedure Jacobi of unit ujacobi.pas.                       *
*                                                                    *
*                                    Pascal version by J-P Moreau    *
*                                          (www.jpmoreau.fr)         *
* ------------------------------------------------------------------ *
* SAMPLE RUN:                                                        *
* Input file (elpro.dat):                                            *
*                                                                    *
*  4                                                                 *
*   4 -2 -1  0                                                       *
*      4  0 -1                                                       *
*         4 -2                                                       *
*            4                                                       *
*                                                                    *
* Output file (elpro.lst):                                           *
*                                                                    *
* ---------------------------------------------------                *
*  Eigenvalues and eigenvectors of a real, symmetric                 *
*     square matrix by Jacobi's iterative method                     *
* ---------------------------------------------------                *
*                                                                    *
*  Dimension of matrix: 4                                            *
*                                                                    *
*  Symmetric matrix is:                                              *
*                                                                    *
*     4.000000    -2.000000    -1.000000     0.000000                *
*    -2.000000     4.000000     0.000000    -1.000000                *
*    -1.000000     0.000000     4.000000    -2.000000                *
*     0.000000    -1.000000    -2.000000     4.000000                *
*                                                                    *
*  Eigenvalues:                                                      *
*     7.00000000                                                     *
*     1.00000000                                                     *
*     3.00000000                                                     *
*     5.00000000                                                     *
*                                                                    *
*  Eigenvectors (in lines):                                          *
*     1.000000    -1.000000    -1.000000     1.000000                *
*     1.000000     1.000000     1.000000     1.000000                *
*     1.000000     1.000000    -1.000000    -1.000000                *
*    -1.000000     1.000000    -1.000000     1.000000                *
*                                                                    *
*  Number of used iterations: 5                                      *
*  Error code: 0                                                     *
* ------------------------------------------------------------------ *
* SHARED VARIABLES WITH UNIT UJACOBI:                                *
*        MAT(N,N)         : given real symmetric square matrix       *
*        Eigenvalues(N)   : eigenvalues stored in unsorted order     *
*        Eigenvectors(N,N): normalized eigenvectors given in lines   *
*********************************************************************}
PROGRAM TEST_UJACOBI;
Uses WinCrt, UJacobi;

Var
    dimen  : integer;        {Dimension of square matrix}
    maxiter: integer;        {Maximum number of Iterations}
    iter   : integer;        {Number of iterations used}
    error  : byte;           {Error code - see unit ujacobi.pas}

    tolerance: float;        {tolerance in answer}

    fp1,fp2: TEXT;           {inpout/output text files}
    i,j    : integer;        {loop variables}
    temp   : float;          {temporary value}

BEGIN
  {initialize pointers declared in unit ujacobi}
  New(MAT); New(Eigenvalues); New(Eigenvectors);
  {open input/output files}
  Assign(fp1,'elpro4.dat'); Reset(fp1);   {read mode}
  Assign(fp2,'elpro4.lst'); Rewrite(fp2); {write mode}
  {read data from input file}
  readln(fp1,dimen);
  for i:=1 to dimen do       {read only upper half triangle}
  begin
    for j:=i to dimen do
    begin
      read(fp1,temp);
      MAT^[i,j]:=temp;
      MAT^[j,i]:=temp
    end;
    readln(fp1)
  end;
  close(fp1);
  {print title and input matrix}
  writeln(fp2,'---------------------------------------------------');
  writeln(fp2,' Eigenvalues and eigenvectors of a real, symmetric');
  writeln(fp2,'    square matrix by Jacobi''s iterative method');
  writeln(fp2,'---------------------------------------------------');
  writeln(fp2);
  writeln(fp2,' Dimension of matrix: ',dimen);
  writeln(fp2);
  writeln(fp2,' Symmetric matrix is:');
  writeln(fp2);
  for i:=1 to dimen do
  begin
    for j:=1 to dimen do
      write(fp2,' ',MAT^[i,j]:12:6);
    writeln(fp2)
  end;
  
  {fix parameters}
  maxiter:=100; tolerance:=1e-8;

  {call jacobi procedure, see unit ujacobi.pas}
  Jacobi(dimen,maxiter,tolerance,iter,error);

  {print results to output file}
  writeln(fp2);
  writeln(fp2,' Eigenvalues:');
  for i:=1 to dimen do
    writeln(fp2,' ',Eigenvalues^[i]:14:8);
  writeln(fp2);
  writeln(fp2,' Eigenvectors (in lines):');
  for i:=1 to dimen do
  begin
    for j:=1 to dimen do
      write(fp2,' ',Eigenvectors^[i,j]:12:6);
    writeln(fp2)
  end;
  writeln(fp2);
  writeln(fp2,' Number of used iterations: ',iter);
  writeln(fp2,' Error code: ',error);
  writeln(fp2,'---------------------------------------------------');
  close(fp2);
  writeln; writeln(' Program done.');
  writeln(' Results in file elpro4.lst.');
  Readkey;
  {free memory used by pointers}
  Dispose(MAT); Dispose(Eigenvalues); Dispose(Eigenvectors);
  {close application}
  DoneWinCrt
END.

{end of file tujacobi.pas}

