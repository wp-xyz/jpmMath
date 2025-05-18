{**********************************************************
  *           Test Gauss method for band matrix             *
  *                                                         *
  *                   Pascal version by J-P Moreau, Paris   *
  *                           (www.jpmoreau.fr)             *
  * ------------------------------------------------------- *
  * Reference:                                              *
  *                                                         *
  * "Numerical Algorithms with C, By Gisela Engeln-Muellges *
  *  and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].    *
  * ------------------------------------------------------- *
  * SAMPLE RUN:                                             *
  *                                                         *
  * Input file: bandmat.dat                                 *
  *                                                         *
  * 5                                                       * 
  * 1                                                       *
  * 1                                                       *
  *  0.0   1.0   -2.0                                       *
  * 10.0   5.0   -7.5                                       *
  *  2.22 -3.25   0.00456                                   *
  *  2.0   4.0   -3.3                                       *
  * -1.0   3.0    0.0                                       *
  *                                                         *
  * Output file: bandmat.lst                                *
  *                                                         *
  * ------------------------------------------------------- *
  * Banded matrix                                           *
  * ------------------------------------------------------- *
  * Dimension     : 5                                       *
  * Subdiagonals  : 1                                       *
  * Superdiagonals: 1                                       *
  *                                                         *
  * Input condensed band matrix:                            *
  *  0.000000  1.000000 -2.000000                           *
  * 10.00000   5.000000 -7.500000                           *
  *  2.220000 -3.250000  0.004560                           *
  *  2.000000  4.000000 -3.300000                           *
  * -1.000000  3.000000  0.000000                           *
  * Input uncondensed band matrix:                          *
  *  1.000000 -2.000000  0.000000  0.000000  0.000000       *
  * 10.00000   5.000000 -7.500000  0.000000  0.000000       *
  *  0.000000  2.220000 -3.250000  0.004560  0.000000       *
  *  0.000000  0.000000  2.000000  4.000000 -3.300000       *
  *  0.000000  0.000000  0.000000 -1.000000  3.000000       *
  * ------------------------------------------------------- *
  * Band without pivot:                                     *
  * ------------------------------------------------------- *
  * Inverse (transposed):                                   *
  * -0.005941 -0.502971 -0.343236  0.236714  0.078905       *
  *  0.100594  0.050297  0.034324 -0.023671 -0.007890       *
  * -0.231916 -0.115958 -0.386526  0.266570  0.088857       *
  *  3.65e-04  1.82e-04  6.08e-04  0.344408  0.114803       *
  *  4.01e-04  2.01e-04  6.69e-04  0.378849  0.459616       *
  * Determinant = -5.62703999998979E+0002                   *
  * ------------------------------------------------------- *
  * Band with pivot                                         *
  * ------------------------------------------------------- *
  * Inverse (transposed):                                   *
  * -0.005941 -0.502971 -0.343236  0.236714  0.078905       *
  *  0.100594  0.050297  0.034324 -0.023671 -0.007890       *
  * -0.231916 -0.115958 -0.386526  0.266570  0.088857       *
  *  3.65e-04  1.82e-04  6.08e-04  0.344408  0.114803       *
  *  4.01e-04  2.01e-04  6.69e-04  0.378849  0.459616       *
  * Determinant =  5.62703999998979E+0002                   *
  * System Solution:                                        *
  * -0.136497 -0.568249 -0.694162  1.202870  0.734290       *
  * ------------------------------------------------------- *
  *                                                         *
  * Uses units WinCrt, Basis_r,fband and fbando.            *
  **********************************************************}

Program Test_band;
Uses WinCrt, Basis_r, fband, fbando;

{-----------------------------------------------------------------
  Nota: Matrices are put in line in a vector of type VECT
        defined in basis_r.pas as: VECT  = Array[0..SIZE] of REAL;
                         1 2 3
        Example : Matrix 4 5 6 is stored as [1,2,3,4,5,6,7,8,9]
                         7 8 9
 -----------------------------------------------------------------}                                                            

VAR                  
  ud, ld, n, i, k, rc, dsign, dim, mode : INTEGER;
  a            : pVECT;     {pointer to real matrix n*n}
  packmat,save : pVECT;     {pointer to real matrix n*dim}
  b            : pVECT;     {pointer to real vector n}
  perm         : IVECT;     {static  integer vector n}
  determ       : REAL;      {matrix determinant}
  solution, matinv : pVECT; {pointers to real matrices n*n}
  fp1, fp2 : TEXT;          {input,output files}


BEGIN {main}
  {dynamic allocation of matrices}
  New(a); New(packmat); New(save); New(b); New(solution); New(matinv);
  if matinv=NIL then LogError('Memory full!');

  {open input and output files}
  Assign(fp1,'bandmat.dat'); Reset(fp1);
  Assign(fp2,'bandmat.lst'); Rewrite(fp2);
           
  WriteHead(fp2,' Banded matrix');

  Readln(fp1,n);
  if n < 1 then LogError('n must be > 0');

  Readln(fp1,ld);
  if ld < 0 then LogError('ld must be > 0');

  Readln(fp1,ud);
  if ud < 0 then LogError('ud must be > 0');

  Writeln(fp2,' Dimension     : ', n );
  Writeln(fp2,' Subdiagonals  : ', ld);
  Writeln(fp2,' Superdiagonals: ', ud);
  Writeln(fp2);

  if ld + ud >= n then LogError('ld + ud must be < n');

  dim := ld + ud + 1 + min (ld, ud);

  SetVec(n*(ld+ud+1),packmat,0.0);
  {read input matrix in condensed form}
  ReadMat (fp1, n, ld + ud + 1, ld + ud + 1, packmat);
  {close input file}
  close(fp1);

  Copy_vector (save,packmat,n*(ld + ud + 1));

  Writeln(fp2,' Input condensed band matrix:');

  WriteMat(fp2, n, ld + ud + 1, ld + ud + 1, packmat);

  SetVec(n*n,a,0.0);
  for i:=0 to n-1 do
    for k:=0 to n-1 do
      if ((ld+k-i>-1) and (ld+k-i<ld+ud+1)) then
        a^[i*n+k]:=packmat^[i*(ld+ud+1)+ld+k-i]
      else
        a^[i*n+k]:=0.0;

  Writeln(fp2,' Input uncondensed band matrix:');
  
  WriteMat (fp2, n, n, n, a);

  WriteHead(fp2,' Band without pivot:');
  Writeln (fp2,' Inverse (transposed): ');

  mode := 0;
  for i := 0 to n-1 do
  begin
    SetVec (n, b, 0.0);
    b^[i] := 1.0;
    {call band method without pivot}
    bando (mode, n, ld, ud, packmat, b, rc);
    if (rc<>0) then LogError ('bando rc<>0');
    WriteVec (fp2, n, n, b);
    mode := 2
  end;

  determ := 1.0;
  for i := 0 to n-1 do determ := determ * packmat^[i*(ld+ud+1)+ld];

  Writeln (fp2,' Determinant = ',determ);
{-------------------------------------------------------------}
  WriteHead(fp2,' Band with pivot');

  Writeln(fp2,' Inverse (transposed):');

  mode := 0;
  for i := 0 to n-1 do
  begin
    SetVec (n, b, 0.0);
    b^[i] := 1.0;
    {call band method with pivot}
    band (mode, n, ld, ud, save, b, perm, dsign, rc);

    if rc<>0 then LogError ('band rc<>0');
    WriteVec(fp2, n, n, b);
    
    if n<=SIZE then
    begin
      solution^[i]:=0.0;
      for k:=0 to n-1 do matinv^[i*n+k]:=b^[k]
    end;  
    mode := 2
  end;

  determ := 1.0*dsign;
  for i := 0 to n-1 do determ := determ * packmat^[i*(ld+ud+1)+ld];

  Writeln(fp2,' Determinant = ',determ);
  
  if n<=SIZE then
  begin
    Writeln(fp2,' System Solution:');
    for k:=0 to n-1 do
      for i:=0 to n-1 do
        solution^[k]:= solution^[k] + matinv^[i*n+k];
    WriteVec(fp2, n, n, solution) 
  end;     
                     
  WriteEnd(fp2); 
  close(fp2);

  Writeln;
  Writeln(' Results in bandmat.lst.');             
  Readkey;
  DoneWinCrt {close program window}
END.
                                           
{end of file tband.pas}