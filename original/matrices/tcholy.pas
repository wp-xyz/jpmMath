{*************************************************************
*   Program to test Cholesky Method (fcholy) for symmetric   *
*              and positive definite matrices                *
*                                                            *
*                     Pascal Version By J-P Moreau, Paris    *
*                              (www.jpmoreau.fr)             *
* ---------------------------------------------------------- *
* Reference:                                                 *
*                                                            *
*   "Numerical Algorithms with C, By Gisela Engeln-Muellges  *
*    and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].     *
* ---------------------------------------------------------- *
* SAMPLE RUN:                                                *
*                                                            *
* Input file: matsym.dat                                     *
*                                                            *
* 4                                                          *
* 5.0 -1.0 -1.0 -1.0  1.0                                    *
*      5.0 -1.0 -1.0  1.0                                    *
*           5.0 -1.0  1.0                                    *
*                5.0  1.0                                    *
*                                                            *
* Output file: matsym.lst                                    *
*                                                            *
* ---------------------------------------------------------- *
*      Solving a linear system by Cholesky Method            *
* ---------------------------------------------------------- *
*  The matrix of left hand coefficients must be a symmetric  *
*  positive definite matrix (else error code = 2)            *
*                                                            *
*  Size of system: 4                                         *
*                                                            *
*  Input Linear System:                                      *
*                                                            *
*   5.000000   -1.000000   -1.000000   -1.000000    1.000000 *
*                                                            *
*  -1.000000    5.000000   -1.000000   -1.000000    1.000000 *
*                                                            *
*  -1.000000   -1.000000    5.000000   -1.000000    1.000000 *
*                                                            *
*  -1.000000   -1.000000   -1.000000    5.000000    1.000000 *
*                                                            *
*  Solution:                                                 *
*                                                            *
*   0.500000    0.500000    0.500000    0.500000             *
* ---------------------------------------------------------- *
*                                                            *
*************************************************************}
PROGRAM Test_Choly; 
Uses WinCrt,Basis_r,Fcholy;

{ Only diagonal and upper-diagonal elements are read }

VAR
  a, org : pMAT;   {types defined in unit fcholy.pas}
  b, x   : pVEC;
  
  tmp    : DOUBLE;
  n, i, j, cas, rc : INTEGER;
  
  fp1, fp2 : TEXT;
  input, output, s : STRING[20];


{main program}
BEGIN
  writeln;  
  write(' Input data file name (without *.dat): '); read(s); 
  
  input := s + '.dat';
  output:= s + '.lst';

  Assign(fp1,input); Reset(fp1);
  Assign(fp2,output); Rewrite(fp2);

  WriteHead(fp2,'     Solving a linear system by Cholesky Method');
  writeln(fp2,' The matrix of left hand coefficients must be a symmetric');
  writeln(fp2,' positive definite matrix (else error code = 2)');
  writeln(fp2);

  read(fp1,n);

  if n < 1 then
  begin
    LogError(' n must be > 0 !');
    rc := 1;
    halt(0)
  end;

  {dynamic allocations}
  New(a); New(org); New(b); New(x);

  for i := 0 to n-1 do
  begin
    {read left hand coefficients of line i}
    for j := i to n-1 do
    begin
      read(fp1,tmp);
      org^[i][j] := tmp;    
      org^[j][i] := tmp
    end;
    {read right hand coefficient of line i}
    readln(fp1,b^[i])
  end;
  close(fp1);                     
     
  writeln(fp2,' Size of system: ', n);
  writeln(fp2);
  writeln(fp2,' Input Linear System:');
  writeln(fp2);

  for i := 0 to n-1 do
  begin
    for j := 0 to n-1 do
    begin
      a^[i][j] := org^[i][j];
      write(fp2,a^[i][j]:12:6)
    end;
    writeln(fp2,b^[i]:12:6);
    writeln(fp2)
  end;

  cas := 0;

  {call Cholesky procedure}
  Choly( cas, n, a, b, x, rc ); 
    
  if rc = 0 then
  begin
    writeln(fp2,' Solution:');
    writeln(fp2);
    for i:=0 to n-1 do write(fp2,x^[i]:12:6)
  end
  else
    writeln(fp2,' ERROR choly: code = ', rc);

  Writeln(fp2);
  WriteEnd(fp2);
  close(fp2);

  Dispose(a); Dispose(org); Dispose(b); Dispose(x);

  writeln;
  writeln(' Results in file ',output);
  readkey; donewincrt
END.

{end of file tcholy.pas}