{********************************************************************
*               Test program for Gauss Seidel method                *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
* (Solve a linear system by Gauss Seidel iterative method).         *
*                                                                   *
* Input file tseidel2.dat contains:                                 *
*                                                                   *
* 16                                                                *
*  4 -1  0  0  -1  0  0  0   0  0  0  0   0  0  0  0                *
* -1  4 -1  0   0 -1  0  0   0  0  0  0   0  0  0  0                *
*  0 -1  4 -1   0  0 -1  0   0  0  0  0   0  0  0  0                *
*  0  0 -1  4   0  0  0 -1   0  0  0  0   0  0  0  0                *
* -1  0  0  0   4 -1  0  0  -1  0  0  0   0  0  0  0                *
*  0 -1  0  0  -1  4 -1  0   0 -1  0  0   0  0  0  0                *
*  0  0 -1  0   0 -1  4 -1   0  0 -1  0   0  0  0  0                *
*  0  0  0 -1   0  0 -1  4   0  0  0 -1   0  0  0  0                *
*  0  0  0  0  -1  0  0  0   4 -1  0  0  -1  0  0  0                *
*  0  0  0  0   0 -1  0  0  -1  4 -1  0   0 -1  0  0                *
*  0  0  0  0   0  0 -1  0   0 -1  4 -1   0  0 -1  0                *
*  0  0  0  0   0  0  0 -1   0  0 -1  4   0  0  0 -1                *
*  0  0  0  0   0  0  0  0  -1  0  0  0   4 -1  0  0                *
*  0  0  0  0   0  0  0  0   0 -1  0  0  -1  4 -1  0                *
*  0  0  0  0   0  0  0  0   0  0 -1  0   0 -1  4 -1                *
*  0  0  0  0   0  0  0  0   0  0  0 -1   0  0 -1  4                *
*  2                                                                *
*  1                                                                *
*  1                                                                *
*  2                                                                *
*  1                                                                *
*  0                                                                *
*  0                                                                *
*  1                                                                *
*  1                                                                *
*  0                                                                *
*  0                                                                *
*  1                                                                *
*  2                                                                *
*  1                                                                *
*  1                                                                *
*  2                                                                *
*                                                                   *
* See full results in output file tseidel2.lst.                     *
*                                                                   *
*  Number of iterations: 30                                         *
*                                                                   *
*  Solution Vector:                                                 *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*   1.000000                                                        *
*                                                                   *
* ----------------------------------------------------------------- *
* Ref.: "Numerical algorithms with C, By Gisela Engeln-Muellges and *
*        Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].            *
*                                                                   *
*                                 TPW Release By J-P Moreau, Paris. *
*                                         (www.jpmoreau.fr)         *
* ----------------------------------------------------------------- *
* Uses file:  Fseidel.pas                                           *
********************************************************************}
PROGRAM TEST_SEIDEL;

Uses WinCrt, Fseidel;

Label fin;

Var
    fp: TEXT;   {I/O files}
    a: pMAT;
    b,x,residu: pVEC;
    omega: Double;
    i,j,n,rc,iter,krit: Integer;

BEGIN

  omega := 1.5;        {relaxation factor}
  krit:=0;          {no special criterion}
  
  {open input file}
  Assign(fp,'tseidel2.dat'); Reset(fp);         

  Readln(fp, n);
  
  if n < 1 then
  begin
    Writeln;
    Writeln(fp,' Dimension must be > 0');
    ReadKey;
    DoneWinCrt
  end;

  New(a); New(b); New(x); New(residu);

  for i:=1 to n do x^[i]:=ZERO;

  For i:=1 to n do
  begin
    for j:=1 to n-1 do Read(fp,a^[i,j]);
    Readln(fp,a^[i,n])
  end;

  For i:=1 to n do Readln(fp,b^[i]);
 
  Close(fp);
                           
  {open output file}                         
  Assign(fp,'tseidel2.lst'); Rewrite(fp);  

  Writeln(fp);
  Writeln(fp,'-----------------------------------------------------------------------');
  Writeln(fp,'    Gauss Seidel method');
  Writeln(fp,'-----------------------------------------------------------------------');
  Writeln(fp,' Dimension of the input matrix = ', n);
  Writeln(fp);
  Writeln(fp,' Input matrix:');
  For i:=1 to n do
  begin
    for j:=1 to n-1 do Write(fp,a^[i,j]:6:2);
    Writeln(fp,a^[i,n]:6:2)
  end;
  Writeln(fp);
  Writeln(fp,' Second member:');
  For i:=1 to n-1 do Write(fp,b^[i]:6:2);
  Writeln(fp,b^[n]:6:2);

  Writeln(fp);
  Writeln(fp,' Solution vector:');

  seidel(krit, n, a, b, omega, x, residu, iter, rc);

  if rc <> 0 then
  begin
    writeln(fp,' Error in seidel: rc <> 0 !');
    goto fin
  end;

  For i:=1 to n-1 do Write(fp,x^[i]:10:6);
  Writeln(fp,x^[n]:10:6);

  Writeln(fp);
  Writeln(fp,' Residual vector (must be near zero):');
  For i:=1 to n-1 do Write(fp,residu^[i]:10:6);
  Writeln(fp,residu^[n]:10:6);

  Writeln(fp);
  Writeln(fp,' Number of iterations: ', iter);
  Writeln(fp,'-----------------------------------------------------------------------');

fin: Writeln;
  Writeln(' Results in tseidel2.lst...');
  Writeln;  

  ReadKey;
  Close(fp);

  Dispose(a); Dispose(b); Dispose(x); Dispose(residu);
  DoneWinCrt

END.              

{ J-P Moreau April, 2005 }