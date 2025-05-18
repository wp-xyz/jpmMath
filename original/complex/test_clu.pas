{*********************************************************************
*    Solving a complex linear system AX = B by LU decomposition      *
*                                                                    *
*                                Pascal version by J-P Moreau, Paris *
*                                         (www.jpmoreau.fr)          *
* ------------------------------------------------------------------ *
* Uses:  units Wincrt, Clu.pas.                                      *
*                                                                    *
* SAMPLE RUN:                                                        *
*                                                                    *  
* Input file (test_clu.dat):                                         *
*                                                                    *
* 15                                                                 *
* 0 0 0 0 6 0 1 0 1 0 4 0 2 0 1 0 2 0 2 0 0 0 3 0 0 0 5 0 0 1  1 1   *
* 6 0 2 0 5 0 2 0 4 0 5 0 5 0 2 0 1 0 2 0 3 0 1 0 5 0 1 0 3 2  2 2   *
* 6 0 2 0 5 0 6 0 3 0 6 0 5 0 0 0 0 0 1 0 3 0 0 0 4 0 0 0 5 3  3 3   *
* 5 0 4 0 3 0 1 0 4 0 4 0 6 0 4 0 6 0 1 0 4 0 2 0 0 0 5 0 3 4  4 4   *
* 3 0 4 0 6 0 4 0 2 0 0 0 4 0 6 0 5 0 5 0 1 0 1 0 6 0 5 0 4 5  5 5   *
* 4 0 2 0 0 0 4 0 1 0 6 0 1 0 1 0 3 0 6 0 1 0 1 0 1 0 1 0 4 6  6 6   *
* 1 0 1 0 4 0 5 0 5 0 0 0 0 0 0 0 0 0 5 0 5 0 2 0 4 0 3 0 5 7  7 7   *
* 2 0 5 0 5 0 0 0 0 0 0 0 0 0 1 0 4 0 1 0 3 0 1 0 2 0 5 0 2 8  8 8   *
* 4 0 4 0 6 0 2 0 3 0 2 0 0 0 5 0 5 0 3 0 5 0 2 0 3 0 4 0 6 9  9 9   *
* 3 0 2 0 0 0 0 0 3 0 4 0 0 0 6 0 5 0 4 0 6 0 1 0 4 0 5 0 3 10 10 10 *
* 1 0 1 0 6 0 2 0 0 0 6 0 0 0 2 0 0 0 0 0 6 0 3 0 4 0 0 0 4 11 11 11 *
* 6 0 3 0 2 0 3 0 6 0 2 0 0 0 2 0 3 0 0 0 5 0 5 0 1 0 6 0 3 12 12 12 *
* 2 0 1 0 5 0 4 0 0 0 4 0 0 0 6 0 3 0 2 0 3 0 1 0 1 0 4 0 6 13 13 13 *
* 3 0 0 0 5 0 1 0 1 0 6 0 0 0 4 0 3 0 6 0 1 0 0 0 1 0 0 0 2 14 14 14 *
* 0 0 3 0 3 0 3 0 4 0 5 0 0 0 3 0 1 0 2 0 4 0 6 0 4 0 1 0 4 15 15 15 *
*                                                                    *
* Output file (test_clu.lst):                                        *
*                                                                    *
*  System solution:                                                  *
*                                                                    *
*  X1= ( -0.057831,  0.092402)                                       *
*  X2= ( -0.885623,  1.415035)                                       *
*  X3= ( -0.106641,  0.170390)                                       *
*  X4= ( -0.421776,  0.673908)                                       *
*  X5= (  0.063331, -0.101190)                                       *
*  X6= ( -0.098265,  0.157007)                                       *
*  X7= (  0.217855, -0.348086)                                       *
*  X8= ( -0.698291,  1.115719)                                       *
*  X9= (  0.869955, -1.390002)                                       *
*  X10= ( -0.338608,  0.541023)                                      *
*  X11= ( -0.463159,  0.740028)                                      *
*  X12= (  0.118877, -0.189941)                                      *
*  X13= (  0.514975, -0.822819)                                      *
*  X14= (  0.013271, -0.021204)                                      *
*  X15= (  0.731170, -1.168252)                                      *
*                                                                    *
**********************************************************************
NOTE: The right-hand side is the two last columns on input file.     }

  Program Test_Clu;

  Uses WinCrt, Clu;

  Var
  A : pCVEC;      { matrix 0:n x 0:n stored in a vector }
  B : pCVEC;      { vecteur 0:n }
  temp : pCVEC;   { vecteur 0:n+1 }
  INDX : pIVEC;   { vecteur 0:n }

  {NOTA : index zero is not used here.}

  d, i, j, n, n2, rc : integer;
  input, output, s : STRING;
  F1, F2 : TEXT;


  Procedure ReadVec1 (VAR fp: TEXT; n, ic: INTEGER; VAR x:pCVEC);
 {====================================================================*
 *                                                                    *
 * Read a complex vector of length x from input file fp. Index starts *
 * at ONE.                                                            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *                                                                    *
 *      fp       pointer to input file.                               *
 *      n        lenght of vector of INTEGER type.                    *
 *      ic       number of items per line (INTEGER).                  *
 *      x        pointer to vector of pVECT type.                     *
 *                                                                    *
 *====================================================================}
  Var
    i : INTEGER;
  Begin
    for i := 1 to n do
    begin
      read(fp,x^[i].R,x^[i].I);
      if (i+1 MOD ic) = 0 then readln(fp)
    end
  End;

  Procedure WriteVec1(VAR fp:TEXT; n, ic:INTEGER; VAR x:pCVEC);
 {====================================================================*
 *                                                                    *
 * Put out vector of length x to  output file. Index starts at ONE.   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *                                                                    *
 *      fp       pointer to output file.                              *
 *      n        lenght of vector of INTEGER type.                    *
 *      ic       number of items per line (INTEGER).                  *
 *      x        pointer to vector of pVECT type.                     *
 *                                                                    *
 *=====================================================================
    ic = number of items per line }
  Var
    i,compte : INTEGER;
  Begin
    compte:=0;
    for i := 1 to n do
    begin
      Inc(compte);
      write(fp,' (',x^[i].R:8:3,',',x^[i].I:8:3,')');
      if compte = ic then
      begin
        writeln(fp);
        compte:=0
      end
    end
  End;


  {main program}
  Begin

  Writeln;
  Write(' Data file name (without .dat): '); read(s);
  input := s + '.dat';
  output := s + '.lst';

  Assign(F1,input); Reset(F1);
  Assign(F2,output); Rewrite(F2);

  readln(F1,n);     { taille du système à résoudre }

  n2:=n+1;

  New(A); New(B); New(temp); New(INDX);

  Writeln(F2);
  Writeln(F2,' COMPLEX LINEAR SYSTEM TO BE SOLVED:');
  writeln(F2,'  N=',n);
  writeln(F2);

  for i:=1 to n do
  begin
    ReadVec1(F1,n2,n2,temp);  {read a line}
    for j:=1 to n do
      A^[i*(n+1)+j] := temp^[j];
    B^[i] := temp^[n+1];
    WriteVec1(F2,n2,n2,temp)  {write a line}
  end;
  close(F1);

  LUDCMP(A,n,INDX,D,rc);

  if rc=0 then LUBKSB(A,n,INDX,B);

  if rc=1 then
    writeln(F2,' The system matrix is singular !')
  else
  begin
    writeln(F2);
    writeln(F2,'  System solution:');
    writeln(F2);  
    for i:=1 to n do
    begin
      write(F2,'  X',i,'=');
      write(F2,' (',B^[i].R:10:6,',',B^[i].I:10:6,')');
      writeln(F2)
    end   
  end;
  Close(F2);
  Writeln;
  Writeln(' Results in file ',output,'.');
  Writeln;
  Readkey;
  DoneWinCrt

  END.

{End of file test_clu.pas}