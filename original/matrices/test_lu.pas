{******************************************************
*  Solving a linear system AX = B by LU decomposition *
*                                                     *
*                 Pascal version by J-P Moreau, Paris *
*                         (www.jpmoreau.fr)           *
* --------------------------------------------------- *
* Uses:  units Wincrt, Basis.pas, Lu.pas              *
*                                                     *
* SAMPLE RUN:                                         *
*                                                     *  
* Input file (test_lu.dat):                           *
*                                                     *
*  4                                                  *
*  8  2    3  12     25.0                             *
*  2  4    7   0.25  13.25                            *
*  3  7    3   5     18.0                             *
* 12  0.25 5   2     19.25                            *
*                                                     *
* Output file (test_lu.lst):                          *
*                                                     *
* --------------------------------------------------- *
*  LINEAR SYSTEM TO BE SOLVED:                        *
* --------------------------------------------------- *
*  N=4                                                *
*                                                     *
*  8.000000  2.000000  3.000000  12.00000  25.00000   *
*  2.000000  4.000000  7.000000  0.250000  13.25000   *
*  3.000000  7.000000  3.000000  5.000000  18.00000   *
*  12.00000  0.250000  5.000000  2.000000  19.25000   *
*                                                     *
*  System solution:                                   *
*                                                     *
*  X1=  1.000000                                      *
*  X2=  1.000000                                      *
*  X3=  1.000000                                      *
*  X4=  1.000000                                      *
* --------------------------------------------------- * 
*                                                     *
******************************************************}
  Program Test_LU;
  Uses WinCrt,Basis,Lu;

  Var
  A : pVECT;      { matrix 0:n x 0:n stored in a vector }
  B : pVECT;      { vecteur 0:n }
  temp : pVECT;   { vecteur 0:n+1 }
  INDX : pIVECT;  { vecteur 0:n }

  { NOTA : zero index is not used here. }

  d, i, j, n, rc : integer;
  input, output, s : STRING;
  F1, F2 : TEXT;



  Begin {main program}

  Writeln;
  Write(' Data file name (without .dat): '); read(s);
  input := s + '.dat';
  output := s + '.lst';

  Assign(F1,input); Reset(F1);
  Assign(F2,output); Rewrite(F2);

  readln(F1,n);     { taille du système à résoudre }

  New(A); New(B); New(temp); New(INDX);

  WriteHead(F2,' LINEAR SYSTEM TO BE SOLVED:');
  writeln(F2,'  N=',n);
  writeln(F2);

  for i:=1 to n do
  begin
    ReadVec1(F1,n+1,n+1,temp);  {read a line}
    for j:=1 to n do
      A^[i*(n+1)+j] := temp^[j];
    B^[i] := temp^[n+1];
    WriteVec1(F2,n+1,n+1,temp)  {write a line}
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
      write(F2,'  X',i,'='); f_aff_reel(F2,B^[i]);
      writeln(F2)
    end   
  end;
  WriteEnd(F2);
  Close(F2);
  Writeln;
  Writeln(' Results in file ',output,'.');
  Writeln;
  Readkey;
  DoneWinCrt

  END.

{End of file test_lu.pas}