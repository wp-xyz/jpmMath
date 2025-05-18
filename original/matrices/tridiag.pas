{************************************************************************
* This program solves a tridiagonal linear set of equations M * U = R   *
* --------------------------------------------------------------------- *
* SAMPLE RUN:                                                           *
*                                                                       *
* Input data file (tridiag.dat):                                        *
*                                                                       *
*  5                                                                    *
*  0.0   1.0   -2.0                                                     *
* 10.0   5.0   -7.5                                                     *
*  2.22 -3.25   0.00456                                                 *
*  2.0   4.0   -3.3                                                     *
* -1.0   3.0    0.0                                                     *
*                                                                       *
* (The first line contains the size of system                           *
* Then the three columns give the three main diagonals of matrix M, the *
* first column begins with a zero and the third column ends with zero)  *
*                                                                       *
* Output file (tridiag.lst):                                            *
*                                                                       *
* --------------------------------------------------------------------  *
*   Tridiagonal system of equations                                     *
* --------------------------------------------------------------------  *
*  upper subdiagonal:                                                   *
* -2.000000 -7.500000  0.004560 -3.300000  0.000000                     *
*  main diagonal:                                                       *
*  1.000000  5.000000 -3.250000  4.000000  3.000000                     *
*  lower subdiagonal:                                                   *
*  0.000000  10.00000  2.220000  2.000000 -1.000000                     *
*  right side member:                                                   *
*  1.000000  1.000000  1.000000  1.000000  1.000000                     *
* --------------------------------------------------------------------  *
*   System solution:                                                    *
*  -0.136497 -0.568249 -0.694162  1.202870  0.734290                    *
* --------------------------------------------------------------------  *
*                                                                       *
* Reference: "Numerical Recipes by W.H. Press, B.P. Flannery, S.A. Teu- *
*             kolsky, W.T. Vetterling, Cambridge University Press, 1987 *
*             [BIBLI 08].                                               *
*                                                                       *
*                      Pascal Version from FORTRAN By J-P Moreau, Paris *
*                                     (www.jpmoreau.fr)                 *
************************************************************************}
  PROGRAM TRIDIAG;
  USES WinCrt, Basis_r;

  VAR
  a, b, c, r, u : pVECT;
  i,n,rc : INTEGER;
  fp1,fp2 : TEXT;
  nom : STRING;

  PROCEDURE ligne;
  Begin
    writeln(fp2)
  End;

  PROCEDURE TRIDAG(A:pVECT;B:pVECT;C:pVECT;R:pVECT;VAR U:pVECT;N:INTEGER;VAR CODE:INTEGER);
  { Solves for a vector U of length N the tridiagonal linear set
  ! M U = R, where A, B and C are the three main diagonals of matrix
  ! M(N,N), the other terms are 0. R is the right side vector. }
  VAR BET:REAL; J:INTEGER;
  GAM : pVECT;
  Begin
  New(GAM);

  IF B^[1]=0 THEN
  begin
    CODE:=1;
    exit
  end;

  BET:=B^[1];
  U^[1]:=R^[1]/BET;
  for J:=2 to N do                  {Decomposition and forward substitution}
  begin
    GAM^[J]:=C^[J-1]/BET;
    BET:=B^[J]-A^[J]*GAM^[J];
    IF BET=0 THEN                   {Algorithm fails}
    begin
      CODE:=2;
      exit
    end;
    U^[J]:=(R^[J]-A^[J]*U^[J-1])/BET
  end;

  for J:=N-1 downto 1 do            {Back substitution}
    U^[J]:=U^[J]-GAM^[J+1]*U^[J+1];
    
  CODE:=0
  End;

  BEGIN {main}
  writeln;
  write(' Input data file name (without .dat): '); readln(nom);

  Assign(fp1,nom+'.dat'); Reset(fp1);
  Assign(fp2,nom+'.lst'); Rewrite(fp2);
  
  read(fp1,n);

  New(a); new(b); new(c); new(r); new(u);

  for i:=1 to n do
  begin
    read(fp1,a^[i],b^[i],c^[i]);
    r^[i]:=1.0
  end; 

  close(fp1);

  writehead(fp2,'  Tridiagonal system of equations');
  writeln(fp2,'  upper subdiagonal:');
  WriteVec1(fp2,n,n,c);
  writeln(fp2,'  main diagonal:');
  WriteVec1(fp2,n,n,b);
  writeln(fp2,'  lower subdiagonal:');
  WriteVec1(fp2,n,n,a);
  writeln(fp2,'  right side member:');
  WriteVec1(fp2,n,n,r);
  writeend(fp2);

  TRIDAG(a,b,c,r,u,n,rc);

  if rc<>0 then
    writeln(fp2,' ERROR IN TRIDAG : RETURN CODE=',rc)
  else
  begin
    writeln(fp2,'  System solution:');
    WriteVec1(fp2,n,n,U);
  end;
  writeend(fp2);

  writeln;
  writeln(' Results in ',nom,'.lst.');
  close(fp2);
  Readkey;
  DoneWinCrt

  END.

{End of file tridiag.pas}