{*******************************************************
*  Solving a symmetric linear system by Gauss method   *
* ---------------------------------------------------- *
* Reference:                                           *
*                                                      *
*   "ALGEBRE Algorithmes et programmes en Pascal       *
*    de Jean-Louis Jardrin - Dunod BO-PRE 1988"        *
*    [BIBLI 10].                                       *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
*                                                      *
* Input file (matsym.dat):                             *
*                                                      *
* 4                                                    *
* 8  2  3  12     25                                   *
*    4  7   0.25  13.25                                *
*       3   5     18                                   *
*           2     19.25                                *
*                                                      *
* Output file (matsym.lst):                            *
*                                                      *
* Size of symmetric system: 4                          *
*                                                      *
* SYSTEM TO BE SOLVED:                                 *
*                                                      *
*  8.0000  2.0000  3.0000 12.0000 25.0000              *
*  2.0000  4.0000  7.0000  0.2500 13.2500              *
*  3.0000  7.0000  3.0000  5.0000 18.0000              *
* 12.0000  0.2500  5.0000  2.0000 19.2500              *
*                                                      *
* System solution:                                     *
*                                                      *
*   x1 =   1.000000                                    *
*   x2 =   1.000000                                    *
*   x3 =   1.000000                                    *
*   x4 =   1.000000                                    *
*                                                      *
*******************************************************}
Program SYSLIN;
uses WinCrt;

CONST NMAX  = 30;
      NMAXP = 31;
TYPE  MAT = ARRAY[1..NMAX,1..NMAXP] of real;
      VEC = ARRAY[1..NMAX] of real;
VAR   i,j,it,n : integer;
      eps    : real;
      A      : MAT;
      X      : VEC;
      F,F1   : TEXT;

{**********************************************************************
*   This procedure solves a symmetric linear system by Gauss method   *
*                     in simple precision                             *
* ------------------------------------------------------------------- *
* The input matrix includes the right hand column, only the upper     *
* triangle terms are given.                                           *
* INPUTS:                                                             *
*           eps:    desired precision (real)                          *
*             n:    size of system (integer)                          *
*             A:    extended matrix of system of type MAT             *
* OUTPUTS:                                                            *
*            it:    error indicator (0=singular system, 1=OK)         *
*             X:    system solution (x1,...,xn) of type VEC           *
**********************************************************************}
Procedure RSLSG(eps:real; n:integer; A:MAT; VAR it:integer; VAR X:VEC);
var i,j,k: integer; q0,s: real;
begin
  it:=1; k:=1;
  Repeat
    if abs(A[k,k]) < eps then it:=0
    else
    begin
      for i:=k+1 to n do
      begin
	q0:=A[k,i]/A[k,k];
	for j:=i to n+1 do A[i,j]:=A[i,j]-q0*A[k,j]
      end;
      Inc(k)
    end
  Until (it=0) OR (k=n);
  if (it=1) AND (abs(A[n,n])>=eps) then
  begin
    X[n]:=A[n,n+1]/A[n,n];
    for i:=n-1 downto 1 do
    begin
      s:=0.0;
      for j:=i+1 to n do s:=s+A[i,j]*X[j];
      X[i]:=(A[i,n+1]-s)/A[i,i]
    end
  end
  else it:=0;
end;


{main program}
BEGIN
  Clrscr;
  Assign(F,'matsym.dat'); Reset(F);
  Assign(F1,'matsym.lst'); Rewrite(F1);

  Read(F,n);

  writeln(F1);
  Writeln(F1,'  Size of symmetric system: ',n);

  writeln(F1);
  writeln(F1,'  SYSTEM TO BE SOLVED:');
  writeln(F1);

  for i:=1 to n do
    for j:=i to n+1 do
    begin
      Read(F,A[i,j]);
      A[j,i]:=A[i,j]
    end;

  Close(F);

  for i:=1 to n do
  begin
    for j:=1 to n+1 do Write(F1,A[i,j]:8:4);
    writeln(F1)
  end;

  eps:=1E-10;

  RSLSG(eps,n,A,it,X);

  case it of
    0: writeln(F1,'  *** ERROR ***');
    1: begin
	 writeln(F1);
	 writeln(F1,'  System solution:');
	 writeln(F1);
	 for i:=1 to n do
	   writeln(F1,'    x',i,' = ',X[i]:10:6)
       end
  end;
  close(F1);
  writeln;
  writeln(' Results in file matsym.lst.');
  Readkey;
  DoneWinCrt
END.

{End of file syslin.pas}