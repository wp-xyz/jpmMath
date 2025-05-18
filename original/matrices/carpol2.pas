{********************************************************
* Calculate the coefficients of the characteristic poly-*
* nomial P(l)=A-l*I of a real matrix A(i,j)             *
* ----------------------------------------------------- *
* Ref.: "Algèbre, Algorithmes et programmes en Pascal   *
*        By Jean-Louis Jardrin, DUNOD Paris, 1988"      *
*        [BIBLI 10].                                    *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
* (Find the caracteristic polynomial of matrix:         *
*           1  0  0  0  0  1                            *
*           1  1  0  0  0 -1                            *
*      A = -1  1  1  0  0  1                            *
*           1 -1  1  1  0 -1                            *
*          -1  1 -1  1  1  1                            *
*           1 -1  1 -1  1 -1  )                         *
*                                                       *
* Input matrix size (max 100): 6                        *
*  Line 1:                                              *
*   Element 1: 1                                        *
*   Element 2: 0                                        *
*   Element 3: 0                                        *
*   Element 4: 0                                        *
*   Element 5: 0                                        *
*   Element 6: 1                                        *
*  Line 2:                                              *
*   Element 1: 1                                        *
*   Element 2: 1                                        *
*   Element 3: 0                                        *
*   Element 4: 0                                        *
*   Element 5: 0                                        *
*   Element 6: -1                                       *
*  Line 3:                                              *
*   Element 1: -1                                       *
*   Element 2: 1                                        *
*   Element 3: 1                                        *
*   Element 4: 0                                        *
*   Element 5: 0                                        *
*   Element 6: 1                                        *
*  Line 4:                                              *
*   Element 1: 1                                        *
*   Element 2: -1                                       *
*   Element 3: 1                                        *
*   Element 4: 1                                        *
*   Element 5: 0                                        *
*   Element 6: -1                                       *
*  Line 5:                                              *
*   Element 1: -1                                       *
*   Element 2: 1                                        *
*   Element 3: -1                                       *
*   Element 4: 1                                        *
*   Element 5: 1                                        *
*   Element 6: 1                                        *
*  Line 6:                                              *
*   Element 1: 1                                        *
*   Element 2: -1                                       *
*   Element 3: 1                                        *
*   Element 4: -1                                       *
*   Element 5: 1                                        *
*   Element 6: -1                                       *
*                                                       *
* Chracteristic polynomial P(l):                        *
* Coefficient for degree 6:  1.0000000000E+000          *
* Coefficient for degree 5: -4.0000000000E+000          *
* Coefficient for degree 4:  0.0000000000E+000          *
* Coefficient for degree 3:  3.0000000000E+001          *
* Coefficient for degree 2: -7.5000000000E+001          *
* Coefficient for degree 1:  7.9000000000E+001          *
* Coefficient for degree 0: -3.2000000000E+001          *
*                                                       *
* (The characteristic polynomial of matrix A is:        *
*   P(l)=l^6-4*l^5+30*l^3-75*l^2+79*l-32 )              *
*                                                       *
*               English version by J-P Moreau, Paris    *
*               (with dynamically allocated vectors).   *
*                         (www.jpmoreau.fr)             *
********************************************************}
Program CARPOL2;
Uses WinCrt;

CONST
      NMAX  = 50;
      NMAXP = 51;

TYPE
      pV  = ^VEC;
      VEC = Array[1..NMAXP] of Double;
      pM  = ^MAT;
      MAT = Array[1..NMAX,1..NMAX] of Double;

VAR
      k,n: Integer;
      A: pM;
      P: pV; 


Procedure Read_data;
Var i,j:Integer;
Begin
  writeln;
  write(' Input matrix size (max ',NMAX,'): '); readln(n);
  for i:=1 to n do
  begin
    writeln; writeln(' Line ',i);
    for j:=1 to n do
    begin
      write('  Element ',j,': '); readln(A^[i,j])
    end
  end
End;

Function TRM(n:integer; C:pM):Double;
Var i:integer; t0:Double;
Begin
  t0:=0.0;
  for i:=1 to n do
    t0:=t0+C^[i,i];
  TRM:=t0
End;

{****************************************************
* The procedure PCMS calculates the coefficients of *
* the characteristic polynomial P(l)=A-l*I of a real*
* square matrix A(i,j).                             *
*                                                   *
* Note: the roots of P(l) are the eigenvalues of    *
*       matrix A(i,j).                              *
* ------------------------------------------------- *
* INPUTS:                                           *
*          n: zize of matrix (integer)              *
*          A: real matrix (n,n)                     *
*             (see main program).                   *
* OUTPUT:                                           *
*          P: vector of coefficients of P(l)        *
*                                                   *
* The types pV (pointer to a real vector) and pM    *
* (pointer to a real matrix) must be declared and   *
* allocated in the calling program.                 *
****************************************************}
Procedure PCMS(n:integer; A:pM; VAR P:pV);
Var i,j,k,l: integer;
    s,t0: Double;
    B,C: pM;
Begin
  New(B); New(C);
  if Odd(n) then P^[1]:=-1.0
            else P^[1]:= 1.0;
  for l:=1 to n do
  begin
    if l=1 then
    begin
      for i:=1 to n do
        for j:=1 to n do
          C^[i,j]:=A^[i,j]
    end
    else
    begin
      for i:=1 to n do
      begin
        for j:=1 to n do
        begin
          s:=0.0;
          for k:=1 to n do
            s:=s+B^[i,k]*A^[k,j];
          C^[i,j]:=s
        end
      end
    end;
    t0:=TRM(n,C)/l;
    P^[l+1]:=-t0*P^[1];
    if l<n then
      for i:=1 to n do
        for j:=1 to n do
          if j=i then B^[i,j]:=C^[i,j]-t0
                 else B^[i,j]:=C^[i,j]
  end;
  Dispose(B); Dispose(C)
End;

{main program}
BEGIN

  New(A); New(P);
  Read_data;
  {call routine to calculate coefficients}
  PCMS(n,A,P);
  {print results}
  writeln;
  writeln(' Characteristic polynomial P(l):');
  For k:=1 to n+1 do
  begin
    write(' Coefficient for degree ',n-k+1,': ');
    writeln(P^[k])
  end;
  {exit section}
  Readkey;
  Dispose(A); Dispose(P);
  DoneWinCrt

END.

{end of file carpol2.pas}

