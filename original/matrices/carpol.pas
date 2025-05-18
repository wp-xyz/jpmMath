{********************************************************
* Calculate the coefficients of the characteristic poly-*
* nomial P(l)=A-l*I of a real tridiagonal matrix A(i,j) *
* ----------------------------------------------------- *
* Ref.: "Algèbre, Algorithmes et programmes en Pascal   *
*        By Jean-Louis Jardrin, DUNOD Paris, 1988"      *
*        [BIBLI 10].                                    *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
* (Find the caracteristic polynomial of matrix:         *
*           4 1 0 0 0 0                                 *
*           2 4 1 0 0 0                                 *
*      A =  0 2 4 1 0 0                                 *
*           0 0 2 4 1 0                                 *
*           0 0 0 2 4 1                                 *
*           0 0 0 0 2 4   )                             *
*                                                       *
* Input matrix size (max 100): 6                        *
* Input 5 elements of subdiagonal:                      *
*   Element 1: 2                                        *
*   Element 2: 2                                        *
*   Element 3: 2                                        *
*   Element 4: 2                                        *
*   Element 5: 2                                        *
* Input 6 elements of main diagonal:                    *
*   Element 1: 4                                        *
*   Element 2: 4                                        *
*   Element 3: 4                                        *
*   Element 4: 4                                        *
*   Element 5: 4                                        *
*   Element 6: 4                                        *
* Input 5 elements of supradiagonal:                    *
*   Element 1: 1                                        *
*   Element 2: 1                                        *
*   Element 3: 1                                        *
*   Element 4: 1                                        *
*   Element 5: 1                                        *
*                                                       *
* Chracteristic polynomial P(l):                        *
* Coefficient for degree 6:  1.0000000000E+000          *
* Coefficient for degree 5: -2.4000000000E+001          *
* Coefficient for degree 4:  2.3000000000E+002          *
* Coefficient for degree 3: -1.1200000000E+003          *
* Coefficient for degree 2:  2.9040000000E+003          *
* Coefficient for degree 1: -3.7760000000E+003          *
* Coefficient for degree 0:  1.9120000000E+003          *
*                                                       *
* (The characteristic polynomial of matrix A is:        *
*   P(l)=l^6-24*l^5+230*l^4-1120*l^3+2904*l^2-3376*l    *
*        +1912   )                                      *
*                                                       *
*               English version by J-P Moreau, Paris    *
*               (with dynamically allocated vectors).   *
*                         (www.jpmoreau.fr)             *
********************************************************} 
Program CARPOL;
Uses WinCrt;

CONST
      NMAX = 100;
      NMAXP = 101;

TYPE
      pV  = ^VEC;
      VEC = Array[1..NMAXP] of Double;

VAR
      k,n: Integer;
      B,C,D,E,P: pV; 


{****************************************************
* The procedure PCTR calculates the coefficients of *
* the characteristic polynomial P(l)=A-l*I of a real*
* square tridiagonal matrix A(i,j).                 *
*                                                   *
* Note: the roots of P(l) are the eigenvalues of    *
*       matrix A(i,j).                              *
* ------------------------------------------------- *
* INPUTS:                                           *
*          n: zize of matrix (integer)              *
*          B: vector of codiagonal terms products   *
*             (see main program).                   *
*          D: vector of main diagonal terms         *
* OUTPUT:                                           *
*          P: vector of coefficients of P(l)        *
*                                                   *
* The type pV (pointer to a real vector) must be    *
* declared and allocated in the calling program.    *
****************************************************} 
Procedure PCTR(n:integer; B,D:pV; VAR P:pV);
Var i,k:integer; T:pV;
Begin
  New(T);
  P^[2]:=D^[1]; T^[1]:=1.0; P^[1]:=-T^[1];
  For k:=2 to n do
  begin
    P^[k+1]:=D^[k]*P^[k]-B^[k-1]*T^[k-1];
    for i:=k downto 3 do
    begin
      T^[i]:=P^[i];
      P^[i]:=-T^[i]+D^[k]*P^[i-1]-B^[k-1]*T^[i-2]
    end;
    T^[2]:=P^[2]; P^[2]:=-T^[2]+D^[k]*P^[1];
    T^[1]:=P^[1]; P^[1]:=-T^[1]
  end;
  Dispose(T)
End;

Procedure Read_data;
Var i:Integer;
Begin
  writeln;
  write(' Input matrix size (max ',NMAX,'): '); readln(n);
  writeln(' Input ',n-1,' elements of subdiagonal:');
  for i:=1 to n-1 do
  begin
    write('   Element ',i,': '); readln(C^[i])
  end;
  writeln(' Input ',n,' elements of main diagonal:');
  for i:=1 to n do
  begin
    write('   Element ',i,': '); readln(D^[i])
  end;
  writeln(' Input ',n-1,' elements of supradiagonal:');
  for i:=1 to n-1 do
  begin
    write('   Element ',i,': '); readln(E^[i])
  end
End;


{main program}
BEGIN

  New(B); New(C); New(D); New(E); New(P);
  Read_data;
  {prepare B vector}
  for k:=1 to n-1 do B^[k]:=C^[k]*E^[k];
  Dispose(C); Dispose(E);
  {call routine to calculate coefficients}
  PCTR(n,B,D,P);
  {print results}
  writeln;
  writeln(' Characteristic polynomial P(l):');
  For k:=1 to n+1 do
    writeln(' Coefficient for degree ',n-k+1,': ',P^[k]);
  {exit section}
  Readkey;
  Dispose(B); Dispose(D); Dispose(P);
  DoneWinCrt

END.

{end of file carpol.pas}