{********************************************************
*   Calculate the coefficients of the characteristic    *
*   polynomial P(l)=A-l*I of a square symmetric real    *
*         matrix A(i,j) by Lanczos's method             *
* ----------------------------------------------------- *
* Ref.: "Algèbre, Algorithmes et programmes en Pascal   *
*        By Jean-Louis Jardrin, DUNOD Paris, 1988"      *
*        [BIBLI 10].                                    *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
* (Find the caracteristic polynomial of matrix:         *
*          10  7   8   7                                *
*      A =  7  5   6   5                                *
*           8  6  10   9                                *
*           7  5   9  10                                *
*                                                       *
* Input matrix size (max 50): 4                         *
*  Line 1:                                              *
*   Element 1: 10                                       *
*   Element 2: 7                                        *
*   Element 3: 8                                        *
*   Element 4: 7                                        *
*  Line 2:                                              *
*   Element 2: 5                                        *
*   Element 3: 6                                        *
*   Element 4: 5                                        *
*  Line 3:                                              *
*   Element 3: 10                                       *
*   Element 4: 9                                        *
*  Line 4:                                              *
*   Element 4: 10                                       *
*                                                       *
* Chracteristic polynomial P(l):                        *
* Coefficient for degree 4:  1.0000000000E+000          *
* Coefficient for degree 3: -3.5000000000E+001          *
* Coefficient for degree 2:  1.4600000000E+002          *
* Coefficient for degree 1: -9.9999999999E+001          *
* Coefficient for degree 0:  1.0000000000E+000          *
*                                                       *
* (The characteristic polynomial of matrix A is:        *
*   P(l)=l^4-35*l^3+146*l^2-100*l+1 )                   *
*                                                       *
*               English version by J-P Moreau, Paris    *
*  (with dynamically allocated vectors and matrices).   *
*                        (www.jpmoreau.fr)              *
********************************************************}
Program CARPOL3;
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
      it,k,n: Integer;
      eps: Double;
      A, Tau: pM;
      D,E,P : pV; 

{read symmetric matrix A}
Procedure Read_data;
Var i,j:Integer;
Begin
  writeln;
  write(' Input matrix size (max ',NMAX,'): '); readln(n);
  for i:=1 to n do
  begin
    writeln; writeln(' Line ',i);
    for j:=i to n do
    begin
      write('  Element ',j,': '); readln(A^[i,j]);
      if j<>i then A^[j,i]:=A^[i,j]
    end
  end
End;

{********************************************************
* The TDSL procedure tridiagonalizes the symmetric real *
* matrix A(n,n) by Lanczos's method. Only the main dia- *
* gonal terms and the supradiagonal terms are calculated*
* ----------------------------------------------------- *
* Inputs:                                               *
*          eps  precision (double)                      *
*          n    size of square A matrix (integer)       *
*          A    square real symmetric matrix (n,n)      *
* Outputs:                                              *
*          it   flag=0 if method failed, else =1        *
*          D    vector of size n+1 storing the n        *
*               elements of the main diagonal.          *
*          E    vector of size n+1 storing the n-1      *
*               elements of the supradiagonal.          *
*          Tau  Utility matrix used by the method.      *
*                                                       *
*   The types pV (pointer to a real vector) and pM      *
*   (pointer to a real matrix) must be declared and     *
*   allocated in the calling program.                   *
********************************************************}
Procedure TDSL(eps:double;n:integer;A:pM;VAR it:integer;
                 VAR D,E:pV; VAR Tau:pM);
Var i,j,k:integer; s,u0,v0:double; W:pV;
Begin
  New(W);   {initialize local vector}
  Tau^[1,1]:=1.0; for i:=2 to n do Tau^[i,1]:=0.0;
  u0:=1.0; it:=1; k:=1;
  Repeat
    for i:=1 to n do
    begin
      s:=0.0;
      for j:=1 to n do s:=s+A^[i,j]*Tau^[j,k];
      W^[i]:=s
    end;
    s:=0.0; for i:=1 to n do s:=s+Tau^[i,k]*W^[i];
    D^[k]:=s/u0;
    if k<>n then
    begin
      for i:=1 to n do
        if k=1 then
          Tau^[i,k+1]:=W^[i]-D^[k]*Tau^[i,k]
        else
          Tau^[i,k+1]:=W^[i]-D^[k]*Tau^[i,k]-E^[k-1]*Tau^[i,k-1];
      v0:=0.0; for i:=1 to n do v0:=v0+SQR(Tau^[i,k+1]);
      if v0<eps then
        it:=0
      else
      begin
        E^[k]:=v0/u0; u0:=v0
      end
    end; 
    k:=k+1
  Until (k>n) or (it=0);
  Dispose(W)
End;
                                    

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
* ------------------------------------------------- *
* (For demo of procedure PCTR alone, see program    *
*  carpol.pas).                                     *
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


{main program}
BEGIN

  New(A); New(Tau); New(D); New(E); New(P);
  Read_data;
  eps:=1e-10;
  {tridiagonalize symmetric A matrix by Lanczos's method}
  TDSL(eps,n,A,it,D,E,TAU);
  {if success, use method for a tridiagonal matrix}
  if it=1 then PCTR(n,E,D,P);
  {print results}
  writeln;
  Case it of
  0: begin
       writeln(' Method failed.')
     end;
  1: begin
       writeln(' Characteristic polynomial P(l):');
       for k:=1 to n+1 do
       begin
         write(' Coefficient for degree ',n-k+1,': ');
         writeln(P^[k])
       end
     end
  end;
  {exit section}
  Readkey;
  Dispose(A); Dispose(Tau); Dispose(D);
  Dispose(E); Dispose(P);
  DoneWinCrt

END.

{end of file carpol3.pas}

