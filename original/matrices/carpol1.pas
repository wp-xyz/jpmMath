{********************************************************
* Calculate the coefficients of the characteristic poly-*
* nomial P(l)=A-l*I of a square complex matrix A(i,j)   *
* ----------------------------------------------------- *
* Ref.: "Algèbre, Algorithmes et programmes en Pascal   *
*        By Jean-Louis Jardrin, DUNOD Paris, 1988"      *
*        [BIBLI 10].                                    *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
* (Find the characteristic polynomial of matrix:        *
*           1 -1  3                                     *
*      A = -1  i -1-2i                                  *
*           i  1 -2+i   )                               *
*                                                       *
* Input matrix size (max 20): 3                         *
*                                                       *
* Line 1                                                *
*  Column 1                                             *
*   Real  part: 1                                       *
*   Imag. part: 0                                       *
*  Column 2                                             *
*   Real  part: -1                                      *
*   Imag. part: 0                                       *
*  Column 3                                             *
*   Real  part: 3                                       *
*   Imag. part: 0                                       *
* Line 2                                                *
*  Column 1                                             *
*   Real  part: -1                                      *
*   Imag. part: 0                                       *
*  Column 2                                             *
*   Real  part: 0                                       *
*   Imag. part: 1                                       *
*  Column 3                                             *
*   Real  part: -1                                      *
*   Imag. part: -2                                      *
* Line 3                                                *
*  Column 1                                             *
*   Real  part: 0                                       *
*   Imag. part: 1                                       *
*  Column 2                                             *
*   Real  part: 1                                       *
*   Imag. part: 0                                       *
*  Column 3                                             *
*   Real  part: -2                                      *
*   Imag. part: 1                                       *
*                                                       *
* Characteristic polynomial P(l):                       *
* Coefficient for degree 3:                             *
* -1.0000000000E+000 + 0.0000000000E+000 i              *
* Coefficient for degree 2:                             *
* -1.0000000000E+000 + 2.0000000000E+000 i              *
* Coefficient for degree 1:                             *
*  3.0000000000E+000 + 1.0000000000E+000 i              *
* Coefficient for degree 0:                             *
*  0.0000000000E+000 + 0.0000000000E+000 i              *
*                                                       *
* (The characteristic polynomial of matrix A is:        *
*   P(l)=-l^3+(-1+2i)l^2+(3+i)l                         *
*                                                       *
*               English version by J-P Moreau, Paris    *
*               (with dynamically allocated vectors).   *
*                         (www.jpmoreau.fr)             *
********************************************************}
Program CARPOL1;
Uses WinCrt, UComplex;

CONST
      NMAX  = 50;
      NMAXP = 51;

TYPE
      pVC  = ^VECC;
      VECC = Array[1..NMAXP] of Complex;
      pMC  = ^MATC;
      MATC = Array[1..NMAX,1..NMAX] of Complex;

VAR
      k,n: Integer;
      A: pMC;
      P: pVC; 


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
      writeln('  Column ',j);
      write('   real  part: '); readln(A^[i,j].r);
      write('   imag. part: '); readln(A^[i,j].i)
    end
  end
End;

Procedure TRC(n:integer; C:pMC; VAR t0:Complex);
Var i:integer; c0:Complex;
Begin
  t0.r:=0.0; t0.i:=0.0;
  for i:=1 to n do
  begin
    c0:=t0;
    CADD(c0,C^[i,i],t0)
  end
End;

{*****************************************************
* The procedure PCCS calculates the complex coeffi-  *
* cients of the characteristic polynomial P(l)=A-l*I *
* of a real square tridiagonal matrix A(i,j).        *
*                                                    *
* Note: the roots of P(l) are the eigenvalues of     *
*       matrix A(i,j).                               *
* -------------------------------------------------- *
* INPUTS:                                            *
*          n: zize of matrix (integer)               *
*          A: complex square matrix (n,n)            *
*             (see main program).                    *
* OUTPUT:                                            *
*          P: vector of complex coefficients of P(l) *
*                                                    *
*  The types pVC (pointer to a complex vector) and   *
*  pMC (pointer to a complex square matrix) must be  *
*  declared and allocated in the calling program.    *
*****************************************************} 
Procedure PCCS(n:integer; A:pMC; VAR P:pVC);
Var i,j,k,l: integer;
    c0,c1,s,t0,z0,z1: Complex;
    B,C: pMC;
Begin
  New(B); New(C);
  z0.r:=0.0; z0.i:=0.0;
  z1.r:=1.0; z1.i:=0.0;
  if Odd(n) then CPRO(-1.0,z1,P^[1])
            else P^[1]:=z1;
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
          s.r:=z0.r; s.i:=z0.i;
          for k:=1 to n do
          begin
            c0.r:=s.r;  c0.i:=s.i;
            CMUL(B^[i,k],A^[k,j],c1);
            CADD(c0,c1,s)
          end;
          C^[i,j].r:=s.r;
          C^[i,j].i:=s.i
        end
      end
    end;
    TRC(n,c,c0);
    CPRO(1.0/l,c0,t0);
    CMUL(t0,P^[1],c0);
    CPRO(-1.0,c0,P^[l+1]);
    if l<n then
      for i:=1 to n do
        for j:=1 to n do
          if j=i then CDIF(C^[i,j],t0,B^[i,j])
                 else B^[i,j]:=C^[i,j]
  end;
  Dispose(B); Dispose(C)
End;

{main program}
BEGIN

  New(A); New(P);
  Read_data;
  {call routine to calculate coefficients}
  PCCS(n,A,P);
  {print results}
  writeln;
  writeln(' Characteristic polynomial P(l):');
  For k:=1 to n+1 do
  begin
    writeln(' Coefficient for degree ',n-k+1,': ');
    write(P^[k].r);
    if P^[K].i>=0 then write(' +')
                  else write(' -');
    writeln(ABS(P^[k].i),' i')               
  end;
  {exit section}
  Readkey;
  Dispose(A); Dispose(P);
  DoneWinCrt

END.

{end of file carpol1.pas}

