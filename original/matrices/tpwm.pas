{**************************************************************
*  Calculate the greatest eigenvalue of a real square matrix  *
*     and the associated eigenvector by the power method      *
* ----------------------------------------------------------- *
* Ref.: "Algèbre - Algorithmes et programmes en Pascal        *
*        By Jean-Louis Jardrin, Dunod Editor, Paris, 1988"    *
*        [BIBLI 10].                                          *
* ----------------------------------------------------------- *
* SAMPLE RUN:                                                 *
*                                                             *
* ----------------------------------------------              *
*  Calculate the greatest eigenvalue of a real                *
*  square matrix and the associated eigenvector               *
*  by the power method.                                       *
* ----------------------------------------------              *
*                                                             *
*  Size of matrix (maximum 25): 4                             *
*                                                             *
*  Line 1                                                     *
*    Element 1: 1                                             *
*    Element 2: 3                                             *
*    Element 3: 0                                             *
*    Element 4: 0                                             *
*                                                             *
*  Line 2                                                     *
*    Element 1: 4                                             *
*    Element 2: 2                                             *
*    Element 3: 0                                             *
*    Element 4: 0                                             *
*                                                             *
*  Line 3                                                     *
*    Element 1: 1                                             *
*    Element 2: -1                                            *
*    Element 3: 5                                             *
*    Element 4: -3                                            *
*                                                             *
*  Line 4                                                     *
*    Element 1: 2                                             *
*    Element 2: 0                                             *
*    Element 3: 4                                             *
*    Element 4: -2                                            *
*                                                             *
*  Precision: 1e-10                                           *
*  Epsilon  : 1e-10                                           *
*  Maximum number of iterations: 27                           *
*                                                             *
*    Eigenvalue: 5.00000000008673E+0000                       *
*                                                             *
*    Eigenvector:                                             *
*    7.49999999991327E-0001                                   *
*    1.00000000000000E+0000                                   *
*   -5.20833333277822E-0001                                   *
*   -8.33333332799904E-0002                                   *            
*                                                             *
*                  English TPW Version By J-P Moreau, Paris.  *
*                             (www.jpmoreau.fr)               *
***************************************************************
  Exact values are: gamma = 5
                    eigenvector = (1/48)(36,48,-25,-4)
--------------------------------------------------------------}                                  
Program Test_PWM;
Uses WinCrt;

CONST
     NMAX = 25;

TYPE
     MAT = Array[1..NMAX,1..NMAX] of Double;
     VEC = Array[1..NMAX] of Double;

VAR
     i,it,m,n: Integer;
     dta,eps,gamma: Double;
     A: MAT;
     X: VEC;


{***********************************************************
* calculate greatest eigenvalue and associated eigenvector *
* by the power method                                      *
* -------------------------------------------------------- *
* INPUTS:                                                  *
*         eps   : smallest number in double precision      *
*         dta   : required precision                       *
*         m     : maximum number of iterations             *
*         n     : size of real square matrix A(n,n)        *
*         A     : real square matrix A(n,n)                *
* OUTPUTS:                                                 *
*         it    : error indicator: -1=no convergence,      *
*                 0=method cannot be applied,              *
*                 1=convergence ok.                        *
*         gamma : greatest eigenvalue (in absolute value)  *
*                 of input matrix A(n,n)                   *
*         X1    : associated eigenvector                   *
***********************************************************}
Procedure PWM(eps,dta:Double; m,n:Integer; VAR A:MAT; VAR it:Integer;
               VAR gamma:double; VAR X1:VEC);
Var  i,j,l: Integer;
     phi,s: Double;
     X0: VEC;
Begin
  for i:=1 to n do X0[i]:=1.0/SQRT(I);
  it:=-1; l:=1;
  While (it=-1) and (l<=m) do
  begin
    gamma:=0.0;
    for i:=1 to n do
    begin
      X1[i]:=0.0;
      for j:=1 to n do X1[i]:=X1[i]+A[i,j]*X0[j];
      if ABS(X1[i])>ABS(gamma) then gamma:=X1[i]
    end;
    if ABS(gamma) < eps then it:=0
    else
    begin
      for i:=1 to n do X1[i]:=X1[i]/gamma;
      phi:=0.0;
      for i:=1 to n do
      begin
        s:=ABS(X1[i]-X0[i]);
        if s>phi then phi:=s
      end;
      if phi<dta then it:=1
      else
      begin
        X0:=X1;
        Inc(l)
      end
    end
  end
End; {of PWM}                                                                                                     

Procedure Read_data;
Var i,j: Integer;
Begin
  writeln;
  write(' Size of matrix (maximum ',NMAX,'): '); readln(n);
  for i:=1 to n do
  begin
    writeln;
    writeln(' Line ',i);
    for j:=1 to n do
    begin
      write('  Element ',j,': '); readln(A[i,j])
    end
  end;
  writeln;
  write(' Precision: '); readln(dta);
  write(' Epsilon  : '); readln(eps);
  write(' Maximum number of iterations: '); readln(m)
End;

{main program}
BEGIN
  writeln(' ----------------------------------------------');
  writeln('  Calculate the greatest eigenvalue of a real  ');
  writeln('  square matrix and the associated eigenvector ');
  writeln('  by the power method.                         ');
  writeln(' ----------------------------------------------');

  Read_Data;

  PWM(eps,dta,m,n,A,it,gamma,X);

  Case it+1 of
    0: writeln('  No convergence !');
    1: writeln('  Method does not apply.');
    2: begin
         writeln;
         writeln('  Eigenvalue: ', gamma);
         writeln;
         writeln('  Eigenvector:');
         for i:=1 to n do writeln('  ',X[i])
       end
  End;
  writeln;
  Readkey; DoneWinCrt
END.

{end of file tpwm.pas}