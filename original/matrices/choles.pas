{*****************************************************************
*  Inversion of a symmetric matrix by Cholesky decomposition.    *
*  The matrix must be positive definite.                         * 
* -------------------------------------------------------------- *
* REFERENCE:                                                     *
*             From a Java Library Created by Vadim Kutsyy,       *
*             "http://www.kutsyy.com".                           *
* -------------------------------------------------------------- * 
* SAMPLE RUN:                                                    *
*                                                                *
* Inversion of a square real symetric matrix by Cholevsky method *
* (The matrix must positive definite).                           *
*                                                                *
* Size = 4                                                       *
*                                                                *
* Matrix A:                                                      *
*   5.000000  -1.000000  -1.000000  -1.000000                    *
*  -1.000000   5.000000  -1.000000  -1.000000                    *
*  -1.000000  -1.000000   5.000000  -1.000000                    *
*  -1.000000  -1.000000  -1.000000   5.000000                    *
*                                                                *
* Determinant = 432.000000                                       *
*                                                                *
* Matrix Inv(A):                                                 *
*   0.250000   0.083333   0.083333   0.083333                    *
*   0.083333   0.250000   0.083333   0.083333                    *
*   0.083333   0.083333   0.250000   0.083333                    *
*   0.083333   0.083333   0.083333   0.250000                    *
*                                                                *
* Do you want a verifivation (y/n)?: y                           *
*                                                                *
* Product Inv(A) * A:                                            *
*   1.000000   0.000000   0.000000   0.000000                    *
*   0.000000   1.000000   0.000000   0.000000                    *
*   0.000000   0.000000   1.000000   0.000000                    *
*   0.000000   0.000000   0.000000   1.000000                    *
*                                                                *
*                  TPW Release 1.1 By Jean-Pierre Moreau, Paris. *
*                                (www.jpmoreau.fr)               *
* -------------------------------------------------------------- *
* Release 1.1 : added verification Inv(A) * A = I.               *
*****************************************************************}
Program Cholesky;
Uses WinCrt;

Const
      SIZE=25;

Type
      pMAT = ^MAT;
      MAT = Array[0..SIZE,0..SIZE] of Double;
      pVEC = ^VEC;
      VEC = Array[0..SIZE] of Double;

Var
      A : pMAT;  {input matrix of size n,n}
      A1: pMAT;  {copy of A for verification}
      B : pMAT;  {Inverse of A of size n,n}
      C : pMAT;  {to store product A1*B=I} 

      i,j,n: Integer;


      Procedure choldc1(n:integer; Var a:pMAT; Var p:pVEC); Forward;


     {  -----------------------------------------------
            Cholesky decomposition.

            input    n  size of matrix
            input    A  Symmetric positive def. matrix
            output  aa  lower deomposed matrix
            uses        choldc1(int,pMAT,pVEC)
        -----------------------------------------------  }
        Procedure choldc(n:integer; A:pMAT; Var aa:pMAT);
	Var i,j:integer; p:pVEC;
        Begin
          New(p);
          for i := 0 to n-1 do 
	    for j := 0 to n-1 do 
	      aa^[i,j] := A^[i,j];

	  choldc1(n, aa, p);

          for i := 0 to n-1 do
          begin
            aa^[i,i] := p^[i];
            for j := i + 1 to n-1 do aa^[i,j] := 0.0
          end;
          Dispose(p)
        End;

      { -----------------------------------------------------
             Inverse of Cholesky decomposition.

             input    n  size of matrix
             input    A  Symmetric positive def. matrix
             output  aa  inverse of lower decomposed matrix
             uses        choldc1(int,pMAT,pVEC)         
        ----------------------------------------------------- }
        Procedure choldcsl(n:integer; A:pMAT; Var aa:pMAT);
	Var i,j,k:integer; sum:double; p:pVEC;
        Begin
          New(p);
          for i := 0 to n-1 do
	    for j := 0 to n-1 do 
	      aa^[i,j] := A^[i,j];

	  choldc1(n, aa, p);

          for i := 0 to n-1 do
          begin
            aa^[i,i] := 1.0 / p^[i];
            for j := i + 1 to n-1 do
            begin
              sum := 0.0;
              for k := i to j-1 do sum := sum - aa^[j,k] * aa^[k,i];
              aa^[j,i] := sum / p^[j];
	    end
	  end;
          Dispose(p)
	End;
 
{  -----------------------------------------------------------------------------
         Computation of Determinant of the matrix using Cholesky decomposition

         input    n  size of matrix
         input    a  Symmetric positive def. matrix
         return      det(a)
         uses        choldc(int,pMAT,pMAT)
   ------------------------------------------------------------------------------ }
        Function choldet(n:integer; a:pMAT): Double;
	Var c:pMAT; d:double; i,j:integer;
        Begin
          New(c);
          d:=1.0;
          choldc(n,a,c);
          for i := 0 to n-1 do  d := d * c^[i,i];
          choldet := d * d;
          Dispose(c)
	End;

 
{  ---------------------------------------------------
         Matrix inverse using Cholesky decomposition

         input    n  size of matrix
         input	  A  Symmetric positive def. matrix
         output  aa  inverse of A
         uses        choldc1(int,pMAT,pVEC)
   --------------------------------------------------- }
        Procedure cholsl(n:integer; A:pMAT; Var aa:pMAT);
	Var i,j,k:integer;
        Begin

          choldcsl(n,A,aa);

          for i := 0 to n-1 do
            for j := i + 1 to n-1 do aa^[i,j] := 0.0;

          for i := 0 to n-1 do
          begin
            aa^[i,i] := aa^[i,i] * aa^[i,i];
            for k := i + 1 to n-1 do
              aa^[i,i] := aa^[i,i] + aa^[k,i] * aa^[k,i];

            for j := i + 1 to n-1 do
              for k := j to n-1 do
                aa^[i,j] := aa^[i,j] + aa^[k,i] * aa^[k,j];
          end;
          for i := 0 to  n-1 do
            for j := 0 to i-1 do
              aa^[i,j] := aa^[j,i]
	End;

{  ------------------------------------------------------
         main method for Cholesky decomposition.

         input         n  size of matrix
         input/output  a  Symmetric positive def. matrix
         output        p  vector of resulting diag of a
         author:       <Vadum Kutsyy, kutsyy@hotmail.com>
   ------------------------------------------------------ }
        Procedure choldc1(n:integer; Var a:pMAT; Var p:pVEC);
        Var i,j,k:integer; sum:double;
        Begin
	  for i := 0 to n-1 do
          begin
            for j := i to n-1 do
            begin
              sum := a^[i,j];
              for k := i - 1 Downto 0 do
                sum := sum - a^[i,k] * a^[j,k];
              if i = j then
              begin
                if sum <= 0 then
                  writeln(' the matrix is not positive definite!');
                  p^[i] := sqrt(sum);
	      end
              else
                a^[j,i] := sum / p^[i]
	    end
          end
	End;

       {print a square real matrix A of size n with caption s
	(n items per line).                                  }
	Procedure MatPrint(s:String;n:integer; A:pMAT);
        Var i,j:integer;
        Begin
          writeln;
          writeln(' ',s);
	  for i:=0 to n-1 do
          begin
	    for j:=0 to n-1 do write(' ',A^[i,j]:10:6);
            writeln
          end
        End;

        Function Check_Matrix(n:integer;A:pMAT):Boolean;
        var i,j,k:integer; sum:double;
        begin
          Check_Matrix:=True;
	  for i := 0 to n-1 do
          begin
            for j := i to n-1 do
            begin
              sum := a^[i,j];
              for k := i - 1 Downto 0 do
                sum := sum - a^[i,k] * a^[j,k];
              if i = j then
                if sum <= 0 then Check_Matrix:=False
	    end
          end
        end;

{******************************************
*    MULTIPLICATION OF TWO SQUARE REAL    *                                     
*    MATRICES                             *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*            B  MATRIX N*N                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*N PRODUCT A*B    *                                     
*                                         *
******************************************}
Procedure MATMULT(n:Integer;A,B:pMAT; VAR C: pMAT);
VAR
    SUM : DOUBLE;
    I,J,K : integer;
BEGIN                                               
  for I:=0 to n-1 do                                                                  
    for J:=0 to n-1 do
    begin                                                                
      SUM:= 0.0;                                                                
      for K:=0 to n-1 do
      SUM:=SUM+A^[I,K]*B^[K,J];                                               
      C^[I,J]:=SUM                                                            
    end                                                                   
END;

{copy MAT A in MAT A1}
Procedure MatCopy(n:Integer; A:pMAT; VAR A1:pMAT);
Var i,j: Integer;
Begin
  For i:=0 to n-1 do
    For j:=0 to n-1 do
      A1^[i,j]:=A^[i,j]
End;

{ main program to demonstrate the use of function cholsl()  }
BEGIN

  New(A); New(A1); New(B);
  writeln(' Inversion of a square real symmetric matrix by Cholesky method');
  writeln(' (The matrix must positive definite).');

  n := 4;
  writeln;
  writeln(' Size = ', n);

  {define lower half of symmetrical matrix}
  A^[0,0]:= 5;
  A^[1,0]:=-1; A^[1,1]:= 5;
  A^[2,0]:=-1; A^[2,1]:=-1; A^[2,2]:= 5;
  A^[3,0]:=-1; A^[3,1]:=-1; A^[3,2]:=-1; A^[3,3]:= 5;

  {define upper half by symmetry}
  For i:=0 to n-1 do
    For j:=i+1 to n-1 do
      A^[i,j]:=A^[j,i];

  MatPrint('Matrix A:',n,A);
  writeln;

  if Check_Matrix(n,A) then
  begin
    MatCopy(n,A,A1);
    writeln(' Determinant = ', choldet(n,A):10:6);
    cholsl(n,A,B);
    MatPrint('Matrix Inv(A):',n,B)
  end
  else
    writeln(' Sorry, this matrix is not positive definite !');

  writeln;
  write(' Do you want a verification (y/n)?: ');
  if Readkey='y' then
  begin
    New(C);
    MatMult(n,A1,B,C);
    Writeln;
    MatPrint('Verification A * Inv(A) = I:',n,C);
    Dispose(C); Readkey
  end;

  Dispose(A); Dispose(A1); Dispose(B);
  DoneWinCrt

END.

{end of file choles.pas}