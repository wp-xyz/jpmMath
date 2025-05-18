{*******************************************************
* Eigenvalues and eigenvectors of a real square matrix *
* by Rutishauser's method and inverse iteration method *
* ---------------------------------------------------- *
* Reference:                                           *
*                                                      *
*   "ALGEBRE Algorithmes et programmes en Pascal       *
*    de Jean-Louis Jardrin - Dunod BO-PRE 1988."       *
*    [BIBLI 10].                                       *
* ---------------------------------------------------- *
* SAMPLE RUN:                                          *
*                                                      *
* Input file: elpro.dat                                *
*                                                      *
* 5                                                    *
*   1   2   3  -7    12                                *
*   2   4   7   3    -1                                *       
*   3   7  10   8     4                                *
*  -7   3   8  -0.75 -9                                *
*  12  -1   4  -9    10                                *
*                                                      *
* Output to screen:                                    *
*                                                      *
*        *** EIGENVALUES AND EIGENVECTORS ***          *
*              OF A REAL SQUARE MATRIX                 *
*              BY RUTISHAUSER'S METHOD                 *
*            AND INVERSE ITERATION METHOD              *
*                                                      *
*  Matrix A:                                           *
*                                                      *
*  1.000000  2.000000  3.000000 -7.000000 12.000000    *
*  2.000000  4.000000  7.000000  3.000000 -1.000000    *
*  3.000000  7.000000 10.000000  8.000000  4.000000    *
* -7.000000  3.000000  8.000000 -0.750000 -9.000000    *
* 12.000000 -1.000000  4.000000 -9.000000 10.000000    *
*                                                      *
*                                                      *
*  Eigenvalue 1: 23.755955                             *
*                                                      *
*  Eigenvector:                                        *
*                                                      *
*  0.705820 -0.012677  0.131486 -0.527499  1.000000    *
*                                                      *
*  Eigenvalue 2:-10.486545                             *
*                                                      *
*  Eigenvector:                                        *
*                                                      *
*  0.467172 -0.007947 -0.507827  1.000000  0.264432    *
*                                                      *
*  Eigenvalue 3: 0.463350                              *
*                                                      *
*  Eigenvector:                                        *
*                                                      *
*  0.208812  1.000000 -0.478027 -0.275000 -0.216915    *
*                                                      *
*  Eigenvalue 4:-7.774580                              *
*                                                      *
*  Eigenvector:                                        *
*                                                      *
*  1.000000 -0.326100  0.209620 -0.147689 -0.815422    *
*                                                      *
*  Eigenvalue 5: 0.463350                              *
*                                                      *
*  Eigenvector:                                        *
*                                                      *
*  0.208812  1.000000 -0.478027 -0.275000 -0.216915    *
*                                                      * 
*******************************************************}
Program ELPRO;
Uses WinCrt;

Const NMAX = 10;

Type  MAT = Array[1..NMAX,1..NMAX] of Double;
      VEC = Array[1..NMAX] of Double;

Var   i,it,j,m,n : integer;
      d1,d2,eps  : Double;
      A, VX      : MAT;
      R          : VEC;
      F          : TEXT;  
      re         : char;

  Procedure Data;
  var i,j: integer;
  begin
    Assign(F,'elpro.dat'); Reset(F);
    Read(F,n);
    Writeln;
    Writeln('  Matrix A:');
    Writeln;
    for i:=1 to n do
    begin
      for j:=1 to n do
      begin
        Read(F,A[i,j]);
        write(A[i,j]:10:6)
      end;
      writeln
    end;
    Writeln;
    Close(F);
    eps:=1E-16;
    d1:=1E-10;
    d2:=1E-8;
    m:=600
  end;

  {**************************************************************
  * Procedure DECCRM determines the lower triangular matrix and *
  * the upper triangukar matrix of Crout's decomposition of a   *
  * given square real matrix, A.                                *
  * ----------------------------------------------------------- *
  * INPUTS:                                                     *
  *         eps: required precision (double)                    *
  *          n : size of matrix A (integer)                     *
  *          A : input matrix (n x n)                           *
  * OUTPUTS:                                                    *
  *          it: flag, =0 if method does not apply              *
  *                    =1 if method is ok.                      *
  *           U: lower triangular matrix.                       *
  *           V: upper triangular matrix.                       *
  **************************************************************}
  Procedure DECCRM(eps:Double; n:integer;A:MAT; VAR it:integer;VAR U,V : MAT);
  VAR  i,j,k : integer;
       s     : Double;
  begin
    if ABS(A[1,1]) < eps then it:=0
    else
    begin
      for i:=1 to n do U[i,1]:=A[i,1];
      V[1,1]:=1.0;
      for j:=2 to n do V[1,j]:=A[1,j]/U[1,1];
      it:=1; k:=2;
      Repeat
	for i:=1 to n do
	  if i < k then U[i,k]:=0.0
	  else
	  begin
	    s:=0.0;
	    for j:=1 to k-1 do s:=s+U[i,j]*V[j,k];
	    U[i,k]:=A[i,k]-s
	  end;
	if ABS(U[k,k]) < eps then it:=0
	else
	begin
	  for j:=1 to n do
	    if j < k then V[k,j]:=0.0
	    else if j=k then V[k,j]:=1.0
            else
	    begin
	      s:=0.0;
	      for i:=1 to k-1 do s:=s+U[k,i]*V[i,j];
	      V[k,j]:=A[k,j]/U[k,k]
	    end;
	  Inc(k)
	end
      Until (it=0) OR (k > n)
    end
  end;

  Procedure MatPrint(title:String; A:MAT; n:integer);
  Var i,j:integer;
  Begin
    Writeln;
    Writeln(title);
    For i:=1 to n do
    begin
      For j:=1 to n do write(' ',A[i,j]:4:8);
      writeln
    end
  End;

  Procedure VecPrint(title:String; X:VEC; n:integer);
  Var i:integer;
  Begin
    Writeln;
    Writeln(title);
    For i:=1 to n do write(' ',X[i]:4:8);
    writeln
  End;


  {********************************************************
  * Calculate the eigenvalues of a real square matrix by  *
  * Rutishauser's Method.                                 *
  * ----------------------------------------------------- *
  * INPUTS:                                               *
  *        eps: absolute precision (double)               *
  *        dta: relative precision (double)               *
  *         m : maximum number of iterations (integer)    *
  *         n : size of matrix A (integer)                *
  *         A : input real square matrix (n x n)          *
  * OUTPUTS:                                              *
  *         it: flag, =-1 if convergence is not obtained  *
  *                   =1 if convergence is ok.            *
  *         R : contains in output the n eigenvalues of A *
  *                                                       *         
  ********************************************************}
  Procedure VAMR(eps,dta:Double; m,n:integer;A:MAT; VAR it:integer;
							 VAR R: VEC);
  VAR  i,j,k,l : integer;
       phi,s,t0: Double;
       U,V     : MAT;
  begin
    t0:=0.0;
    l:=1;
    Repeat
      for i:=1 to n do R[i]:=A[i,i];

      DECCRM(eps,n,A,it,U,V);

      if it=0 then
      begin
	for i:=1 to n do A[i,i]:=A[i,i]+1.0;
	t0:=t0+1.0
      end
      else
      begin
	for i:=1 to n do
	  for j:=1 to n do
	  begin
	    s:=0.0;
	    for k:=1 to n do s:=s+V[i,k]*U[k,j];
	    A[i,j]:=s
	  end;
	phi:=0.0;
        for i:=1 to n do
	begin
	  s:=ABS(A[i,i]-R[i]);
	  if s > phi then phi:=s
	end;
	if phi < dta then
	  for i:=1 to n do R[i]:=A[i,i]-t0
	else
	begin
	  Inc(l);
	  it:=-1
	end
      end
    Until (l > m) OR (it=1)
  end;

  {***********************************************************
  * Procedure IIM calculates a real eigenvalue and the asso- *
  * ciated eigenvector of a real square matrix the inverse   *
  * iteration method.                                        *
  * -------------------------------------------------------- *
  * INPUTS:                                                  *
  *         eps : absolute precision (double)                *
  *         dta : relative precision (double)                *
  *          m  : maximum number of iterations (integer)     *
  *          n  : size of matrix A                           *
  *          A  : input real square matrix (n x n)           *
  * OUTPUTS:                                                 *
  *          it : flag, =-1 if convergence is not obtained   *
  *                     =1 if convergence is ok.             *
  *        Gamma: starting value for the eigenvalue as input *
  *               approximation of the eigenvalue with preci-*
  *               sion dta in output.                        *
  *          X1 : contains in output the associated eigen-   *
  *               vector.                                    *
  *                                                          *
  ***********************************************************}
  Procedure IIM(eps,dta:Double; m,n:integer; A:MAT; VAR it:integer;
					   VAR gamma:Double; VAR X1:VEC);
  VAR  i,j,k,l,l0 : integer;
       p0,phi,s,t0: Double;
       W,X0       : VEC;
       LP : Array[1..NMAX] of integer;
  begin
    for i:=1 to n do A[i,i]:=A[i,i]-gamma;
    for k:=1 to n-1 do
    begin
      p0:=A[k,k]; l0:=k;
      for i:=k+1 to n do
	if ABS(A[i,k]) > ABS(p0) then
        begin
	  p0:=A[i,k]; l0:=i
	end;
      LP[k]:=l0;
      if ABS(p0) < eps then
      begin
	p0:=eps; A[l0,k]:=eps
      end;
      if l0 <> k then
	for j:=k to n do
	begin
	  t0:=A[k,j]; A[k,j]:=A[l0,j]; A[l0,j]:=t0
	end;
      for i:=k+1 to n do
      begin
	A[i,k]:=A[i,k]/p0;
	for j:=k+1 to n do
	  A[i,j]:=A[i,j]-A[i,k]*A[k,j]
      end
    end; {end k loop}

    if ABS(A[n,n]) < eps then A[n,n]:=eps;
    for i:=1 to n do X0[i]:=1.0/SQRT(i);
    it:=-1; l:=1;
    While (it=-1) AND (l<=m) do
    begin
      for i:=1 to n do W[i]:=X0[i];
      for k:=1 to n-1 do
      begin
        l0:=LP[k];
        if l0<>k then
        begin
          t0:=W[k]; W[k]:=W[l0]; W[l0]:=t0
        end;
        for i:=k+1 to n do W[i]:=W[i]-A[i,k]*W[k]
      end;
      X1[n]:=W[n]/A[n,n];
      for i:=n-1 downto 1 do
      begin
        s:=0.0;
        for j:=i+1 to n do s:=s+A[i,j]*X1[j];
	X1[i]:=(W[i]-s)/A[i,i]
      end;
      p0:=0.0;
      for i:=1 to n do
        if ABS(X1[i]) > ABS(p0) then p0:=X1[i];
      for i:=1 to n do
        X1[i]:=X1[i]/p0;
      phi:=0.0;
      for i:=1 to n do
      begin
        s:=ABS(X1[i]-X0[i]);
        if s > phi then phi:=s
      end;

      if phi < dta then
      begin
        gamma:=gamma+1.0/p0;
        it:=1
      end
      else
      begin
        for i:=1 to n do X0[i]:=X1[i];
        Inc(l)
      end
    end { while }
  end;    	  						    

  {*****************************************************
  * INPUTS:                                            *
  * EPS : precision (Double)                           *
  * D1  : precision d1 (Double)                        *
  * D2  : precision d2 (Double)                        *
  * M   : maximum number of iterations (integer)       *
  * N   : order of matrix A (integer)                  *
  * A   : input matrix to study (of MAT type)          *
  * -------------------------------------------------- *
  * OUTPUTS:                                           *
  * IT  : -1 if no convergence is obtained (integer)   *
  * R   : table of eigenvalues (of VEC type)           *
  * VX  : table of eigenvectors (of MAT type)          *
  *****************************************************}
  Procedure EPMRI(eps,d1,d2:Double; m,n:integer; A:MAT; VAR it:integer;
						 VAR R:VEC; VAR VX:MAT);
  var  i,j,k: integer; X: VEC;
  begin

    VAMR(eps,d2,m,n,A,it,R);

    j:=1;
    While (it=1) AND (j<=n) do
    begin
      IIM(eps,d1,m,n,A,it,R[j],X);
      if it=1 then
      begin
	for i:=1 to n do VX[i,j]:=X[i];
	Inc(j)
      end
    end
  end;

  {main program}
  BEGIN
    Writeln;
    Writeln(' ':13,'*** EIGENVALUES AND EIGENVECTORS ***');
    Writeln(' ':19,'OF A REAL SQUARE MATRIX');
    Writeln(' ':19,'BY RUTISHAUSER''S METHOD');
    Writeln(' ':17,'AND INVERSE ITERATION METHOD');

    Data;    {read data from input text file}

    EPMRI(eps,d1,d2,m,n,A,it,R,VX);

    if it=-1 then Writeln(CHR(10),'  No convergence !')
    else
      for j:=1 to n do
      begin
	writeln;
	write(CHR(10),'  Eigenvalue ',j,':');
        writeln(R[j]:10:6);
        writeln;
	writeln('  Eigenvector:');
	writeln;
	for i:=1 to n do write(VX[i,j]:10:6);
        writeln; writeln;
	writeln(' < press space bar to continue...>');
        Repeat
	  re:=readkey
	Until re=' ';
	gotoxy(wherex,wherey-1); clreol
      end;
    Readkey;
    DoneWinCrt
  END.

{End of file elpro.pas}