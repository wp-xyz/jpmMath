{************************************************************************
* SOLVE A LINEAR SYSTEM AX=B BY THE SINGULAR VALUE DECOMPOSITION METHOD *
*                                                                       *
* This method is recommanded when the matrix A is ill-conditionned or   *
* near-singular, or there are fewer equations than unknowns.            *
* --------------------------------------------------------------------- *
* SAMPLE RUN:                                                           *
*  1. Solve normal linear system of size 5 x 5:                         *
*        Matrix A:          Vector B:                                   *
*    1  2  0  0      0         3                                        *
*    2  4  7  0      0        13                                        *
*    0  7 10  8      0        25                                        *
*    0  0  8 -0.75  -9        -1.75                                     *
*    0  0  0 -9     10         1                                        *
*                                                                       *
* Input file contains:                                                  *
* 5                                                                     *
* 5                                                                     *
* 5                                                                     *
*   1.000000   2.000000   0.000000   0.000000   0.000000                *
*   2.000000   4.000000   7.000000   0.000000   0.000000                *
*   0.000000   7.000000  10.000000   8.000000   0.000000                * 
*   0.000000   0.000000   8.000000  -0.750000  -9.000000                *
*   0.000000   0.000000   0.000000  -9.000000  10.000000                *
*   3.000000                                                            *
*  13.000000                                                            *
*  25.000000                                                            *
*  -1.750000                                                            *
*   1.000000                                                            *
*                                                                       *
* Output file contains:                                                 *
* M = 5                                                                 *
* N = 5                                                                 *
* Items per line: 5                                                     *
* Matrix A:                                                             *
*      1.00      2.00      0.00      0.00      0.00                     *
*      2.00      4.00      7.00      0.00      0.00                     *
*      0.00      7.00     10.00      8.00      0.00                     *
*      0.00      0.00      8.00     -0.75     -9.00                     *
*      0.00      0.00      0.00     -9.00     10.00                     *
*                                                                       *
* Vector B:                                                             *
*      3.00     13.00     25.00     -1.75      1.00                     *
*                                                                       *
*  Matrix U:                                                            *
* -0.034422  0.690358 -0.717675  0.071746 -0.044894                     *
* -0.311812 -0.614665 -0.549190  0.413973  0.227984                     *
* -0.663543  0.222546  0.320343  0.484542 -0.415673                     *
* -0.483350  0.237931  0.181439 -0.208391  0.795873                     *
*  0.477150  0.198631  0.218615  0.738425  0.373911                     *
*                                                                       *
*  Vector W:                                                            *
*     19.12      0.78      2.53     12.54      9.16                     *
*                                                                       *
*  Matrix V:                                                            *
* -0.034422 -0.690358 -0.717675  0.071746  0.044894                     *
* -0.311812  0.614665 -0.549190  0.413973 -0.227984                     *
* -0.663543 -0.222546  0.320343  0.484542  0.415673                     *
* -0.483350 -0.237931  0.181439 -0.208391 -0.795873                     *
*  0.477150 -0.198631  0.218615  0.738425 -0.373911                     *
*                                                                       *
*  Solution:                                                            *
*  1.000000  1.000000  1.000000  1.000000  1.000000                     *
*                                                                       *
*  2. Solve linear system of size 4 x 5 (m < n):                        *
*        Matrix A:          Vector B:                                   *
*    1  2  0  0      0         3                                        *
*    2  4  7  0      0        13                                        *
*    0  7 10  8      0        25                                        *
*    0  0  8 -0.75  -9        -1.75                                     *
*                                                                       *
* Input file contains:                                                  *
* 4                                                                     *
* 5                                                                     *
* 5                                                                     *
*   1.000000   2.000000   0.000000   0.000000   0.000000                *
*   2.000000   4.000000   7.000000   0.000000   0.000000                *
*   0.000000   7.000000  10.000000   8.000000   0.000000                * 
*   0.000000   0.000000   8.000000  -0.750000  -9.000000                *
*   3.000000                                                            *
*  13.000000                                                            *
*  25.000000                                                            *
*  -1.750000                                                            *
*                                                                       *
* Output file contains:                                                 *
* M = 4                                                                 *
* N = 5                                                                 *
* Items per line: 5                                                     *
* Matrix A:                                                             *
*      1.00      2.00      0.00      0.00      0.00                     *
*      2.00      4.00      7.00      0.00      0.00                     *
*      0.00      7.00     10.00      8.00      0.00                     *
*      0.00      0.00      8.00     -0.75     -9.00                     *
*                                                                       *
* Vector B:                                                             *
*      3.00     13.00     25.00     -1.75                               *
*  Matrix U:                                                            *
* -0.048526  0.080513 -0.252041  0.000000  0.963140                     *
* -0.420052  0.055591 -0.869634  0.000000 -0.253383                     *
* -0.769221  0.500339  0.396761 -0.000000  0.023246                     *
* -0.479062 -0.860284  0.150972 -0.000000  0.087286                     *
*                                                                       *
*  Vector W:                                                            *
*     17.71      9.95      4.18      0.00      1.65                     *
*                                                                       *
*  Matrix V:                                                            *
* -0.050190  0.019275 -0.476492  0.832543  0.277376                     *
* -0.404496  0.390720 -0.288408 -0.416271  0.653651                     *
* -0.816982 -0.149796 -0.218240 -0.000000 -0.512322                     *
* -0.327269  0.467357  0.732420  0.364238  0.073239                     *
*  0.243515  0.778527 -0.325129 -0.030353 -0.477457                     *
*                                                                       *
*  Solution:                                                            *
*  0.375463  1.312268  1.000000  0.726765  1.022770                     *
*                                                                       *
*  3. Solve linear system of size 5 x 5 (near-singular):                *
*        Matrix A:           Vector B:                                  *
*    1  2  0     0      0        3                                      *
*    2  4  1e-6  0      0        6.000001                               *
*    0  7 10     8      0       25                                      *
*    0  0  8    -0.75  -9       -1.75                                   *
*    0  0  0    -9     10        1                                      *
*                                                                       *
* Input file contains:                                                  *
* 5                                                                     *
* 5                                                                     *
* 5                                                                     *
*   1.000000   2.000000   0.000000   0.000000   0.000000                *
*   2.000000   4.000000   0.000001   0.000000   0.000000                *
*   0.000000   7.000000  10.000000   8.000000   0.000000                * 
*   0.000000   0.000000   8.000000  -0.750000  -9.000000                *
*   0.000000   0.000000   0.000000  -9.000000  10.000000                *
*   3.000000                                                            *
*   6.000001                                                            *
*  25.000000                                                            *
*  -1.750000                                                            *
*   1.000000                                                            *
*                                                                       *
* Output file contains:                                                 *
* M = 5                                                                 *
* N = 5                                                                 *
* Items per line: 5                                                     *
* Matrix A:                                                             *
*      1.00      2.00      0.00      0.00      0.00                     *
*      2.00      4.00      0.00      0.00      0.00                     *
*      0.00      7.00     10.00      8.00      0.00                     *
*      0.00      0.00      8.00     -0.75     -9.00                     *
*      0.00      0.00      0.00     -9.00     10.00                     *
*                                                                       *
* Vector B:                                                             *
*      3.00      6.00     25.00     -1.75      1.00                     *
*                                                                       *
*  Matrix U:                                                            *
* -0.029197  0.894427 -0.435067  0.093070 -0.034671                     *
* -0.058394 -0.447214 -0.870134  0.186141 -0.069341                     *
* -0.649209  0.000000  0.207861  0.720418 -0.127747                     *
* -0.500300  0.000000 -0.091238 -0.280151  0.814181                     *
*  0.569179  0.000000  0.045304  0.599335  0.561053                     *
*                                                                       *
*  Vector W:                                                            *
*     18.34      0.00      4.28     11.55      8.75                     *
*                                                                       *
*  Matrix V:                                                            *
* -0.007961 -0.859939 -0.508355  0.040295 -0.019809                     *
* -0.263732  0.429970 -0.676684  0.517265 -0.141801                     *
* -0.572267 -0.174666  0.315179  0.429750  0.598313                     *
* -0.542088 -0.157890  0.309309  0.050175 -0.763559                     *
*  0.555908 -0.142101  0.297764  0.737300 -0.196212                     *
*                                                                       *
*  Solution:                                                            *
*  0.222075  1.388963  0.841992  0.857168  0.871451                     *
*                                                                       *
* --------------------------------------------------------------------- *
* Reference:  "Numerical Recipes By W.H. Press, B. P. Flannery,         *
*              S.A. Teukolsky and W.T. Vetterling, Cambridge            *
*              University Press, 1986" [BIBLI 08].                      *
*                                                                       *
*                               TPW Release 2.0 By J-P Moreau, Paris    *
*                                        (www.jpmoreau.fr)              *
* --------------------------------------------------------------------- *
* Release 2.0: Added dynamic allocations, I/O files and two more        *
*              examples.                                                *
************************************************************************}
PROGRAM TEST_SVBKSB;
Uses WinCrt;

Const
      MMAX=100;
      NMAX=100;

Type
     pMAT = ^MAT;
     MAT = Array[1..MMAX,1..NMAX] of REAL;
     pVEC = ^VEC;
     VEC = Array[1..MMAX] of REAL;

Var
     i,j,k,l,m,n,nx: Integer;
     A,U,V: pMAT;
     B,W,X: pVEC;
     WMIN,WMAX: REAL;

     fp:TEXT;
     filename: String;


Procedure svbksb(var u:pMAT; var w:pVEC; var v:pMAT; m,n:Integer; b:pVEC; var x:pVEC);
{------------------------------------------------------------------------------------------- 
! Solves A · X = B for a vector X, where A is specified by the arrays u, w, v as returned by 
! svdcmp. m and n are the dimensions of a, and will be equal for square matrices. b(1:m) is 
! the input right-hand side. x(1:n) is the output solution vector. No input quantities are 
! destroyed, so the routine may be called sequentially with different b’s. 
!------------------------------------------------------------------------------------------}
Var i,j,jj: Integer; 
    s: REAL;
    tmp: pVEC;
Begin
  New(tmp);
  for j:=1 to n do  {Calculate UTB}
  begin
    s:=0.0;
    if w^[j] <> 0.0 then  {Nonzero result only if wj is nonzero}
    begin 
      for i:=1 to m do s:=s+u^[i,j]*b^[i];
      s:=s/w^[j]          {This is the divide by wj }
    end; 
    tmp^[j]:=s
  end; 
  for j:=1 to n do       {Matrix multiply by V to get answer }
  begin
    s:=0.0;
    for jj:=1 to n do s:=s+v^[j,jj]*tmp^[jj];
    x^[j]:=s
  end;
  Dispose(tmp) 
End; 

{Note that a typical use of svdcmp and svbksb superficially resembles the 
 typical use of ludcmp and lubksb: In both cases, you decompose the left-hand 
 matrix A just once, and then can use the decomposition either once or many times 
 with different right-hand sides. The crucial difference is the "editing" of the singular
 values before using SVDCMP. }

Function IMin(a,b:Integer): Integer;
Begin
  if a<=b then IMin:=a
  else IMin:=b
End;

Function IMax(a,b:Integer): Integer;
Begin
  if a>=b then IMax:=a
  else IMax:=b
End;

Function Max(a,b:REAL): REAL;
Begin
  if a>=b then Max:=a
  else Max:=b
End;

Function pythag(a,b:REAL): REAL;
{Computes sqrt(a*a + b*b) without destructive underflow or overflow}
Var absa,absb: REAL;
Begin 
  absa:=abs(a);
  absb:=abs(b);
  if absa > absb then 
    pythag:=absa*sqrt(1.0+(absb/absa)*(absb/absa))
  else 
    if absb = 0.0 then 
      pythag:=0.0
    else
      pythag:=absb*sqrt(1.0+(absa/absb)*(absa/absb))
End;

Function Sign(a,b : REAL) : REAL;
Begin
  if (b <0.0) then Sign := - Abs(a)
              else Sign :=   Abs(a)
End;

Procedure svdcmp(var a:pMAT; m,n:Integer; var w:pVEC; var v:pMAT); 
{-------------------------------------------------------------------------------------- 
  Given a matrix a(1:m,1:n), this routine computes its singular value decomposition, 
  A := U · W · Vt. The matrix U replaces a on output. The diagonal matrix of singular
  values W is output as a vector w(1:n). The matrix V (not the transpose Vt) is output 
  as v(1:n,1:n). 
--------------------------------------------------------------------------------------}
Label 1,2,3;
Var i,its,j,jj,k,l,nm: Integer; 
    anorm,c,f,g,h,s,scale,x,y,z: REAL;
    rv1: pVEC;
Begin
  New(rv1);
  g:=0.0;  {Householder reduction to bidiagonal form }
  scale:=0.0;
  anorm:=0.0;
for i:=1 to n do
begin
  l:=i+1;
  rv1^[i]:=scale*g;
  g:=0.0;
  s:=0.0;
  scale:=0.0;
  if i <= m then
  begin 
    for k:=i to m do scale:=scale+abs(a^[k,i]);
    if scale <> 0.0 then
    begin 
      for k:=i to m do
      begin
        a^[k,i]:=a^[k,i]/scale;
        s:=s+a^[k,i]*a^[k,i]
      end; 
      f:=a^[i,i];
      g:=-Sign(sqrt(s),f);
      h:=f*g-s;
      a^[i,i]:=f-g;
      for j:=l to n do
      begin
        s:=0.0;
        for k:=i to m do s:=s+a^[k,i]*a^[k,j];
        f:=s/h;
        for k:=i to m do a^[k,j]:=a^[k,j]+f*a^[k,i]
      end; 
      for k:=i to m do a^[k,i]:=scale*a^[k,i]
    end 
  end; 
  w^[i]:=scale*g;
  g:=0.0;
  s:=0.0;
  scale:=0.0;
  if (i <= m) and (i <> n) then
  begin 
    for k:=l to n do scale:=scale+abs(a^[i,k]);
    if scale <> 0.0 then
    begin 
      for k:=l to n do
      begin
        a^[i,k]:=a^[i,k]/scale;
        s:=s+a^[i,k]*a^[i,k]
      end; 
      f:=a^[i,l];
      g:=-Sign(sqrt(s),f);
      h:=f*g-s;
      a^[i,l]:=f-g;
      for k:=l to n do rv1^[k]:=a^[i,k]/h;
      for j:=l to m do
      begin
        s:=0.0;
        for k:=l to n do s:=s+a^[j,k]*a^[i,k];
        for k:=l to n do a^[j,k]:=a^[j,k]+s*rv1^[k]
      end; 
      for k:=l to n do a^[i,k]:=scale*a^[i,k]
    end 
  end; 
  anorm:=Max(anorm,(abs(w^[i])+abs(rv1^[i])))
end; {for i:=1 to n}
 
for i:=n downto 1 do    {Accumulation of right-hand transformations }
begin
  if i < n then
  begin 
    if g <> 0.0 then
    begin 
      for j:=l to n do  {Double division to avoid possible underflow }
        v^[j,i]:=(a^[i,j]/a^[i,l])/g;
      for j:=l to n do
      begin
        s:=0.0;
        for k:=l to n do s:=s+a^[i,k]*v^[k,j];
        for k:=l to n do v^[k,j]:=v^[k,j]+s*v^[k,i]
      end 
    end; 
    for j:=l to n do
    begin
      v^[i,j]:=0.0;
      v^[j,i]:=0.0
    end 
  end; 
  v^[i,i]:=1.0;
  g:=rv1^[i];
  l:=i
end; 

for i:=IMin(m,n) Downto 1 do  {Accumulation of left-hand transformations }
begin
  l:=i+1;
  g:=w^[i];
  for j:=l to n do a^[i,j]:=0.0;
  if g <> 0.0 then
  begin 
    g:=1.0/g;
    for j:=l to n do
    begin
      s:=0.0;
      for k:=l to m do s:=s+a^[k,i]*a^[k,j];
      f:=(s/a^[i,i])*g;
      for k:=i to m do a^[k,j]:=a^[k,j]+f*a^[k,i]
    end; 
    for j:=i to m do a^[j,i]:=a^[j,i]*g
  end   
  else
    for j:= i to m do a^[j,i]:=0.0;
  a^[i,i]:=a^[i,i]+1.0
end; 

for k:=n Downto 1 do  {Diagonalization of the bidiagonal form: Loop over}
begin                 {singular values, and over allowed iterations } 
for its:=1 to 30 do
begin
for l:=k Downto 1 do  {Test for splitting }
begin
  nm:=l-1;  {Note that rv1(1) is always zero }
  if (abs(rv1^[l])+anorm) = anorm then goto 2;
  if (abs(w^[nm])+anorm) = anorm then goto 1
end; 
1: c:=0.0;  {Cancellation of rv1(l), if l > 1 }
s:=1.0;
for i:=l to k do
begin
  f:=s*rv1^[i];
  rv1^[i]:=c*rv1^[i];
  if (abs(f)+anorm) = anorm then goto 2; 
  g:=w^[i];
  h:=pythag(f,g);
  w^[i]:=h;
  h:=1.0/h;
  c:=g*h;
  s:=-(f*h);
  for j:=1 to m do
  begin
    y:=a^[j,nm];
    z:=a^[j,i];
    a^[j,nm]:=(y*c)+(z*s);
    a^[j,i]:=-(y*s)+(z*c)
  end 
end; 
2: z:=w^[k];
if l = k then       {Convergence }
begin
  if z < 0.0 then   {Singular value is made nonnegative }
  begin 
    w^[k]:=-z;
    for j:=1 to n do v^[j,k]:=-v^[j,k]
  end; 
  goto 3 
end; 
if its = 30 then  writeln(' No convergence in svdcmp'); 
x:=w^[l];            {Shift from bottom 2-by-2 minor }
nm:=k-1;
y:=w^[nm];
g:=rv1^[nm];
h:=rv1^[k];
f:=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
g:=pythag(f,1.0);
f:=((x-z)*(x+z)+h*((y/(f+Sign(g,f)))-h))/x;
c:=1.0;             {Next QR transformation: }
s:=1.0;
for j:=l to nm do
begin
  i:=j+1;
  g:=rv1^[i];
  y:=w^[i];
  h:=s*g;
  g:=c*g;
  z:=pythag(f,h);
  rv1^[j]:=z;
  c:=f/z;
  s:=h/z;
  f:= (x*c)+(g*s);
  g:=-(x*s)+(g*c);
  h:=y*s;
  y:=y*c;
  for jj:=1 to n do
  begin
    x:=v^[jj,j];
    z:=v^[jj,i];
    v^[jj,j]:= (x*c)+(z*s);
    v^[jj,i]:=-(x*s)+(z*c)
  end; 
  z:=pythag(f,h);
  w^[j]:=z;   {Rotation can be arbitrary if z = 0 }
  if z <> 0.0 then
  begin 
    z:=1.0/z;
    c:=f*z;
    s:=h*z
  end; 
  f:= (c*g)+(s*y);
  x:=-(s*g)+(c*y);
  for jj:=1 to m do
  begin
    y:=a^[jj,j];
    z:=a^[jj,i];
    a^[jj,j]:= (y*c)+(z*s);
    a^[jj,i]:=-(y*s)+(z*c)
  end 
end;  {for j:=l to nm}
rv1^[l]:=0.0;
rv1^[k]:=f;
w^[k]:=x

end;   {for its:=1 to 30}
 
3:end; {for k:=n downto 1}

Dispose(rv1)
End;
 
Procedure WriteMat(Title:String;A:pMAT;m,n,p,q:Integer);
Var i,j,l: Integer;
Begin
  writeln(fp,title);
  l:=0;
  For i:=1 to m do
    for j:=1 to n do
    begin
      Inc(l);
      if l mod k <> 0 then write(fp,A^[i,j]:p:q)
      else writeln(fp,A^[i,j]:p:q)
    end;
  writeln(fp)
End;

Procedure WriteVec(Title:String;V:pVEC;n,p,q:Integer);
Var i,l: Integer;
Begin
  writeln(fp,title);
  l:=0;
  For i:=1 to n do
  begin
    Inc(l);
    if l mod k <> 0 then Write(fp,V^[i]:p:q)
    else writeln(fp,V^[i]:p:q)
  end;
  writeln(fp)
End;


{main program}
BEGIN

  New(A); New(U); New(V); New(B); New(W); New(X);

  Writeln;
  Write(' Input data file name: '); Readln(filename); 
  Assign(fp,filename); Reset(fp);
{size of linear system}
  Readln(fp, m);
  Readln(fp, n);
{read number of items per line}
  Readln(fp, k);
{Ser matrix A to zero}
  nx:=IMax(m,n);
  For i:=1 to nx do
    For j:=1 to nx do
      A^[i,j]:=0.0;
{Read matrix A }
  l:=0;
  For i:=1 to m do
    For j:=1 to n do
    begin
      Inc(l);
      if l mod k <>0 then read(fp,A^[i,j])
      else Readln(fp,A^[i,j])
    end;
{read right-hand side B}
  For i:=1 to nx do B^[i]:=0.0;
  For i:=1 to m do Readln(fp,B^[i]);
  Close (fp);

  Assign(fp,'tsvbksb.lst'); Rewrite(fp);
  Writeln(fp,' M = ', m);
  Writeln(fp,' N = ', n);
  Writeln(fp,' Items per line: ',k);
  WriteMat(' Matrix A:',A,m,n,10,2);
  WriteVec(' Vector B:',B,m,10,2);
 
{ Save A in U }
  For i:=1 to m do
    For j:=1 to n do
      U^[i,j]:=A^[i,j];

{call singular value decomposition subroutine }
  SVDCMP(U,m,n,W,V);
  WriteMat('  Matrix U:',U,m,n,10,6);
  WriteVec('  Vector W:',W,n,10,2);
  WriteMat('  Matrix V:',V,n,n,10,6);
{seek highest value of W's and set near-zero
 values to exactly zero (for near-singular cases) }    
  WMAX:=0.0;
  For j:=1 to n do
    if W^[j] > WMAX then WMAX:=W^[j];
  WMIN:=WMAX*1e-6;
  For j:=1 to n do
    if W^[j] < WMIN then W^[j]:=0.0;
{call solver for SVD matrix }
  SVBKSB(U,W,V,m,n,B,X);
{print solution }
  WriteVec('  Solution:',X,n,10,6);
{exit section}
  Close(fp);
  Writeln;
  Writeln(' Results in file tsvbksb.lst.');
  ReadKey;

  Dispose(A); Dispose(U); Dispose(V); Dispose(B); Dispose(W); Dispose(X);
  DoneWinCrt

END.

{end of file tsvbksb.pas}