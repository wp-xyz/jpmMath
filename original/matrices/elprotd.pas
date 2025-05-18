{***************************************************
*  Calculate the eigenvalues and the eigenvectors  *
*      of a real symmetric tridiagonal matrix      *
* ------------------------------------------------ *
* Ref.: "NUMERICAL RECIPES, Cambridge University   *
*        Press, 1986" [BIBLI 08].                  *
*                                                  *
* SAMPLE RUN:                                      *
* (find eigenvalues and eigenvectors of matrix:    *
*          1  2  0  0     0                        *
*          2  4  7  0     0                        *
*     A =  0  7 10  8     0                        *
*          0  0  8 -0.75 -9                        *
*          0  0  0 -9    10  )                     *
*                                                  *
* Eigenvalues 1:   -0.78071442                     *
* Eigenvector:                                     *
* 1.000000 -0.890357  0.322363  0.344649  0.287721 *
*                                                  *
* Eigenvalues 2:    2.53046878                     *
* Eigenvector:                                     *
* 1.000000  0.765234 -0.446362 -0.252815 -0.304616 *
*                                                  *
* Eigenvalues 3:   -9.15659229                     *
* Eigenvector:                                     *
*-0.056408  0.286458 -0.522285  1.000000  0.469812 *
*                                                  *
* Eigenvalues 4:   12.53989066                     *
* Eigenvector:                                     *
* 0.097161  0.560616  0.656182 -0.282210  1.000000 *
*                                                  *
* Eigenvalues 5:   19.11694726                     *
* Eigenvector:                                     *
* 0.051876  0.469920  1.000000  0.728439 -0.719095 *
*                                                  *
*             Pascal version by J-P Moreau, Paris  *
*                     (www.jpmoreau.fr)            *
***************************************************}
PROGRAM TEST_TQLI;
Uses WinCrt;

Const
      NMAX = 30;

Type
      pM  = ^MAT;
      MAT = Array[1..NMAX,1..NMAX] of real;
      pV  = ^VEC;
      VEC = Array[1..NMAX] of real;

Var
      D, E : pV;
      Z    : pM;
      i,j,n: integer;
      max  : real;


Function Sign(a,b:real):real;
begin
  if b>=0 then Sign:=ABS(a)
          else Sign:=-ABS(a)
end;

{read vectors D and E from screen}
Procedure Read_data(VAR n:integer; VAR D,E:pV);
var i:integer;
Begin
  writeln;
  write(' Input size of matrix: '); readln(n);
  writeln(' Input ',n-1,' elements of subdiagonal:');
  E^[1]:=1.0;   {arbitrary value}
  for i:=1 to n-1 do
  begin
    write('   Element ',i,': '); readln(E^[i+1])
  end;
  writeln(' Input ',n,' elements of main diagonal:');
  for i:=1 to n do
  begin
    write('   Element ',i,': '); readln(D^[i])
  end
End;

{***************************************************
* This subroutine implements the QL algorithm with *
* implicit shifts to determine the eigenvalues and *
* eigenvectors of a real symmetric tridiagonal     *
* matrix.                                          *
* ------------------------------------------------ *
* INPUTS:                                          *
*         D    Elements of main diagonal           *
*         E    Elements of subdiagonal (from E(2)  *
*              to E(n))                            *
*         n    size of matrix                      *
*         Z    identity matrix                     *
* OUTPUTS                                          *
*         D    vector storing the n eigenvalues    *
*         Z    matrix storing the n eigenvectors   *
*              in columns.                         *
***************************************************}
Procedure TQLI(Var D,E:pV;n:integer; VAR Z:pM);
Label 1,2,fin;
Var b,c,dd,f,g,p,r,s: real;
    i,iter,k,l,m:integer;
Begin
 if n>1 then
 begin
  for i:=2 to n do E^[i-1]:=E^[i];
  E^[n]:=0.0;
  for l:=1 to n do
  begin
    iter:=0;
1:  for m:=l to n-1 do
    begin
      dd:=ABS(D^[m])+ABS(D^[m+1]);
      if ABS(E^[m])+dd=dd then goto 2
    end;
    m:=n;
2:  if m<>l then
    begin
      if iter=30  then goto fin;
      iter:=iter+1;
      g:=(D^[l+1]-D^[l])/(2.0*E^[l]);
      r:=SQRT(g*g+1.0);
      g:=D^[m]-D^[l]+E^[l]/(g+SIGN(r,g));
      s:=1.0; c:=1.0; p:=0.0;
      for i:=m-1 Downto l do
      begin
        f:=s*E^[i]; b:=c*E^[i];
	if ABS(f)>=ABS(g) then
        begin
	  c:=g/f;
	  r:=SQRT(c*c+1.0);
	  E^[i+1]:=f*r;
	  s:=1.0/r;
	  c:=c*s
        end
	else
        begin
	  s:=f/g;
	  r:=SQRT(s*s+1.0);
	  E^[i+1]:=g*r;
	  c:=1.0/r;
	  s:=s*c
        end;
	g:=D^[i+1]-p;
	r:=(D^[i]-g)*s+2.0*c*b;
	p:=s*r;
	D^[i+1]:=g+p;
	g:=c*r-b;
	for k:=1 to n do
        begin
          f:=Z^[k,i+1];
	  Z^[k,i+1]:=s*Z^[k,i]+c*f;
	  Z^[k,i]:=c*Z^[k,i]-s*f
        end
      end; 
      D^[l]:=D^[l]-p;
      E^[l]:=g; E^[m]:=0.0;
      goto 1
    end
  end; {of l loop}
fin:if iter>0 then
    writeln(' ',iter,' iterations.');
 end
end;

{main program}
BEGIN

  {allocate vectors and matrix}
  New(D); New(E); New(Z);

  read_data(n,D,E);

  {initialize Z to identity matrix}
  for i:=1 to n do
    for j:=1 to n do
      if i<>j then
        Z^[i,j] := 0.0
      else   
        Z^[i,j] := 1.0;

  {calculate eigenvalues and eigenvectors}
  TQLI(D,E,n,Z);

  Dispose(E);

  {print results}
  for  j:=1 to n do
  begin
    writeln;
    writeln(' Eigenvalue ',j,': ', D^[j]:12:8);
    writeln(' Eigenvector:');
    max:=ABS(Z^[1,j]);
    for i:=2 to n do
      if ABS(Z^[i,j])>ABS(max) then max:=Z^[i,j];
    for i:=1 to n do
      write(' ',(Z^[i,j]/max):10:6);
    writeln
  end;

  writeln;
  readkey;
  Dispose(D); Dispose(Z);
  DoneWinCrt

END.

{end of file elprotd.pas}

