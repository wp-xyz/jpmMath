{********************************************************
*    Approximation of a discrete real function F(x) by  *
*    least squares                                      *
* ----------------------------------------------------- *
* Ref.: "Méthodes de calcul numérique, Tome 2 By Claude *
*        Nowakowski, PSI Edition, 1984" [BIBLI 04].     *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
*                                                       *
* Number of points    : 11                              *
*                                                       *
* Degree of polynomial: 3                               *
*                                                       *
* Function to approximate:                              *
*  X(1), Y(1) = 0 0                                     *
*  X(2), Y(2) = 0.1 0.2                                 *
*  X(3), Y(3) = 0.2 0.4                                 *
*  X(4), Y(4) = 0.3 0.6                                 *
*  X(5), Y(5) = 0.4 0.8                                 *
*  X(6), Y(6) = 0.5 1                                   *
*  X(7), Y(7) = 0.6 0.8                                 *
*  X(8), Y(8) = 0.7 0.6                                 *
*  X(9), Y(9) = 0.8 0.4                                 *
*  X(10), Y(10) = 0.9 0.2                               *
*  X(11), Y(11) = 1 0                                   *
*                                                       *
* Polynomial approximation of degree 3 (11 points)      *
* Coefficients of polynomial:                           *
*  A(0) =    -0.069930070                               *
*  A(1) =     3.496503496                               *
*  A(2) =    -3.496503496                               *                      
*  A(3) =     0.000000000                               *
*                                                       *
* Approximated function:                                *
*        X           Y                                  *
*    0.000000   -0.069930                               *
*    0.100000    0.244755                               *
*    0.200000    0.489510                               *
*    0.300000    0.664336                               *
*    0.400000    0.769231                               *
*    0.500000    0.804196                               *
*    0.600000    0.769231                               *
*    0.700000    0.664336                               *
*    0.800000    0.489510                               *
*    0.900000    0.244755                               *
*    1.000000   -0.069930                               *
*                                                       *
*                    TPW version by J-P Moreau, Paris.  *
*                           (www.jpmoreau.fr)           *
********************************************************}
Program Approx;
Uses WinCrt;

Const  SIZE = 25;

Var    i,ij,j,k,n,n1,m,m1,m2: Integer;
       C: Array[0..SIZE,0..SIZE] of REAL;
       A,B,X,Xc,Y,Yx: Array[0..SIZE] of REAL;
       p,xx,s,yc: REAL;

       Function IntPower(x:REAL;k:Integer): REAL;
       begin
         if x=0 then IntPower:=0
                else IntPower:=Exp(k*Ln(x))
       end;

{main program}
BEGIN
  writeln;
  write(' Number of points    : '); readln(n);
  n:=n-1;
  write(' Degree of polynomial: '); readln(m);
  n1:=n+1; m1:=m+1; m2:=m+2;
  writeln;
  Writeln(' Function to approximate:');
  for i:=1 to n1 do
  begin
    write('  X(',i,'), Y(',i,') = '); readln(X[i],Y[i])
  end;
  for k:=1 to m2 do
  begin
    Xc[k]:=0;
    for i:=1 to n1 do Xc[k]:=Xc[k]+IntPower(X[i],k)
  end;
  yc:=0;
  for i:=1 to n1 do yc:=yc+Y[i];
  for k:=1 to m do
  begin
    Yx[k]:=0;
    for i:=1 to n1 do Yx[k]:=Yx[k]+Y[i]*IntPower(X[i],k)
  end;
  for i:=1 to m1 do
    for j:=1 to m1 do
    begin
      ij:=i+j-2;
      if (i=1) and (j=1) then C[1,1]:=n1
                         else C[i,j]:=Xc[ij]
    end;
  B[1]:=yc; for i:=2 to m1 do B[i]:=Yx[i-1];
  for k:=1 to m do
    for i:=k+1 to m1 do
    begin
      B[i]:=B[i]-C[i,k]/C[k,k]*B[k];
      for j:=k+1 to m1 do
        C[i,j]:=C[i,j]-C[i,k]/C[k,k]*C[k,j]
    end;
  A[m1]:=B[m1]/C[m1,m1];
  for i:=m downto 1 do
  begin
    s:=0;
    for k:=i+1 to m1 do s:=s+C[i,k]*A[k];
    A[i]:=(B[i]-s)/C[i,i]
  end;
  writeln;
  writeln(' Polynomial approximation of degree ',m,' (',n+1,' points)');
  writeln(' Coefficients of polynomial:');
  for i:=1 to m1 do writeln('  A(',i-1,') = ',A[i]:15:9);
  writeln;
  writeln(' Approximated function:');
  writeln(' ':8,'X',' ':11,'Y');
  for i:=1 to n1 do
  begin
    xx:=X[i]; p:=0;
    for k:=1 to m1 do p:=p*xx+A[m1+1-k];
    writeln(xx:12:6,p:12:6)
  end;
  writeln; Readkey;
  Donewincrt
END.

{end of file approx.pas}