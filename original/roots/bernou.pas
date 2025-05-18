{***********************************************
* Roots of a polynomial by Bernouilli's method *
* -------------------------------------------- *
* Reference: "Methodes de calcul numerique By  *
*             Claude Nowakowski, PSI Editions, *
*             France, 1981" [BIBLI 04].        *
*                                              *
*             Pascal version by J-P Moreau     *
*                  (www.jpmoreau.fr)           *
* -------------------------------------------- *
* SAMPLE RUN:                                  *
*                                              *
* Roots of a polynomial by Bernouilli's method *
*                                              *
* Input order of polynomial: 4                 *
*                                              *
* Input coefficients of the polynomial:        *
*                                              *
*   A( 0 ) = 1                                 *
*   A( 1 ) = 2                                 *
*   A( 2 ) = -13                               *
*   A( 3 ) = -14                               *
*   A( 4 ) = 24                                *
*                                              *
* ===== R O O T S =====                        *
*                                              *
*   X1 =   -4.0000                             *
*   X2 =    3.0000                             *
*   X3 =   -2.0001                             *
*   X4 =    1.0001                             *
*                                              *
***********************************************}
Program Test_Bernouilli;
Uses WinCrt;

Var
        k,nd : integer;
        A,Y  : Array[0..25] of double;


{ Bernouilli's subroutine }
Procedure Bernouilli;
Label 5,50,100;
Var k,l,n,nroot : integer;
    er,s,x,x1 : double;
Begin
  writeln;
  writeln(' ===== R O O T S =====');
  writeln;
  n:=nd;
  if ABS(A[0])<1e-12 then
  begin
    writeln(' Error: A(0) must be <> 0 !');
    exit
  end;
  for k:=0 to n do A[k]:=A[k]/A[0];
  nroot:=0;
5: x1:=0.0; nroot:=nroot+1;
  for k:=2 to n do Y[k]:=0.0;
  Y[1]:=1.0;
  for l:=1 to 50 do
  begin
    s:=0.0;
    for k:=1 to n do s:=s+Y[k]*A[k];
    Y[0]:=-s;
    x:=Y[0]/Y[1];
    er:=ABS(x1-x);
    if l<n then goto 50;
    if er<5e-5 then goto 100;
50: x1:=x;
    for k:=n downto 1 do Y[k]:=Y[k-1]
  end; {l loop}
  writeln(' Error: no convergence (error=',er,')');
  exit;
100: writeln('   X',nroot,' = ',x:9:4);
  n:=n-1;
  if n=0 then exit;
  A[1]:=A[1]+x;
  for k:=2 to n do A[k]:=A[k]+x*A[k-1];
  goto 5
End;


{main program}
BEGIN
  clrscr;
  writeln;
  writeln(' Roots of a polynomial by Bernouilli''s method');
  writeln;
  write(' Input order of polynomial: '); read(nd);
  writeln;

  writeln(' Input coefficients of the polynomial:');
  writeln;
  for k:=0 to nd do
  begin
    write('   A(',k,') = '); read(A[k])
  end;

  {call Bernouilli procedure and writeln results}
  Bernouilli;

  Readkey; DoneWinCrt

END.

{end of file bernou.pas}