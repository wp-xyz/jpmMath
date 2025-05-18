{********************************************************************
* This program solves a Toeplitz linear system:                     *
*                N                                                  *
*              S     R      x  = y   (i=1,...,N)                    *
*               j=1   N+i-j  j    i                                 *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
*   N=10                                                            *
*   R = (-10,-12,3,4,5,6,7,8,9,5,11,12,13,14,15,16,17,18,19)        *
*   Y = (1,1,1,1,1,1,1,1,30,1)                                      *
*                                                                   *
* The solution vector is:                                           *
*  1   4.145258                                                     *
*  2  -7.202165                                                     *
*  3   0.136163                                                     *
*  4   0.473885                                                     *
*  5   0.811607                                                     *
*  6   1.149329                                                     *
*  7   1.487051                                                     *
*  8   1.824773                                                     *
*  9  -3.637505                                                     *
*  10   2.500217                                                    *
*                                                                   *
* ----------------------------------------------------------------- *
* Reference: "Numerical Recipes by W.H. Press, B.P. Flannery, S.A.  *
*             Teukolsky, W.T. Vetterling, Cambridge University      *
*             Press, 1987"                                          *
*                                                                   *
*                             Pascal Version By J-P Moreau, Paris   *
*                                     (www.jpmoreau.fr)             *
********************************************************************}
Program Test_Toepltz;

  Type

  VEC1 = Array[1..19] of double;
  VEC  = Array[1..10] of double;


  Const

  R:VEC1 = (-10.0, -12.0,  3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 5.0,
             11.0,  12.0, 13.0,14.0,15.0,16.0,17.0,18.0,19.0);
  Y:VEC  = (1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,30.0, 1.0);

  Var

  i, N: Integer;

  X: VEC;



Procedure TOEPLZ(R:VEC1; Var X:VEC; Y:VEC; N:Integer);
 {--------------------------------------------------------------------
 * This program solves a Toeplitz linear system:                     *
 *                N                                                  *
 *              S     R      x  = y   (i=1,...,N)                    *
 *               j=1   N+i-j  j    i                                 *
 *                                                                   *
 *  Input consists of vectors R of size 2N-1 and Y of size N.        *
 *  Output is solution vector X of size N.                           *
 *  N is the size of the linear system.                              *
 --------------------------------------------------------------------}
 Label 99, 100;
 Const
  NMAX=100; ZERO=0.0;
  {NMAX is the maximum expected value of N. }
Type
  VEC2 = Array[1..NMAX] of double;
Var
  G, H: VEC2;
  PP, PT1, PT2, QQ, QT1, QT2, SD, SGD, SGN, SHN, SXN: Double;
  J, K, M, M1, M2: Integer;

Begin
if R[N] = ZERO then goto 99;
X[1]:=Y[1]/R[N];             {initialize for the recursion. }
if N = 1 then goto 100;
G[1]:=R[N-1]/R[N];
H[1]:=R[N+1]/R[N];

for M:=1 to N do             {main loop over the recursion. }
begin
  M1:=M+1;
  SXN:=-Y[M1];
  SD:=-R[N];
  for J:=1 to M do
  begin
    SXN:=SXN+R[N+M1-J]*X[J];
    SD:=SD + R[N+M1-J]*G[M-J+1]
  end;
  if SD = ZERO then goto 99;
  X[M1]:=SXN/SD;             {whence x. }
  for J:=1 to M do
    X[J] := X[J] - X[M1]*G[M-J+1];
  if M1 = N then goto 100;   {normal exit}
  SGN:=-R[N-M1];             {compute numerator and denominator for G and H.}
  SHN:=-R[N+M1];
  SGD:=-R[N];
  for J:=1 to M do
  begin
    SGN := SGN + R[N+J-M1]*G[J];
    SHN := SHN + R[N-J+M1]*H[J];
    SGD := SGD + R[N+J-M1]*H[M-J+1]
  end;

  if (SD = ZERO) OR (SGD = ZERO) then goto 99;
  G[M1]:=SGN/SGD;             {whence G and H }
  H[M1]:=SHN/SD;
  K:=M;
  M2:=(M+1) Div 2;
  PP:=G[M1];
  QQ:=H[M1];
  for J:=1 to M2 do
  begin
    PT1:=G[J];
    PT2:=G[K];
    QT1:=H[J];
    QT2:=H[K];
    G[J]:=PT1-PP*QT2;
    G[K]:=PT2-PP*QT1;
    H[J]:=QT1-QQ*PT2;
    H[K]:=QT2-QQ*PT1;
    Dec(K)
  end
end; {of M loop}

99:Writeln(' Levinson method fails: singular principal minor.');

100: End; {TOEPLZ}


{main program}
BEGIN

  N := 10;

  TOEPLZ(R,X,Y,N);

  Writeln;
  Writeln(' The solution vector is:');
  for i:=1 to N do
    Writeln('  ',i,'  ',X[i]:9:6);
  Readln

END.

{end of file toeplitz.pas}