{********************************************************************
* This program solves a Vandermonde linear system:                  *
*               N      k-1                                          *
*              S     xi    wi = qk  (k=1,...,N)                     *
*               i=1                                                 *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
*   N=10                                                            *
*   X(I)=I (I=1,10)                                                 *
*   Q(I)=1, except Q(9)=30                                          *
*                                                                   *
* The solution vector is:                                           *
*    1       1.004315                                               *
*    2      -0.038120                                               *
*    3       0.149603                                               *
*    4      -0.342631                                               *
*    5       0.503472                                               *
*    6      -0.493403                                               *
*    7       0.322222                                               *
*    8      -0.135218                                               *
*    9       0.033085                                               *
*   10      -0.003596                                               *
*                                                                   *
* ----------------------------------------------------------------- *
* Reference: "Numerical Recipes by W.H. Press, B.P. Flannery, S.A.  *
*             Teukolsky, W.T. Vetterling, Cambridge University      *
*             Press, 1987"                                          *
*                                                                   *
*                             Pascal Version By J-P Moreau, Paris   *
*                                    (www.jpmoreau.fr)              *
********************************************************************}
Program Test_Vander;
Type
    VEC = Array[1..10] of double;
Var
    i,N: integer;
    X, Q, W: VEC;


Procedure Vander(X:VEC; Var W:VEC; Q:VEC; N:Integer);
{---------------------------------------------------------------------
!                                        N     k-1
! Solves the Vandermonde linear system S     xi    wi = qk (k=1,...,N)
!                                        i=1
! Input consists of the vectors X and Q, each of length N, the vector
! W is output.
!---------------------------------------------------------------------}
Const NMAX=100; ZERO=0.0; ONE=1.0;
{ NMAX is the maximum expected value of N.}
Var B, S, T, XX: Double;
    i,j,k,k1: Integer;
    C: VEC;

Begin
if N=1 then
  W[1]:=Q[1]
else
begin
  For i:=1 to N do   {initialize array C }
    C[I]:=ZERO;
  C[N]:=-X[1];       {coefficients of the master polynomial}
                     {are found by recursion. }
  For i:=2 to N do
  begin
    XX:=-X[i];
    For j:=N+1-i to N-1 do
      C[j]:=C[j]+XX*C[j+1];
   C[N]:=C[N]+XX;
  end;
  For i:=1 to N do   {!each subfactor in turn}
  begin
    XX:=X[i];
    T:=ONE;
    B:=ONE;
    S:=Q[N];
    k:=N;
    For j:=2 to N do   {is synthetically divided, }
    begin
      k1:=k-1;
      B:=C[k]+XX*B;
      S:=S+Q[k1]*B;    {matrix-multiplied by the right-hand side, }
      T:=XX*T+B;
      k:=k1
    end;
    W[i]:=S/T          {and suppliedwith a denominator.  }
  end
end

End;


{main program}
BEGIN


  N:=10;

  For i:=1 to N do X[i]:=1.0*i;
  For i:=1 to N do Q[i]:=1.0;
  Q[9]:=30.0;

  Vander(X,W,Q,N);

  writeln;
  writeln(' The solution vector is:');
  for i:=1 to N do
    writeln(' ',i:4, '  ',W[i]:13:6);
  Readln;

END.

{end of file tvander.pas}