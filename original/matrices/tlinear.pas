{********************************************************
*   SOLVE A LINEAR SYSTEM BY TRIANGULARIZATION METHOD   *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
*                                                       *
* Linear system AX = B                                  *
*                                                       *
* Matrix A:                                             *
*   2.00   1.00  -2.00                                  *
*   3.00  -1.00   1.00                                  *
*   7.00   5.00  -3.00                                  *
*                                                       *
* Right side:                                           *
*   B(1) =   1.00                                       *
*   B(2) =   0.00                                       *
*   B(3) =   0.00                                       *
*                                                       *
* Solution of AX=B:                                     *
*   X(1) =   0.0625                                     *
*   X(2) =  -0.5000                                     *
*   X(3) =  -0.6875                                     *
*                                                       *
* ----------------------------------------------------- *
* Reference: "Méthodes de calcul numérique - Tome 1 By  *
*             Claude Nowakowski, PS1 1981" [BIBLI 07].  *
*                                                       *
*                    TPW Release By J-P Moreau, Paris.  *
*                           (www.jpmoreau.fr)           *
********************************************************}
Program Tlinear;

Uses Wincrt1;

Const SIZE = 25;

Type
      MAT = Array[1..SIZE,1..SIZE] of Double;
      VEC = Array[1..SIZE] of Double;

Var
      A: MAT;
      B, X: VEC;

      I, J, K, N: Integer;
      S: Double;

BEGIN

  N := 3;

  A[1,1] := 2.0; A[1,2] :=  1.0; A[1,3] := -2.0;
  A[2,1] := 3.0; A[2,2] := -1.0; A[2,3] :=  1.0;
  A[3,1] := 7.0; A[3,2] :=  5.0; A[3,3] := -3.0;

  B[1] := 1.0;
  B[2] := 0.0;
  B[3] := 0.0;

  Clrscr;

  Writeln;
  Writeln(' Linear system AX = B');
  Writeln;
  Writeln(' Matrix A:');
  FOR I := 1 TO N do
  begin
    FOR J := 1 TO N do Write(' ', A[I, J]:6:2);
    Writeln
  end;

  Writeln;
  Writeln(' Right Side:');
  FOR I := 1 TO N do Writeln('   B(', I, ') = ', B[I]:6:2);

  {Transform A into triangular matrix}
  FOR K := 1 TO N - 1 do
  begin
    FOR I := K + 1 TO N do
    begin
      B[I] := B[I] - A[I, K] / A[K, K] * B[K];
      FOR J := K + 1 TO N do
        A[I, J] := A[I, J] - A[I, K] / A[K, K] * A[K, J]
    end
  end;

  {Solve triangular system}
  X[N] := B[N] / A[N, N];
  FOR I := N - 1 DOWNTO 1 DO
  begin
    S := 0.0;
    FOR K := I + 1 TO N DO  S := S + A[I, K] * X[K];
    X[I] := (B[I] - S) / A[I, I]
  end;

  {Print results: }
  Writeln;
  Writeln(' Solution of AX=B:');
  FOR I := 1 TO N DO Writeln('   X(', I, ') = ', X[I]:8:4);

  ReadKey;
  DoneWinCrt

END.

{DATA 3
 DATA 2,1,-2,3,-1,1,7,5,-3
 DATA 1,0,0

 end of file Tlinear.pas }