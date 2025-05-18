{************************************************************
*    Integration of a discrete function by the weighting    *
*    coefficients Method                                    *
* --------------------------------------------------------- *
* SAMPLE RUN:                                               *
*                                                           *
* Number of points: 7                                       *
*                                                           *
* Input the 7 points, example: 1.5  0.666666                *
*                                                           *
* 0  ?  1 1                                                 *
* 1  ?  1.5 0.666666                                        *
* 2  ?  2 0.5                                               *
* 3  ?  2.5 0.4                                             *
* 4  ?  3 0.333333                                          *
* 5  ?  3.5 0.285714                                        *
* 6  ?  4 0.25                                              *
*                                                           *
* Coefficients A(I):                                        *
* A( 0) =  1.46428571428525E-0001                           *
* A( 1) =  7.71428571428824E-0001                           * 
* A( 2) =  9.64285714279978E-0002                           *
* A( 3) =  9.71428571429267E-0001                           * 
* A( 4) =  9.64285714280966E-0002                           *
* A( 5) =  7.71428571428745E-0001                           *
* A( 6) =  1.46428571428545E-0001                           *
*                                                           *
* Integral of F(X) from  1.00 to  4.00:                     *
*   I =     1.386657                                        *
*                                                           *
* --------------------------------------------------------- *
* Reference: "Méthodes de calcul numérique - Tome 1         *
*             By Claude Nowakowski, PS1 1981" [BIBLI07].    *
*                                                           *
*                       TPW Release By J-P Moreau, Paris.   *
*                              (www.jpmoreau.fr)            *
*************************************************************
See explanations in file dinteg.txt
-----------------------------------  }
Program weighting_Coeffs;

Uses WinCrt1;

Const SIZE = 25;

Type  MAT = Array[0..SIZE,0..SIZE] of Double;
      VEC = Array[0..SIZE] of Double;

Var   X,Y,A,B: VEC;
      T: MAT;
      I,J,K,N: Integer;
      AA,BB,S: Double;

BEGIN

  Writeln;
  Write(' Number of points: '); Readln(N);

  Writeln;
  Writeln(' Input the ', N, ' points, example: 1.5  0.666666');
  Writeln;

  N := N - 1;

  FOR K := 0 TO N do
  begin
    Write(' ', K,'  ?  '); Readln(X[K], Y[K])
  end;

{ Calculate the coefficients of the linear system }
  FOR J := 0 TO N do
  begin
    T[0, J] := 1.0;
    FOR I := 1 TO N do
      T[I, J] := T[I - 1, J] * X[J]
  end;

{ right side coefficients }
  AA := 1.0; BB := 1.0;
  FOR I := 0 TO N do
  begin
    AA := AA * X[0]; BB := BB * X[N];
    B[I] := (BB - AA) / (1.0*(I + 1))
  end;

{ Calculate coefficients A(I) and triangularize }
  FOR K := 0 TO N - 1 do
  begin
    FOR I := K + 1 TO N do
    begin
      B[I] := B[I] - T[I, K] / T[K, K] * B[K];
      FOR J := K + 1 TO N do
        T[I, J] := T[I, J] - T[I, K] / T[K, K] * T[K, J]
    end
  end;

{ Solve triangular linear system T * A = B }
  A[N] := B[N] / T[N, N];
  FOR I := N - 1 DownTo 0 do
  begin
    S := 0.0;
    FOR K := I + 1 TO N do
      S := S + T[I, K] * A[K];
    A[I] := (B[I] - S) / T[I, I]
  end;

  Writeln;
  Writeln(' Coefficients A(I):');
  FOR I := 0 TO N do
    Writeln(' A(', I:2, ') = ', A[I]);

{ Calculate Integral of discrete F(X) }
  S := 0.0;
  FOR I := 0 TO N do S := S + A[I] * Y[I];

{ Print final result }
  Writeln;
  Writeln(' Integral of F(X) from ', X[0]:5:2, ' to ', X[N]:5:2, ':');
  Writeln('   I = ', S:12:6);
  Writeln;

  ReadKey;
  DoneWinCrt

END.