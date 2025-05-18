{***********************************************************************
*                                                                      *
*           Interpolate a function by trigonometric polynoms           *
*                                                                      *
* -------------------------------------------------------------------- *
*   SAMPLE RUN:                                                        *
*                                                                      *
*   Number of points: 3                                                *
*                                                                      *
*   Number of harmonics: 3                                             *
*                                                                      *
*                                                                      *
*   Function to interpolate:                                           *
*       X           Y                                                  *
*     0.000000   3.141593                                              *
*     0.897598   0.897598                                              *
*     1.795196   1.795196                                              *
*     2.692794   2.692794                                              *
*     3.590392   3.590392                                              *
*     4.487990   4.487990                                              *
*     5.385587   5.385587                                              *
*     6.283185   3.141593                                              *
*                                                                      *
*   Fourier coefficients:                                              *
*   A(0) =    3.141593  B(0) =    0.000000                             *
*   A(1) =   -0.000000  B(1) =   -1.863881                             *
*   A(2) =   -0.000000  B(2) =   -0.715810                             *
*   A(3) =   -0.000000  B(3) =   -0.204871                             *
*                                                                      *
*   Interpolation:                                                     *
*       X           Y                                                  *
*     0.000000   3.141593                                              *
*     0.448799   1.573507                                              *
*     0.897598   0.897598                                              *
*     1.346397   1.174039                                              *
*     1.795196   1.795196                                              *
*     2.243995   2.293325                                              *
*     2.692794   2.692794                                              *
*     3.141593   3.141593                                              *
*     3.590392   3.590392                                              *
*     4.039191   3.989860                                              *
*     4.487990   4.487990                                              *
*     4.936788   5.109147                                              *
*     5.385587   5.385587                                              *
*     5.834386   4.709678                                              *
*     6.283185   3.141593                                              *
*                                                                      *
* -------------------------------------------------------------------- *
*    "Méthodes de calcul numérique, Tome 2 By Claude Nowakowski,       *
*     PSI Edition, 1984" [BIBLI 04].                                   *
*                                                                      *
*                            Pascal Release 1.0 By J-P Moreau, Paris.  *
*                                      (www.jpmoreau.fr)               *
***********************************************************************}
Program TrigPoly;

Uses WinCrt;

Label 10, 90, 100, 200;

Type  Vec = Array[0..50] of Double;

Var
    IP, K, L, M, N, NH, NN: Integer;
    A, B, TF: Vec;
    C1, CF, CP, CT, F0, S1, SP, X: Double;
    Q, S, U0, U1, U2: Double;

BEGIN

    ClrScr;
10: Writeln;
    Write(' Number of points   : '); Readln(N);
    Writeln;
    Write(' Number of harmonics: '); Readln(NH);
    Writeln;

    IF N < NH THEN
    begin
      Writeln(' Error, N must be >= NH.');
      GOTO 10
    end;

    NN := 2 * N + 1;

{ TF stores function values }
    FOR K := 0 TO NN do
      TF[K] := K * 2.0 * PI / NN;

    TF[0] := (TF[0] + TF[NN]) / 2.0;
    TF[NN] := TF[0];

    CF := 2.0 / NN;

    Writeln;
    Writeln(' Function to interpolate:');
    Writeln('      X          Y');
    FOR K := 0 TO NN do
    begin
      X := K * 2.0 * PI / NN;
      Writeln(X:11:6, TF[K]:11:6)
    end;

    Writeln; ReadKey;

    CT := PI * CF;

{ Initialize }
    S1 := SIN(CT); C1 := COS(CT);
    CP := 1.0; SP := 0.0; IP := 0; F0 := TF[0];

{ Calculate Fourier coefficients A, B, for each harmonic }
90: U2 := 0.0; U1 := 0.0; M := NN - 1;
100:U0 := TF[M] + 2.0 * CP * U1 - U2;
    U2 := U1; U1 := U0; M := M - 1;
    IF M > 0 THEN GOTO 100;

    A[IP] := CF * (F0 + CP * U1 - U2);
    B[IP] := CF * SP * U1;

    IF IP > NH THEN GOTO 200;
    Q := C1 * CP - S1 * SP;
    SP := C1 * SP + S1 * CP;
    CP := Q;
    IP := IP + 1;
    GOTO 90;

200:A[0] := 0.5 * A[0];
    Writeln;
    Writeln(' Fourier coefficients:');
    FOR K := 0 TO NH do
      Writeln(' A(',K,') = ',A[K]:11:6,'  B(',K,') = ', B[K]:11:6);

{ Interpolate function }
    Writeln;
    Writeln(' Interpolation:');
    Writeln('      X          Y');
    FOR K := 0 TO NN * 2 do
    begin
      X := K * 2.0 * PI / (2 * NN);
      S := A[0];
      FOR L := 1 TO NH do
        S := S + A[L] * COS(X * L) + B[L] * SIN(X * L);
      Writeln(X:11:6, S:11:6)
    end;

    ReadKey;
    DoneWinCrt

END. {of main program

End of file trigpoly.pas}