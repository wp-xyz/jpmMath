{************************************************************
*    Interpolate a function F(x) by continuous fractions    *
* --------------------------------------------------------- *
* SAMPLE RUN:                                               *
* (Interpolate function e(x) between x=0 and x=2)           *
*                                                           *
* Number of points: 3                                       *
* X(0), Y(0): 0 1                                           *
* X(1), Y(1): 1 2.71828                                     *
* X(2), Y(2): 2 7.38906                                     *
*                                                           *
* Coefficients D(K):                                        *
* D(0) =      1.000000                                      *
* D(1) =      0.581977                                      *
* D(2) =     -3.718271                                      *
*                                                           *
* X = 1.5                                                   *
*                                                           *
* For X = 1.5    Y =    4.351909                            *
*                                                           *
* --------------------------------------------------------- *
* Ref.: "Méthodes de calcul numérique, Tome 2 By Claude     *
*        Nowakowski, PSI Edition, 1984" [BIBLI 04].         *
*                                                           *
*                      Pascal Release By J-P Moreau, Paris. *
*                              (www.jpmoreau.fr)            *
************************************************************}
Program Confract;

Uses WinCrt;

Const SIZE = 50;

Var
    K,L,M,N,N1: Integer;
    X,Y,D: Array[0..SIZE] of double;
    DD, DL, S, XX: Double;

BEGIN

  Writeln;
  Write(' Number of points: '); Readln(N1);

  N := N1 - 1;
  M := N;

{Read data from screen}
  FOR K := 0 TO N do
  begin
    Write(' X(',K,'), Y(',K,'): '); Readln(X[K], Y[K])
  end;

{Calculate coefficients D(K) }
  FOR K := 0 TO N do D[K] := Y[K];

  FOR L := 1 TO M do
  begin
    FOR K := L TO N do
    begin
      DD := (X[K] - X[L - 1]) / (D[K] - D[L - 1]);
      IF K <> L THEN
        D[K] := DD
      ELSE
        DL := DD
    end;
    D[L] := DL
  end;

{print coefficients}
  Writeln;
  Writeln(' Coefficients D(K):');
  FOR K := 0 TO N do
    Writeln(' D(',K,') = ', D[K]:10:6);

{Interpolate for X=XX }
  Writeln;
  Write(' X= '); Readln(XX);

{Evaluate continuous fraction}
  S := (XX - X[N - 1]) / D[N];
  FOR K := N - 1 DOWNTO 1 DO
    S := (XX - X[K - 1]) / (D[K] + S);

  S := S + D[0];

  Writeln(' For X= ', XX:10:6,' Y= ', S:10:6);

  ReadKey;
  DoneWinCrt

END.

{end of file confract.pas}