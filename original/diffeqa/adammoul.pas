{*******************************************************************
*   Solve Y' = F(X,Y) with initial conditions using the Adams-     *
*   Moulton Prediction-Correction Method                           *
* ---------------------------------------------------------------- *
* SAMPLE RUN                                                       *
* (Integrate Y' = -Y + X/(1+X)^2 from X=0 to X=1 with initial      *
*  condition Y(0) = 1 )                                            *
*                                                                  *
*     X          Y         Y True   |Y-Y True|                     *
* ---------------------------------------------                    *
*  0.000000   1.000000   1.000000    0.000000                      *
*  0.050000   0.952381   0.952381    0.000000                      *
*  0.100000   0.909091   0.909091    0.000000                      *
*  0.150000   0.869569   0.869565    0.000004  1                   *
*  0.200000   0.833340   0.833333    0.000006  1                   *
*  0.250000   0.800009   0.800000    0.000009  1                   *
*  0.300000   0.769241   0.769231    0.000010  1                   *
*  0.350000   0.740752   0.740741    0.000011  1                   *
*  0.400000   0.714298   0.714286    0.000012  1                   *
*  0.450000   0.689668   0.689655    0.000012  1                   *
*  0.500000   0.666679   0.666667    0.000013  1                   *
*  0.550000   0.645174   0.645161    0.000013  1                   *
*  0.600000   0.625013   0.625000    0.000013  1                   *
*  0.650000   0.606073   0.606061    0.000013  1                   *
*  0.700000   0.588248   0.588235    0.000013  1                   *
*  0.750000   0.571441   0.571429    0.000013  1                   *
*  0.800000   0.555568   0.555556    0.000012  1                   *
*  0.850000   0.540553   0.540541    0.000012  1                   *
*  0.900000   0.526327   0.526316    0.000012  1                   *
*  0.950000   0.512832   0.512821    0.000011  1                   *
*  1.000000   0.500011   0.500000    0.000011  1                   *
* ---------------------------------------------------------------- *
* REFERENCE: "Méthode de calcul numérique- Tome 2 - Programmes en  *
*             Basic et en Pascal By Claude Nowakowski, Edition du  *
*             P.S.I., 1984".                                       *
*                                                                  *
*                              TPW Release By J-P Moreau, Paris.   *
*                                       (www.jpmoreau.fr)          *
********************************************************************
(See explanation file adambash.txt)
 ---------------------------------                      }
Program Adam_Moulton;

Uses WinCrt1;

Label 100, 200;

Var
    X, Y: Array[0..3] of Double;
    C1,C2,C3,C4, EC, ER, F, H, XX, YY, YC, YP: Double;
    K, L:Integer;

Procedure P1000(XX,YY:Double; Var C:Double);
Begin
  C := -YY + XX / (Sqr(1.0 + XX))
End;

Procedure P2000(XX:Double; Var F:Double);
Begin
  F := 1.0 / (1.0 + XX)
End;


{mainprogram}
BEGIN
  H := 0.05;                 {integration step}
  X[0] := 0.0; Y[0] := 1.0;  {Initial conditions}
  EC := 0.000001;            {Precision}

  ClrScr;
  Writeln;
  Writeln('     X          Y         Y True   |Y-Y True| ');
  Writeln(' ---------------------------------------------');
  Writeln(' ',X[0]:9:6, '  ',Y[0]:9:6, '  ',Y[0]:9:6, '   ',X[0]:9:6);

  {Start with Runge-Kutta}
  FOR K := 0 TO 1 do
  begin
    XX := X[K]; YY := Y[K]; P1000(XX,YY,C1);
    XX := X[K] + H / 2.0; YY := Y[K] + H / 2.0 * C1; P1000(XX,YY,C2);
    YY := Y[K] + H / 2.0 * C2; P1000(XX,YY,C3);
    X[K + 1] := X[K] + H;
    XX := X[K + 1]; YY := Y[K] + H * C3; P1000(XX,YY,C4);
    Y[K + 1] := Y[K] + H * (C1 + 2 * C2 + 2 * C3 + C4) / 6.0;
    XX := X[K + 1]; P2000(XX,F); ER := ABS(F - Y[K + 1]);

    Writeln(' ',X[K+1]:9:6, '  ',Y[K+1]:9:6, '  ',F:9:6, '   ',ER:9:6)

  end;

100: K := 2;
  XX := X[K]; YY := Y[K]; P1000(XX,YY,C1);
  XX := X[K - 1]; YY := Y[K - 1]; P1000(XX,YY,C2);
  XX := X[K - 2]; YY := Y[K - 2]; P1000(XX,YY,C3);
  X[K + 1] := X[K] + H;
  YP := Y[K] + H / 12.0 * (23.0 * C1 - 16.0 * C2 + 5.0 * C3);

  L := 0;

200: XX := X[K + 1]; YY := YP; P1000(XX,YY,C1);
  XX := X[K]; YY := Y[K]; P1000(XX,YY,C2);
  XX := X[K - 1]; YY := Y[K - 1]; P1000(XX,YY,C3);
  YC := Y[K] + H / 12.0 * (5.0 * C1 + 8.0 * C2 - C3);

  IF ABS(YP - YC) > EC THEN
  begin
    YP := YC; L := L + 1; GOTO 200
  end;

  Y[K + 1] := YC; XX := X[K + 1]; P2000(XX,F); ER := ABS(F - Y[K + 1]);

  Writeln(' ',X[K+1]:9:6, '  ',Y[K+1]:9:6, '  ',F:9:6, '   ',ER:9:6,'  ',L);

  FOR K := 0 TO 2 do
  begin
    X[K] := X[K + 1]; Y[K] := Y[K + 1]
  end;

  IF X[2] < 1.0 THEN GOTO 100;

  ReadKey;
  DoneWinCrt

END.

{end of file Adammoul.pas}