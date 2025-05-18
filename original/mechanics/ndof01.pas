{**********************************************************************************
*  Response of a N d.o.f. Mass-Spring System with Damping to a sinusoidal Force   *
*                         By Transfer Matrices Method                             *
* ------------------------------------------------------------------------------- *
* EXPLANATION:                                                                    *
* Find response of Mass 2m in 3 d.o.f. undamped system below:                     *
*        Nodes: 1          2    3          4    5          6    7                 * 
*                                                                                 *
*                          --> x1          --> x2          --> x3                 *
*              \|    k     ------    2k    ------    k     ------                 *
*              \|--/\/\/\--| 2m |--/\/\/\--| m  |--/\/\/\--| 3m |                 *
*              \|          ------          ------          ------                 *
*                          --> F1(t)                                              *
*                                                                                 *
* The differential equations of the system motions is:                            *
*               | 2m  0  0 | | x1" |   | 3k  -2k  0 | | x1 |   | 0 |              *
*               |  0  m  0 | | x2" | + |-2k   3k -k | | x2 | = | 0 |              *
*               |  0  0 3m | | x3" |   |  0   -k  k | | x3 |   | 0 |              *
*                                                                                 *
*             | X1 |  jwt                                                         *
* Noting  X = | X2 | e      the matrix form of equations becomes:                 *
*             | X3 |                                                              *
*                                                                                 *
*             | 3k - 2mw2    -2k        0     | | X1 |                            *
*             |   -2k      3k - mw2    -k     | | X2 | =  0                       *
*             |     0         -k     k - 3mw2 | | X3 |                            *
*                                                                                 *
* The matrix determinant is null for w1 = 0.10517 sqr(k/m), w2 = 0.80861 sqr(k/m) *
* and w3 = 3.9195 sqr(k/m), hence the resonance fresuencies are:                  *
*             F1 = w1/2pi = 0.05161 sqr(k/m)                                      *
*             F2 = w2/2pi = 0.1431 sqr(k/m)                                       *
*             F3 = w3/2pi = 0.3151 sqr(k/m)                                       *
* with the corresponding modal vectors:                                           *
*                    |   1   |         |   1    |         |    1    |             *
*             phi1 = | 1.395 |  phi2 = | 0.6914 |  phi3 = | -2.420  |             *
*                    | 2.038 |         |-0.4849 |         |  0.2249 |             *
*                                                                                 *
* ------------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                     *
*                                                                                 *
*  Kind of elements                                                               *
*  ----------------                                                               *
*   1: Spring                                                                     *
*   2: Mass                                                                       *
*   3: Viscous Damper                                                             *
*   4: Spring + Viscous Damper in parallel                                        *
*   5: Spring with structural Damping                                             *
*   6: Sinusoidal Force                                                           *
*                                                                                 *
*  Number of elements: 7                                                          *
*                                                                                 *
*  Kind of element 1: 1                                                           *
*  Kind of element 2: 2                                                           *
*  Kind of element 3: 1                                                           *
*  Kind of element 4: 2                                                           *
*  Kind of element 5: 1                                                           *
*  Kind of element 6: 2                                                           *
*  Kind of element 7: 6                                                           *
*                                                                                 *
*  Mass #1 = 2                                                                    *
*  Mass #2 = 1                                                                    *
*  Mass #3 = 3                                                                    *
*                                                                                 *
*  Spring #1 = 1                                                                  *
*  Spring #2 = 2                                                                  *
*  Spring #3 = 1                                                                  *
*                                                                                 *
*  Excitation Force #1 F.COS(PHI) = 1000                                          *
*  Excitation Force #1 F.SIN(PHI) = 0                                             *
*                                                                                 *
*  For which node number do you want the response: 3                              *
*                                                                                 *
*  Fixed-Fixed System, code=1.                                                    *
*  Fixed-Free System,  code=2.                                                    *
*  Free-Fixed System,  code=3.                                                    *
*  Free-Free System,   code=4.                                                    *
*                                                                                 *
*  Code = 2                                                                       *
*                                                                                 *
*  Frequency Sweep                                                                *
*  ---------------                                                                *
*  Starting Frequency: 0.30                                                       *
*  Ending Frequency..: 0.32                                                       *
*  Frequency Step....: 0.001                                                      *
*                                                                                 *
*  Freq=  0.300  Displacement=     96.1  Phase=  0.0 Deg.                         *
*  Freq=  0.301  Displacement=    101.2  Phase=  0.0 Deg.                         *
*  Freq=  0.302  Displacement=    107.1  Phase=  0.0 Deg.                         *
*  Freq=  0.303  Displacement=    114.0  Phase=  0.0 Deg.                         *
*  Freq=  0.304  Displacement=    122.2  Phase=  0.0 Deg.                         *
*  Freq=  0.305  Displacement=    132.1  Phase=  0.0 Deg.                         *
*  Freq=  0.306  Displacement=    144.1  Phase=  0.0 Deg.                         *
*  Freq=  0.307  Displacement=    159.3  Phase=  0.0 Deg.                         *
*  Freq=  0.308  Displacement=    178.8  Phase=  0.0 Deg.                         *
*  Freq=  0.309  Displacement=    204.7  Phase=  0.0 Deg.                         *
*  Freq=  0.310  Displacement=    240.9  Phase=  0.0 Deg.                         *
*  Freq=  0.311  Displacement=    294.8  Phase=  0.0 Deg.                         *
*  Freq=  0.312  Displacement=    383.9  Phase=  0.0 Deg.                         *
*  Freq=  0.313  Displacement=    558.2  Phase=  0.0 Deg.                         *
*  Freq=  0.314  Displacement=   1051.8  Phase=  0.0 Deg.                         *
*  Freq=  0.315  Displacement=  12214.5  Phase=  0.0 Deg. <-- 3rd resonance       *
*  Freq=  0.316  Displacement=   1226.3  Phase=180.0 Deg.                         *
*  Freq=  0.317  Displacement=    574.1  Phase=180.0 Deg.                         *
*  Freq=  0.318  Displacement=    370.7  Phase=180.0 Deg.                         *
*  Freq=  0.319  Displacement=    271.5  Phase=180.0 Deg.                         *
*  Freq=  0.320  Displacement=    212.8  Phase=180.0 Deg.                         *
*                                                                                 *
* Example #2:                                                                     *
*        Nodes: 1          2    3          4    5          6    7                 *
*                          --> x1          --> x2    k     --> x3                 *
*              \|    k     *----*    2k    *----*--/\/\/\--*----*                 *
*              \|--/\/\/\--| 2m |--/\/\/\--| m  |  ---*    | 3m |                 *
*              \|          *----*          *----*---| |----*----*                 *
*                          --> F1(t)               ---*                           *
*                                                   E                             *
*                                                                                 *
*  Number of elements: 7                                                          *
*                                                                                 *
*  Kind of element 1: 1                                                           *
*  Kind of element 2: 2                                                           *
*  Kind of element 3: 1                                                           *
*  Kind of element 4: 2                                                           *
*  Kind of element 5: 5                                                           *
*  Kind of element 6: 2                                                           *
*  Kind of element 7: 6                                                           *
*                                                                                 *
*  Mass #1 = 2                                                                    *
*  Mass #2 = 1                                                                    *
*  Mass #3 = 3                                                                    *
*                                                                                 *
*  Spring #1 = 1                                                                  *
*  Spring #2 = 2                                                                  *
*                                                                                 *
*  Spring with structural Damping #1  K = 1                                       *
*  Spring with structural Damping #1  E = 0.05                                    *
*                                                                                 *
*  Excitation Force #1 F.COS(PHI) = 1000                                          *
*  Excitation Force #1 F.SIN(PHI) = 0                                             *
*                                                                                 *
*  For which node number do you want the response: 3                              *
*                                                                                 *
*  Fixed-Fixed System, code=1.                                                    *
*  Fixed-Free System,  code=2.                                                    *
*  Free-Fixed System,  code=3.                                                    *
*  Free-Free System,   code=4.                                                    *
*                                                                                 *
*  Code = 2                                                                       *
*                                                                                 *
*  Frequency Sweep                                                                *
*  ---------------                                                                *
*  Starting Frequency: 0.30                                                       *
*  Ending Frequency..: 0.32                                                       *
*  Frequency Step....: 0.001                                                      *
*                                                                                 *
*  Freq=  0.300  Displacement=     95.7  Phase=  3.5 Deg.                         *
*  Freq=  0.301  Displacement=    100.6  Phase=  3.9 Deg.                         *
*  Freq=  0.302  Displacement=    106.3  Phase=  4.5 Deg.                         *
*  Freq=  0.303  Displacement=    113.0  Phase=  5.1 Deg.                         *
*  Freq=  0.304  Displacement=    120.9  Phase=  5.8 Deg.                         *
*  Freq=  0.305  Displacement=    130.4  Phase=  6.7 Deg.                         *
*  Freq=  0.306  Displacement=    141.8  Phase=  7.8 Deg.                         *
*  Freq=  0.307  Displacement=    156.0  Phase=  9.1 Deg.                         *
*  Freq=  0.308  Displacement=    173.9  Phase= 10.8 Deg.                         *
*  Freq=  0.309  Displacement=    197.2  Phase= 13.0 Deg.                         *
*  Freq=  0.310  Displacement=    228.3  Phase= 15.9 Deg.                         *
*  Freq=  0.311  Displacement=    271.8  Phase= 20.1 Deg.                         *
*  Freq=  0.312  Displacement=    334.9  Phase= 26.5 Deg.                         *
*  Freq=  0.313  Displacement=    429.1  Phase= 37.0 Deg.                         *
*  Freq=  0.314  Displacement=    557.6  Phase= 55.2 Deg.                         *
*  Freq=  0.315  Displacement=    644.2  Phase= 84.1 Deg. <-- 3rd resonance       *
*  Freq=  0.316  Displacement=    562.8  Phase=-65.6 Deg.                         *
*  Freq=  0.317  Displacement=    422.0  Phase=-45.6 Deg.                         *
*  Freq=  0.318  Displacement=    317.1  Phase=-34.2 Deg.                         *
*  Freq=  0.319  Displacement=    247.5  Phase=-27.3 Deg.                         *
*  Freq=  0.320  Displacement=    200.3  Phase=-22.8 Deg.                         *
*                                                                                 *
* Note: with a moderate structural damping added to spring #3 (eta=5%), the 3rd   *
*       resonance frequency is practically unchanged, but the response is much    *
*       lower and the phases are completely different.                            *
* ------------------------------------------------------------------------------- *
* REFERENCE: "Mécanique des vibrations linéaires By M. Lalanne, P. Berthier,      *
*             J. Der Hagopian, Masson, Paris 1980" [BIBLI 16].                    *
*                                                                                 *
*                                           Pascal Release By J-P Moreau, Paris.  *
*                                                    (www.jpmoreau.fr)            *
**********************************************************************************}
Program Ndof01;

Uses WinCrt;

Label 200,680,770,840,910,1000,1090,1180,1200,1300;

Const
      NMAX = 25;
      TINY = 1E-12;

Type
      MAT5 = Array[1..5,1..5] of Double;
      VEC  = Array[1..NMAX] of Double;
      VEC5 = Array[1..5] of Double;

Var

    N,N1, I,I1, J,J1, K, I9, M9, K9, C3, C6, C8, C9, V9, F9: Integer;
    C4,D,F,F3,F4,F5, M2, S4, Sum, T3,W,W2: Double;

    C: Array[1..NMAX] of Integer;
    T,T1,T2,A: MAT5;
    V1,V2: VEC5;
    B: Array[1..2,1..2] of Double;
    E1,E2: Array[1..2] of Double;
    C1,C2,E3,K1,K2,K3,M1: VEC;
    F1,F2: VEC;

    Ans: String[3];


Procedure S2450;  {Mass Matrix}
Begin
  M9 := M9 + 1;
  A[1, 2] := -M1[M9] * W2;
  A[3, 4] := A[1, 2]
End;

Procedure S2500;  {Stiffness Matrix}
Begin
  K9 := K9 + 1;
  A[2, 1] := 1.0 / K1[K9];
  A[4, 3] := A[2, 1]
End;

Procedure S2550;  {Viscous damping Matrix}
Begin
  C9 := C9 + 1;
  A[4, 1] := -1.0 / C1[C9] / W;
  A[2, 3] := -A[4, 1]
End;

Procedure S2600;  {Spring + viscous damper in parallel Matrix}
Begin
  V9 := V9 + 1;
  A[2, 1] := K2[V9] / (Sqr(K2[V9]) + Sqr(C2[V9]) * W2);
  A[4, 3] := A[2, 1];
  A[4, 1] := -C2[V9] * W / (Sqr(K2[V9]) + Sqr(C2[V9]) * W2);
  A[2, 3] := -A[4, 1]
End;

Procedure S2650;  {Spring with structural damping Matrix}
Begin
  I9 := I9 + 1;
  A[2, 1] := 1.0 / (K3[I9] * (1.0 + Sqr(E3[I9])));
  A[4, 3] := A[2, 1];
  A[4, 1] := -E3[I9] / (K3[I9] * (1.0 + Sqr(E3[I9])));
  A[2, 3] := -A[4, 1]
End;

Procedure S2700;  {Force Matrix}
Begin
  F9 := F9 + 1;
  A[1, 5] := -F1[F9];
  A[3, 5] := -F2[F9]
End;


{main program}
BEGIN

M9 := 0; K9 := 0; C9 := 0; V9 := 0; I9 := 0; F9 := 0;

200: Writeln;
Writeln(' Kind of elements');
Writeln(' ----------------');
Writeln('  1; Spring');
Writeln('  2; Mass');
Writeln('  3; Viscous Damper');
Writeln('  4; Spring + Viscous Damper in parallel');
Writeln('  5; Spring with structural Damping');
Writeln('  6; Sinusoidal Force');
Writeln;
Write(' Number of elements: '); Readln(N);

Writeln;
FOR I := 1 TO N DO
begin
  Write(' Kind of element ', I, ': '); Readln(C[I]);
  IF C[I] = 1 THEN
    K9 := K9 + 1
  ELSE IF C[I] = 2 THEN
    M9 := M9 + 1
  ELSE IF C[I] = 3 THEN
    C9 := C9 + 1
  ELSE IF C[I] = 4 THEN
    V9 := V9 + 1
  ELSE IF C[I] = 5 THEN
    I9 := I9 + 1
  ELSE IF C[I] = 6 THEN
    F9 := F9 + 1
  ELSE
  begin
    Writeln(' Unknown Element.');
    Halt
  end
end;

IF F9 <> 0 THEN GOTO 680;
Writeln;
Writeln(' There is no excitation Force, Redo!');
GOTO 200;

680: Writeln;
IF M9 = 0 THEN GOTO 770;
{input M9 masses}
FOR I := 1 TO M9 DO
begin
  Write(' Mass #', I,' = '); Readln(M1[I])
end;
Writeln;

770: IF K9 = 0 THEN GOTO 840;
{input K9 springs}
FOR I := 1 TO K9 DO
begin
  Write(' Spring #', I, ' = '); Readln(K1[I])
end;
Writeln;

840: IF C9 = 0 THEN GOTO 910;
{input C9 viscous Dampers}
FOR I := 1 TO C9 DO
begin
  Write(' Viscous Damper #', I, ' = '); Readln(C1[I])
end;
Writeln;

910: IF V9 = 0 THEN GOTO 1000;
{input V9 springs + viscous Dampers in parallel}
FOR I := 1 TO V9 DO
begin
  Write(' Spring + Viscous Damper #', I, '  K = '); Readln(K2[I]);
  Write(' Spring + Viscous Damper #', I, '  C = '); Readln(C2[I])
end;
Writeln;

1000: IF I9 = 0 THEN GOTO 1090;
{input I9 springs with structural Damping}
FOR I := 1 TO I9 DO
begin
  Write(' Spring with structural Damping #', I, '  K = '); Readln(K3[I]);
  Write(' Spring with structural Damping #', I, '  E = '); Readln(E3[I])
end;
Writeln;

1090: {input F9 sinusoidal forces}
FOR I := 1 TO F9 DO
begin
  Write(' Excitation Force #', I, '  F.COS(PHI) = '); Readln(F1[I]);
  Write(' Excitation Force #', I, '  F.SIN(PHI) = '); Readln(F2[I])
end;

1180: Writeln;
Write(' For which node number do you want the response: '); Readln(N1);

1200: Writeln;
Writeln(' Fixed-Fixed System, code:=1.');
Writeln(' Fixed-Free System,  code:=2.');
Writeln(' Free-Fixed System,  code:=3.');
Writeln(' Free-Free System,   code:=4.');
Writeln;
Write(' Code = '); Readln(C8);
1300: Writeln;
Writeln(' Frequency Sweep');
Writeln(' ---------------');
Write(' Starting Frequency; '); Readln(F3);
Write(' Ending Frequency..; '); Readln(F4);
Write(' Frequency Step....; '); Readln(F5);
Writeln;
{end of data section}

{Main frequency loop}
F:=F3;
Repeat
  M9 := 0; K9 := 0; C9 := 0; V9 := 0; I9 := 0; F9 := 0; C6 := 1;
  IF C6 = N1 THEN
    {T2=identity matrix}
    FOR I := 1 TO 5 DO
      FOR J := 1 TO 5 DO
        IF J = I THEN
          T2[I, J] := 1.0
        ELSE
          T2[I, J] := 0.0;
  W := 2 * PI * F;
  W2 := W * W;
  {T=0}
  FOR I := 1 TO 5 DO
    FOR J := 1 TO 5 DO
      T[I, J] := 0.0;
  FOR I := 1 TO N DO
  begin
    {A:=identity matrix}
    FOR K := 1 TO 5 DO
      FOR J := 1 TO 5 DO
        IF J = K THEN
          A[K, J] := 1.0
        ELSE
          A[K, J] := 0.0;
    C3 := C[I];
    {Select according to kind of element}
    IF C3 = 1 THEN
      S2500
    ELSE IF C3 = 2 THEN
      S2450
    ELSE IF C3 = 3 THEN
      S2550
    ELSE IF C3 = 4 THEN
      S2600
    ELSE IF C3 = 5 THEN
      S2650
    ELSE IF C3 = 6 THEN
      S2700;

    IF I = 1 THEN
      {T=A}
      FOR K := 1 TO 5 DO
        FOR J := 1 TO 5 DO
          T[K, J] := A[K, J]
    ELSE
    begin
      {T1=A MPY T}
      FOR I1 := 1 TO 5 DO
        FOR J := 1 TO 5 DO
        begin
          Sum := 0.0;
          FOR K := 1 TO 5 DO
          begin
            Sum := Sum + A[I1, K] * T[K, J];
            T1[I1, J] := Sum
          end
        end;
      {T=T1}
      FOR K := 1 TO 5 DO
        FOR J := 1 TO 5 DO
          T[K, J] := T1[K, J]
    end;
    C6 := C6 + 1;
    IF C6 = N1 THEN
      {T2=T}
      FOR K := 1 TO 5 DO
        FOR J := 1 TO 5 DO
          T2[K, J] := T[K, J]
  end; {I loop}

  IF C8 = 1 THEN
  begin
    {Fixed-Fixed}
    I1 := 2; J1 := 1
  end
  ELSE IF C8 = 2 THEN
  begin
    {Fixed-Free}
    I1 := 1; J1 := 1
  end
  ELSE IF C8 = 3 THEN
  begin
    {Free-Free}
    I1 := 1; J1 := 2
  end
  ELSE IF C8 = 4 THEN
  begin
    {Free-Fixed}
    I1 := 2; J1 := 2
  end;

  B[1, 1] := T[I1, J1];
  B[1, 2] := T[I1, J1 + 2];
  B[2, 1] := T[I1 + 2, J1];
  B[2, 2] := T[I1 + 2, J1 + 2];
  E1[1] := -T[I1, 5];
  E1[2] := -T[I1 + 2, 5];
  D := B[1, 1] * B[2, 2] - B[1, 2] * B[2, 1];

  {Zero Divide Protection}
  IF ABS(D) < TINY THEN
  begin
    Writeln(' CAUTION, SYSTEM IS SINGULAR!');
    GOTO 1200   {Change limit conditions}
  end;

  E2[1] := (B[2, 2] * E1[1] - B[1, 2] * E1[2]) / D;
  E2[2] := (B[1, 1] * E1[2] - B[2, 1] * E1[1]) / D;

  {V1=0}
  FOR I := 1 TO 5 DO V1[I] := 0.0;
  V1[5] := 1.0;
  IF J1 <> 1 THEN
  begin
    V1[2] := E2[1]; V1[4] := E2[2]
  end
  ELSE
  begin
    V1[1] := E2[1]; V1[3] := E2[2]
  end;

  {V2:=T2 MPY V1}
  FOR I := 1 TO 5 DO
  begin
    Sum := 0.0;
    FOR J := 1 TO 5 DO
      Sum := Sum + T2[I, J] * V1[J];
    V2[I] := Sum
  end;
  M2 := SQRT(Sqr(V2[2]) + Sqr(V2[4]));
  IF ABS(M2) < TINY THEN
    T3 := 0.0
  ELSE
  begin
    C4 := V2[2] / M2;
    S4 := -V2[4] / M2;
    IF C4 = -1.0 THEN
      T3 := PI
    ELSE IF C4 = 1.0 THEN
      T3 := 0.0
    ELSE
      IF S4 >= 0.0 THEN
        {T3:=ACOS(C4)}
        T3 := ArcTan(SQRT(1.0 - C4*C4) / C4)
      ELSE
        {T3:=2*PI-ACOS(C4)}
        T3 := 2.0 * PI - ArcTan(SQRT(1.0 - C4*C4) / C4);
  end;
  {Convert phase in degrees}
  T3 := T3 / PI * 180.0;
  {Writeln current results}
  Writeln(' Freq=', F:7:3,'  Displacement=', M2:7:1,'  Phase= ',T3:5:1,' Deg.');
  F:=F+F5
Until F>F4;
{end frequency loop}

Writeln;
Writeln(' Change limit conditions; LIM');
Writeln(' Change frequency sweep ; SWE');
Writeln(' Change response node...; NOD');
Writeln(' Exit program...........; EXI');
Writeln;
Write(' Your answer: '); Readln(Ans);

IF Ans = 'LIM' THEN GOTO 1200;
IF Ans = 'SWE' THEN GOTO 1300;
IF Ans = 'NOD' THEN GOTO 1180;

DoneWinCrt

END.

{end of file ndof01.pas}