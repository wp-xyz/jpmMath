{**********************************************************************************
* Frequencies and Modes of Mass-spring Systems without Damping By Transfer Method *
* ------------------------------------------------------------------------------- *
* EXPLANATION:                                                                    *
* Find resonance frequencies and modes of 2 d.o.f. system below:                  *
*                                                                                 *
*                          --> x1          --> x2                                 *
*              \|    k     ------    k     ------    k     |/                     *
*              \|--/\/\/\--| 3m |--/\/\/\--| m  |--/\/\/\--|/                     *
*              \|          -----           -----           |/                     *
*                          --> F1(t)       --> F2(t)                              *
*                                                                                 *
* composed of two masses, 3m and m, and three springs of stiffness k, with fixed- *
* fixed limit conditions (damping and friction here neglected).                   *
*                                                                                 *
* The differential system of mass motions, x1(t) and x2(t), is given by:          *
*                   3mx1" + 2kx1 - kx2 = F1(t)                                    *
*                    mx2" +2kx2  - kx1 = F2(t)                                    *
* where F1(t) anf F2(t) are respectively the forces acting on masses 3m and m.    *
*                                                                                 *
* Or with matrix notation:              M x" + K x = F(t)                         *
*                                                                                 *
*                      | 3m  0 |                           | 2k  -k |             *
* with mass matrix M = |       |  and stiffness matrix K = |        |             *
*                      | 0   m |                           | -k  2k |             *
*                                                                                 *
*                         | x1 |                      | F1(t) |                   *
* displacement vector X = |    |  force vector F(t) = |       |                   *
*                         | x2 |                      | F2(t) |                   *
*                                                                                 *
* For F1(t)=F2(t)=0 (free motion), the exact solution is:                         *
*                     f1 = 0.10693 sqrt(k/m) Hz                                   *
*                     f2 = 0.23686 sqrt(k/m) Hz                                   *
* with corresponding modal vectors:                                               *
*                | X1 |   |  1   |           | X1 |   |   1  |                    *
*         phi1 = |    | = |      |    phi2 = |    | = |      |                    *
*                | X2 |   |0.6457|           | X2 |   |-4.646|                    *
* ------------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                     *
*                                                                                 *
*  Frequencies and Modes of Mass-spring Systems                                   *
*  By Transfer Method.                                                            *
*                                                                                 *
*  Kind of elements                                                               *
*  ----------------                                                               *
*  C(I)=1: Spring                                                                 *
*  C(I)=2: Mass                                                                   *
*                                                                                 *
*  Number of Elements: 5                                                          *
*                                                                                 *
*  Kind of element 1 C(1) = 1                                                     *
*  Kind of element 2 C(2) = 2                                                     *
*  Kind of element 3 C(3) = 1                                                     *
*  Kind of element 4 C(4) = 2                                                     *
*  Kind of element 5 C(1) = 1                                                     *
*                                                                                 *
*  Fixed-Fixed System, code=1.                                                    *
*  Fixed-Free System,  code=2.                                                    *
*  Free-Fixed System,  code=3.                                                    *
*  Free-Free System,   code=4.                                                    *
*                                                                                 *
*  Code = 1                                                                       *
*                                                                                 *
*  Mass #1 = 3                                                                    *
*  Mass #2 = 1                                                                    *
*                                                                                 *
*  Spring #1 = 1                                                                  *
*  Spring #2 = 1                                                                  *
*  Spring #3 = 1                                                                  *
*                                                                                 *
*  Frequency Sweep: (S)WEEP                                                       *
*  Calculate Mode.: (M)ODE                                                        *
*  Exit Program...: (E)XIT                                                        *
*                                                                                 *
*  Answer = S                                                                     *
*                                                                                 *
*  Frequency Sweep                                                                *
*  ---------------                                                                *
*  Starting Frequency: 0                                                          *
*  Ending Frequency..: 1                                                          *
*  Frequency Step....: 0.1                                                        *
*                                                                                 *
* Frequency= 0.00000000000000E+0000 Hz  Determinant= 3.00000000000000E+0000       *
* Frequency= 1.00000000000000E-0001 Hz  Determinant= 3.09290228614616E-0001       *
* Frequency= 2.00000000000000E-0001 Hz  Determinant=-2.15207544198299E+0000       *
* Frequency= 3.00000000000000E-0001 Hz  Determinant= 1.24481939188828E+0001       *
* Frequency= 4.00000000000000E-0001 Hz  Determinant= 7.21639165290047E+0001       *
* Frequency= 5.00000000000000E-0001 Hz  Determinant= 2.16270437893292E+0002       *
* Frequency= 6.00000000000000E-0001 Hz  Determinant= 4.95264630803773E+0002       *
* Frequency= 7.00000000000000E-0001 Hz  Determinant= 9.70864895339590E+0002       *
* Frequency= 8.00000000000000E-0001 Hz  Determinant= 1.71601115886700E+0003       *
* Frequency= 9.00000000000000E-0001 Hz  Determinant= 2.81486487603940E+0003       *
* Frequency= 1.00000000000000E+0000 Hz  Determinant= 4.36280902879725E+0003       *
*                                                                                 *
* Note: A minimum value of matrix determinant (here T[2,1]) occurs near a         *
*       resonance frequency, here f=0.1 hz. Refine sweep to improve solution.     *
*                                                                                 *
*  Frequency Sweep: (S)WEEP                                                       *
*  Calculate Mode.: (M)ODE                                                        *
*  Exit Program...: (E)XIT                                                        *
*                                                                                 *
*  Answer = S                                                                     *
*                                                                                 *
*  Frequency Sweep                                                                *
*  ---------------                                                                *
*  Starting Frequency: 0                                                          *
*  Ending Frequency..: 0.2                                                        *
*  Frequency Step....: 0.01                                                       *
*                                                                                 *
* Frequency= 0.00000000000000E+0000 Hz  Determinant= 3.00000000000000E+0000       *
* Frequency= 1.00000000000000E-0002 Hz  Determinant= 2.96846402228021E+0000       *
* Frequency= 2.00000000000000E-0002 Hz  Determinant= 2.87441716548520E+0000       *
* Frequency= 3.00000000000000E-0002 Hz  Determinant= 2.71954265870803E+0000       *
* Frequency= 4.00000000000000E-0002 Hz  Determinant= 2.50664588377048E+0000       *
* Frequency= 5.00000000000000E-0002 Hz  Determinant= 2.23965437522305E+0000       *
* Frequency= 6.00000000000000E-0002 Hz  Determinant= 1.92361782034494E+0000       *
* Frequency= 7.00000000000000E-0002 Hz  Determinant= 1.56470805914406E+0000       *
* Frequency= 8.00000000000000E-0002 Hz  Determinant= 1.17021908435703E+0000       *
* Frequency= 9.00000000000000E-0002 Hz  Determinant= 7.48567041449202E-0001       *
* Frequency= 1.00000000000000E-0001 Hz  Determinant= 3.09290228614617E-0001       *
* Frequency= 1.10000000000000E-0001 Hz  Determinant=-1.36950903223961E-0001       *
* Frequency= 1.20000000000000E-0001 Hz  Determinant=-5.78373750415060E-0001       *
* Frequency= 1.30000000000000E-0001 Hz  Determinant=-1.00207355657850E+0000       *
* Frequency= 1.40000000000000E-0001 Hz  Determinant=-1.39402341260537E+0000       *
* Frequency= 1.50000000000000E-0001 Hz  Determinant=-1.73907425665808E+0000       *
* Frequency= 1.60000000000000E-0001 Hz  Determinant=-2.02095487417030E+0000       *
* Frequency= 1.70000000000000E-0001 Hz  Determinant=-2.22227189784700E+0000       *
* Frequency= 1.80000000000000E-0001 Hz  Determinant=-2.32450980766444E+0000       *
* Frequency= 1.90000000000000E-0001 Hz  Determinant=-2.30803093087016E+0000       *
*                                                                                 *
* Note: Now the minimum absolute value is for f = 0.11 hz. Refine again value.    *
*                                                                                 *
*  Frequency Sweep: (S)WEEP                                                       *
*  Calculate Mode.: (M)ODE                                                        *
*  Exit Program...: (E)XIT                                                        *
*                                                                                 *
*  Answer = S                                                                     *
*                                                                                 *
*  Frequency Sweep                                                                *
*  ---------------                                                                *
*  Starting Frequency: 0.10                                                       *
*  Ending Frequency..: 0.12                                                       *
*  Frequency Step....: 0.0005                                                     *
*                                                                                 *
* Frequency= 1.00000000000000E-0001 Hz  Determinant= 3.09290228614616E-0001       *
* Frequency= 1.00500000000000E-0001 Hz  Determinant= 2.87050179054777E-0001       *
* Frequency= 1.01000000000000E-0001 Hz  Determinant= 2.64793891697701E-0001       *
* Frequency= 1.01500000000000E-0001 Hz  Determinant= 2.42522779754481E-0001       *
* Frequency= 1.02000000000000E-0001 Hz  Determinant= 2.20238263449666E-0001       *
* Frequency= 1.02500000000000E-0001 Hz  Determinant= 1.97941770021258E-0001       *
* Frequency= 1.03000000000000E-0001 Hz  Determinant= 1.75634733720712E-0001       *
* Frequency= 1.03500000000000E-0001 Hz  Determinant= 1.53318595812939E-0001       *
* Frequency= 1.04000000000000E-0001 Hz  Determinant= 1.30994804576305E-0001       *
* Frequency= 1.04500000000000E-0001 Hz  Determinant= 1.08664815302630E-0001       *
* Frequency= 1.05000000000000E-0001 Hz  Determinant= 8.63300902971892E-0002       *
* Frequency= 1.05500000000000E-0001 Hz  Determinant= 6.39920988787099E-0002       *
* Frequency= 1.06000000000000E-0001 Hz  Determinant= 4.16523173793777E-0002       *
* Frequency= 1.06500000000000E-0001 Hz  Determinant= 1.93122291448298E-0002       *
* Frequency= 1.07000000000000E-0001 Hz  Determinant=-3.02667546584057E-0003       *
* Frequency= 1.07500000000000E-0001 Hz  Determinant=-2.53628990800852E-0002       *
* Frequency= 1.08000000000000E-0001 Hz  Determinant=-4.76949373119040E-0002       *
* Frequency= 1.08500000000000E-0001 Hz  Determinant=-7.00212787618394E-0002       *
* Frequency= 1.09000000000000E-0001 Hz  Determinant=-9.23404050169805E-0002       *
* Frequency= 1.09500000000000E-0001 Hz  Determinant=-1.14650790650961E-0001       *
* Frequency= 1.10000000000000E-0001 Hz  Determinant=-1.36950903223962E-0001       *
* Frequency= 1.10500000000000E-0001 Hz  Determinant=-1.59239203282708E-0001       *
* Frequency= 1.11000000000000E-0001 Hz  Determinant=-1.81514144360470E-0001       *
* Frequency= 1.11500000000000E-0001 Hz  Determinant=-2.03774172977065E-0001       *
* Frequency= 1.12000000000000E-0001 Hz  Determinant=-2.26017728638852E-0001       *
* Frequency= 1.12500000000000E-0001 Hz  Determinant=-2.48243243838741E-0001       *
* Frequency= 1.13000000000000E-0001 Hz  Determinant=-2.70449144056184E-0001       *
* Frequency= 1.13500000000000E-0001 Hz  Determinant=-2.92633847757178E-0001       *
* Frequency= 1.14000000000000E-0001 Hz  Determinant=-3.14795766394267E-0001       *
* Frequency= 1.14500000000000E-0001 Hz  Determinant=-3.36933304406539E-0001       *
* Frequency= 1.15000000000000E-0001 Hz  Determinant=-3.59044859219629E-0001       *
* Frequency= 1.15500000000000E-0001 Hz  Determinant=-3.81128821245717E-0001       *
* Frequency= 1.16000000000000E-0001 Hz  Determinant=-4.03183573883528E-0001       *
* Frequency= 1.16500000000000E-0001 Hz  Determinant=-4.25207493518332E-0001       *
* Frequency= 1.17000000000000E-0001 Hz  Determinant=-4.47198949521946E-0001       *
* Frequency= 1.17500000000000E-0001 Hz  Determinant=-4.69156304252732E-0001       *
* Frequency= 1.18000000000000E-0001 Hz  Determinant=-4.91077913055596E-0001       *
* Frequency= 1.18500000000000E-0001 Hz  Determinant=-5.12962124261990E-0001       *
* Frequency= 1.19000000000000E-0001 Hz  Determinant=-5.34807279189912E-0001       *
* Frequency= 1.19500000000000E-0001 Hz  Determinant=-5.56611712143907E-0001       *
*                                                                                 *
* Note: now our best first resonance frequency is: 0.107 hz (exact value 0.10693).*
*       We can now try to calculate the first modal vector:                       *
*                                                                                 *
*  Frequency Sweep: (S)WEEP                                                       *
*  Calculate Mode.: (M)ODE                                                        *
*  Exit Program...: (E)XIT                                                        *
*                                                                                 *
*  Answer = M                                                                     *
*                                                                                 *
*  Calculate a Mode                                                               *
*  ----------------                                                               *
*  Resonance Frequency: 0.107                                                     *
*                                                                                 *
* Node #1 Force= 1.00000000000000E+0000 Displacement= 0.00000000000000E+0000      *
* Node #2 Force= 1.00000000000000E+0000 Displacement= 1.00000000000000E+0000      *
* Node #3 Force=-3.55965209456865E-0001 Displacement= 1.00000000000000E+0000      *
* Node #4 Force=-3.55965209456865E-0001 Displacement= 6.44034790543135E-0001      *
* Node #5 Force=-6.47061466008975E-0001 Displacement= 6.44034790543135E-0001      *
* Node #6 Force=-6.47061466008975E-0001 Displacement=-3.02667546583979E-0003      *
*                                                                                 *
* Note: here the modal vector found is [1.000 0.644], exact value is [1 0.6457].  *
*       A small error occurs because the frequency is still approximate.          *
*       Proceed in the same way to approximate 2nd frequency and modal vector.    *
*       This program is only intended to demonstrate the transfer method, it is   *
*       not suitable to solve, as is, real problems.                              *
* ------------------------------------------------------------------------------- *
* REFERENCE: "Mécanique des vibrations linéaires By M. Lalanne, P. Berthier,      *
*             J. Der Hagopian, Masson, Paris 1980" [BIBLI 16].                    *
*                                                                                 *
*                                           Pascal Release By J-P Moreau, Paris.  *
*                                                    (www.jpmoreau.fr)            *
**********************************************************************************}
Program Modes;

Uses WinCrt;

Label 680, 1000, 1280, Fin;

Const NMAX = 25;

Type
     MAT2 = Array[1..2,1..2] of Double;
     VEC = Array[1..NMAX] of Double;

Var
    I,N,C1,C6,C8,K9,M9: Integer;
    A,T,T1: MAT2;
    K1,M1: VEC;
    F,F2,F3,F4,W,W2,X: Double;
    C: Array[1..NMAX] of Integer;
    V1,V2: Array[1..2] of Double;
    Ans:String[4];

  Procedure Mass;
  Begin
    M9 := M9 + 1;
    A[1, 1] := 1.0; A[1, 2] := 0.0;
    A[2, 1] := 0.0; A[2, 2] := 1.0;
    A[1, 2] := -M1[M9] * W2
  end;

  Procedure Stiffness;
  Begin
    K9 := K9 + 1;
    A[1, 1] := 1.0; A[1, 2] := 0.0;
    A[2, 1] := 0.0; A[2, 2] := 1.0;
    A[2, 1] := 1.0 / K1[K9]
  End;

{main program}
BEGIN

  ClrScr;
  Writeln;
  Writeln(' Frequencies and Modes of Mass-spring Systems');
  Writeln(' By Transfer Method.');
  Writeln;
  Writeln(' Kind of elements');
  Writeln('-----------------');
  Writeln('  C(I)=1: Spring');
  Writeln('  C(I)=2: Mass');
  Writeln;
  Write(' Number of Elements: '); Readln(N);
  M9 := 0; K9 := 0;

  Writeln;
  FOR I := 1 TO N do
  begin
    Write(' Kind of element', I, ' C(', I, ') = ');
    Readln(C[I]);
    IF C[I] = 1 THEN K9 := K9 + 1;
    IF C[I] = 2 THEN M9 := M9 + 1
  end;

  Writeln;
  Writeln(' Fixed-Fixed System, code=1.');
  Writeln(' Fixed-Free System,  code=2.');
  Writeln(' Free-Fixed System,  code=3.');
  Writeln(' Free-Free System,   code=4.');
  Writeln;
  Write(' Code = '); Readln(C8);
  Writeln;
  IF M9 <> 0 THEN  {There are masses}
    FOR I := 1 TO M9 do
    begin
      Write(' Mass #', I,' = ');
      Readln(M1[I])
    end;
  Writeln;
  IF K9 <> 0 THEN  {There are springs}
    FOR I := 1 TO K9 do
    begin
      Write(' Spring #', I,' = ');
      Readln(K1[I])
    end;

680: Writeln;
  Writeln(' Frequency Sweep: (S)WEEP');
  Writeln(' Calculate Mode.: (M)ODE');
  Writeln(' Exit Program...: (E)XIT');
  Writeln;
  Write(' Answer = '); Readln(Ans);

  IF Ans[1] = 'M' THEN GOTO 1000;
  IF Ans[1] = 'E' THEN GOTO Fin;

  Writeln;
  Writeln(' Frequency Sweep');
  Writeln(' ---------------');
  Write(' Starting Frequency: '); Readln(F2);
  Write(' Ending Frequency..: '); Readln(F3);
  Write(' Frequency Step....: '); Readln(F4);
  GOTO 1280;

1000: Writeln(' Calculate a Mode');
  Writeln(' ----------------');
  Write(' Resonance Frequency: '); Readln(F);
  M9 := 0; K9 := 0;
  W := 2 * PI * F;
  W2 := W * W;
  V1[1] := 0.0; V1[2] := 0.0; C6 := 1;
  IF (C8 = 3) OR (C8 = 4) THEN   {system is free-fixed or free-free}
    V1[2] := 1.0
  ELSE
    V1[1] := 1.0;
  Writeln(' Node #', C6,' Force=', V1[1],' Displacement=', V1[2]);
  FOR I := 1 TO N do
  begin
    C6 := C6 + 1;
    C1 := C[I];
    IF C1 = 1 THEN
      Stiffness
    ELSE
      Mass;
    {V2 := A MPY V1}
    V2[1] := A[1, 1] * V1[1] + A[1, 2] * V1[2];
    V2[2] := A[2, 1] * V1[1] + A[2, 2] * V1[2];
    V1[1] := V2[1]; V1[2] := V2[2];
    Writeln(' Node #', C6,' Force=', V1[1],' Displacement=', V1[2])
  end;
  GOTO 680;

1280: {Frequency sweep}
  F:=F2;
  Repeat
    M9 := 0; K9 := 0;
    W := 2.0 * PI * F; W2 := W * W;
    T[1,1]:= 0.0; T[1,2]:=0.0; T[2,1]:= 0.0; T[2,2]:=0.0;
    FOR I := 1 TO N do
    begin
      C1 := C[I];
      IF C1 = 1 THEN
        Stiffness
      ELSE
        Mass;
      IF I = 1 THEN
      begin
        {T=A}
        T[1, 1] := A[1, 1]; T[1, 2] := A[1, 2];
        T[2, 1] := A[2, 1]; T[2, 2] := A[2, 2]
      end
      ELSE
      begin
        {T1=A MPY T}
        T1[1, 1] := A[1, 1] * T[1, 1] + A[1, 2] * T[2, 1];
        T1[1, 2] := A[1, 1] * T[1, 2] + A[1, 2] * T[2, 2];
        T1[2, 1] := A[2, 1] * T[1, 1] + A[2, 2] * T[2, 1];
        T1[2, 2] := A[2, 1] * T[1, 2] + A[2, 2] * T[2, 2];
        {T=T1}
        T[1, 1] := T1[1, 1]; T[1, 2] := T1[1, 2];
        T[2, 1] := T1[2, 1]; T[2, 2] := T1[2, 2]
      end
    end;
    IF C8 = 1 THEN
      X := T[2, 1]
    ELSE IF C8 = 2 THEN
      X := T[1, 1]
    ELSE IF C8 = 3 THEN
      X := T[2, 2]
    ELSE
      X := T[1, 2];
    Writeln(' Frequency=', F,' Hz  Determinant=', X);
    F:=F+F4
  Until F >= F3;
  Writeln;
  GOTO 680;

Fin: DoneWinCrt

END.
 
{end of file modes.pas}