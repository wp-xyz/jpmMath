{**********************************************************************
*   Response of a 1dof Mass-Spring System with viscous damping to a   *
*   periodic input force:  M X" + C X' + K X = F(t)                   *
* ------------------------------------------------------------------- *
* Main Variables:                                                     *
*                  M: Mass                                            *
*                  K: Stiffness                                       *
*                  C: Viscous Damping Coefficient                     *
*                  F(t): Periodic Force of Frequency F8               *
*                                                                     *
*                         N                                           *
* Note:    F(t) = A0/2 + Sum [A(P) COS(PWt) + B(P) SIN(PWt)]          *
*                        P=1                                          *
* ------------------------------------------------------------------- *
* SAMPLE RUN:                                                         *
*                                                                     *
*  Mass M = 50                                                        *
*  Stiffness K = 1e6                                                  *
*  Viscous damping coefficient: 0.1                                   *
*                                                                     *
*                   N                                                 *
*    F(t) = A0/2 + Sum [A(P) COS(PWt) + B(P) SIN(PWt)]                *
*                  P=1                                                *
*                                                                     *
*  N = 1                                                              *
*                                                                     *
*  Excitation frequency = 20                                          *
*                                                                     *
*  A0 = 0                                                             *
*                                                                     *
*  A(1) = 0                                                           *
*                                                                     *
*  B(1) = 10000                                                       *
*                                                                     *
*  Time Scanning                                                      *
*                                                                     *
*  Excitation Period = 0.050                                          *
*                                                                     *
*  Starting time = 0                                                  *
*  Ending time...= 0.05                                               *
*  Time step.....= 0.0025                                             *
*                                                                     *
*   Time     Displacement                                             *
* ------------------------                                            *
*  0.0000    -0.000003                                                *
*  0.0025     0.014682                                                *
*  0.0050     0.027930                                                *
*  0.0075     0.038444                                                *
*  0.0100     0.045195                                                *
*  0.0125     0.047521                                                *
*  0.0150     0.045196                                                *
*  0.0175     0.038447                                                *
*  0.0200     0.027935                                                *
*  0.0225     0.014688                                                *
*  0.0250     0.000003                                                *
*  0.0275    -0.014682                                                *
*  0.0300    -0.027930                                                *
*  0.0325    -0.038444                                                *
*  0.0350    -0.045195                                                *
*  0.0375    -0.047521                                                *
*  0.0400    -0.045196                                                *
*  0.0425    -0.038447                                                *
*  0.0450    -0.027935                                                *
*  0.0475    -0.014688                                                *
*  0.0500    -0.000003                                                *
* ------------------------------------------------------------------- *
* REFERENCE: "Mécanique des vibrations linéaires By M. Lalanne,       *
*             P. Berthier, J. Der Hagopian, Masson, Paris 1980"       *
*             [BIBLI 16].                                             *
*                                                                     *
*                               Pascal Release By J-P Moreau, Paris.  *
*                                        (www.jpmoreau.fr)            *
**********************************************************************}
Program ONEDOF02;

Uses WinCrt;

Const NMAX=25;

Var
    A0,A1,C,C0,C1,C2,D,F8,M,K,S0,T,T0,T1,T2,T3,W,X3,Y1,Y2: Double;
    A,B,T4,X: Array[1..NMAX] of Double;   {coefficients of periodic force}
    I,J,N,N1: Integer;

BEGIN

  Writeln;
  Write(' Mass M = '); Readln(M);
  Write(' Stiffness K = '); Readln(K);
  Write(' Viscous damping coefficient: '); Readln(C);

  Writeln;
  Writeln('                N');
  Writeln(' F(t) = A0/2 + Sum [A(P) COS(PWt) + B(P) SIN(PWt)]');
  Writeln('               P=1');
  Writeln;
  Write(' N = '); Readln(N);
  Writeln;
  Write(' Excitation frequency = '); Readln(F8);
  W := 2.0 * PI * F8;
  Writeln;
  Write(' A0 = '); Readln(A0);
  Writeln;
  FOR I := 1 TO N do
    Write(' A(', I,')= '); Readln(A[I]);
  Writeln;
  FOR I := 1 TO N do
    Write(' B(', I,')= '); Readln(B[I]);
  Writeln;
  Writeln(' Time Scanning');
  Writeln;
  T0 := 2.0 * PI / W;
  Writeln(' Excitation Period = ', T0:7:3);
  Writeln;
  Write(' Starting time = '); Readln(T1);
  Write(' Ending time...= '); Readln(T2);
  Write(' Time step.....= '); Readln(T3);

  N1 := Round(1 + (T2 - T1) / T3);  {Number of response points}
  J := 0;
  Y1 := 1E+20;
  Y2 := 0.0;

  {Main loop}
  T:=T1;
  Repeat
    J := J + 1;
    X[J] := A0 / K / 2.0;
    T4[J] := T;
    FOR I := 1 TO N DO
    Begin
      A1 := I * W * T;
      C0 := COS(A1);
      S0 := SIN(A1);
      C1 := K - M * I * I * W * W;
      C2 := C * I * W;
      D := C1 * C1 + C2 * C2;
      IF D = 0.0 THEN
      begin
        Writeln('STOP - ZERO DIVIDE ERROR!');
        Halt
      end
      ELSE
        X[J] := X[J] + A[I] * (C1 * C0 + C2 * S0) / D + B[I] * (C1 * S0 - C2 * C0) / D
    End;
    X3 := X[J];
    IF Y1 > X3 THEN Y1 := X3;
    IF Y2 < X3 THEN Y2 := X3;
    T:=T+T3;
  Until T >= T2+T3;

  J := 0;
  ClrScr;
  Writeln;
  Writeln('   Time     Displacement ');
  Writeln(' ------------------------');
  T:=T1;
  Repeat
    J := J + 1;
    Writeln(T:8:4,'  ',X[J]:11:6);
    T:=T+T3
  Until T >= T2+T3;

  ReadKey;
  DoneWinCrt

END.

{end of file 1dof02.pas}