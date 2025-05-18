{**********************************************************************
* Response of a 1dof Mass-Spring System with damping to a sinusoidal  *
* input force                                                         *
* ------------------------------------------------------------------- *
* Main Variables:                                                     *
*                  M: Mass                                            *
*                  K: Stiffness                                       *
*                  C: Viscous Damping Coefficient                     *
*                  E: Structural Damping Coefficient                  *
*                 F0: Input Force                                     *
* ------------------------------------------------------------------- *
* SAMPLE RUN:                                                         *
*                                                                     *
*  Mass M..... = 50                                                   *
*  Stiffness K = 1e6                                                  *  
*  Force F0... = 10000                                                *
*                                                                     *
*  (V)iscous or (S)tructural Damping: S                               *
*                                                                     *
*  Structural damping factor: 0.05                                    *
*                                                                     *
*  Frequency Scanning                                                 *
*                                                                     *
*  Resonance Frequency of undamped system = 22.508                    *
*                                                                     *
*  Starting frequency = 20                                            *
*  Ending frequency...= 25                                            *
*  Frequency step.....= 0.25                                          *
*                                                                     *
* Frequency   Displacement   Phase (deg)                              *
* --------------------------------------                              * 
*  20.00      0.046234       13.366013                                *
*  20.25      0.050756       14.701446                                *
*  20.50      0.056293       16.347699                                *
*  20.75      0.063206       18.423074                                *
*  21.00      0.072037       21.111461                                *
*  21.25      0.083609       24.711299                                *
*  21.50      0.099181       29.729354                                *
*  21.75      0.120525       37.058260                                *
*  22.00      0.149218       48.252814                                *
*  22.25      0.181993       65.500492                                *
*  22.50      0.199980       89.194985                                *
*  22.75      0.183564      -66.609176                                *
*  23.00      0.149839      -48.520651                                *
*  23.25      0.119585      -36.721514                                *
*  23.50      0.097048      -29.028193                                *
*  23.75      0.080680      -23.790753                                *
*  24.00      0.068578      -20.053143                                *
*  24.25      0.059388      -17.273961                                *
*  24.50      0.052222      -15.136029                                *
*  24.75      0.046502      -13.444959                                *
*  25.00      0.041843      -12.076311                                *
*                                                                     *
* ------------------------------------------------------------------- *
* REFERENCE: "Mécanique des vibrations linéaires By M. Lalanne,       *
*             P. Berthier, J. Der Hagopian, Masson, Paris 1980"       *
*             [BIBLI 16].                                             *
*                                                                     *
*                               Pascal Release By J-P Moreau, Paris.  *
*                                        (www.jpmoreau.fr)            *
**********************************************************************}
Program SIN1DOF;

Uses WinCrt;

Var
    C,D,E,F,F0,F1,F2,F3,F4,M, K, I1,M1,O1,T1,T3,TMP,W: REAL;
    C2: Integer;
    B: Char;

BEGIN

  Writeln;
  Write(' Mass M..... = '); readln(M);
  Write(' Stiffness K = '); readln(K);
  Write(' Force F0... = '); readln(F0);
  Writeln;

  Repeat
    Write(' (V)iscous or (S)tructural Damping: '); Readln(B)
  Until B IN ['V','S'];

  Writeln;
  if B='V' then
  begin
    Write(' Viscous damping coefficient C = '); Readln(C);
    C2 := 1
  end
  else
  begin
    Write(' Structural damping factor: '); Readln(E);
    C2 := 2
  end;

  Writeln;
  Writeln(' Frequency Scanning');
  Writeln;
  F4 := 1.0 / 2.0 / PI * SQRT(K / M);
  Writeln(' Resonance Frequency of undamped system = ', F4:7:3,' Hz.');
  Writeln;
  Write(' Starting frequency = '); Readln(F1);
  Write(' Ending frequency...= '); Readln(F2);
  Write(' Frequency step.....= '); Readln(F3);

  F := F1;

  ClrScr;
  Writeln(' Frequency   Displacement   Phase (deg) ');
  Writeln(' -------------------------------------- ');

  While F<=F2 do
  begin
    W := 2.0 * PI * F;
    TMP:=K-M*W*W;
    IF C2 = 1 THEN {Viscous damping}
    begin
      D  := TMP*TMP + C*W*C*W;
      O1 := TMP / SQRT(D);
      I1 := C * W / SQRT(D);
      M1 := F0 / SQRT(D)
    end
    else   {Structural damping}
    begin
      D  := TMP*TMP + E*K*E*K;
      O1 := TMP / SQRT(D);
      I1 := E * K / SQRT(D);
      M1 := F0 / SQRT(D)
    end;

    IF M1 = 0.0 THEN
      T1 := 0.0
    else
    begin
      IF O1 = -1 THEN
        T1 := PI
      ELSE IF O1 = 1 THEN
        T1 := 0.0
      else
      begin
        IF I1 > 0 THEN
          {T1 = ACOS(O1) }
          T1 := ArcTan(SQRT(1.0 - O1*O1) / O1)
        else
          T1 := 2.0 * PI - ArcTan(SQRT(1.0 - O1*O1) / O1)
      end
    end;

    T3 := T1 / PI * 180;

    Writeln(F:7:2, M1:14:6,'  ',T3:14:6); 
    F := F + F3
  end;

  ReadKey;
  DoneWinCrt

END.

{end of file 1dof01.pas}