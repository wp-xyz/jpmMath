{************************************************************
*                STUDENT  T-PROBABILITY  LAW                *
* --------------------------------------------------------- *
* SAMPLE RUN:                                               *
* (Calculate Student T-probability (lower tail and upper    *
*  tail) for T=0.257).                                      *
*                                                           *
*  X= 0.2570000                                             *
*  PROB1= 0.6000294                                         *
*  ERROR=0                                                  *
*  X= 0.2570000                                             *
*  PROB2= 0.3999705                                         *
*  ERROR=0                                                  *
*  PROB1+PROB2= 0.9999998                                   *
*                                                           *
* --------------------------------------------------------- *
* Ref.:"JOURNAL OF APPLIED STATISTICS (1968) VOL.17, P.189, *
*       & VOL.19, NO.1".                                    *
*                                                           *
*                     Pascal Release By J-P Moreau, Paris.  *
*                               (www.jpmoreau.fr)           *
************************************************************}
PROGRAM student_test;

  {Calculate x power n}
  function Power(x:REAL; n:Integer):REAL;
  var
    i,m : integer;
  begin
    Result := 1.0;
    if n=0 then
        exit;
    m :=  n;
    if n < 0 then m := -n;
    for i:=1 to m do
      Result := x*Result;
    if n < 0 then Result := 1.0 / Result;
  end;


  function ProbST(T: REAL; IDF: Integer; Var ErrCode: Integer): REAL;
  { ---------------------------------------------------------------------
  !        ALGORITHM AS 3  APPL. STATIST. (1968) VOL.17, P.189
  !
  !        STUDENT T PROBABILITY (LOWER TAIL)
  ! --------------------------------------------------------------------}
  const
    G1:REAL = 1.0 / PI;
  var
    A, B, C, F, S, FK: REAL;
    IM2,IOE,K,KS: Integer;
  begin
    ErrCode := 1;
    Result := 0.0;

    if IDF < 1 then
      Exit;

    ErrCode := 0;
    F := 1.0 * IDF;
    A := T / SQRT(F);
    B := F / (F + T * T);
    IM2 := IDF - 2;
    IOE := IDF mod 2;
    S := 1.0;
    C := 1.0;
    F := 1.0;
    KS := 2 + IOE;
    FK := KS;
    if IM2 >= 2 then
    begin
      K:=KS;
      while K <= IM2 do
      begin
        C := C * B * (FK - 1.0) / FK;
        S := S + C;
        if S = F then break;
        F := S;
        FK := FK + 2.0;
        Inc(K, 2)
      end;
    end;
    if IOE <> 1 then begin
      Result := 0.5 + 0.5 * A * sqrt(B) * S;
      Exit;
    end;
    if IDF = 1 then S := 0.0;
    Result := 0.5 + (A * B * S + arctan(A)) * G1;
  end;

  function StudNT (T, DOFF:REAL; Var ErrCode:Integer): REAL;
  { ----------------------------------------------------------------
  !     ALGORITHM AS 27  APPL. STATIST. VOL.19, NO.1
  !
  !     Calculate the upper tail area under Student's t-distribution
  !
  !     Translated from Algol by Alan Miller
  ! ---------------------------------------------------------------}
  var
    V, X, TT: REAL;
    A1, A2, A3, A4, A5, B1, B2, C1, C2, C3, C4, C5, D1, D2: REAL;
    E1, E2, E3, E4, E5, F1, F2, G1, G2, G3, G4, G5, H1, H2: REAL;
    I1, I2, I3, I4, I5, J1, J2: REAL;
    Positive: Boolean;
  begin
    A1:=0.09979441; A2:=-0.581821; A3:=1.390993; A4:=-1.222452; A5:=2.151185;
    B1:=5.537409; B2:=11.42343;
    C1:=0.04431742; C2:=-0.2206018; C3:=-0.03317253; C4:=5.679969; C5:=-12.96519;
    D1:=5.166733; D2:=13.49862;
    E1:=0.009694901; E2:=-0.1408854; E3:=1.88993; E4:=-12.75532; E5:=25.77532;
    F1:=4.233736; F2:=14.3963;
    G1:=-9.187228E-5; G2:=0.03789901; G3:=-1.280346; G4:=9.249528; G5:=-19.08115;
    H1:=2.777816; H2:=16.46132;
    I1:=5.79602E-4; I2:=-0.02763334; I3:=0.4517029; I4:=-2.657697; I5:=5.127212;
    J1:=0.5657187; J2:=21.83269;

  { Check that number of degrees of freedom > 4 }
    if DOFF < 2.0 then
    begin
      ErrCode := 1;
      Result := -1.0;
      exit;
    end;

    if DOFF <= 4.0 then
      ErrCode := Round(DOFF)
    else
      ErrCode := 0;

    { Evaluate series }
    V := 1.0 / DOFF;
    Positive := (T >= 0.0);
    TT := abs(T);
    X := (1.0 + TT * (((A1 + V * (A2 + V * (A3 + V * (A4 + V * A5)))) / (1.0 - V * (B1 - V * B2))) +
         TT * (((C1 + V * (C2 + V * (C3 + V * (C4 + V * C5)))) / (1.0 - V * (D1 - V * D2))) +
         TT * (((E1 + V * (E2 + V * (E3 + V * (E4 + V * E5)))) / (1.0 - V * (F1 - V * F2))) +
         TT * (((G1 + V * (G2 + V * (G3 + V * (G4 + V * G5)))) / (1.0 - V * (H1 - V * H2))) +
         TT * ((I1 + V * (I2 + V * (I3 + V * (I4 + V * I5)))) /
         (1.0 - V * (J1 - V * J2))))))));
    X := 0.5 * Power(X,-8);
    if Positive then
      Result := X
    else
      Result := 1.0 - X;
  end;

{main program}

var
  X, XNU, PROB1, PROB2: REAL;
  NU, ERROR: INTEGER;

begin
  X:=0.257;
  NU:=19;

  Prob1 := ProbST(X, NU, ERROR);
  writeln;
  writeln('  X=', X:10:7);
  writeln('  PROB1=', PROB1:10:7);
  writeln('  ERROR=', ERROR);
  WriteLn;

  X:=0.257;
  XNU:=19.0;

  Prob2 := StudNT(X, XNU, ERROR);
  writeln('  X=', X:10:7);
  writeln('  PROB2=', Prob2:10:7);
  writeln('  ERROR=', ERROR);

  writeln('  PROB1+PROB2=', Prob1 + Prob2:10:7);
  writeln;

  Write('Press ENTER to close...');
  ReadLn;
end.

{end of file Student.pas}
