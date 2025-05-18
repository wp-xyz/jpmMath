{*********************************************************************
*         NORMAL AND INVERSE NORMAL PROBABILITY FUNCTIONS            *
* ------------------------------------------------------------------ *
* SAMPLE RUN:                                                        *
*                                                                    *
* Value of variable: 0.2                                             *
*                                                                    *
* Probability =  3.91042693975456E-0001                              *
*                                                                    *
* Verify:                                                            *
*   X =  2.00000286102295E-0001                                      *
*                                                                    *
*                             Pascal Version By J-P Moreau, Paris.   *
*                                      (www.jpmoreau.fr)             *
*********************************************************************}
Program Normal_Law;

Uses WinCrt;

Var
    P, X, Y: Double;


FUNCTION PHI (u: DOUBLE): Double;
{Standard Normal Probability Function}
Begin
  PHI := (1.0/Sqrt(2.0 * PI)) * Exp(-0.5*u*u)
End;

FUNCTION NORMAL (P: DOUBLE): Double;
{Inverse Standard Normal Function}
Label fin;
Var
    chimax,chisqval,epsilon,minchisq,maxchisq: Double;
Begin
  epsilon := 1e-6;
  chimax := 1e6;
  minchisq := 0.0;
  maxchisq := chimax;
  IF P <= 0.0 THEN
  begin
    NORMAL := maxchisq;
    goto fin
  end
  ELSE
  begin
    IF P >= 1.0 THEN
    begin
      NORMAL := 0.0;
      goto fin
    end
  end;
  chisqval := 0.5;
  WHILE ABS(maxchisq - minchisq) > epsilon DO
  begin
    IF PHI(chisqval) < P THEN
      maxchisq := chisqval
    ELSE
      minchisq := chisqval;
    chisqval := (maxchisq + minchisq) * 0.5
  end;
  NORMAL := chisqval;
fin: End;


{main program}
BEGIN

  Writeln;
  Write(' Value of variable: '); Readln(X);
  Writeln;
  P := PHI(X);
  Writeln(' Probability = ', P);
  Writeln;
  Writeln(' Verify:');
  Y := NORMAL(P);
  Writeln('   X = ', Y);
  Writeln;
  ReadKey;
  DoneWinCrt

END.

{end of file normal.pas}