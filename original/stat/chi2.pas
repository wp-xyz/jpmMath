{*************************************************************************
*                  CHI2 AND INVERSE CHI2 FUNCTIONS                       *
* ---------------------------------------------------------------------- *
*  Reference:                                                            *
*                                                                        *
*  http://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html        *
*                                                                        *
* The following JavaScript functions for calculating chi-square probabi- *
* lities and critical values were adapted by John Walker from C imple-   *
* mentations written by Gary Perlman of Wang Institute, Tyngsboro.       *
* Both the original C code and the JavaScript edition are in the public  *
* domain.                                                                *
* ---------------------------------------------------------------------- *
* SAMPLE RUN:                                                            *
*                                                                        *
*    CHI2 Law:                                                           *
*    ========                                                            *
*  Value of variable: 100                                                *
*  Degree of freedom: 100                                                *
*                                                                        *
*  Probability =  4.81191684527957E-0001                                 *
*                                                                        *
*  Verify Inverse Law:                                                   *
*    X =  1.00000000148592E+0002                                         *
*                                                                        *
*                                   TPW Version By J-P Moreau, Paris.    *
*                                           (www.jpmoreau.fr)            *
*************************************************************************}
Program Stats01;
Uses WinCrt;

Var
    P,X,Y: Double;
    idf: Integer;


FUNCTION poz (z:DOUBLE): Double;
{------------------------------------------------------------------------
 POZ  --  probability of standard normal z value

 Adapted from a polynomial approximation in:
 Ibbetson D, Algorithm 209
 Collected Algorithms of the CACM 1963 p. 616

 Note:
 This routine has six digits accuracy, so it is only useful for absolute
 z values < 6.  For z values >=  6.0,  poz() returns 0.0.
 -----------------------------------------------------------------------}
Var
    Y, X, w, zmax: DOUBLE;
Begin
    zmax := 6.0;
    IF z = 0.0 THEN
      X := 0.0
    ELSE
    begin
      Y := 0.5 * ABS(z);
      IF Y >= zmax * 0.5 THEN
        X := 1.0
      ELSE IF Y < 1.0 THEN
      begin
        w := Y * Y;
        X := ((((((((0.000124818987*w -0.001075204047)*w + 0.005198775019)*w -0.019198292004)*w + 0.059054035642)*w
             -0.151968751364)*w + 0.319152932694) * w - 0.5319230073)*w + 0.7978845605929999)*Y*2
      end
      ELSE
      begin
        Y := Y - 2;
        X := (((((((((((((-0.000045255659*Y + 0.00015252929)*Y -0.000019538132)*Y - 6.769049860000001E-04)*Y
             +0.001390604284)*Y-0.00079462082)*Y -0.002034254874)*Y + 0.006549791214)*Y - 0.010557625006)*Y
             +0.011630447319)*Y-9.279453341000001E-03)*Y + 0.005353579108)*Y - 0.002141268741)*Y
             +0.000535310849)*Y + 0.999936657524
      end
    end;
    IF z > 0.0 THEN
      poz := ((X + 1) * 0.5)
    ELSE
      poz := ((1 - X) * 0.5)
    
End;


{chi2 function}
FUNCTION chi2 (X: DOUBLE; df: INTEGER): Double;
{--------------------------------------------------------
 Adapted From:
 Hill, I. D. and Pike, M. C. Algorithm 299
 Collected Algorithms for the CACM 1967 p. 243
 Updated for rounding errors based on remark in
 ACM TOMS June 1985, page 185
 -------------------------------------------------------}
Var
    even: Boolean;   {parity of df}
    a,bigx,c,e,lnpi, ipi,s,Y,z: Double;
Begin
    bigx := 200.0;
    lnpi := Ln(Sqrt(PI));
    ipi := 1.0 / Ln(PI);
    IF (X <= 0.0) OR (df < 1) THEN chi2 := 1.0;
    a := 0.5 * X;
    IF 2*(df Div 2) = df THEN
      even := True
    ELSE
      even := False;
    IF df > 1 THEN  Y := EXP(-a);
    IF even = True THEN
      s := Y
    ELSE
      s := 2.0 * poz(-Sqrt((X)));
    IF df > 2 THEN
    begin
        X := 0.5 * (df - 1);
        IF even THEN
          z := 1.0
        ELSE
          z := 0.5;
        IF (a > bigx) THEN
        begin
            IF even THEN
              e := 0.0
            ELSE
              e := lnpi;
            c := Ln(a);
            WHILE z <= X do
            begin
              e := Ln(z) + e;
              s := s + EXP(c * z - a - e);
              z := z + 1.0
            end;
            chi2 := s
        end
        ELSE
        begin
            IF even THEN
              e := 1.0
            ELSE
              e := ipi / Sqrt(a);
            c := 0.0;
            WHILE z <= X do
            begin
              e := e * (a / z);
              c := c + e;
              z := z + 1.0
            end;
            chi2 := c * Y + s
        end
    end
    ELSE
      chi2 := s

End;

{Inverse chi2 function}
FUNCTION invchi2 (P: DOUBLE; idf: INTEGER): Double;
Var chimax,chisqval,epsilon,minchisq,maxchisq: Double;
Begin
    epsilon := 1e-6;
    chimax  := 1e6;
    minchisq := 0.0;
    maxchisq := chimax;
    IF P <= 0.0 THEN
      invchi2 := maxchisq
    ELSE
      IF P >= 1.0 THEN invchi2 := 0.0;
    chisqval := 1.0 * idf / Sqrt(P);
    WHILE ABS(maxchisq - minchisq) > epsilon DO
    begin
      IF chi2(chisqval, idf) < P THEN
        maxchisq := chisqval
      ELSE
        minchisq := chisqval;
      chisqval := (maxchisq + minchisq) * 0.5
    end;
    invchi2 := chisqval
End;

{main program}
BEGIN

  Writeln;
  Writeln('  CHI2 Law:');
  Writeln('  ========');
  Write(' Value of variable: '); Readln(X);
  Write(' Degree of freedom: '); Readln(idf);
  Writeln;
  P := chi2(X, idf);
  Writeln(' Probability = ', P);
  Writeln;
  Writeln(' Verify Inverse Law:');
  Y := invchi2(P, idf);
  Writeln('   X = ', Y);
  Writeln;

  ReadKey;
  DoneWinCrt

END.

{end of file chi2.pas}