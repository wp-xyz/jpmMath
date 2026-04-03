{*************************************************************
*       Calculate Function Gamma with a complex argument     *
* ---------------------------------------------------------- *
* EXPLANATION:                                               *
* Purpose: This program computes the gamma function G(z)     *
*          or Ln[G(z)] for a complex argument using          *
*          subroutine CGAMA                                  *
* Input :  x  --- Real part of z                             *
*          y  --- Imaginary part of z                        *
*          KF --- Function code                              *
*          KF=0 for Ln[G(z)]                                 *
*          KF=1 for G(z)                                     *
* Output:  GR --- Real part of Ln[G(z)] or G(z)              *
*          GI --- Imaginary part of Ln[G(z)] or G(z)         *
* Examples:                                                  *
*    x         y           Re[G(z)]           Im[G(z)]       *
*  --------------------------------------------------------  *
*   2.50      5.00     .2267360319E-01    -.1172284404E-01   *
*   5.00     10.00     .1327696517E-01     .3639011746E-02   *
*   2.50     -5.00     .2267360319E-01     .1172284404E-01   *
*   5.00    -10.00     .1327696517E-01    -.3639011746E-02   *
*                                                            *
*    x         y          Re[LnG(z)]         Im[LnG(z)]      *
*  --------------------------------------------------------  *
*   2.50      5.00    -.3668103262E+01     .5806009801E+01   *
*   5.00     10.00    -.4285507444E+01     .1911707090E+02   *
*   2.50     -5.00    -.3668103262E+01    -.5806009801E+01   *
*   5.00    -10.00    -.4285507444E+01    -.1911707090E+02   *
* ---------------------------------------------------------- *
* SAMPLE RUNS:                                               *
*                                                            *
* Please enter KF, x and y: 1 2.5 5.0                        *
*                                                            *
*    x         y           Re[G(z)]           Im[G(z)]       *
*  --------------------------------------------------------  *
*   2.50      5.00       0.0226736032      -0.0117228440     *
*                                                            *
* Please enter KF, zx and zy: 0 2.5 5.0                      *
*                                                            *
*    x         y          Re[LnG(z)]         Im[LnG(z)]      *
*  --------------------------------------------------------  *
*   2.50      5.00      -3.6681032624      5.8060098006      *
*                                                            *
* ---------------------------------------------------------- *
* REFERENCE:                                                 *
*    "Fortran Routines for Computation of Special Functions, *
*     jin.ece.uiuc.edu/routines/routines.html".              *
*                                                            *
*                       Pascal Release By J-P Moreau, Paris. *
*                                (www.jpmoreau.fr)           *
*************************************************************}
program CGammaLaz;

type
  Float = Double; //Extended;

{ Auxiliary functions }

function sinh(x: Float): Float;
var
  expx: Float;
begin
  expx := exp(x);
  Result := 0.5 * (expx - 1.0/expx);
end;

function cosh(x: Float): Float;
var
  expx: Float;
begin
  expx := exp(x);
  Result := 0.5 * (expx + 1.0/expx)
end;

{ Calculate x power n }
function Power(x: Float; n: Integer): Float;
var
  i, m: integer;
begin
  Result := 1.0;
  if n = 0 then
    exit;

  m :=  n;
  if n < 0 then m := -n;
  for i := 1 to m do
    Result := x * Result;
  if n < 0 then
    Result := 1.0 / Result;
end;


{===============================================================================
  Purpose: Compute the gamma function G(z) or Ln[G(z)] for a complex argument
  Input:   x  --- Real part of z
           y  --- Imaginary part of z
           KF --- Function code
                     KF=0 for Ln[G(z)]
                     KF=1 for G(z)
  Output:  GR --- Real part of Ln[G(z)] or G(z)
           GI --- Imaginary part of Ln[G(z)] or G(z)
===============================================================================}
procedure CGamma(X, Y: Float; KF: Integer; out GR, GI:Float);
const
  A: Array[1..10] of Float = (
    8.333333333333333E-02, -2.777777777777778E-03,
    7.936507936507937E-04, -5.952380952380952E-04,
    8.417508417508418E-04, -1.917526917526918E-03,
    6.410256410256410E-03, -2.955065359477124E-02,
    1.796443723688307E-01, -1.39243221690590
  );
var
  G0,GR1,GI1,SR,SI,T,TH,TH1,TH2,X0,X1,Y1,Z1,Z2: Float;
  J,K,NA: Integer;
begin
  if (Y = 0.0) and (X = Int(X)) and (X <= 0.0) then
  begin
    GR := 1E+300;  {arbitrary big number}
    GI := 0.0;
    Exit;
  end;

  if X < 0.0 then
  begin
    X1 := X;
    Y1 := Y;
    X := -X;
    Y := -Y
  end;
  X0 := X;

  if X <= 7.0 then
  begin
    NA := Round(7.0 - X);
    X0 := X + NA
  end;

  Z1 := sqrt(X0*X0 + Y*Y);
  TH := arctan(Y / X0);
  GR := (X0 - 0.5)*Ln(Z1) - TH*Y - X0 + 0.5*Ln(2.0*PI);
  GI := TH*(X0 - 0.5) + Y*Ln(Z1) - Y;
  for K:=1 to 10 do
  begin
    T := Power(Z1, (1-2*K));
    GR := GR + A[K] * T * cos((2.0*K - 1.0) * TH);
    GI := GI - A[K] * T * sin((2.0*K - 1.0) * TH);
  end;

  if X <= 7.0 then
  begin
    GR1 := 0.0;
    GI1 := 0.0;
    for J := 0 to NA - 1 do
    begin
      GR1 := GR1 + 0.5 * Ln(Sqr(X+J) + Y*Y);
      GI1 := GI1 + ArcTan(Y / (X+J));
    end;
    GR := GR - GR1;
    GI := GI - GI1
  end;

  if X1 < 0.0 then
  begin
    Z1 := sqrt(X*X + Y*Y);
    TH1 := ArcTan(Y / X);
    SR := -sin(PI*X) * cosh(PI * Y);
    SI := -cos(PI*X) * sinh(PI * Y);
    Z2 := sqrt(SR*SR + SI*SI);
    TH2 := ArcTan(SI / SR);
    if SR < 0.0 then TH2 := PI + TH2;
    GR := Ln(PI / (Z1*Z2)) - GR;
    GI := -TH1 - TH2 - GI;
    X := X1;
    Y := Y1
  end;

  if KF = 1 then
  begin
    G0 := exp(GR);
    GR := G0 * cos(GI);
    GI := G0 * sin(GI);
  end;
end;

{ Main program }

var
  X, Y, GR, GI: Float;
  KF: Integer;
  w, d: Integer;

begin
  WriteLn;
  Write('  Please enter KF, Re(z) and Im(z): ');
  ReadLn(KF, X, Y);
  WriteLn;
  if KF = 1 then
    WriteLn('      Re[z]     Im[z]            Re[G(z)]                 Im[G(z)]')
  else
  if KF = 0 then
    WriteLn('      Re[z]     Im[z]           Re[LnG(z)]               Im[LnG(z)]')
  else
  begin
    WriteLn('KF is only allowed to be 0 (for LnGamma) or 1 (for Gamma)');
    WriteLn;
    Write('Press ENTER to close...');
    ReadLn;
    Halt;
  end;
  WriteLn  ('    ---------------------------------------------------------------------');

  CGamma(X, Y, KF, GR, GI);

  if SizeOf(Float) = SizeOf(Double) then
  begin
    d := 15;
    w := 25;
  end else
  begin
    d := 18;
    w := 28;
  end;

  WriteLn(' ', X:10:2, Y:10:2, GR:w:d, GI:w:d);
  WriteLn;

  Write('Press ENTER to close...');
  ReadLn;
end.

{end of file mcgama.pas}
