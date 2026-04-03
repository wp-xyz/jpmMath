unit jpmSpecial;
{$mode objfpc}{$H+}

{-------------------------------------------------------------------------------
                              Unit jpmSpecial
     Special mathematical functions: Gamma, Beta, Incomplete Beta/Gamma, Erf.
     Algorithms from Numerical Recipes (Press et al.) and J-P Moreau originals.
-------------------------------------------------------------------------------}

interface

uses
  SysUtils, Math, jpmtypes;

function GammaFunc(x: Float): Float;
function LnGamma(x: Float): Float;
function BetaFunc(a, b: Float): Float;
function LnBeta(a, b: Float): Float;
function IncompleteBeta(x, a, b: Float): Float;
function IncompleteGamma(a, x: Float): Float;
function IncompleteGammaQ(a, x: Float): Float;
function ErfFunc(x: Float): Float;
function ErfcFunc(x: Float): Float;

procedure self_test;

implementation

{ --------------------------------------------------------------------------- }
{ LnGamma — ln(Gamma(x)) via 6-term Lanczos approximation, x > 0             }
{ --------------------------------------------------------------------------- }
function LnGamma(x: Float): Float;
var
  cof: array[1..6] of Float;
  stp, xx, tmp, ser: Float;
  j: integer;
begin
  cof[1] :=  76.18009173;
  cof[2] := -86.50532033;
  cof[3] :=  24.01409822;
  cof[4] :=  -1.231739516;
  cof[5] :=   0.120858003e-2;
  cof[6] :=  -0.536382e-5;
  stp := 2.50662827465;

  xx  := x - 1.0;
  tmp := xx + 5.5;
  tmp := (xx + 0.5) * Ln(tmp) - tmp;
  ser := 1.0;
  for j := 1 to 6 do
  begin
    xx  := xx + 1.0;
    ser := ser + cof[j] / xx
  end;
  result := tmp + Ln(stp * ser)
end;

{ --------------------------------------------------------------------------- }
{ GammaFunc — Gamma(x) for x > 0                                              }
{ --------------------------------------------------------------------------- }
function GammaFunc(x: Float): Float;
begin
  result := Exp(LnGamma(x))
end;

{ --------------------------------------------------------------------------- }
{ LnBeta — ln(Beta(a,b))                                                      }
{ --------------------------------------------------------------------------- }
function LnBeta(a, b: Float): Float;
begin
  result := LnGamma(a) + LnGamma(b) - LnGamma(a + b)
end;

{ --------------------------------------------------------------------------- }
{ BetaFunc — Beta(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b)                         }
{ --------------------------------------------------------------------------- }
function BetaFunc(a, b: Float): Float;
begin
  result := Exp(LnBeta(a, b))
end;

{ --------------------------------------------------------------------------- }
{ BetaCF — continued fraction for incomplete beta (used by IncompleteBeta)    }
{ Lentz-style modified Steed algorithm; replaces original goto-based BETACF.  }
{ --------------------------------------------------------------------------- }
function BetaCF(a, b, x: Float): Float;
const
  ITMAX = 200;
  EPS   = 3e-7;
  FPMIN = 1e-300;
var
  am, bm, az, qab, qap, qam, bz, em, tem, d, ap, bp, app, bpp, aold: Float;
  m: integer;
  done: boolean;
begin
  am  := 1.0;
  bm  := 1.0;
  az  := 1.0;
  qab := a + b;
  qap := a + 1.0;
  qam := a - 1.0;
  bz  := 1.0 - qab * x / qap;
  done := false;
  m    := 1;
  while (m <= ITMAX) and not done do
  begin
    em  := m;
    tem := em + em;
    d   := em * (b - em) * x / ((qam + tem) * (a + tem));
    ap  := az + d * am;
    bp  := bz + d * bm;
    d   := -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem));
    app := ap + d * az;
    bpp := bp + d * bz;
    aold := az;
    if Abs(bpp) < FPMIN then bpp := FPMIN;
    am  := ap  / bpp;
    bm  := bp  / bpp;
    az  := app / bpp;
    bz  := 1.0;
    if Abs(az - aold) < EPS * Abs(az) then
      done := true;
    m := m + 1
  end;
  result := az
end;

{ --------------------------------------------------------------------------- }
{ IncompleteBeta — regularized incomplete beta I_x(a,b), result in [0,1]      }
{ --------------------------------------------------------------------------- }
function IncompleteBeta(x, a, b: Float): Float;
var
  bt: Float;
begin
  if (x < 0.0) or (x > 1.0) then
  begin
    result := 0.0;
    exit
  end;

  if (x = 0.0) or (x = 1.0) then
    bt := 0.0
  else
    bt := Exp(LnGamma(a + b) - LnGamma(a) - LnGamma(b)
              + a * Ln(x) + b * Ln(1.0 - x));

  if x < (a + 1.0) / (a + b + 2.0) then
    result := bt * BetaCF(a, b, x) / a
  else
    result := 1.0 - bt * BetaCF(b, a, 1.0 - x) / b
end;

{ --------------------------------------------------------------------------- }
{ GammaSeries — series expansion for P(a,x); used when x < a+1               }
{ --------------------------------------------------------------------------- }
function GammaSeries(a, x: Float): Float;
const
  ITMAX = 200;
  EPS   = 3e-7;
var
  ap, del, sum: Float;
  n: integer;
begin
  if x <= 0.0 then
  begin
    result := 0.0;
    exit
  end;
  ap  := a;
  del := 1.0 / a;
  sum := del;
  n   := 1;
  while n <= ITMAX do
  begin
    ap  := ap + 1.0;
    del := del * x / ap;
    sum := sum + del;
    if Abs(del) < Abs(sum) * EPS then
    begin
      result := sum * Exp(-x + a * Ln(x) - LnGamma(a));
      exit
    end;
    n := n + 1
  end;
  result := sum * Exp(-x + a * Ln(x) - LnGamma(a))
end;

{ --------------------------------------------------------------------------- }
{ GammaCF — continued fraction for Q(a,x) via Lentz method; used when x>=a+1 }
{ --------------------------------------------------------------------------- }
function GammaCF(a, x: Float): Float;
const
  ITMAX = 200;
  EPS   = 3e-7;
  FPMIN = 1e-300;
var
  b, c, d, h, del, an: Float;
  i: integer;
begin
  b := x + 1.0 - a;
  c := 1.0 / FPMIN;
  d := 1.0 / b;
  if Abs(d) < FPMIN then d := FPMIN;
  h := d;
  i := 1;
  while i <= ITMAX do
  begin
    an := -i * (i - a);
    b  := b + 2.0;
    d  := an * d + b;
    if Abs(d) < FPMIN then d := FPMIN;
    c  := b + an / c;
    if Abs(c) < FPMIN then c := FPMIN;
    d   := 1.0 / d;
    del := d * c;
    h   := h * del;
    if Abs(del - 1.0) < EPS then
    begin
      result := Exp(-x + a * Ln(x) - LnGamma(a)) * h;
      exit
    end;
    i := i + 1
  end;
  result := Exp(-x + a * Ln(x) - LnGamma(a)) * h
end;

{ --------------------------------------------------------------------------- }
{ IncompleteGamma — regularized lower incomplete gamma P(a,x)                 }
{ --------------------------------------------------------------------------- }
function IncompleteGamma(a, x: Float): Float;
begin
  if x < 0.0 then
  begin
    result := 0.0;
    exit
  end;
  if x = 0.0 then
  begin
    result := 0.0;
    exit
  end;
  if x < a + 1.0 then
    result := GammaSeries(a, x)
  else
    result := 1.0 - GammaCF(a, x)
end;

{ --------------------------------------------------------------------------- }
{ IncompleteGammaQ — upper incomplete gamma Q(a,x) = 1 - P(a,x)              }
{ --------------------------------------------------------------------------- }
function IncompleteGammaQ(a, x: Float): Float;
begin
  result := 1.0 - IncompleteGamma(a, x)
end;

{ --------------------------------------------------------------------------- }
{ ErfcFunc — complementary error function erfc(x) = 1 - erf(x)               }
{ Chebyshev rational approximation (Numerical Recipes), max error ~1.2e-7     }
{ --------------------------------------------------------------------------- }
function ErfcFunc(x: Float): Float;
var
  t, z, ans: Float;
begin
  z := Abs(x);
  t := 1.0 / (1.0 + 0.5 * z);
  ans := t * Exp(-z * z - 1.26551223
         + t * ( 1.00002368
         + t * ( 0.37409196
         + t * ( 0.09678418
         + t * (-0.18628806
         + t * ( 0.27886807
         + t * (-1.13520398
         + t * ( 1.48851587
         + t * (-0.82215223
         + t *   0.17087294)))))))));
  if x >= 0.0 then
    result := ans
  else
    result := 2.0 - ans
end;

{ --------------------------------------------------------------------------- }
{ ErfFunc — error function erf(x)                                             }
{ --------------------------------------------------------------------------- }
function ErfFunc(x: Float): Float;
begin
  if x >= 0.0 then
    result := 1.0 - ErfcFunc(x)
  else
    result := ErfcFunc(-x) - 1.0
end;

{ --------------------------------------------------------------------------- }
{ self_test — verifies all functions against known values                      }
{ --------------------------------------------------------------------------- }
procedure self_test;
const
  TOL = 1e-4;

  procedure check(const name: string; computed, expected: Float);
  var
    ok: string;
  begin
    if Abs(computed - expected) < TOL then
      ok := 'OK'
    else
      ok := '*** FAIL ***';
    WriteLn(Format('  %-30s computed=%10.6f  expected=%10.6f  %s',
                   [name, computed, expected, ok]))
  end;

begin
  WriteLn('=== jpmSpecial self_test ===');
  WriteLn;

  { Gamma }
  check('GammaFunc(0.5)',  GammaFunc(0.5),   1.772454);
  check('GammaFunc(1.0)',  GammaFunc(1.0),   1.000000);
  check('GammaFunc(5.0)',  GammaFunc(5.0),  24.000000);

  { LnGamma }
  check('LnGamma(10.0)',   LnGamma(10.0),   12.801827);

  { Beta }
  check('BetaFunc(2,3)',   BetaFunc(2, 3),   0.083333);

  { IncompleteBeta }
  check('IncompleteBeta(0.5,2,3)', IncompleteBeta(0.5, 2, 3), 0.6875);

  { IncompleteGamma }
  check('IncompleteGamma(1,1)',    IncompleteGamma(1.0, 1.0),  0.632121);

  { Erf / Erfc }
  check('ErfFunc(1.0)',   ErfFunc(1.0),   0.842701);
  check('ErfcFunc(1.0)',  ErfcFunc(1.0),  0.157299);

  WriteLn;
  WriteLn('=== done ===');
end;

end.
