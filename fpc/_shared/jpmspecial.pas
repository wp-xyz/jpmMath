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

{ CAT 9 — Orthogonal polynomials and special functions }
function LegendrePn(n: integer; x: Float): Float;   { P_n(x) }
function LegendreQn(n: integer; x: Float): Float;   { Q_n(x), |x|<1 }
function HermiteHn(n: integer; x: Float): Float;    { physicists H_n(x) }
function LaguerreLn(n: integer; x: Float): Float;   { L_n(x) }
function AiryAi(x: Float): Float;                   { Airy Ai(x) }
function AiryBi(x: Float): Float;                   { Airy Bi(x) }
function EllipticK(k: Float): Float;                { complete elliptic integral K(k) }
function EllipticE(k: Float): Float;                { complete elliptic integral E(k) }
function Hypergeometric2F1(a, b, c, x: Float): Float; { Gauss 2F1, |x|<1 }

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
                   [name, computed, expected, ok]));
    if ok = '*** FAIL ***' then
      SelfTestFail(name + ': computed=' + FloatToStr(computed) + ' expected=' + FloatToStr(expected));
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

  { CAT 9 — Orthogonal polynomials }
  check('LegendrePn(2,0.5)',    LegendrePn(2, 0.5),   -0.125000);
  check('LegendrePn(3,0.5)',    LegendrePn(3, 0.5),   -0.437500);
  check('LegendreQn(1,0.5)',    LegendreQn(1, 0.5),   -0.725347);
  check('HermiteHn(3,1.0)',     HermiteHn(3, 1.0),    -4.000000);
  check('LaguerreLn(3,1.0)',    LaguerreLn(3, 1.0),   -0.666667);
  check('AiryAi(0.0)',          AiryAi(0.0),            0.355028);
  check('AiryBi(0.0)',          AiryBi(0.0),            0.614927);
  check('EllipticK(0.0)',       EllipticK(0.0),         1.570796);
  check('EllipticE(0.0)',       EllipticE(0.0),         1.570796);
  check('EllipticK(0.5)',       EllipticK(0.5),         1.685750);
  check('Hypergeometric2F1(0.5,0.5,1,0.25)', Hypergeometric2F1(0.5,0.5,1,0.25), 1.073182);

  WriteLn;
  WriteLn('=== done ===');
end;


{ ------------------------------------------------------------------ }
{ CAT 9 — Orthogonal polynomials and special functions               }
{ ------------------------------------------------------------------ }

function LegendrePn(n: integer; x: Float): Float;
var
  k: integer;
  p0, p1, p2: Float;
begin
  if n = 0 then begin result := 1; exit end;
  if n = 1 then begin result := x; exit end;
  p0 := 1; p1 := x;
  for k := 2 to n do
  begin
    p2 := ((2*k - 1)*x*p1 - (k - 1)*p0) / k;
    p0 := p1; p1 := p2
  end;
  result := p1
end;

function LegendreQn(n: integer; x: Float): Float;
var
  k: integer;
  q0, q1, q2: Float;
begin
  { Q0(x) = 0.5*ln((1+x)/(1-x)), Q1(x) = x*Q0(x)-1 }
  if Abs(x) >= 1 then begin result := 0; exit end;
  q0 := 0.5 * Ln((1 + x) / (1 - x));
  if n = 0 then begin result := q0; exit end;
  q1 := x * q0 - 1;
  if n = 1 then begin result := q1; exit end;
  for k := 2 to n do
  begin
    q2 := ((2*k - 1)*x*q1 - (k - 1)*q0) / k;
    q0 := q1; q1 := q2
  end;
  result := q1
end;

function HermiteHn(n: integer; x: Float): Float;
var
  k: integer;
  h0, h1, h2: Float;
begin
  { physicists: H_{n+1} = 2x H_n - 2n H_{n-1} }
  if n = 0 then begin result := 1; exit end;
  if n = 1 then begin result := 2*x; exit end;
  h0 := 1; h1 := 2*x;
  for k := 1 to n - 1 do
  begin
    h2 := 2*x*h1 - 2*k*h0;
    h0 := h1; h1 := h2
  end;
  result := h1
end;

function LaguerreLn(n: integer; x: Float): Float;
var
  k: integer;
  l0, l1, l2: Float;
begin
  if n = 0 then begin result := 1; exit end;
  if n = 1 then begin result := 1 - x; exit end;
  l0 := 1; l1 := 1 - x;
  for k := 1 to n - 1 do
  begin
    l2 := ((2*k + 1 - x)*l1 - k*l0) / (k + 1);
    l0 := l1; l1 := l2
  end;
  result := l1
end;

function AiryAi(x: Float): Float;
{ Power series for small |x|; asymptotic for large x }
var
  z, s, t: Float;
  k: integer;
begin
  if x > 4 then
  begin
    z := (2/3) * Power(x, 1.5);
    result := Exp(-z) / (2 * Sqrt(Pi) * Power(x, 0.25))
  end
  else if x < -4 then
  begin
    z := (2/3) * Power(-x, 1.5);
    result := (Sin(z + Pi/4)) / (Sqrt(Pi) * Power(-x, 0.25))
  end
  else
  begin
    { power series: Ai(x) = c1*f(x) - c2*g(x)
      f = sum x^{3k}/(3k)!*prod(3j-2), g = sum x^{3k+1}/(3k+1)!*prod(3j-1) }
    s := 0; t := 1;
    for k := 0 to 30 do
    begin
      s := s + t;
      t := t * x * x * x / ((3*k+2)*(3*k+3));
      if Abs(t) < 1e-15 * Abs(s) then break
    end;
    { c1 = Ai(0) = 1/(3^(2/3)*Gamma(2/3)), approximate = 0.3550280539 }
    result := 0.3550280539 * s;
    s := 0; t := x;
    for k := 0 to 30 do
    begin
      s := s + t;
      t := t * x * x * x / ((3*k+3)*(3*k+4));
      if Abs(t) < 1e-15 * Abs(s) then break
    end;
    { c2 = -Ai'(0) = 1/(3^(1/3)*Gamma(1/3)), approximate = 0.2588194037 }
    result := result - 0.2588194037 * s
  end
end;

function AiryBi(x: Float): Float;
var
  z, s, t: Float;
  k: integer;
begin
  if x > 4 then
  begin
    z := (2/3) * Power(x, 1.5);
    result := Exp(z) / (Sqrt(Pi) * Power(x, 0.25))
  end
  else if x < -4 then
  begin
    z := (2/3) * Power(-x, 1.5);
    result := Cos(z + Pi/4) / (Sqrt(Pi) * Power(-x, 0.25))
  end
  else
  begin
    s := 0; t := 1;
    for k := 0 to 30 do
    begin
      s := s + t;
      t := t * x * x * x / ((3*k+2)*(3*k+3));
      if Abs(t) < 1e-15 * Abs(s) then break
    end;
    result := 0.6149266274 * s;  { Bi(0) = 3^(1/6)/Gamma(2/3) }
    s := 0; t := x;
    for k := 0 to 30 do
    begin
      s := s + t;
      t := t * x * x * x / ((3*k+3)*(3*k+4));
      if Abs(t) < 1e-15 * Abs(s) then break
    end;
    result := result + 0.4482883574 * s  { Bi'(0) = 3^(5/6)/Gamma(1/3) }
  end
end;

function EllipticK(k: Float): Float;
{ AGM method: K(k) = Pi / (2 * AGM(1, sqrt(1-k^2))) }
var
  a, b, tmp: Float;
  i: integer;
begin
  if Abs(k) >= 1 then begin result := 1e308; exit end;
  a := 1;
  b := Sqrt(1 - k*k);
  for i := 1 to 30 do
  begin
    tmp := (a + b) / 2;
    b   := Sqrt(a * b);
    a   := tmp;
    if Abs(a - b) < 1e-14 * a then break
  end;
  result := Pi / (2 * a)
end;

function EllipticE(k: Float): Float;
{ AGM-based algorithm for E(k) }
var
  a, b, c, tmp, pwr: Float;
  i: integer;
begin
  if Abs(k) >= 1 then begin result := 1; exit end;
  a := 1; b := Sqrt(1 - k*k);
  c := k*k / 2;
  pwr := 1;
  for i := 1 to 30 do
  begin
    tmp := (a + b) / 2;
    b   := Sqrt(a * b);
    pwr := pwr * 2;
    c   := c - pwr * (tmp - b) * (tmp - b);
    a   := tmp;
    if Abs(a - b) < 1e-14 * a then break
  end;
  result := Pi * (1 - c) / (2 * a)
end;

function Hypergeometric2F1(a, b, c, x: Float): Float;
{ Gauss hypergeometric series 2F1(a,b;c;x), |x| < 1 }
var
  s, t: Float;
  n: integer;
begin
  if Abs(x) >= 1 then begin result := 0; exit end;
  s := 1; t := 1;
  for n := 0 to 99 do
  begin
    t := t * (a + n) * (b + n) * x / ((c + n) * (n + 1));
    s := s + t;
    if Abs(t) < 1e-14 * Abs(s) then break
  end;
  result := s
end;

end.
