unit jpmBessel;
{$mode objfpc}{$H+}

{-------------------------------------------------------------------------------
                               Unit jpmBessel
     Bessel functions of real argument: J0, J1, Y0, Y1, I0, I1, K0, K1,
     Jn (integer order), and zeros of J0/J1.
     Polynomial approximations from Abramowitz & Stegun and Numerical Recipes.
-------------------------------------------------------------------------------}

interface

uses
  SysUtils, Math, jpmtypes;

function BesselJ0(x: Float): Float;
function BesselJ1(x: Float): Float;
function BesselY0(x: Float): Float;
function BesselY1(x: Float): Float;
function BesselI0(x: Float): Float;
function BesselI1(x: Float): Float;
function BesselK0(x: Float): Float;
function BesselK1(x: Float): Float;
function BesselJn(n: integer; x: Float): Float;
function BesselJ0Zero(k: integer): Float;
function BesselJ1Zero(k: integer): Float;

procedure self_test;

implementation

{ --------------------------------------------------------------------------- }
{ BesselJ0 — J_0(x), Bessel function first kind order 0                       }
{ A&S 9.4.1 / 9.4.3, Numerical Recipes approach, error < 1.6e-7              }
{ --------------------------------------------------------------------------- }
function BesselJ0(x: Float): Float;
var
  ax, z, xx, y, ans1, ans2: Float;
begin
  ax := Abs(x);
  if ax < 8.0 then
  begin
    y    := x * x;
    ans1 := 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
            + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
    ans2 := 57568490411.0 + y * (1029532985.0 + y * (9494680.718
            + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
    result := ans1 / ans2
  end
  else
  begin
    z    := 8.0 / ax;
    y    := z * z;
    xx   := ax - 0.785398164;
    ans1 := 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
            + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
    ans2 := -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5
            + y * (0.7621095161e-6 + y * (-0.934935152e-7))));
    result := Sqrt(0.636619772 / ax) * (Cos(xx) * ans1 - z * Sin(xx) * ans2)
  end
end;

{ --------------------------------------------------------------------------- }
{ BesselJ1 — J_1(x), Bessel function first kind order 1                       }
{ A&S 9.4.2 / 9.4.4, error < 1.3e-7                                          }
{ --------------------------------------------------------------------------- }
function BesselJ1(x: Float): Float;
var
  ax, z, xx, y, ans1, ans2: Float;
begin
  ax := Abs(x);
  if ax < 8.0 then
  begin
    y    := x * x;
    ans1 := x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
            + y * (-2972611.439 + y * (15704.48260 + y * (-30.16116360))))));
    ans2 := 144725228442.0 + y * (2300535178.0 + y * (18583304.74
            + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
    result := ans1 / ans2
  end
  else
  begin
    z    := 8.0 / ax;
    y    := z * z;
    xx   := ax - 2.356194491;
    ans1 := 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
            + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
    ans2 := 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5
            + y * (-0.88228987e-6 + y * 0.105787412e-6)));
    result := Sqrt(0.636619772 / ax) * (Cos(xx) * ans1 - z * Sin(xx) * ans2);
    if x < 0.0 then
      result := -result
  end
end;

{ --------------------------------------------------------------------------- }
{ BesselY0 — Y_0(x), Bessel function second kind order 0, x > 0              }
{ A&S 9.4.1 / 9.4.3, error < 1.4e-7                                          }
{ --------------------------------------------------------------------------- }
function BesselY0(x: Float): Float;
var
  z, xx, y, ans1, ans2: Float;
begin
  if x < 8.0 then
  begin
    y    := x * x;
    ans1 := -2957821389.0 + y * (7062834065.0 + y * (-512359803.6
            + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));
    ans2 := 40076544269.0 + y * (745249964.8 + y * (7189466.438
            + y * (47447.26470 + y * (226.1030244 + y * 1.0))));
    result := (ans1 / ans2) + 0.636619772 * BesselJ0(x) * Ln(x)
  end
  else
  begin
    z    := 8.0 / x;
    y    := z * z;
    xx   := x - 0.785398164;
    ans1 := 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
            + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
    ans2 := -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5
            + y * (0.7621095161e-6 + y * (-0.934935152e-7))));
    result := Sqrt(0.636619772 / x) * (Sin(xx) * ans1 + z * Cos(xx) * ans2)
  end
end;

{ --------------------------------------------------------------------------- }
{ BesselY1 — Y_1(x), Bessel function second kind order 1, x > 0              }
{ A&S 9.4.2 / 9.4.4, error < 1.4e-7                                          }
{ --------------------------------------------------------------------------- }
function BesselY1(x: Float): Float;
var
  z, xx, y, ans1, ans2: Float;
begin
  if x < 8.0 then
  begin
    y    := x * x;
    ans1 := -0.4900604943e13 + y * (0.1275274390e13 + y * (-0.5153438139e11
            + y * (0.7349264551e9 + y * (-0.4237922726e7 + y * 0.8511937935e4))));
    ans2 := 0.2499580570e14 + y * (0.4244419664e12 + y * (0.3733650367e10
            + y * (0.2245904002e8 + y * (0.1020426050e6 + y * (0.3549632885e3 + y)))));
    result := x * (ans1 / ans2) + 0.636619772 * (BesselJ1(x) * Ln(x) - 1.0 / x)
  end
  else
  begin
    z    := 8.0 / x;
    y    := z * z;
    xx   := x - 2.356194491;
    ans1 := 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
            + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
    ans2 := 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5
            + y * (-0.88228987e-6 + y * 0.105787412e-6)));
    result := Sqrt(0.636619772 / x) * (Sin(xx) * ans1 + z * Cos(xx) * ans2)
  end
end;

{ --------------------------------------------------------------------------- }
{ BesselI0 — I_0(x), modified Bessel function first kind order 0              }
{ A&S 9.8.1 / 9.8.2  (Numerical Recipes), error < 1.6e-7 for |x|<=3.75       }
{                                          error < 1.9e-7 for |x|>3.75        }
{ --------------------------------------------------------------------------- }
function BesselI0(x: Float): Float;
var
  ax, y, ans: Float;
begin
  ax := Abs(x);
  if ax <= 3.75 then
  begin
    y   := (x / 3.75);
    y   := y * y;
    ans := 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
           + y * (0.2659732 + y * (0.0360768 + y * 0.0045813)))));
    result := ans
  end
  else
  begin
    y   := 3.75 / ax;
    ans := (Exp(ax) / Sqrt(ax))
           * (0.39894228 + y * (0.01328592 + y * (0.00225319
           + y * (-0.00157565 + y * (0.00916281 + y * (-0.02057706
           + y * (0.02635537 + y * (-0.01647633 + y * 0.00392377))))))));
    result := ans
  end
end;

{ --------------------------------------------------------------------------- }
{ BesselI1 — I_1(x), modified Bessel function first kind order 1              }
{ A&S 9.8.3 / 9.8.4, error < 8e-9                                             }
{ --------------------------------------------------------------------------- }
function BesselI1(x: Float): Float;
var
  ax, y, ans: Float;
begin
  ax := Abs(x);
  if ax <= 3.75 then
  begin
    y   := (x / 3.75);
    y   := y * y;
    ans := ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
           + y * (0.02658733 + y * (0.00301532 + y * 0.00032411))))));
    result := ans
  end
  else
  begin
    y   := 3.75 / ax;
    ans := (Exp(ax) / Sqrt(ax))
           * (0.39894228 + y * (-0.03988024 + y * (-0.00362018
           + y * (0.00163801 + y * (-0.01031555 + y * (0.02282967
           + y * (-0.02895312 + y * (0.01787654 + y * (-0.00420059)))))))));
    if x < 0.0 then
      result := -ans
    else
      result := ans
  end
end;

{ --------------------------------------------------------------------------- }
{ BesselK0 — K_0(x), modified Bessel function second kind order 0, x > 0     }
{ A&S 9.8.5 / 9.8.6, error < 1e-8                                             }
{ --------------------------------------------------------------------------- }
function BesselK0(x: Float): Float;
var
  y, ans: Float;
begin
  if x <= 2.0 then
  begin
    y   := x * x / 4.0;
    ans := (-Ln(x / 2.0) * BesselI0(x))
           + (-0.57721566 + y * (0.42278420 + y * (0.23069756
           + y * (0.03488590 + y * (0.00262698 + y * (0.00010750 + y * 0.0000074))))));
    result := ans
  end
  else
  begin
    y   := 2.0 / x;
    ans := (Exp(-x) / Sqrt(x))
           * (1.25331414 + y * (-0.07832358 + y * (0.02189568
           + y * (-0.01062446 + y * (0.00587872 + y * (-0.00251540 + y * 0.00053208))))));
    result := ans
  end
end;

{ --------------------------------------------------------------------------- }
{ BesselK1 — K_1(x), modified Bessel function second kind order 1, x > 0     }
{ A&S 9.8.7 / 9.8.8, error < 8e-9                                             }
{ --------------------------------------------------------------------------- }
function BesselK1(x: Float): Float;
var
  y, ans: Float;
begin
  if x <= 2.0 then
  begin
    y   := x * x / 4.0;
    ans := (Ln(x / 2.0) * BesselI1(x))
           + (1.0 / x) * (1.0 + y * (0.15443144 + y * (-0.67278579
           + y * (-0.18156897 + y * (-0.01919402 + y * (-0.00110404 + y * (-0.00004686)))))));
    result := ans
  end
  else
  begin
    y   := 2.0 / x;
    ans := (Exp(-x) / Sqrt(x))
           * (1.25331414 + y * (0.23498619 + y * (-0.03655620
           + y * (0.01504268 + y * (-0.00780353 + y * (0.00325614 + y * (-0.00068245)))))));
    result := ans
  end
end;

{ --------------------------------------------------------------------------- }
{ BesselJn — J_n(x) for integer n >= 0, via Miller backward recurrence        }
{ Stable downward recurrence from large m down to 0, then normalise with J0.  }
{ --------------------------------------------------------------------------- }
function BesselJn(n: integer; x: Float): Float;
const
  ACC  = 160.0;
  BIGNO  = 1.0e10;
  BIGNI  = 1.0e-10;
var
  j, jsum, m: integer;
  ax, bj, bjm, bjp, sum, tox, ans: Float;
  sign: boolean;
begin
  ax := Abs(x);
  if n = 0 then
  begin
    result := BesselJ0(x);
    exit
  end;
  if n = 1 then
  begin
    result := BesselJ1(x);
    exit
  end;
  if ax = 0.0 then
  begin
    result := 0.0;
    exit
  end;

  if ax > Float(n) then
  begin
    { Forward recurrence — safe when x > n }
    tox := 2.0 / ax;
    bjm := BesselJ0(ax);
    bj  := BesselJ1(ax);
    j   := 1;
    while j < n do
    begin
      bjp := Float(j) * tox * bj - bjm;
      bjm := bj;
      bj  := bjp;
      j   := j + 1
    end;
    ans := bj
  end
  else
  begin
    { Miller backward recurrence }
    tox  := 2.0 / ax;
    m    := 2 * ((n + Trunc(Sqrt(ACC * Float(n)))) div 2);
    jsum := 0;
    bjp  := 0.0;
    sum  := 0.0;
    ans  := 0.0;
    bj   := 1.0;
    j    := m;
    while j > 0 do
    begin
      bjm := Float(j) * tox * bj - bjp;
      bjp := bj;
      bj  := bjm;
      if Abs(bj) > BIGNO then
      begin
        bj  := bj  * BIGNI;
        bjp := bjp * BIGNI;
        ans := ans * BIGNI;
        sum := sum * BIGNI
      end;
      if jsum <> 0 then
        sum := sum + bj;
      jsum := 1 - jsum;
      if j = n then
        ans := bjp;
      j := j - 1
    end;
    sum  := 2.0 * sum - bj;
    ans  := ans / sum
  end;

  sign := (x < 0.0) and ((n mod 2) = 1);
  if sign then
    result := -ans
  else
    result := ans
end;

{ --------------------------------------------------------------------------- }
{ BesselJ0Zero — k-th positive zero of J0, k=1,2,...                          }
{ McMahon: j_0k ~ beta + 1/(8*beta) - 31*4/(3*(8*beta)^3) + ...              }
{ then one Newton step for higher accuracy                                     }
{ --------------------------------------------------------------------------- }
function BesselJ0Zero(k: integer): Float;
const
  PI = 3.14159265358979323846;
var
  beta, b8, r: Float;
begin
  beta := (Float(k) - 0.25) * PI;
  b8   := 8.0 * beta;
  { For J0: mu=0, correction = +(1/b8)*(1 - 31/(6*b8^2)*(...)) }
  r    := beta + (1.0 / b8) * (1.0 - (31.0 / (6.0 * b8 * b8))
          * (1.0 - (3779.0 / (40.0 * b8 * b8))));
  { Newton step: r <- r - J0(r)/J0'(r) = r - J0(r)/(-J1(r)) = r + J0(r)/J1(r) }
  r    := r - BesselJ0(r) / (-BesselJ1(r));
  result := r
end;

{ --------------------------------------------------------------------------- }
{ BesselJ1Zero — k-th positive zero of J1, k=1,2,...                          }
{ McMahon: j_1k ~ beta - 3/(8*beta) - ...                                     }
{ then one Newton step                                                         }
{ --------------------------------------------------------------------------- }
function BesselJ1Zero(k: integer): Float;
const
  PI = 3.14159265358979323846;
var
  beta, b8, r: Float;
begin
  beta := (Float(k) + 0.25) * PI;
  b8   := 8.0 * beta;
  r    := beta - (3.0 / b8) * (1.0 - (11.0 / (4.0 * b8 * b8))
          * (1.0 - (147.0 / (8.0 * b8 * b8))));
  { Newton step using J1'(r) = J0(r) - J1(r)/r }
  r    := r - BesselJ1(r) / (BesselJ0(r) - BesselJ1(r) / r);
  result := r
end;

{ --------------------------------------------------------------------------- }
{ self_test — verify against known mathematical table values                   }
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
    WriteLn(Format('  %-32s computed=%12.6f  expected=%12.6f  %s',
                   [name, computed, expected, ok]));
    if ok = '*** FAIL ***' then
      SelfTestFail(name + ': computed=' + FloatToStr(computed) + ' expected=' + FloatToStr(expected));
  end;

begin
  WriteLn('=== jpmBessel self_test ===');
  WriteLn;

  { J0 }
  check('BesselJ0(0.0)',        BesselJ0(0.0),      1.000000);
  check('BesselJ0(2.4048)',     BesselJ0(2.4048),   0.000000);
  check('BesselJ0(5.0)',        BesselJ0(5.0),     -0.177597);

  { J1 }
  check('BesselJ1(0.0)',        BesselJ1(0.0),      0.000000);
  check('BesselJ1(3.8317)',     BesselJ1(3.8317),   0.000000);
  check('BesselJ1(5.0)',        BesselJ1(5.0),     -0.327579);

  { Y0, Y1 — just smoke-test finite values }
  check('BesselY0(1.0)',        BesselY0(1.0),      0.088257);
  check('BesselY1(2.0)',        BesselY1(2.0),     -0.107032);

  { I0 }
  check('BesselI0(0.0)',        BesselI0(0.0),      1.000000);
  check('BesselI0(2.0)',        BesselI0(2.0),      2.279585);

  { I1 }
  check('BesselI1(0.0)',        BesselI1(0.0),      0.000000);
  check('BesselI1(2.0)',        BesselI1(2.0),      1.590636);

  { K0, K1 — smoke-test finite values }
  check('BesselK0(1.0)',        BesselK0(1.0),      0.421024);
  check('BesselK1(1.0)',        BesselK1(1.0),      0.601907);

  { Jn }
  check('BesselJn(2,5.0)',      BesselJn(2, 5.0),   0.046565);

  { Zeros }
  check('BesselJ0Zero(1)',      BesselJ0Zero(1),    2.404826);
  check('BesselJ1Zero(1)',      BesselJ1Zero(1),    3.831706);

  WriteLn;
  WriteLn('=== done ===');
end;

end.
