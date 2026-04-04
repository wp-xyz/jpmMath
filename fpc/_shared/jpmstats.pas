{-------------------------------------------------------------------------------
                                 Unit jpmStats                        
         Statistical function based on Jean-Pierre Moraus's code. 

This unit collects generally available statistical calculations from
Jean-Pierre Moraus's sample projects.
-------------------------------------------------------------------------------}

unit jpmStats;

{$mode ObjFPC}{$H+}

interface

uses
  SysUtils, Math, jpmTypes, jpmSpecial;

{ --- existing --- }
function NormalDist(u: Float): Float;
function InvNormalDist(P: Float): Float;

{ --- CAT 1.1: Chi-Square --- }
function ChiSquareDist(chi2, df: Float): Float;
function InvChiSquare(p, df: Float): Float;

{ --- CAT 1.2: F-Distribution --- }
function FDist(f, df1, df2: Float): Float;
function InvFDist(p, df1, df2: Float): Float;

{ --- CAT 1.4: Discrete Distributions --- }
function BinomialProb(n, k: integer; p: Float): Float;
function BinomialCDF(n, k: integer; p: Float): Float;
function PoissonProb(k: integer; lambda: Float): Float;
function PoissonCDF(k: integer; lambda: Float): Float;

{ --- CAT 1.5: Moments --- }
function StatMean(var data: TFloatArray; n: integer): Float;
function StatVariance(var data: TFloatArray; n: integer): Float;
function StatStdDev(var data: TFloatArray; n: integer): Float;
function StatSkewness(var data: TFloatArray; n: integer): Float;
function StatKurtosis(var data: TFloatArray; n: integer): Float;
function StatMedian(var data: TFloatArray; n: integer): Float;

procedure self_test;

implementation

{ Standard Normal Probability Function (Gaussian curve)
  Taken from the "normal.pas" sample project where it is named "Phi()". }
function NormalDist(u: Float): Float;
begin
  Result := 1.0/Sqrt(2.0*pi) * exp(-0.5*u*u)
end;

{ Inverse Standard Normal Function 
  Taken from the "normal.pas" sample project where it is named "Normal()". }
function InvNormalDist(P: Float): Float;
var
  chimax, chisqval, epsilon, minchisq, maxchisq: Float;
begin
  epsilon := 1e-6;
  chimax := 1e6;
  minchisq := 0.0;
  maxchisq := chimax;

  if P <= 0.0 then
  begin
    Result := maxchisq;
    exit;
  end;

  if P >= 1.0 then
  begin
    Result := 0.0;
    exit;
  end;

  chisqval := 0.5;
  while abs(maxchisq - minchisq) > epsilon do
  begin
    IF NormalDist(chisqval) < P then
      maxchisq := chisqval
    else
      minchisq := chisqval;
    chisqval := (maxchisq + minchisq) * 0.5
  end;
  Result := chisqval;
end;

(* Alternative (direct) method of calculation
function InvNormalDist(P: Float): Float;
var
  a: Double;
begin
  if P <= 0 then
    Result := NaN
  else
  if P > 37 then  // to be more exact: replace this by InvNormalDist(MIN_FLOAT)
    Result := 0
  else
  begin
    a := -2.0 * (ln(sqrt(2.0*pi)) + ln(P));
    if a >= 0 then
      Result := sqrt(a)
    else
      Result := 0.0; //NaN;
  end;
end;               *)

{ ---------------------------------------------------------------------------
  CAT 1.1 — Chi-Square Distribution
  --------------------------------------------------------------------------- }

{ CDF of chi-square distribution: P(X <= chi2) with df degrees of freedom.
  Uses the regularized lower incomplete gamma: P(df/2, chi2/2). }
function ChiSquareDist(chi2, df: Float): Float;
begin
  if chi2 <= 0.0 then
  begin
    result := 0.0;
    exit
  end;
  result := IncompleteGamma(df / 2.0, chi2 / 2.0)
end;

{ Inverse chi-square: find chi2 such that ChiSquareDist(chi2, df) = p (bisection). }
function InvChiSquare(p, df: Float): Float;
const
  EPSILON = 1e-7;
  MAXITER = 200;
var
  lo, hi, mid, fmid: Float;
  iter: integer;
begin
  if p <= 0.0 then
  begin
    result := 0.0;
    exit
  end;
  if p >= 1.0 then
  begin
    result := 1e15;
    exit
  end;
  lo := 0.0;
  hi := 1000.0;
  { expand hi until CDF(hi) >= p }
  while ChiSquareDist(hi, df) < p do
    hi := hi * 2.0;
  iter := 0;
  mid  := (lo + hi) / 2.0;
  fmid := ChiSquareDist(mid, df);
  while (Abs(hi - lo) > EPSILON) and (iter < MAXITER) do
  begin
    if fmid < p then
      lo := mid
    else
      hi := mid;
    mid  := (lo + hi) / 2.0;
    fmid := ChiSquareDist(mid, df);
    iter := iter + 1
  end;
  result := mid
end;

{ ---------------------------------------------------------------------------
  CAT 1.2 — F-Distribution
  --------------------------------------------------------------------------- }

{ CDF of F-distribution: P(X <= f) with df1, df2 degrees of freedom.
  Uses complement form of the regularized incomplete beta. }
function FDist(f, df1, df2: Float): Float;
var
  x: Float;
begin
  if f <= 0.0 then
  begin
    result := 0.0;
    exit
  end;
  x := df2 / (df2 + df1 * f);
  result := 1.0 - IncompleteBeta(x, df2 / 2.0, df1 / 2.0)
end;

{ Inverse F-distribution: find f such that FDist(f, df1, df2) = p (bisection). }
function InvFDist(p, df1, df2: Float): Float;
const
  EPSILON = 1e-7;
  MAXITER = 300;
var
  lo, hi, mid, fmid: Float;
  iter: integer;
begin
  if p <= 0.0 then
  begin
    result := 0.0;
    exit
  end;
  if p >= 1.0 then
  begin
    result := 1e15;
    exit
  end;
  lo := 0.0;
  hi := 1000.0;
  while FDist(hi, df1, df2) < p do
    hi := hi * 2.0;
  iter := 0;
  mid  := (lo + hi) / 2.0;
  fmid := FDist(mid, df1, df2);
  while (Abs(hi - lo) > EPSILON) and (iter < MAXITER) do
  begin
    if fmid < p then
      lo := mid
    else
      hi := mid;
    mid  := (lo + hi) / 2.0;
    fmid := FDist(mid, df1, df2);
    iter := iter + 1
  end;
  result := mid
end;

{ ---------------------------------------------------------------------------
  CAT 1.4 — Discrete Distributions
  --------------------------------------------------------------------------- }

{ Binomial probability P(X=k) for Binomial(n,p).
  Computed in log space to avoid overflow: exp(lnC(n,k) + k*ln(p) + (n-k)*ln(1-p)). }
function BinomialProb(n, k: integer; p: Float): Float;
var
  lnProb: Float;
begin
  if (k < 0) or (k > n) then
  begin
    result := 0.0;
    exit
  end;
  if p = 0.0 then
  begin
    if k = 0 then result := 1.0 else result := 0.0;
    exit
  end;
  if p = 1.0 then
  begin
    if k = n then result := 1.0 else result := 0.0;
    exit
  end;
  lnProb := LnGamma(n + 1) - LnGamma(k + 1) - LnGamma(n - k + 1)
            + k * Ln(p) + (n - k) * Ln(1.0 - p);
  result := Exp(lnProb)
end;

{ Binomial CDF: P(X <= k) = sum_{i=0}^{k} BinomialProb(n,i,p). }
function BinomialCDF(n, k: integer; p: Float): Float;
var
  i: integer;
  s: Float;
begin
  if k < 0 then
  begin
    result := 0.0;
    exit
  end;
  if k >= n then
  begin
    result := 1.0;
    exit
  end;
  s := 0.0;
  for i := 0 to k do
    s := s + BinomialProb(n, i, p);
  result := s
end;

{ Poisson probability P(X=k) for Poisson(lambda).
  Computed in log space: exp(-lambda + k*ln(lambda) - ln(k!)). }
function PoissonProb(k: integer; lambda: Float): Float;
begin
  if k < 0 then
  begin
    result := 0.0;
    exit
  end;
  if lambda = 0.0 then
  begin
    if k = 0 then result := 1.0 else result := 0.0;
    exit
  end;
  result := Exp(-lambda + k * Ln(lambda) - LnGamma(k + 1))
end;

{ Poisson CDF: P(X <= k) = sum_{i=0}^{k} PoissonProb(i, lambda). }
function PoissonCDF(k: integer; lambda: Float): Float;
var
  i: integer;
  s: Float;
begin
  if k < 0 then
  begin
    result := 0.0;
    exit
  end;
  s := 0.0;
  for i := 0 to k do
    s := s + PoissonProb(i, lambda);
  result := s
end;

{ ---------------------------------------------------------------------------
  CAT 1.5 — Moments
  --------------------------------------------------------------------------- }

function StatMean(var data: TFloatArray; n: integer): Float;
var
  i: integer;
  s: Float;
begin
  s := 0.0;
  for i := 0 to n - 1 do
    s := s + data[i];
  result := s / n
end;

{ Population variance (divides by n). }
function StatVariance(var data: TFloatArray; n: integer): Float;
var
  i: integer;
  m, s, d: Float;
begin
  m := StatMean(data, n);
  s := 0.0;
  for i := 0 to n - 1 do
  begin
    d := data[i] - m;
    s := s + d * d
  end;
  result := s / n
end;

function StatStdDev(var data: TFloatArray; n: integer): Float;
begin
  result := Sqrt(StatVariance(data, n))
end;

function StatSkewness(var data: TFloatArray; n: integer): Float;
var
  i: integer;
  m, sd, s, d: Float;
begin
  m  := StatMean(data, n);
  sd := StatStdDev(data, n);
  s  := 0.0;
  for i := 0 to n - 1 do
  begin
    d := (data[i] - m) / sd;
    s := s + d * d * d
  end;
  result := s / n
end;

{ Excess kurtosis (Fisher's definition): kurt - 3. }
function StatKurtosis(var data: TFloatArray; n: integer): Float;
var
  i: integer;
  m, sd, s, d: Float;
begin
  m  := StatMean(data, n);
  sd := StatStdDev(data, n);
  s  := 0.0;
  for i := 0 to n - 1 do
  begin
    d := (data[i] - m) / sd;
    s := s + d * d * d * d
  end;
  result := s / n - 3.0
end;

{ Median: copies data, sorts the copy, returns the middle value. }
function StatMedian(var data: TFloatArray; n: integer): Float;
var
  tmp: TFloatArray;
  i, j: integer;
  t: Float;
begin
  SetLength(tmp, n);
  for i := 0 to n - 1 do
    tmp[i] := data[i];
  { simple insertion sort }
  for i := 1 to n - 1 do
  begin
    t := tmp[i];
    j := i - 1;
    while (j >= 0) and (tmp[j] > t) do
    begin
      tmp[j + 1] := tmp[j];
      j := j - 1
    end;
    tmp[j + 1] := t
  end;
  if n mod 2 = 1 then
    result := tmp[n div 2]
  else
    result := (tmp[n div 2 - 1] + tmp[n div 2]) / 2.0
end;

{ ---------------------------------------------------------------------------
  self_test
  --------------------------------------------------------------------------- }

procedure self_test;
const
  TOL = 1e-3;

  procedure check(const name: string; computed, expected: Float);
  var
    ok: string;
  begin
    if Abs(computed - expected) < TOL then
      ok := 'OK'
    else
      ok := '*** FAIL ***';
    WriteLn(Format('  %-38s computed=%10.6f  expected=%10.6f  %s',
                   [name, computed, expected, ok]));
    if ok = '*** FAIL ***' then
      SelfTestFail(name + ': computed=' + FloatToStr(computed) + ' expected=' + FloatToStr(expected));
  end;

var
  d1, d2: TFloatArray;
begin
  WriteLn('=== jpmStats self_test ===');
  WriteLn;

  { Chi-Square }
  check('ChiSquareDist(3.84,1)',    ChiSquareDist(3.84, 1.0),  0.9500);
  check('InvChiSquare(0.95,1)',     InvChiSquare(0.95, 1.0),   3.8415);

  { F-Distribution }
  check('FDist(4.0,1,30)',          FDist(4.0, 1.0, 30.0),     0.9450);

  { Binomial }
  check('BinomialProb(10,3,0.5)',   BinomialProb(10, 3, 0.5),  0.117188);

  { Poisson }
  check('PoissonProb(2,3.0)',       PoissonProb(2, 3.0),       0.224042);

  { Moments — [1,2,3,4,5] }
  SetLength(d1, 5);
  d1[0] := 1.0;  d1[1] := 2.0;  d1[2] := 3.0;  d1[3] := 4.0;  d1[4] := 5.0;
  check('StatMean([1..5])',         StatMean(d1, 5),            3.0);

  { Moments — [2,4,4,4,5,5,7,9] }
  SetLength(d2, 8);
  d2[0] := 2.0;  d2[1] := 4.0;  d2[2] := 4.0;  d2[3] := 4.0;
  d2[4] := 5.0;  d2[5] := 5.0;  d2[6] := 7.0;  d2[7] := 9.0;
  check('StatStdDev([2,4,4,4,5,5,7,9])', StatStdDev(d2, 8),    2.0);

  { Median — [3,1,4,1,5] }
  SetLength(d1, 5);
  d1[0] := 3.0;  d1[1] := 1.0;  d1[2] := 4.0;  d1[3] := 1.0;  d1[4] := 5.0;
  check('StatMedian([3,1,4,1,5])', StatMedian(d1, 5),           3.0);

  WriteLn;
  WriteLn('=== done ===');
end;

end.

