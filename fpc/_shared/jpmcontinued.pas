
unit jpmContinued;
{$mode objfpc}{$H+}

interface

uses SysUtils, Math, jpmtypes;

{ Maximum number of continued fraction terms }
const
  CF_MAX = 64;

type
  TCFArray = array[0..CF_MAX-1] of int64;

{ Convert a float to continued fraction coefficients [a0; a1, a2, ...]
  n = number of terms to compute }
procedure FloatToCF(x: Float; var a: TCFArray; var n: integer; maxterms: integer);

{ Convert continued fraction coefficients back to float }
function CFToFloat(var a: TCFArray; n: integer): Float;

{ Compute convergents p[k]/q[k] for k=0..n-1 }
procedure CFConvergents(var a: TCFArray; n: integer;
  var p, q: array of int64);

{ Evaluate a generalized continued fraction:
  b0 + a1/(b1 + a2/(b2 + ...))
  using Lentz method, maxiter iterations }
function GenCFEval(var acoeff, bcoeff: TFloatArray; maxiter: integer): Float;

{ Classic continued fraction for sqrt(N) - returns period length }
function SqrtCF(n: int64; var a: TCFArray; var period: integer): boolean;

{ Best rational approximation to x with denominator <= maxden }
procedure BestRational(x: Float; maxden: int64; var num, den: int64);

procedure self_test;

implementation

procedure FloatToCF(x: Float; var a: TCFArray; var n: integer; maxterms: integer);
var
  xi, frac: Float;
  k: integer;
begin
  if maxterms > CF_MAX then maxterms := CF_MAX;
  n := 0;
  xi := x;
  for k := 0 to maxterms - 1 do
  begin
    a[k] := Trunc(xi);
    if xi < 0 then
      if xi - a[k] < 0 then dec(a[k]);
    frac := xi - a[k];
    inc(n);
    if Abs(frac) < 1e-10 then break;
    xi := 1.0 / frac
  end
end;

function CFToFloat(var a: TCFArray; n: integer): Float;
var
  k: integer;
  r: Float;
begin
  if n <= 0 then begin result := 0; exit; end;
  r := a[n-1];
  for k := n-2 downto 0 do
    r := a[k] + 1.0 / r;
  result := r
end;

procedure CFConvergents(var a: TCFArray; n: integer;
  var p, q: array of int64);
var
  k: integer;
begin
  if n <= 0 then exit;
  p[0] := a[0]; q[0] := 1;
  if n = 1 then exit;
  p[1] := a[1]*a[0] + 1; q[1] := a[1];
  for k := 2 to n-1 do
  begin
    p[k] := a[k]*p[k-1] + p[k-2];
    q[k] := a[k]*q[k-1] + q[k-2]
  end
end;

function GenCFEval(var acoeff, bcoeff: TFloatArray; maxiter: integer): Float;
var
  f, c, d, delta, ai, bi: Float;
  i: integer;
begin
  { Lentz method }
  f := bcoeff[0];
  if Abs(f) < 1e-300 then f := 1e-300;
  c := f;
  d := 0;
  for i := 1 to maxiter - 1 do
  begin
    if i <= high(acoeff) then ai := acoeff[i] else break;
    if i <= high(bcoeff) then bi := bcoeff[i] else break;
    d := bi + ai * d;
    if Abs(d) < 1e-300 then d := 1e-300;
    c := bi + ai / c;
    if Abs(c) < 1e-300 then c := 1e-300;
    d := 1.0 / d;
    delta := c * d;
    f := f * delta;
    if Abs(delta - 1.0) < 1e-14 then break
  end;
  result := f
end;

function SqrtCF(n: int64; var a: TCFArray; var period: integer): boolean;
var
  m, d, a0, ak: int64;
  k: integer;
begin
  result := false;
  period := 0;
  a0 := Trunc(Sqrt(n));
  if a0 * a0 = n then exit; { perfect square }
  a[0] := a0;
  m := 0; d := 1; ak := a0;
  k := 1;
  repeat
    m := d * ak - m;
    d := (n - m * m) div d;
    if d = 0 then break;
    ak := (a0 + m) div d;
    if k < CF_MAX then
    begin
      a[k] := ak;
      inc(k)
    end;
    inc(period)
  until (ak = 2 * a0) or (k >= CF_MAX);
  result := true
end;

procedure BestRational(x: Float; maxden: int64; var num, den: int64);
var
  a: TCFArray;
  n, k: integer;
  p0, p1, q0, q1, pk, qk: int64;
begin
  FloatToCF(x, a, n, 20);
  p0 := 1; q0 := 0;
  p1 := a[0]; q1 := 1;
  for k := 1 to n-1 do
  begin
    pk := a[k]*p1 + p0;
    qk := a[k]*q1 + q0;
    if qk > maxden then break;
    p0 := p1; q0 := q1;
    p1 := pk; q1 := qk
  end;
  num := p1; den := q1
end;

procedure self_test;
var
  a: TCFArray;
  n, period: integer;
  x, approx: Float;
  num, den: int64;
  p, q: array[0..CF_MAX-1] of int64;
begin
  writeln('=== jpmContinued self_test ===');
  writeln;

  { Test 1: FloatToCF and CFToFloat for Pi }
  x := Pi;
  FloatToCF(x, a, n, 10);
  write('CF(Pi,10 terms): [');
  write(a[0]);
  for n := 1 to 9 do write(';',a[n]);
  writeln(']  expected [3;7,15,1,292,1,1,1,2,1]');
  approx := CFToFloat(a, 10);
  writeln('CFToFloat back: ', approx:12:8, '  expected ', Pi:12:8);
  if Abs(approx - Pi) < 1e-10 then writeln('PASS') else writeln('FAIL');
  writeln;

  { Test 2: CFConvergents for Pi - check 22/7 and 355/113 }
  FloatToCF(Pi, a, n, 6);
  CFConvergents(a, 6, p, q);
  writeln('Pi convergents:');
  writeln('  p1/q1 = ', p[1], '/', q[1], '  expected 22/7');
  writeln('  p3/q3 = ', p[3], '/', q[3], '  expected 333/106');
  if (p[1] = 22) and (q[1] = 7) then writeln('  22/7 PASS') else writeln('  22/7 FAIL');
  writeln;

  { Test 3: sqrt(2) CF = [1;2,2,2,...] }
  SqrtCF(2, a, period);
  writeln('CF(sqrt(2)): a[0]=', a[0], ' period=', period,
          '  expected a[0]=1, period=1, a[1]=2');
  writeln('  a[1]=', a[1], '  expected 2');
  if (a[0]=1) and (a[1]=2) and (period=1) then writeln('PASS') else writeln('FAIL');
  writeln;

  { Test 4: BestRational for Pi, max den 1000 → 355/113 }
  BestRational(Pi, 1000, num, den);
  writeln('BestRational(Pi, 1000): ', num, '/', den,
          '  expected 355/113');
  if (num=355) and (den=113) then writeln('PASS') else writeln('FAIL');
  writeln;

  { Test 5: BestRational for sqrt(2), max den 100 → 99/70 }
  BestRational(Sqrt(2), 100, num, den);
  writeln('BestRational(sqrt2, 100): ', num, '/', den,
          '  expected 99/70');
  if (num=99) and (den=70) then writeln('PASS') else writeln('FAIL');
  writeln;

  writeln('=== done ===')
end;

end.
