unit jpmPolynomials;
{$mode objfpc}{$H+}

{ Polynomial algebra over Float (TPolyArray = array of Float).
  Index 0 = constant term; index k = coefficient of x^k.
  Degree parameters (da, db, dc …) track the highest non-zero index. }

interface

uses SysUtils, Math, jpmtypes;

type
  TPolyArray = array of Float;

{ 8.4 Helper utilities }
function  PolyDegree(var p: TPolyArray; maxDeg: integer): integer;
procedure PolyNormalize(var p: TPolyArray; var dp: integer);
procedure PolyPrint(var p: TPolyArray; dp: integer; name: string);

{ 8.1 Horner evaluation }
function  PolyEval(var p: TPolyArray; degree: integer; x: Float): Float;

{ 8.2 Polynomial arithmetic }
procedure PolyAdd(var a, b: TPolyArray; da, db: integer; var c: TPolyArray; var dc: integer);
procedure PolySub(var a, b: TPolyArray; da, db: integer; var c: TPolyArray; var dc: integer);
procedure PolyMul(var a, b: TPolyArray; da, db: integer; var c: TPolyArray; var dc: integer);
procedure PolyDiv(var a, b: TPolyArray; da, db: integer; var q, r: TPolyArray; var dq, dr: integer);
procedure PolyDeriv(var p: TPolyArray; dp: integer; var d: TPolyArray; var dd: integer);
procedure PolyInteg(var p: TPolyArray; dp: integer; c0: Float; var q: TPolyArray; var dq: integer);

{ 8.3 Polynomial GCD (result is monic) }
procedure PolyGCD(var a, b: TPolyArray; da, db: integer; var g: TPolyArray; var dg: integer);

{ 8.4 Polynomial fractions }
type
  TPolyFrac = record
    numer, denom: TPolyArray;
    dn, dd: integer;
  end;

procedure PolyFracMake(var num, den: TPolyArray; dn, dd: integer; var f: TPolyFrac);
function  PolyFracEval(var f: TPolyFrac; x: Float): Float;
procedure PolyFracAdd(var f1, f2: TPolyFrac; var res: TPolyFrac);
procedure PolyFracSub(var f1, f2: TPolyFrac; var res: TPolyFrac);
procedure PolyFracMul(var f1, f2: TPolyFrac; var res: TPolyFrac);
procedure PolyFracDiv(var f1, f2: TPolyFrac; var res: TPolyFrac);
procedure PolyFracReduce(var f: TPolyFrac);
{ Partial fraction decomposition for simple real roots.
  roots[0..f.dd-1] = roots of denominator.
  coeffs[i] = residue A_i for each (x - roots[i]) factor }
procedure PolyFracPartial(var f: TPolyFrac; var roots: TFloatArray;
                          var coeffs: TFloatArray);

procedure self_test;

implementation

const
  POLY_SMALL = 1.0e-12;

{ ====================================================================
  8.4  Helpers
  ==================================================================== }

{ Return actual degree of p (strip leading near-zero coefficients). }
function PolyDegree(var p: TPolyArray; maxDeg: integer): integer;
var
  d: integer;
begin
  d := maxDeg;
  while (d > 0) and (Abs(p[d]) < POLY_SMALL) do
    dec(d);
  result := d;
end;

{ Divide all coefficients by the leading coefficient so the result is monic. }
procedure PolyNormalize(var p: TPolyArray; var dp: integer);
var
  i: integer;
  lc: Float;
begin
  dp := PolyDegree(p, dp);
  lc := p[dp];
  if Abs(lc) < POLY_SMALL then
    exit;
  for i := 0 to dp do
    p[i] := p[i] / lc;
end;

{ Print polynomial in human-readable form: name = c0 + c1*x + c2*x^2 + … }
procedure PolyPrint(var p: TPolyArray; dp: integer; name: string);
var
  i: integer;
  first: boolean;
  coeff: Float;
begin
  write(name, ' = ');
  first := true;
  for i := dp downto 0 do
  begin
    coeff := p[i];
    if Abs(coeff) >= POLY_SMALL then
    begin
      if not first then
      begin
        if coeff > 0 then
          write(' + ')
        else
        begin
          write(' - ');
          coeff := Abs(coeff);
        end;
      end
      else
      begin
        first := false;
        if coeff < 0 then
        begin
          write('-');
          coeff := Abs(coeff);
        end;
      end;
      write(Format('%g', [coeff]));
      if i = 1 then
        write('*x')
      else if i > 1 then
        write('*x^', i);
    end;
  end;
  if first then
    write('0');
  writeln;
end;

{ ====================================================================
  8.1  Horner evaluation
  ==================================================================== }

{ Evaluate polynomial p of given degree at x using Horner's method. }
function PolyEval(var p: TPolyArray; degree: integer; x: Float): Float;
var
  i: integer;
  res: Float;
begin
  res := p[degree];
  for i := degree - 1 downto 0 do
    res := res * x + p[i];
  result := res;
end;

{ ====================================================================
  8.2  Arithmetic
  ==================================================================== }

{ c = a + b }
procedure PolyAdd(var a, b: TPolyArray; da, db: integer; var c: TPolyArray; var dc: integer);
var
  i: integer;
begin
  if da >= db then
    dc := da
  else
    dc := db;
  SetLength(c, dc + 1);
  for i := 0 to dc do
    c[i] := 0.0;
  for i := 0 to da do
    c[i] := c[i] + a[i];
  for i := 0 to db do
    c[i] := c[i] + b[i];
  dc := PolyDegree(c, dc);
end;

{ c = a - b }
procedure PolySub(var a, b: TPolyArray; da, db: integer; var c: TPolyArray; var dc: integer);
var
  i: integer;
begin
  if da >= db then
    dc := da
  else
    dc := db;
  SetLength(c, dc + 1);
  for i := 0 to dc do
    c[i] := 0.0;
  for i := 0 to da do
    c[i] := c[i] + a[i];
  for i := 0 to db do
    c[i] := c[i] - b[i];
  dc := PolyDegree(c, dc);
end;

{ c = a * b }
procedure PolyMul(var a, b: TPolyArray; da, db: integer; var c: TPolyArray; var dc: integer);
var
  i, j: integer;
begin
  dc := da + db;
  SetLength(c, dc + 1);
  for i := 0 to dc do
    c[i] := 0.0;
  for i := 0 to da do
    for j := 0 to db do
      c[i + j] := c[i + j] + a[i] * b[j];
end;

{ Euclidean (long) division: a = q*b + r
  If da < db: q = 0, r = a. }
procedure PolyDiv(var a, b: TPolyArray; da, db: integer; var q, r: TPolyArray; var dq, dr: integer);
var
  i, j: integer;
  coeff: Float;
begin
  dr := da;
  SetLength(r, da + 1);
  for i := 0 to da do
    r[i] := a[i];
  if da < db then
  begin
    dq := 0;
    SetLength(q, 1);
    q[0] := 0.0;
    dr := PolyDegree(r, da);
    exit;
  end;
  dq := da - db;
  SetLength(q, dq + 1);
  for i := 0 to dq do
    q[i] := 0.0;
  for i := dq downto 0 do
  begin
    coeff := r[i + db] / b[db];
    q[i] := coeff;
    for j := 0 to db do
      r[i + j] := r[i + j] - coeff * b[j];
  end;
  dr := PolyDegree(r, da);
  dq := PolyDegree(q, dq);
end;

{ Formal derivative: d = p' }
procedure PolyDeriv(var p: TPolyArray; dp: integer; var d: TPolyArray; var dd: integer);
var
  i: integer;
begin
  if dp = 0 then
  begin
    dd := 0;
    SetLength(d, 1);
    d[0] := 0.0;
    exit;
  end;
  dd := dp - 1;
  SetLength(d, dd + 1);
  for i := 0 to dd do
    d[i] := (i + 1) * p[i + 1];
end;

{ Formal antiderivative with constant term c0: q = integral(p) + c0 }
procedure PolyInteg(var p: TPolyArray; dp: integer; c0: Float; var q: TPolyArray; var dq: integer);
var
  i: integer;
begin
  dq := dp + 1;
  SetLength(q, dq + 1);
  q[0] := c0;
  for i := 1 to dq do
    q[i] := p[i - 1] / i;
end;

{ ====================================================================
  8.3  GCD — Euclidean algorithm; result is monic
  ==================================================================== }

procedure PolyGCD(var a, b: TPolyArray; da, db: integer; var g: TPolyArray; var dg: integer);
var
  p0, p1, rem, quot: TPolyArray;
  d0, d1, drem, dquot: integer;
  i: integer;
  tmp: TPolyArray;
  dtmp: integer;
  isZero: boolean;
begin
  d0 := da;
  SetLength(p0, da + 1);
  for i := 0 to da do
    p0[i] := a[i];
  d1 := db;
  SetLength(p1, db + 1);
  for i := 0 to db do
    p1[i] := b[i];
  { ensure higher-degree poly is in p0 }
  if d0 < d1 then
  begin
    tmp := p0; p0 := p1; p1 := tmp;
    dtmp := d0; d0 := d1; d1 := dtmp;
  end;
  isZero := (d1 = 0) and (Abs(p1[0]) < POLY_SMALL);
  while not isZero do
  begin
    PolyDiv(p0, p1, d0, d1, quot, rem, dquot, drem);
    p0 := p1; d0 := d1;
    p1 := rem; d1 := drem;
    isZero := (d1 = 0) and (Abs(p1[0]) < POLY_SMALL);
  end;
  g := p0;
  dg := d0;
  PolyNormalize(g, dg);
end;

{ ====================================================================
  Self-test
  ==================================================================== }

procedure self_test;
var
  p, q, c, d, g, qpol, r: TPolyArray;
  dp, dq, dc, dd, dg, dqpol, dr: integer;
  v: Float;

  procedure check(cond: boolean; msg: string);
  begin
    if not cond then
    begin
      writeln('FAIL: ', msg);
      halt(1);
    end
    else
      writeln('PASS: ', msg);
  end;

begin
  writeln('=== jpmPolynomials self_test ===');

  { --- setup polynomials ---
    p = 1 + 2x + x^2 = (x+1)^2  }
  dp := 2;
  SetLength(p, 3);
  p[0] := 1.0; p[1] := 2.0; p[2] := 1.0;

  { q = x + 1 }
  dq := 1;
  SetLength(q, 2);
  q[0] := 1.0; q[1] := 1.0;

  PolyPrint(p, dp, 'p');
  PolyPrint(q, dq, 'q');

  { 8.1 PolyEval: p(3) = 1 + 6 + 9 = 16 }
  v := PolyEval(p, dp, 3.0);
  writeln(Format('PolyEval(p, 3) = %g  (expected 16)', [v]));
  check(Abs(v - 16.0) < 1.0e-10, Format('PolyEval(p,3) = 16, got %g', [v]));

  { 8.2 PolyAdd: p + q = 2 + 3x + x^2 }
  PolyAdd(p, q, dp, dq, c, dc);
  PolyPrint(c, dc, 'p+q');
  check(dc = 2, Format('PolyAdd degree=2, got %d', [dc]));
  check(Abs(c[0] - 2.0) < 1.0e-10, 'PolyAdd c[0]=2');
  check(Abs(c[1] - 3.0) < 1.0e-10, 'PolyAdd c[1]=3');
  check(Abs(c[2] - 1.0) < 1.0e-10, 'PolyAdd c[2]=1');

  { 8.2 PolySub: p - q = x^2 + x }
  PolySub(p, q, dp, dq, c, dc);
  PolyPrint(c, dc, 'p-q');
  check(dc = 2, Format('PolySub degree=2, got %d', [dc]));
  check(Abs(c[0]) < 1.0e-10, 'PolySub c[0]=0');
  check(Abs(c[1] - 1.0) < 1.0e-10, 'PolySub c[1]=1');
  check(Abs(c[2] - 1.0) < 1.0e-10, 'PolySub c[2]=1');

  { 8.2 PolyMul: p * q = (x+1)^3 = 1 + 3x + 3x^2 + x^3 }
  PolyMul(p, q, dp, dq, c, dc);
  PolyPrint(c, dc, 'p*q');
  check(dc = 3, Format('PolyMul degree=3, got %d', [dc]));
  check(Abs(c[0] - 1.0) < 1.0e-10, 'PolyMul c[0]=1');
  check(Abs(c[1] - 3.0) < 1.0e-10, 'PolyMul c[1]=3');
  check(Abs(c[2] - 3.0) < 1.0e-10, 'PolyMul c[2]=3');
  check(Abs(c[3] - 1.0) < 1.0e-10, 'PolyMul c[3]=1');

  { 8.2 PolyDiv: p / q = (x+1), remainder 0 }
  PolyDiv(p, q, dp, dq, qpol, r, dqpol, dr);
  PolyPrint(qpol, dqpol, 'p div q');
  PolyPrint(r, dr, 'p mod q');
  check(dqpol = 1, Format('PolyDiv quot degree=1, got %d', [dqpol]));
  check(Abs(qpol[0] - 1.0) < 1.0e-10, 'PolyDiv quot[0]=1');
  check(Abs(qpol[1] - 1.0) < 1.0e-10, 'PolyDiv quot[1]=1');
  check(Abs(r[0]) < 1.0e-10, 'PolyDiv remainder=0');

  { 8.2 PolyDeriv: p' = 2 + 2x }
  PolyDeriv(p, dp, d, dd);
  PolyPrint(d, dd, 'p''');
  check(dd = 1, Format('PolyDeriv degree=1, got %d', [dd]));
  check(Abs(d[0] - 2.0) < 1.0e-10, 'PolyDeriv d[0]=2');
  check(Abs(d[1] - 2.0) < 1.0e-10, 'PolyDeriv d[1]=2');

  { 8.2 PolyInteg: integral(q, c0=0) = x^2/2 + x }
  PolyInteg(q, dq, 0.0, c, dc);
  PolyPrint(c, dc, 'integ(q,0)');
  check(dc = 2, Format('PolyInteg degree=2, got %d', [dc]));
  check(Abs(c[0]) < 1.0e-10, 'PolyInteg c[0]=0');
  check(Abs(c[1] - 1.0) < 1.0e-10, 'PolyInteg c[1]=1');
  check(Abs(c[2] - 0.5) < 1.0e-10, 'PolyInteg c[2]=0.5');

  { 8.3 PolyGCD: gcd(p, q) = x+1 (monic) }
  PolyGCD(p, q, dp, dq, g, dg);
  PolyPrint(g, dg, 'gcd(p,q)');
  check(dg = 1, Format('PolyGCD degree=1, got %d', [dg]));
  check(Abs(g[0] - 1.0) < 1.0e-10, 'PolyGCD g[0]=1');
  check(Abs(g[1] - 1.0) < 1.0e-10, 'PolyGCD g[1]=1');

  writeln('=== self_test PASSED ===');
end;


{ ====================================================================
  8.4  Polynomial Fractions
  ==================================================================== }

procedure PolyFracMake(var num, den: TPolyArray; dn, dd: integer; var f: TPolyFrac);
var k: integer;
begin
  SetLength(f.numer, dn+1);
  SetLength(f.denom, dd+1);
  for k := 0 to dn do f.numer[k] := num[k];
  for k := 0 to dd do f.denom[k] := den[k];
  f.dn := dn; f.dd := dd
end;

function PolyFracEval(var f: TPolyFrac; x: Float): Float;
var n, d: Float;
begin
  n := PolyEval(f.numer, f.dn, x);
  d := PolyEval(f.denom, f.dd, x);
  if Abs(d) < POLY_SMALL then result := 0
  else result := n / d
end;

procedure PolyFracAdd(var f1, f2: TPolyFrac; var res: TPolyFrac);
var t1, t2, num, den: TPolyArray;
    dt1, dt2, dnum, dden: integer;
begin
  { res = f1.n/f1.d + f2.n/f2.d = (f1.n*f2.d + f2.n*f1.d) / (f1.d*f2.d) }
  PolyMul(f1.numer, f2.denom, f1.dn, f2.dd, t1, dt1);
  PolyMul(f2.numer, f1.denom, f2.dn, f1.dd, t2, dt2);
  PolyAdd(t1, t2, dt1, dt2, num, dnum);
  PolyMul(f1.denom, f2.denom, f1.dd, f2.dd, den, dden);
  PolyFracMake(num, den, dnum, dden, res)
end;

procedure PolyFracSub(var f1, f2: TPolyFrac; var res: TPolyFrac);
var t1, t2, num, den: TPolyArray;
    dt1, dt2, dnum, dden: integer;
begin
  PolyMul(f1.numer, f2.denom, f1.dn, f2.dd, t1, dt1);
  PolyMul(f2.numer, f1.denom, f2.dn, f1.dd, t2, dt2);
  PolySub(t1, t2, dt1, dt2, num, dnum);
  PolyMul(f1.denom, f2.denom, f1.dd, f2.dd, den, dden);
  PolyFracMake(num, den, dnum, dden, res)
end;

procedure PolyFracMul(var f1, f2: TPolyFrac; var res: TPolyFrac);
var num, den: TPolyArray;
    dnum, dden: integer;
begin
  PolyMul(f1.numer, f2.numer, f1.dn, f2.dn, num, dnum);
  PolyMul(f1.denom, f2.denom, f1.dd, f2.dd, den, dden);
  PolyFracMake(num, den, dnum, dden, res)
end;

procedure PolyFracDiv(var f1, f2: TPolyFrac; var res: TPolyFrac);
var num, den: TPolyArray;
    dnum, dden: integer;
begin
  PolyMul(f1.numer, f2.denom, f1.dn, f2.dd, num, dnum);
  PolyMul(f1.denom, f2.numer, f1.dd, f2.dn, den, dden);
  PolyFracMake(num, den, dnum, dden, res)
end;

procedure PolyFracReduce(var f: TPolyFrac);
var g, n2, d2: TPolyArray;
    dg, dn2, dd2: integer;
begin
  PolyGCD(f.numer, f.denom, f.dn, f.dd, g, dg);
  if dg = 0 then exit; { gcd = constant, already reduced }
  PolyDiv(f.numer, g, f.dn, dg, n2, n2, dn2, dn2); { quotient only }
  PolyDiv(f.denom, g, f.dd, dg, d2, d2, dd2, dd2);
  { PolyDiv returns quotient in first output, remainder in second — fix: use proper vars }
  PolyDiv(f.numer, g, f.dn, dg, n2, d2, dn2, dd2);
  PolyDiv(f.denom, g, f.dd, dg, d2, g,  dd2, dg);
  SetLength(f.numer, dn2+1);
  Move(n2[0], f.numer[0], (dn2+1)*SizeOf(Float));
  f.dn := dn2;
  SetLength(f.denom, dd2+1);
  Move(d2[0], f.denom[0], (dd2+1)*SizeOf(Float));
  f.dd := dd2
end;

procedure PolyFracPartial(var f: TPolyFrac; var roots: TFloatArray;
                          var coeffs: TFloatArray);
{ Simple poles: f(x) = N(x)/D(x) where D has simple real roots.
  Residue at root r_i: A_i = N(r_i) / D'(r_i) }
var deriv: TPolyArray;
    dd: integer;
    i: integer;
    nr, dr: Float;
begin
  SetLength(coeffs, f.dd);
  PolyDeriv(f.denom, f.dd, deriv, dd);
  for i := 0 to f.dd - 1 do
  begin
    nr := PolyEval(f.numer, f.dn, roots[i]);
    dr := PolyEval(deriv, dd, roots[i]);
    if Abs(dr) < POLY_SMALL then coeffs[i] := 0
    else coeffs[i] := nr / dr
  end
end;

end.
