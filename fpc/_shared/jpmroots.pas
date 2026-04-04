
{ Root finding routines ported from Jean-Pierre Moreau's programs }

unit jpmRoots;

{$mode objfpc}{$H+}

interface

uses
  SysUtils, Math,
  jpmTypes;

type
  TComplexRoot = record
    Re, Im: Float;
  end;
  TComplexRootArray = array of TComplexRoot;

{ Bisection method: finds root of AFunc in [a,b] }
function BisectionRoot(AFunc: TFunction1; a, b, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;

{ Newton-Raphson: finds root given function and its derivative }
function NewtonRoot(AFunc, ADerivFunc: TFunction1; x0, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;

{ Secant method: finds root using two initial guesses }
function SecantRoot(AFunc: TFunction1; x0, x1, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;

{ Regula Falsi (false position): bracketed secant method }
function RegulaFalsiRoot(AFunc: TFunction1; a, b, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;

{ Brent's method: robust bracketed root finder }
function BrentRoot(AFunc: TFunction1; a, b, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;

{ Mueller's method: uses three points, handles complex roots }
function MuellerRoot(AFunc: TFunction1; x0, x1, x2, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;

{ Steffensen's method: fixed-point iteration with Aitken acceleration }
function SteffensenRoot(AFunc: TFunction1; x0, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;

{ Aitken's delta-squared method }
function AitkenRoot(AFunc: TFunction1; x0, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;

{ Bairstow's method: finds all roots of a polynomial (real + complex) }
{ Coeffs[0] is constant term, Coeffs[n] is leading coefficient }
function BairstowRoots(const Coeffs: array of Float; Tol: Float;
  MaxIter: integer): TComplexRootArray;

procedure self_test;

implementation

{***********************************************************************
* BisectionRoot                                                        *
* Finds a root of AFunc(x)=0 in interval [a,b] by bisection.         *
* Requires AFunc(a)*AFunc(b) < 0 (sign change in interval).          *
***********************************************************************}
function BisectionRoot(AFunc: TFunction1; a, b, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;
var
  fa, fb, fc, c: Float;
begin
  fa := AFunc(a);
  fb := AFunc(b);
  Iterations := 0;
  Result := a;
  if fa * fb > 0 then
    exit; { no sign change - return a as error indicator }
  repeat
    c := (a + b) / 2.0;
    fc := AFunc(c);
    Inc(Iterations);
    if fa * fc <= 0 then
    begin
      b := c;
      fb := fc;
    end
    else
    begin
      a := c;
      fa := fc;
    end;
  until (Abs(b - a) < Tol) or (Iterations >= MaxIter);
  Result := (a + b) / 2.0;
end;

{***********************************************************************
* NewtonRoot                                                           *
* Newton-Raphson iteration: x(n+1) = x(n) - f(x)/f'(x)              *
***********************************************************************}
function NewtonRoot(AFunc, ADerivFunc: TFunction1; x0, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;
var
  x, fx, fpx, dx: Float;
begin
  x := x0;
  Iterations := 0;
  repeat
    fx  := AFunc(x);
    fpx := ADerivFunc(x);
    if Abs(fpx) < 1e-15 then
      break; { avoid division by zero }
    dx := fx / fpx;
    x  := x - dx;
    Inc(Iterations);
  until (Abs(dx) < Tol) or (Iterations >= MaxIter);
  Result := x;
end;

{***********************************************************************
* SecantRoot                                                           *
* Secant method: chord through (x0,f(x0)) and (x1,f(x1))            *
***********************************************************************}
function SecantRoot(AFunc: TFunction1; x0, x1, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;
var
  y0, y1, x2, denom: Float;
begin
  Iterations := 0;
  y0 := AFunc(x0);
  y1 := AFunc(x1);
  repeat
    denom := y1 - y0;
    if Abs(denom) < 1e-15 then
      denom := 1e-15; { guard against zero denominator }
    x2 := (x0 * y1 - x1 * y0) / denom;
    x0 := x1;  y0 := y1;
    x1 := x2;  y1 := AFunc(x1);
    Inc(Iterations);
  until (Abs(x1 - x0) < Tol) or (Iterations >= MaxIter);
  Result := x1;
end;

{***********************************************************************
* RegulaFalsiRoot                                                      *
* False position method: bracketed secant, maintains sign-change.     *
***********************************************************************}
function RegulaFalsiRoot(AFunc: TFunction1; a, b, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;
var
  fa, fb, fc, c: Float;
begin
  fa := AFunc(a);
  fb := AFunc(b);
  Iterations := 0;
  Result := a;
  if fa * fb > 0 then
    exit;
  repeat
    c  := (a * fb - b * fa) / (fb - fa);
    fc := AFunc(c);
    Inc(Iterations);
    if fa * fc <= 0 then
    begin
      b := c;
      fb := fc;
    end
    else
    begin
      a := c;
      fa := fc;
    end;
  until (Abs(fc) < Tol) or (Abs(b - a) < Tol) or (Iterations >= MaxIter);
  Result := c;
end;

{***********************************************************************
* BrentRoot                                                            *
* Brent's method: combines bisection, secant and inverse quadratic.   *
* Reference: Brent (1973) "Algorithms for Minimization without Deriv" *
***********************************************************************}
function BrentRoot(AFunc: TFunction1; a, b, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;
var
  fa, fb, fc, c, d, e, s, p, q, r, tol1, xm: Float;

begin
  fa := AFunc(a);
  fb := AFunc(b);
  Iterations := 0;
  Result := b;
  if fa * fb > 0 then
    exit;
  c  := a;
  fc := fa;
  d  := b - a;
  e  := d;
  repeat
    if fb * fc > 0 then
    begin
      c  := a;  fc := fa;
      d  := b - a;
      e  := d;
    end;
    if Abs(fc) < Abs(fb) then
    begin
      a := b;  fa := fb;
      b := c;  fb := fc;
      c := a;  fc := fa;
    end;
    tol1 := 2.0 * 1e-15 * Abs(b) + 0.5 * Tol;
    xm   := 0.5 * (c - b);
    if (Abs(xm) <= tol1) or (fb = 0) then
    begin
      Result := b;
      exit;
    end;
    { Try inverse quadratic interpolation or secant }
    if (Abs(e) >= tol1) and (Abs(fa) > Abs(fb)) then
    begin
      s := fb / fa;
      if a = c then
      begin
        p := 2.0 * xm * s;
        q := 1.0 - s;
      end
      else
      begin
        q := fa / fc;
        r := fb / fc;
        p := s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q := (q - 1.0) * (r - 1.0) * (s - 1.0);
      end;
      if p > 0 then
        q := -q
      else
        p := -p;
      if (2.0 * p < (3.0 * xm * q - Abs(tol1 * q))) and
         (2.0 * p < Abs(e * q)) then
      begin
        e := d;
        d := p / q;
        { use interpolation }
      end
      else
      begin
        d := xm;
        e := d;
      end;
    end
    else
    begin
      d := xm;
      e := d;
    end;
    a  := b;
    fa := fb;
    if Abs(d) > tol1 then
      b := b + d
    else if xm > 0 then
      b := b + tol1
    else
      b := b - tol1;
    fb := AFunc(b);
    Inc(Iterations);
  until Iterations >= MaxIter;
  Result := b;
end;

{***********************************************************************
* MuellerRoot                                                          *
* Mueller's method: quadratic interpolation through three points.     *
* Finds real roots (complex roots not returned here).                 *
***********************************************************************}
function MuellerRoot(AFunc: TFunction1; x0, x1, x2, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;
var
  f0, f1, f2, h1, h2, d1, d2, a, b, c, rad, den, dx, xnew: Float;
begin
  Iterations := 0;
  f0 := AFunc(x0);
  f1 := AFunc(x1);
  f2 := AFunc(x2);
  repeat
    h1 := x1 - x0;
    h2 := x2 - x1;
    d1 := (f1 - f0) / h1;
    d2 := (f2 - f1) / h2;
    a  := (d2 - d1) / (h2 + h1);
    b  := a * h2 + d2;
    c  := f2;
    rad := b * b - 4.0 * a * c;
    if rad < 0 then
      rad := 0; { take real part only }
    rad := Sqrt(rad);
    if Abs(b + rad) > Abs(b - rad) then
      den := b + rad
    else
      den := b - rad;
    if Abs(den) < 1e-15 then
      den := 1e-15;
    dx   := -2.0 * c / den;
    xnew := x2 + dx;
    x0 := x1;  f0 := f1;
    x1 := x2;  f1 := f2;
    x2 := xnew;
    f2 := AFunc(x2);
    Inc(Iterations);
  until (Abs(dx) < Tol) or (Iterations >= MaxIter);
  Result := x2;
end;

{***********************************************************************
* SteffensenRoot                                                       *
* Steffensen's method: g(x) = x - f(x), accelerated with Aitken's.  *
***********************************************************************}
function SteffensenRoot(AFunc: TFunction1; x0, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;
var
  x1, x2, denom, dx: Float;
begin
  Iterations := 0;
  repeat
    x1 := x0 + AFunc(x0);
    x2 := x1 + AFunc(x1);
    denom := x2 - 2.0 * x1 + x0;
    if Abs(denom) < 1e-15 then
      denom := 1e-15;
    dx := (x2 - x1) * (x2 - x1) / denom;
    x0 := x2 - dx;
    Inc(Iterations);
  until (Abs(dx) < Tol) or (Iterations >= MaxIter);
  Result := x0;
end;

{***********************************************************************
* AitkenRoot                                                           *
* Aitken's delta-squared acceleration of fixed-point iteration.       *
* Treats g(x) = x - f(x) as fixed-point map.                         *
***********************************************************************}
function AitkenRoot(AFunc: TFunction1; x0, Tol: Float;
  MaxIter: integer; var Iterations: integer): Float;
var
  x1, x2, denom, xnew: Float;
begin
  Iterations := 0;
  repeat
    x1 := x0 - AFunc(x0);
    x2 := x1 - AFunc(x1);
    denom := x2 - 2.0 * x1 + x0;
    if Abs(denom) < 1e-15 then
    begin
      Result := x2;
      exit;
    end;
    xnew := x0 - (x1 - x0) * (x1 - x0) / denom;
    x0   := xnew;
    Inc(Iterations);
  until (Abs(AFunc(x0)) < Tol) or (Iterations >= MaxIter);
  Result := x0;
end;

{***********************************************************************
* BairstowRoots                                                        *
* Bairstow's method for finding all roots of a polynomial.            *
* Coeffs[0]=constant, Coeffs[deg]=leading coeff.                     *
* Returns array of TComplexRoot with all roots (real and complex).    *
***********************************************************************}
function BairstowRoots(const Coeffs: array of Float; Tol: Float;
  MaxIter: integer): TComplexRootArray;
var
  n, i, iter, rootIdx: integer;
  r, s, dr, ds: Float;
  a, b2, c2: array of Float;
  disc, sqrtDisc: Float;
  done: boolean;
begin
  n := High(Coeffs); { degree of polynomial }
  Result := nil;
  SetLength(Result, 0);
  if n < 1 then
    exit;
  SetLength(a, n + 1);
  for i := 0 to n do
    a[i] := Coeffs[i];
  rootIdx := 0;
  SetLength(Result, n);
  { Divide out roots until degree < 2 }
  while n >= 2 do
  begin
    { Initial estimates for quadratic factor x^2 - r*x - s }
    r := 0.5;
    s := 0.5;
    SetLength(b2, n + 1);
    SetLength(c2, n + 1);
    iter := 0;
    done := false;
    repeat
      b2[n] := a[n];
      if n >= 1 then
        b2[n-1] := a[n-1] + r * b2[n];
      for i := n-2 downto 0 do
        b2[i] := a[i] + r * b2[i+1] + s * b2[i+2];
      c2[n] := b2[n];
      if n >= 1 then
        c2[n-1] := b2[n-1] + r * c2[n];
      for i := n-2 downto 0 do
        c2[i] := b2[i] + r * c2[i+1] + s * c2[i+2];
      { Solve 2x2 system for dr, ds }
      { Jacobian: [c2[2] c2[3]; c2[1] c2[2]] }
      { det = c2[2]^2 - c2[3]*c2[1] }
      disc := c2[2] * c2[2] - c2[3] * c2[1];
      if Abs(disc) < 1e-20 then
        disc := 1e-20;
      dr := (-b2[1] * c2[2] + b2[0] * c2[3]) / disc;
      ds := (-b2[0] * c2[2] + b2[1] * c2[1]) / disc;
      r  := r + dr;
      s  := s + ds;
      Inc(iter);
      if (Abs(dr) < Tol) and (Abs(ds) < Tol) then
        done := true;
    until done or (iter >= MaxIter);
    { Extract roots of quadratic x^2 - r*x - s = 0 }
    disc := r * r + 4.0 * s;
    if disc >= 0 then
    begin
      sqrtDisc := Sqrt(disc);
      Result[rootIdx].Re := (r + sqrtDisc) / 2.0;
      Result[rootIdx].Im := 0.0;
      Inc(rootIdx);
      Result[rootIdx].Re := (r - sqrtDisc) / 2.0;
      Result[rootIdx].Im := 0.0;
      Inc(rootIdx);
    end
    else
    begin
      sqrtDisc := Sqrt(-disc);
      Result[rootIdx].Re   := r / 2.0;
      Result[rootIdx].Im   := sqrtDisc / 2.0;
      Inc(rootIdx);
      Result[rootIdx].Re   := r / 2.0;
      Result[rootIdx].Im   := -sqrtDisc / 2.0;
      Inc(rootIdx);
    end;
    { Deflate polynomial by the quadratic factor }
    n := n - 2;
    SetLength(a, n + 1);
    for i := 0 to n do
      a[i] := b2[i + 2];
  end;
  { Handle remaining linear or constant factor }
  if n = 1 then
  begin
    Result[rootIdx].Re := -a[0] / a[1];
    Result[rootIdx].Im := 0.0;
    Inc(rootIdx);
  end;
  SetLength(Result, rootIdx);
end;

{***********************************************************************
* self_test                                                            *
* Tests all root-finding methods on f(x) = x^3 - x - 2               *
* Exact root: ~1.521379706804                                          *
* Also tests Bairstow on x^3 - 6x^2 + 11x - 6 (roots: 1, 2, 3)     *
***********************************************************************}

function TestFunc(x: Float): Float;
begin
  Result := x * x * x - x - 2.0;
end;

function TestFuncDeriv(x: Float): Float;
begin
  Result := 3.0 * x * x - 1.0;
end;

procedure self_test;
var
  root: Float;
  iters: integer;
  roots: TComplexRootArray;
  i: integer;
begin
  WriteLn('=== jpmRoots Self Test ===');
  WriteLn('Target function: f(x) = x^3 - x - 2');
  WriteLn('Expected root: ~1.521379706804');
  WriteLn;

  iters := 0;
  root := BisectionRoot(@TestFunc, 1.0, 2.0, 1e-10, 100, iters);
  WriteLn('Bisection:      root=', root:18:12, '  iters=', iters);
  SelfTestCheck(Abs(root - 1.521379706804) < 1e-8, 'BisectionRoot');

  iters := 0;
  root := NewtonRoot(@TestFunc, @TestFuncDeriv, 1.5, 1e-10, 100, iters);
  WriteLn('Newton-Raphson: root=', root:18:12, '  iters=', iters);
  SelfTestCheck(Abs(root - 1.521379706804) < 1e-8, 'NewtonRoot');

  iters := 0;
  root := SecantRoot(@TestFunc, 1.0, 2.0, 1e-10, 100, iters);
  WriteLn('Secant:         root=', root:18:12, '  iters=', iters);
  SelfTestCheck(Abs(root - 1.521379706804) < 1e-8, 'SecantRoot');

  iters := 0;
  root := RegulaFalsiRoot(@TestFunc, 1.0, 2.0, 1e-10, 100, iters);
  WriteLn('Regula Falsi:   root=', root:18:12, '  iters=', iters);
  SelfTestCheck(Abs(root - 1.521379706804) < 1e-8, 'RegulaFalsiRoot');

  iters := 0;
  root := BrentRoot(@TestFunc, 1.0, 2.0, 1e-10, 100, iters);
  WriteLn('Brent:          root=', root:18:12, '  iters=', iters);
  SelfTestCheck(Abs(root - 1.521379706804) < 1e-8, 'BrentRoot');

  iters := 0;
  root := MuellerRoot(@TestFunc, 1.0, 1.5, 2.0, 1e-10, 100, iters);
  WriteLn('Mueller:        root=', root:18:12, '  iters=', iters);
  SelfTestCheck(Abs(root - 1.521379706804) < 1e-8, 'MuellerRoot');

  iters := 0;
  root := SteffensenRoot(@TestFunc, 1.5, 1e-10, 100, iters);
  WriteLn('Steffensen:     root=', root:18:12, '  iters=', iters);
  SelfTestCheck(Abs(root - 1.521379706804) < 1e-8, 'SteffensenRoot');

  iters := 0;
  root := AitkenRoot(@TestFunc, 1.5, 1e-10, 100, iters);
  WriteLn('Aitken:         root=', root:18:12, '  iters=', iters);
  SelfTestCheck(Abs(root - 1.521379706804) < 1e-8, 'AitkenRoot');

  WriteLn;
  WriteLn('Bairstow test: x^3 - 6x^2 + 11x - 6 (roots: 1, 2, 3)');
  { Coeffs: index 0=constant, index n=leading. Polynomial: -6 + 11x - 6x^2 + x^3 }
  roots := BairstowRoots([-6.0, 11.0, -6.0, 1.0], 1e-10, 200);
  for i := 0 to High(roots) do
  begin
    if Abs(roots[i].Im) < 1e-8 then
      WriteLn('  root[', i, '] = ', roots[i].Re:12:8)
    else
      WriteLn('  root[', i, '] = ', roots[i].Re:12:8, ' + ', roots[i].Im:12:8, 'i');
  end;
  WriteLn('=== End Self Test ===');
end;

end.
