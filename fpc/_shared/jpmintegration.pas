{ Integration routines extracted from Jean-Pierre Moreaus's programs }

unit jpmIntegration;

{$mode objFPC}{$H+}

interface

uses
  Classes, SysUtils, Math,
  jpmTypes;

function RombergIntegral(AFunc: TFunction1; a, b, Prec: Float; var ObtPrec: Float;
  var Iterations: integer; MinIterations, MaxIterations: integer): Float;

function SimpsonIntegral(f: TFunction1; a, b: Float; n: integer): Float;
function AdaptiveSimpson(f: TFunction1; a, b, tol: Float): Float;
function GaussLegendre(f: TFunction1; a, b: Float; n: integer): Float;
procedure GetGLNodesWeights(n: integer; var nodes, weights: TFloatArray);
function TrapezoidIntegral(f: TFunction1; a, b: Float; n: integer): Float;
procedure self_test;

implementation

{*******************************************************
* Integral of a function AFunc(X) by Romberg's method  *
* ---------------------------------------------------- *
* INPUTS:                                              *
*              a  begin value of x variable            *
*              b  end value of x variable              *
*           Prec  desired precision                    *
*  MinIterations  minimum number of iterations         *
*  MaxIterations  maximum number of iterations         *
*                                                      *
* OUTPUTS:                                             *
*        ObtPrec  obtained precision for integral      *
*     Iterations  number of iterations done            *
*                                                      *
* RETURNED VALUE  the integral of AFunc(X) from a to b *
*                                                      *
*******************************************************}
function RombergIntegral(AFunc: TFunction1; a, b, Prec: Float; var ObtPrec: Float;
  var Iterations: integer; MinIterations, MaxIterations: integer): Float;
const
  MAXITER = 15;
var
  i,j        : integer;
  pas,r,s,ta : Float;
  t          : array of array of Float;
begin
  if MaxIterations > MAXITER then
    MaxIterations := MAXITER;
  r := AFunc(a);
  ta := (r + AFunc(b)) / 2;
  Iterations := 0;
  pas := b - a;
  SetLength(t, MaxIterations+1, MaxIterations+1);
  t[0,0] := ta*pas;
  repeat
    Inc(Iterations);
    pas := pas/2;
    s := ta;
    for i := 1 to pred(1 shl Iterations) do   {2^n-1}
      s := s + AFunc(a + pas*i);
    t[0, Iterations] := s*pas;
    r := 1;
    for i := 1 to Iterations do
    begin
      r := r*4;
      j := Iterations-i;
      t[i,j] := (r*t[i-1, j+1] - t[i-1, j]) / (r-1)
    end;
    obtprec := abs(t[Iterations, 0] - t[Iterations-1, 0])
  until (Iterations >= MaxIterations) or ((ObtPrec < Prec) and (Iterations >= MinIterations));
  Result := t[Iterations, 0]
end;

{*******************************************************
* Composite Simpson's Rule                             *
* n must be even; if odd, n is incremented by 1        *
*******************************************************}
function SimpsonIntegral(f: TFunction1; a, b: Float; n: integer): Float;
var
  h, s: Float;
  i: integer;
begin
  if odd(n) then
    inc(n);
  h := (b - a) / n;
  s := f(a) + f(b);
  i := 1;
  while i <= n - 1 do
  begin
    if odd(i) then
      s := s + 4 * f(a + i * h)
    else
      s := s + 2 * f(a + i * h);
    inc(i)
  end;
  result := h / 3 * s
end;

{*******************************************************
* Adaptive Simpson's Rule (recursive)                  *
* Splits interval when error > tol, max depth = 50     *
*******************************************************}
function AdaptiveSimpson(f: TFunction1; a, b, tol: Float): Float;

  function SimpsonRule(x1, x2: Float): Float;
  var
    mid: Float;
  begin
    mid := (x1 + x2) / 2;
    result := (x2 - x1) / 6 * (f(x1) + 4 * f(mid) + f(x2))
  end;

  function Recurse(x1, x2, tol0, whole: Float; depth: integer): Float;
  var
    mid, left, right: Float;
  begin
    mid := (x1 + x2) / 2;
    left := SimpsonRule(x1, mid);
    right := SimpsonRule(mid, x2);
    if (depth >= 50) or (abs(left + right - whole) <= 15 * tol0) then
      result := left + right + (left + right - whole) / 15
    else
      result := Recurse(x1, mid, tol0 / 2, left, depth + 1) +
                Recurse(mid, x2, tol0 / 2, right, depth + 1)
  end;

var
  whole: Float;
begin
  whole := SimpsonRule(a, b);
  result := Recurse(a, b, tol, whole, 0)
end;

{*******************************************************
* Fill Gauss-Legendre nodes and weights for order n    *
* Supported: n = 2, 3, 4, 5, 6, 8, 10                 *
*******************************************************}
procedure GetGLNodesWeights(n: integer; var nodes, weights: TFloatArray);
begin
  case n of
    2:
      begin
        setlength(nodes, 2);
        setlength(weights, 2);
        nodes[0]   := -0.5773502691896257;
        nodes[1]   :=  0.5773502691896257;
        weights[0] :=  1.0;
        weights[1] :=  1.0
      end;
    3:
      begin
        setlength(nodes, 3);
        setlength(weights, 3);
        nodes[0]   := -0.7745966692414834;
        nodes[1]   :=  0.0;
        nodes[2]   :=  0.7745966692414834;
        weights[0] :=  0.5555555555555556;
        weights[1] :=  0.8888888888888889;
        weights[2] :=  0.5555555555555556
      end;
    4:
      begin
        setlength(nodes, 4);
        setlength(weights, 4);
        nodes[0]   := -0.8611363115940526;
        nodes[1]   := -0.3399810435848563;
        nodes[2]   :=  0.3399810435848563;
        nodes[3]   :=  0.8611363115940526;
        weights[0] :=  0.3478548451374538;
        weights[1] :=  0.6521451548625461;
        weights[2] :=  0.6521451548625461;
        weights[3] :=  0.3478548451374538
      end;
    5:
      begin
        setlength(nodes, 5);
        setlength(weights, 5);
        nodes[0]   := -0.9061798459366170;
        nodes[1]   := -0.5384693101056831;
        nodes[2]   :=  0.0;
        nodes[3]   :=  0.5384693101056831;
        nodes[4]   :=  0.9061798459366170;
        weights[0] :=  0.2369268850561891;
        weights[1] :=  0.4786286704993665;
        weights[2] :=  0.5688888888888889;
        weights[3] :=  0.4786286704993665;
        weights[4] :=  0.2369268850561891
      end;
    6:
      begin
        setlength(nodes, 6);
        setlength(weights, 6);
        nodes[0]   := -0.9324695142031521;
        nodes[1]   := -0.6612093864662645;
        nodes[2]   := -0.2386191860831969;
        nodes[3]   :=  0.2386191860831969;
        nodes[4]   :=  0.6612093864662645;
        nodes[5]   :=  0.9324695142031521;
        weights[0] :=  0.1713244923791704;
        weights[1] :=  0.3607615730481386;
        weights[2] :=  0.4679139345726910;
        weights[3] :=  0.4679139345726910;
        weights[4] :=  0.3607615730481386;
        weights[5] :=  0.1713244923791704
      end;
    8:
      begin
        setlength(nodes, 8);
        setlength(weights, 8);
        nodes[0]   := -0.9602898564975363;
        nodes[1]   := -0.7966664774136267;
        nodes[2]   := -0.5255324099163290;
        nodes[3]   := -0.1834346424956498;
        nodes[4]   :=  0.1834346424956498;
        nodes[5]   :=  0.5255324099163290;
        nodes[6]   :=  0.7966664774136267;
        nodes[7]   :=  0.9602898564975363;
        weights[0] :=  0.1012285362903763;
        weights[1] :=  0.2223810344533745;
        weights[2] :=  0.3137066458778873;
        weights[3] :=  0.3626837833783620;
        weights[4] :=  0.3626837833783620;
        weights[5] :=  0.3137066458778873;
        weights[6] :=  0.2223810344533745;
        weights[7] :=  0.1012285362903763
      end;
    10:
      begin
        setlength(nodes, 10);
        setlength(weights, 10);
        nodes[0]   := -0.9739065285171717;
        nodes[1]   := -0.8650633666889845;
        nodes[2]   := -0.6794095682990244;
        nodes[3]   := -0.4333953941292472;
        nodes[4]   := -0.1488743389816312;
        nodes[5]   :=  0.1488743389816312;
        nodes[6]   :=  0.4333953941292472;
        nodes[7]   :=  0.6794095682990244;
        nodes[8]   :=  0.8650633666889845;
        nodes[9]   :=  0.9739065285171717;
        weights[0] :=  0.0666713443086881;
        weights[1] :=  0.1494513491505806;
        weights[2] :=  0.2190863625159820;
        weights[3] :=  0.2692667193099963;
        weights[4] :=  0.2955242247147529;
        weights[5] :=  0.2955242247147529;
        weights[6] :=  0.2692667193099963;
        weights[7] :=  0.2190863625159820;
        weights[8] :=  0.1494513491505806;
        weights[9] :=  0.0666713443086881
      end
    else
      begin
        setlength(nodes, 0);
        setlength(weights, 0)
      end
  end
end;

{*******************************************************
* Gauss-Legendre Quadrature                            *
* Supports n = 2, 3, 4, 5, 6, 8, 10                   *
* Transform [a,b] -> [-1,1]: x = (b+a)/2 + (b-a)/2*t  *
*******************************************************}
function GaussLegendre(f: TFunction1; a, b: Float; n: integer): Float;
var
  nodes, weights: TFloatArray;
  mid, half, s: Float;
  i: integer;
begin
  GetGLNodesWeights(n, nodes, weights);
  mid  := (b + a) / 2;
  half := (b - a) / 2;
  s := 0;
  for i := 0 to length(nodes) - 1 do
    s := s + weights[i] * f(mid + half * nodes[i]);
  result := half * s
end;

{*******************************************************
* Composite Trapezoidal Rule with n intervals          *
*******************************************************}
function TrapezoidIntegral(f: TFunction1; a, b: Float; n: integer): Float;
var
  h, s: Float;
  i: integer;
begin
  h := (b - a) / n;
  s := (f(a) + f(b)) / 2;
  for i := 1 to n - 1 do
    s := s + f(a + i * h);
  result := h * s
end;

{*******************************************************
* Self-test: integrate sin(x) from 0 to Pi = 2.0       *
*******************************************************}

function TestSin(x: Float): Float;
begin
  result := sin(x)
end;

function TestExp(x: Float): Float;
begin
  result := exp(x)
end;

procedure self_test;
var
  res, err, obtPrec: Float;
  iters: integer;
begin
  writeln('=== jpmIntegration self_test ===');
  writeln;

  { SimpsonIntegral }
  res := SimpsonIntegral(@TestSin, 0, Pi, 100);
  err := abs(res - 2.0);
  writeln('SimpsonIntegral(sin, 0, Pi, 100)    = ', res:20:15, '  error=', FormatFloat('0.000E+00', err));
  if err < 2e-8 then
    writeln('  PASS')
  else
    writeln('  FAIL (expected error < 2e-8)');

  { AdaptiveSimpson }
  res := AdaptiveSimpson(@TestSin, 0, Pi, 1e-8);
  err := abs(res - 2.0);
  writeln('AdaptiveSimpson(sin, 0, Pi, 1e-8)   = ', res:20:15, '  error=', FormatFloat('0.000E+00', err));
  if err < 1e-8 then
    writeln('  PASS')
  else
    writeln('  FAIL (expected error < 1e-8)');

  { GaussLegendre n=5 }
  res := GaussLegendre(@TestSin, 0, Pi, 5);
  err := abs(res - 2.0);
  writeln('GaussLegendre(sin, 0, Pi, 5)        = ', res:20:15, '  error=', FormatFloat('0.000E+00', err));
  if err < 1e-6 then
    writeln('  PASS')
  else
    writeln('  FAIL (expected error < 1e-6)');

  { TrapezoidIntegral }
  res := TrapezoidIntegral(@TestSin, 0, Pi, 1000);
  err := abs(res - 2.0);
  writeln('TrapezoidIntegral(sin, 0, Pi, 1000) = ', res:20:15, '  error=', FormatFloat('0.000E+00', err));
  if err < 1e-5 then
    writeln('  PASS')
  else
    writeln('  FAIL (expected error < 1e-5)');

  { RombergIntegral regression }
  res := RombergIntegral(@TestSin, 0, Pi, 1e-10, obtPrec, iters, 3, 15);
  err := abs(res - 2.0);
  writeln('RombergIntegral(sin, 0, Pi)         = ', res:20:15, '  error=', FormatFloat('0.000E+00', err), '  iters=', iters);
  if err < 1e-8 then
    writeln('  PASS')
  else
    writeln('  FAIL (expected error < 1e-8)');

  { GaussLegendre exp(x) from 0..1 = e-1 }
  res := GaussLegendre(@TestExp, 0, 1, 4);
  err := abs(res - (exp(1) - 1));
  writeln('GaussLegendre(exp, 0, 1, 4)         = ', res:20:15, '  error=', FormatFloat('0.000E+00', err));
  if err < 1e-8 then
    writeln('  PASS')
  else
    writeln('  FAIL (expected error < 1e-8)');

  writeln;
  writeln('=== self_test done ===')
end;

end.

