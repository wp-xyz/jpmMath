{ Integration routines extracted from Jean-Pierre Moreaus's programs }

unit jpmIntegration;

{$mode objFPC}{$H+}

interface

uses
  Classes, SysUtils, Math,
  jpmTypes;

function GaussIntegral(AFunc: TFunction1; Order: Integer;
  x1, x2: Float): Float;
function GaussIntegral(AFunc: TFunction2; Order: Integer;
  x1, x2, y1, y2: Float): Float;
function GaussIntegral(AFunc: TFunction3; Order: Integer;
  x1, x2, y1, y2, z1, z2: Float): Float;

function RombergIntegral(AFunc: TFunction1; a, b, Prec: Float; var ObtPrec: Float;
  var Iterations: integer; MinIterations, MaxIterations: integer): Float;

function SimpsonIntegral(f: TFunction1; a, b: Float; n: integer): Float;
function AdaptiveSimpson(f: TFunction1; a, b, tol: Float): Float;
function GaussLegendre(f: TFunction1; a, b: Float; n: integer): Float;
procedure GetGLNodesWeights(n: integer; var nodes, weights: TFloatArray);
function TrapezoidIntegral(f: TFunction1; a, b: Float; n: integer): Float;
procedure self_test;

implementation

{**************************************************************
* Integration by Gauss method of a real function F=F(X) or    *
* F=F(X,Y) or F=F(X,Y,Z). The integral is calculated by using *
* from 2 to 10 Gauss points.                                  *
* ----------------------------------------------------------- *
* INPUTS:                                                     *
*         AFunc  Function of 1, 2 or 3 real variables         *
*         Order  Number of Gauss points, order of polynomial  *
'        x1, x2  begin and end values of x variable           *
'        y1, y2  begin and end values of y variable           *
'        z1, z2  begin and end values of z variable           *
*                                                             *
* RETURNED VALUE: the integral of F(X) from x1 to x2 (or      *
*   F(x,y) from x1 to x2 and y1 to y2, or F(x, y, z) from     *
*   x1 to x2 and y1 to y2 and z1 to z2.                       *
* ----------------------------------------------------------- *
* Ref.: "Mécanique des vibrations linéaires By  M. Lalanne,   *
*       P. Berthier and J. Der Hagopian, Masson, Paris, 1980" *
*       [BIBLI 16].                                           *
* ----------------------------------------------------------- *
*                          TPW Version By J-P Moreau, Paris.  *
*                                          (www.jpmoreau.fr)  *
*                                                             *
* Extracted as stand-alone procedure by W. Pamler             *
**************************************************************}
procedure CalcGaussCoeffs(n: Integer; var A, H: TFloatArray);
var
  i, j: Integer;
begin
  SetLength(A, n);
  SetLength(H, n);
  case n of
    2: begin
         A[0] := -0.5773502691896257;  H[0] := 1.0;
         A[1] :=  0.5773502691896257;  H[1] := 1.0;
       end;
    3: begin
         A[0] := -0.7745966692414834;  H[0] := 0.5555555555555556;
         A[1] :=  0.0;                 H[1] := 0.8888888888888888;
         A[2] :=  0.7745966692414834;  H[2] := 0.5555555555555556;
       end;
    4: begin
         A[0] := -0.8611363115940526;  H[0] := 0.3478548451374538;
         A[1] := -0.3399810435848563;  H[1] := 0.6521451548625461;
         A[2] :=  0.3399810435848563;  H[2] := 0.6521451548625461;
         A[3] :=  0.8611363115940526;  H[3] := 0.3478548451374538;
       end;
    5: begin
         A[0] := -0.9061798459386640;  H[0] := 0.2369268850561891;
         A[1] := -0.5384693101056831;  H[1] := 0.4786286704993665;
         A[2] :=  0.0;                 H[2] := 0.5688888888888889;
         A[3] :=  0.5384693101056831;  H[3] := 0.4786286704993665;
         A[4] :=  0.9061798459386640;  H[4] := 0.2369268850561891;
       end;
    6: begin
         A[0] := -0.9324695142031521;  H[0] := 0.1713244923791704;
         A[1] := -0.6612093864662645;  H[1] := 0.3607615730481386;
         A[2] := -0.2386191860831969;  H[2] := 0.4679139345726910;
         A[3] :=  0.2386191860831969;  H[3] := 0.4679139345726910;
         A[4] :=  0.6612093864662645;  H[4] := 0.3607615730481386;
         A[5] :=  0.9324695142031521;  H[5] := 0.1713244923791704;
       end;
    7: begin
         A[0] := -0.9491079123427585;  H[0] := 0.1294849661688697;
         A[1] := -0.7415311855993945;  H[1] := 0.2797053914892766;
         A[2] := -0.4058451513773972;  H[2] := 0.3818300505051189;
         A[3] :=  0.0;                 H[3] := 0.4179591836734694;
         A[4] :=  0.4058451513773972;  H[4] := 0.3818300505051189;
         A[5] :=  0.7415311855993945;  H[5] := 0.2797053914892766;
         A[6] :=  0.9491079123427585;  H[6] := 0.1294849661688697;
       end;
    8: begin
         A[0] := -0.9602898564975363;  H[0] := 0.1012285362903763;
         A[1] := -0.7966664774136267;  H[1] := 0.2223810344533745;
         A[2] := -0.5255324099163290;  H[2] := 0.3137066458778873;
         A[3] := -0.1834346424956498;  H[3] := 0.3626837833783620;
         A[4] :=  0.1834346424956498;  H[4] := 0.3626837833783620;
         A[5] :=  0.5255324099163290;  H[5] := 0.3137066458778873;
         A[6] :=  0.7966664774136267;  H[6] := 0.2223810344533745;
         A[7] :=  0.9602898564975363;  H[7] := 0.1012285362903763;
       end;
    9: begin
         A[0] := -0.9681602395076261;  H[0] := 0.0812743883615744;
         A[1] := -0.8360311073266358;  H[1] := 0.1806481606948574;
         A[2] := -0.6133714327005904;  H[2] := 0.2606106964029354;
         A[3] := -0.3242534234038089;  H[3] := 0.3123470770400029;
         A[4] :=  0.0;                 H[4] := 0.3302393550012598;
         A[5] :=  0.3242534234038089;  H[5] := 0.3123470770400029;
         A[6] :=  0.6133714327005904;  H[6] := 0.2606106964029354;
         A[7] :=  0.8360311073266358;  H[7] := 0.1806481606948574;
         A[8] :=  0.9681602395076261;  H[8] := 0.0812743883615744;
       end;
   10: begin
         A[0] := -0.9739065285171717;  H[0] := 0.0666713443086881;
         A[1] := -0.8650633666889845;  H[1] := 0.1494513491505806;
         A[2] := -0.6794095682990244;  H[2] := 0.2190863625159820;
         A[3] := -0.4333953941292472;  H[3] := 0.2692667193099963;
         A[4] := -0.1488743389816312;  H[4] := 0.2955242247147529;
         A[5] :=  0.1488743389816312;  H[5] := 0.2955242247147529;
         A[6] :=  0.4333953941292472;  H[6] := 0.2692667193099963;
         A[7] :=  0.6794095682990244;  H[7] := 0.2190863625159820;
         A[8] :=  0.8650633666889845;  H[8] := 0.1494513491505806;
         A[9] :=  0.9739065285171717;  H[9] := 0.0666713443086881;
       end;
  end;
end;

function GaussIntegral(AFunc: TFunction1; Order:Integer;
  x1, x2: Float): Float;
var
  A: TFloatArray = nil;
  H: TFloatArray = nil;
  x: Float;
  a1, b1: Float;
  i: Integer;
begin
  Result := 0.0;
  CalcGaussCoeffs(Order, A, H);
  a1 := -(x1 - x2) / 2;
  b1 :=  (x1 + x2) / 2;
  for i := 0 to Order-1 do
  begin
    x := a1 * A[i] + b1;
    Result := Result + AFunc(x) * H[i] * a1;
  end;
end;

function GaussIntegral(AFunc: TFunction2; Order: Integer;
  x1, x2, y1, y2: Float): Float;
var
  A: Array of Float = nil;
  H: Array of Float = nil;
  x, y: Float;
  a1, a2, b1, b2: Float;
  i, j: Integer;
begin
  Result := 0.0;
  CalcGaussCoeffs(Order, A, H);
  a1 := -(x1 - x2) / 2;
  b1 :=  (x1 + x2) / 2;
  a2 := -(y1 - y2) / 2;
  b2 :=  (y1 + y2) / 2;
  for i := 0 to Order - 1 do
  begin
    x := a1 * A[i] + b1;
    for j := 0 to Order - 1 do
    begin
      y := a2 * A[j] + b2;
      Result := Result + AFunc(x,y) * H[i] * H[j] * a1 * a2;
    end;
  end;
end;

function GaussIntegral(AFunc: TFunction3; Order: Integer;
  x1, x2, y1, y2, z1, z2: Float): Float;
var
  A: Array of Float = nil;
  H: Array of Float = nil;
  x, y, z: Float;
  a1, a2, a3, b1, b2, b3: Float;
  i, j, k: Integer;
begin
  Result := 0.0;
  CalcGaussCoeffs(Order, A, H);
  a1 := -(x1 - x2) / 2;
  b1 :=  (x1 + x2) / 2;
  a2 := -(y1 - y2) / 2;
  b2 :=  (y1 + y2) / 2;
  a3 := -(z1 - z2) / 2;
  b3 :=  (z1 + z2) / 2;
  for i := 0 to Order-1 do
  begin
    x := a1 * A[i] + b1;
    for j := 0 to Order-1 do
    begin
      y := a2 * A[j] + b2;
      for k := 0 to Order-1 do
      begin
        z := a3 * A[k] + b3;
        Result := Result + AFunc(x,y,z) * H[i] * H[j] * H[k] * a1 * a2 * a3;
      end;
    end;
  end;
end;

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
  i, j: integer;
  pas, r, s, ta: Float;
  t: array of array of Float;
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

