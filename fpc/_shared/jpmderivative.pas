unit jpmDerivative;
{$mode objfpc}{$H+}

interface

uses SysUtils, Math, jpmtypes;

type
  TFloatArray  = array of Float;
  TFuncND      = function(var x: TFloatArray; n: integer): Float;
  TVectorFunc  = function(var x: TFloatArray; n: integer;
                           var y: TFloatArray; m: integer): boolean;
  TFloatMatrix = array of TFloatArray;

{ First derivative by central differences: (f(x+h) - f(x-h)) / (2h) }
function Deriv1(f: TFunction1; x, h: Float): Float;

{ Second derivative: (f(x+h) - 2f(x) + f(x-h)) / h^2 }
function Deriv2(f: TFunction1; x, h: Float): Float;

{ Third derivative: (f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)) / (2h^3) }
function Deriv3(f: TFunction1; x, h: Float): Float;

{ Fourth derivative: (f(x+2h) - 4f(x+h) + 6f(x) - 4f(x-h) + f(x-2h)) / h^4 }
function Deriv4(f: TFunction1; x, h: Float): Float;

{ Improved derivative via Richardson extrapolation; order = 1 or 2 }
function DerivRichardson(f: TFunction1; x, h: Float; order: integer): Float;

{ Romberg differentiation for first derivative; error estimate returned in err }
function DerivRomberg(f: TFunction1; x, h: Float; var err: Float): Float;

{ Central difference gradient of an N-dimensional function }
procedure Gradient(f: TFuncND; var x: TFloatArray; n: integer; h: Float;
  var grad: TFloatArray);

{ Numerical Jacobian J[i][j] = d(y_i)/d(x_j) via central differences }
procedure Jacobian(f: TVectorFunc; var x: TFloatArray; n: integer;
  var y: TFloatArray; m: integer; h: Float; var J: TFloatMatrix);

procedure self_test;

implementation

const
  ROMBERG_SIZE = 5;

type
  TRombergRow = array[0..ROMBERG_SIZE - 1] of Float;

{ -----------------------------------------------------------------------
  Basic finite-difference formulas
  ----------------------------------------------------------------------- }

function Deriv1(f: TFunction1; x, h: Float): Float;
begin
  result := (f(x + h) - f(x - h)) / (2.0 * h)
end;

function Deriv2(f: TFunction1; x, h: Float): Float;
begin
  result := (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)
end;

function Deriv3(f: TFunction1; x, h: Float): Float;
begin
  result := (f(x + 2.0 * h) - 2.0 * f(x + h)
             + 2.0 * f(x - h) - f(x - 2.0 * h)) / (2.0 * h * h * h)
end;

function Deriv4(f: TFunction1; x, h: Float): Float;
begin
  result := (f(x + 2.0 * h) - 4.0 * f(x + h) + 6.0 * f(x)
             - 4.0 * f(x - h) + f(x - 2.0 * h)) / (h * h * h * h)
end;

{ helper: dispatch to first or second central difference }
function InternalDeriv(f: TFunction1; x, h: Float; order: integer): Float;
begin
  case order of
    1: result := Deriv1(f, x, h);
    2: result := Deriv2(f, x, h)
  else
    result := Deriv1(f, x, h)
  end
end;

{ -----------------------------------------------------------------------
  Richardson extrapolation: D(h/2) and D(h) combined as (4*D(h/2)-D(h))/3
  ----------------------------------------------------------------------- }
function DerivRichardson(f: TFunction1; x, h: Float; order: integer): Float;
var
  d1, d2: Float;
begin
  d1 := InternalDeriv(f, x, h, order);
  d2 := InternalDeriv(f, x, h * 0.5, order);
  result := (4.0 * d2 - d1) / 3.0
end;

{ -----------------------------------------------------------------------
  Romberg differentiation (5-step table, first derivative)
  Based on: Engeln-Muellges & Uhlig, "Numerical Algorithms with C", 1996
  ----------------------------------------------------------------------- }
function DerivRomberg(f: TFunction1; x, h: Float; var err: Float): Float;
var
  d           : TRombergRow;
  i, j, m     : integer;
  h2, d1, d2  : Float;
begin
  for i := 0 to ROMBERG_SIZE - 1 do
    d[i] := 0.0;

  h2   := 2.0 * h;
  d[0] := (f(x + h) - f(x - h)) / h2;

  for j := 1 to ROMBERG_SIZE - 1 do
  begin
    d1   := d[0];
    h2   := h;
    h    := h * 0.5;
    d[0] := (f(x + h) - f(x - h)) / h2;
    m    := 4;
    for i := 1 to j do
    begin
      d2   := d[i];
      d[i] := (m * d[i - 1] - d1) / (m - 1);
      d1   := d2;
      m    := m * 4
    end
  end;

  err    := Abs(d[ROMBERG_SIZE - 1] - d[ROMBERG_SIZE - 2]);
  result := d[ROMBERG_SIZE - 1]
end;

{ -----------------------------------------------------------------------
  Gradient: grad[i] = (f(x + h*e_i) - f(x - h*e_i)) / (2h)
  ----------------------------------------------------------------------- }
procedure Gradient(f: TFuncND; var x: TFloatArray; n: integer; h: Float;
  var grad: TFloatArray);
var
  i, k     : integer;
  xp, xm   : TFloatArray;
  fp, fm   : Float;
begin
  SetLength(grad, n);
  SetLength(xp,   n);
  SetLength(xm,   n);
  for i := 0 to n - 1 do
  begin
    for k := 0 to n - 1 do
    begin
      xp[k] := x[k];
      xm[k] := x[k]
    end;
    xp[i]   := x[i] + h;
    xm[i]   := x[i] - h;
    fp       := f(xp, n);
    fm       := f(xm, n);
    grad[i]  := (fp - fm) / (2.0 * h)
  end
end;

{ -----------------------------------------------------------------------
  Jacobian: J[i][j] = d(y_i)/d(x_j)  (central differences)
  ----------------------------------------------------------------------- }
procedure Jacobian(f: TVectorFunc; var x: TFloatArray; n: integer;
  var y: TFloatArray; m: integer; h: Float; var J: TFloatMatrix);
var
  row, col, k  : integer;
  xp, xm       : TFloatArray;
  yp, ym       : TFloatArray;
begin
  SetLength(J, m);
  for row := 0 to m - 1 do
    SetLength(J[row], n);
  SetLength(xp, n);
  SetLength(xm, n);
  SetLength(yp, m);
  SetLength(ym, m);
  for col := 0 to n - 1 do
  begin
    for k := 0 to n - 1 do
    begin
      xp[k] := x[k];
      xm[k] := x[k]
    end;
    xp[col] := x[col] + h;
    xm[col] := x[col] - h;
    f(xp, n, yp, m);
    f(xm, n, ym, m);
    for row := 0 to m - 1 do
      J[row][col] := (yp[row] - ym[row]) / (2.0 * h)
  end;
  { keep y in sync with current x }
  f(x, n, y, m)
end;

{ -----------------------------------------------------------------------
  Test functions (global so they satisfy TFunction1 / TFuncND signatures)
  ----------------------------------------------------------------------- }

function TSin(x: Float): Float;
begin
  result := Sin(x)
end;

function TCubic(x: Float): Float;
begin
  result := x * x * x
end;

function TQuartic(x: Float): Float;
begin
  result := x * x * x * x
end;

{ f(x,y) = x^2 + y^2  (used for gradient test) }
function TGradXY(var x: TFloatArray; n: integer): Float;
begin
  { n implicitly defines the valid index range; only x[0..1] used }
  result := x[0] * x[0] + x[1] * x[1];
  if n < 0 then n := 0   { suppress unused-parameter hint }
end;

{ f([x,y]) = [x^2 + y^2, x*y]  (used for Jacobian test) }
function TVec2(var x: TFloatArray; n: integer;
               var y: TFloatArray; m: integer): boolean;
begin
  y[0]   := x[0] * x[0] + x[1] * x[1];
  y[1]   := x[0] * x[1];
  result := true;
  if (n < 0) or (m < 0) then ;   { suppress unused-parameter hints }
end;

{ -----------------------------------------------------------------------
  self_test  — no ReadLn; hard-coded inputs; results written to stdout
  ----------------------------------------------------------------------- }
procedure self_test;
var
  res, err    : Float;
  pt, grad    : TFloatArray;
  jx, jy      : TFloatArray;
  J           : TFloatMatrix;
  i           : integer;
begin
  WriteLn('=== jpmDerivative self_test ===');
  WriteLn;

  res := Deriv1(@TSin, Pi / 4.0, 1e-5);
  WriteLn('Deriv1(sin, Pi/4, 1e-5)         = ', res:12:8,
          '  expected  0.70710678');

  res := Deriv2(@TSin, Pi / 4.0, 1e-5);
  WriteLn('Deriv2(sin, Pi/4, 1e-5)         = ', res:12:8,
          '  expected -0.70710678');

  res := Deriv3(@TCubic, 2.0, 1e-4);
  WriteLn('Deriv3(x^3,  x=2, h=1e-4)       = ', res:12:8,
          '  expected  6.00000000');

  res := Deriv4(@TQuartic, 1.0, 1e-3);
  WriteLn('Deriv4(x^4,  x=1, h=1e-3)       = ', res:12:8,
          '  expected 24.00000000');

  res := DerivRichardson(@TSin, Pi / 4.0, 1e-3, 1);
  WriteLn('DerivRichardson(sin,Pi/4,1e-3,1)= ', res:12:8,
          '  expected  0.70710678');

  res := DerivRichardson(@TSin, Pi / 4.0, 1e-3, 2);
  WriteLn('DerivRichardson(sin,Pi/4,1e-3,2)= ', res:12:8,
          '  expected -0.70710678');

  err := 0.0;
  res := DerivRomberg(@TSin, Pi / 4.0, 0.1, err);
  WriteLn('DerivRomberg(sin, Pi/4, h=0.1)  = ', res:12:8,
          '  err = ', err, '  expected 0.70710678, err < 1e-10');

  WriteLn;
  SetLength(pt, 2);
  pt[0] := 1.0;
  pt[1] := 2.0;
  Gradient(@TGradXY, pt, 2, 1e-5, grad);
  Write('Gradient(x^2+y^2, [1,2], 1e-5)  = [');
  for i := 0 to 1 do
  begin
    if i > 0 then
      Write(', ');
    Write(grad[i]:10:6)
  end;
  WriteLn(']  expected [2.0, 4.0]');

  WriteLn;
  SetLength(jx, 2);
  SetLength(jy, 2);
  jx[0] := 1.0;
  jx[1] := 2.0;
  Jacobian(@TVec2, jx, 2, jy, 2, 1e-5, J);
  WriteLn('Jacobian([x^2+y^2, x*y], x=[1,2], h=1e-5):');
  WriteLn('  J[0] = [', J[0][0]:8:4, ', ', J[0][1]:8:4,
          ']  expected [2.0, 4.0]');
  WriteLn('  J[1] = [', J[1][0]:8:4, ', ', J[1][1]:8:4,
          ']  expected [2.0, 1.0]');

  WriteLn;
  WriteLn('=== self_test done ===')
end;

end.
