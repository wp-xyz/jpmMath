unit jpmChebyshev;
{$mode objfpc}{$H+}

interface

uses SysUtils, Math, jpmtypes;

const
  CHEB_MAX = 50;

type
  TChebCoeffs = array[1..CHEB_MAX] of Float;

{ Compute n Chebyshev coefficients for f on [a,b] }
procedure ChebFit(f: TFunction1; a, b: Float; var c: TChebCoeffs; n: integer);

{ Evaluate Chebyshev approximation at x using m terms }
function ChebEval(a, b: Float; var c: TChebCoeffs; m: integer; x: Float): Float;

{ Compute Chebyshev coefficients of the integral (integral vanishes at a) }
procedure ChebInteg(a, b: Float; var c: TChebCoeffs; var cint: TChebCoeffs; n: integer);

{ Compute Chebyshev coefficients of the derivative }
procedure ChebDeriv(a, b: Float; var c: TChebCoeffs; var cder: TChebCoeffs; n: integer);

procedure self_test;

implementation

procedure ChebFit(f: TFunction1; a, b: Float; var c: TChebCoeffs; n: integer);
var
  bma, bpa, fac, y, sum: Float;
  fvals: TChebCoeffs;
  j, k: integer;
begin
  bma := 0.5 * (b - a);
  bpa := 0.5 * (b + a);
  for k := 1 to n do
  begin
    y       := cos(Pi * (k - 0.5) / n);
    fvals[k] := f(y * bma + bpa)
  end;
  fac := 2.0 / n;
  for j := 1 to n do
  begin
    sum := 0.0;
    for k := 1 to n do
      sum := sum + fvals[k] * cos(Pi * (j - 1) * (k - 0.5) / n);
    c[j] := fac * sum
  end
end;

function ChebEval(a, b: Float; var c: TChebCoeffs; m: integer; x: Float): Float;
var
  d, dd, sv, y, y2: Float;
  j: integer;
begin
  d  := 0.0;
  dd := 0.0;
  y  := (2.0 * x - a - b) / (b - a);
  y2 := 2.0 * y;
  for j := m downto 2 do
  begin
    sv := d;
    d  := y2 * d - dd + c[j];
    dd := sv
  end;
  result := y * d - dd + 0.5 * c[1]
end;

procedure ChebInteg(a, b: Float; var c: TChebCoeffs; var cint: TChebCoeffs; n: integer);
var
  con, fac, sum: Float;
  j: integer;
begin
  con := 0.25 * (b - a);
  sum := 0.0;
  fac := 1.0;
  for j := 2 to n - 1 do
  begin
    cint[j] := con * (c[j - 1] - c[j + 1]) / (j - 1);
    sum      := sum + fac * cint[j];
    fac      := -fac
  end;
  cint[n] := con * c[n - 1] / (n - 1);
  sum      := sum + fac * cint[n];
  cint[1]  := 2.0 * sum
end;

procedure ChebDeriv(a, b: Float; var c: TChebCoeffs; var cder: TChebCoeffs; n: integer);
var
  con: Float;
  j: integer;
begin
  cder[n]     := 0.0;
  cder[n - 1] := 2.0 * (n - 1) * c[n];
  if n >= 3 then
    for j := n - 2 downto 1 do
      cder[j] := cder[j + 2] + 2.0 * j * c[j + 1];
  con := 2.0 / (b - a);
  for j := 1 to n do
    cder[j] := cder[j] * con
end;

{ -----------------------------------------------------------------------
  Test functions (global scope to satisfy TFunction1 signature)
  ----------------------------------------------------------------------- }

function TSin(x: Float): Float;
begin
  result := sin(x)
end;

{ -----------------------------------------------------------------------
  self_test — no ReadLn; hard-coded inputs; results written to stdout
  ----------------------------------------------------------------------- }
procedure self_test;
const
  N   = 20;
  A   = -1.0;
  B   =  1.0;
var
  c, cder, cint: TChebCoeffs;
  val, expected, diff: Float;
  j: integer;
  coeffDecayOk: boolean;
begin
  WriteLn('=== jpmChebyshev self_test ===');
  WriteLn;

  { --- ChebFit --- }
  ChebFit(@TSin, A, B, c, N);

  { ChebEval at 0.5 }
  val      := ChebEval(A, B, c, N, 0.5);
  expected := sin(0.5);
  WriteLn('ChebEval(0.5)          = ', val:12:8,
          '  expected ', expected:12:8, '  diff ', Abs(val - expected):10:2, 'e');

  { ChebEval at -1.0 }
  val      := ChebEval(A, B, c, N, -1.0);
  expected := sin(-1.0);
  WriteLn('ChebEval(-1.0)         = ', val:12:8,
          '  expected ', expected:12:8, '  diff ', Abs(val - expected):10:2, 'e');

  { ChebEval at 0.0 }
  val      := ChebEval(A, B, c, N, 0.0);
  expected := sin(0.0);
  WriteLn('ChebEval(0.0)          = ', val:12:8,
          '  expected ', expected:12:8, '  diff ', Abs(val - expected):10:2, 'e');

  { Coefficient decay check }
  coeffDecayOk := true;
  for j := 15 to N do
    if Abs(c[j]) >= 1.0e-10 then
      coeffDecayOk := false;
  if coeffDecayOk then
    WriteLn('c[15..20] all < 1e-10  : PASS')
  else
  begin
    WriteLn('c[15..20] decay check  : FAIL');
    for j := 15 to N do
      WriteLn('  c[', j, '] = ', c[j])
  end;

  WriteLn;

  { --- ChebDeriv --- }
  ChebDeriv(A, B, c, cder, N);
  val      := ChebEval(A, B, cder, N, 0.5);
  expected := cos(0.5);
  WriteLn('ChebDeriv eval(0.5)    = ', val:12:8,
          '  expected ', expected:12:8, '  diff ', Abs(val - expected):10:2, 'e');

  WriteLn;

  { --- ChebInteg: integral of sin from 0 to 0.5 = -cos(0.5)+cos(0) --- }
  ChebInteg(A, B, c, cint, N);
  val      := ChebEval(A, B, cint, N, 0.5) - ChebEval(A, B, cint, N, 0.0);
  expected := -cos(0.5) + cos(0.0);
  diff     := Abs(val - expected);
  WriteLn('ChebInteg [0..0.5]     = ', val:12:8,
          '  expected ', expected:12:8, '  diff ', diff:10:2, 'e');

  WriteLn;
  WriteLn('=== self_test done ===')
end;

end.
