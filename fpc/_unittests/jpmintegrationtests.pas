unit jpmIntegrationTests;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, fpcunit, testutils, testregistry;

type

  TTestIntegrations = class(TTestCase)
  published
    procedure TestRombergIntegration;
  end;

implementation

uses
  jpmIntegration;

function SinFunc(x: Double): Double;
begin
  Result := sin(x);
end;

function SinExpFunc(x: double): Double;
begin
  Result := sin(x) * exp(-x);
end;

procedure TTestIntegrations.TestRombergIntegration;
const
  EPS = 1E-10;
var
  actual, expected: Double;
  prec: Double;
  nIt: Integer;
begin
  // Integral of sin(x) from 0 to pi
  expected := 2.0;
  actual := RombergIntegral(@SinFunc, 0, pi, EPS, prec, nIt, 0, 50);
  CheckEquals(expected, actual, EPS, 'Romberg Integration Test #1: Integral sin(x) mismatch');

  // Integral of sin(x)*exp(-x) from 0 to 10
  expected := 0.5 - (sin(10.0) + cos(10.0)) / (2.0 * exp(10));   // Result from Wolfram Alpha
  actual := RombergIntegral(@SinExpFunc, 0, 10, EPS, prec, nIt, 0, 50);
  CheckEquals(expected, actual, EPS, 'Romberg Integration Test #2: Integral sin(x)*exp(-x) mismatch');
end;

initialization
  RegisterTest(TTestIntegrations);

end.

