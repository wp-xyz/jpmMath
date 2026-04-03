unit jpmIntegrationTests;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, fpcunit, testutils, testregistry;

type

  TTestIntegrations = class(TTestCase)
  published
    procedure TestGaussIntegration;
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

function XFUnc(x: Double): Double;
begin
  Result := x;
end;

function XYFunc(x, y: Double): Double;
begin
  Result := x * y;
end;

function XYZFunc(x, y, z: Double): Double;
begin
  Result := x * y * z;
end;

function SinXYZFunc(x, y, z: Double): Double;
begin
  Result := sin(x * y * z);
end;

procedure TTestIntegrations.TestGaussIntegration;
const
  EPS = 1E-10;
var
  n: Integer;
  actual, expected: Double;
begin
  // Integral of x dx between x=0..1
  expected := 0.5;
  for n := 2 to 10 do
  begin
    actual := GaussIntegral(@XFunc, n, 0,1);
    CheckEquals(expected, actual, EPS, 'Gauss Integral f(x)=x, order ' + IntToStr(n) + ': Result mismatch');
  end;

  // Integral x*y dxdy between x=0..1 and y=0..1
  expected := 0.25;
  for n := 2 to 10 do
  begin
    actual := GaussIntegral(@XYFunc, n , 0,1, 0,1);
    CheckEquals(expected, actual, EPS, 'Gauss Integral f(x,y)=x*y, order ' + IntToStr(n) + ': Result mismatch');
  end;

  // Integral x*y*z dxdydz between x=0..1 and y=0..1 and z=0..1
  expected := 0.125;
  for n := 2 to 10 do
  begin
    actual := GaussIntegral(@XYZFunc, n , 0,1, 0,1, 0,1);
    CheckEquals(expected, actual, EPS, 'Gauss Integral f(x,y,z)=x*y*z, order ' + IntToStr(n) + ': Result mismatch');
  end;

  // Integral sin(x*y*z) dxdydz between x=0..1 and y=0..1 and z=0..1
  expected := 0.1224345;  // result provided by Wolfram Alpha
  for n := 2 to 10 do
  begin
    actual := GaussIntegral(@SinXYZFunc, n , 0,1, 0,1, 0,1);
    CheckEquals(expected, actual, 1E-5, 'Gauss Integral f(x,y,z)=sin(x*y*z), order ' + IntToStr(n) + ': Result mismatch');
  end;

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

