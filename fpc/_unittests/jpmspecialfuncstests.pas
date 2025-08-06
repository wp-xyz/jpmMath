unit jpmSpecialFuncsTests;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, fpcunit, testutils, testregistry,
  jpmTypes, jpmSpecialFunc;

type

  TTestGamma = class(TTestCase)
  private
    Expected: array[1..11] of Float;
    ExpectedInt: array[1..10] of Float;
  protected
    procedure SetUp; override;
    procedure TearDown; override;
  published
    procedure TestIntegerArguments;
    procedure TestFloatArguments;
  end;

implementation

procedure TTestGamma.SetUp;
begin
  // Values calculated by Wolfram Alpha
  expected[ 1] := 9.5135076986687318362924871772654021925505786260883773430500007704342;  // x = 0.1
  expected[ 2] := 4.5908437119988030532047582759291520034341099982934030177888531362300;  // x = 0.2
  expected[ 3] := 2.9915689876875906283125165159049177911128060249217151127441196509563;  // x = 0.3
  expected[ 4] := 2.2181595437576882230590540219076794507705665017714695822419777526461;  // x = 0.4
  expected[ 5] := 1.7724538509055160272981674833411451827975494561223871282138077898529;  // x = 0.5
  expected[ 6] := 1.4891922488128171023943333883213422813205990387599247353386795640450;  // x = 0.6
  expected[ 7] := 1.2980553326475577856811711791528116177841411705539462479216453882541;  // x = 0.7
  expected[ 8] := 1.1642297137253033736363209382684586931419617688911877529848944678618;  // x = 0.8
  expected[ 9] := 1.0686287021193193548973053356944807781698387850609731790493706839815;  // x = 0.9
  expected[10] := 1.0;                                                                    // x = 1.0
  expected[11] := 0.9513507698668731836292487177265402192550578626088377343050000770434;  // x = 1.1

  expectedInt[1] := 1.0;   // x = 1.0
  expectedInt[2] := 1.0;
  expectedInt[3] := 2.0;
  expectedInt[4] := 6.0;
  expectedInt[5] := 24.0;
  expectedInt[6] := 120.0;
  expectedInt[7] := 720.0;
  expectedInt[8] := 5040.0;
  expectedInt[9] := 40320.0;
  expectedInt[10] := 362880.0;
end;

procedure TTestGamma.TearDown;
begin
  //
end;

procedure TTestGamma.TestIntegerArguments;
var
  x, y: Double;
  i: Integer;
begin
  for i := Low(ExpectedInt) to High(ExpectedInt) do
  begin
    x := 1.0 * i;
    y := Gamma(x);
    CheckEquals(ExpectedInt[i], y, 'Gamma(x) calculation mismatch for x = ' + FloatToStr(x) +
      ', difference ' + Format('%.3e', [y - ExpectedInt[i]]));
  end;
end;

procedure TTestGamma.TestFloatArguments;
const
  EPS = 1E-15;
var
  x, y: Double;
  i: Integer;
  dx: Float;
begin
  dx := 0.1;
  for i := Low(Expected) to High(Expected) do
  begin
    x := dx * i;
    y := Gamma(x);
    CheckEquals(Expected[i], y, EPS, 'Gamma(x) calculation mismatch for x = ' + FloatToStr(x) +
      ', difference ' + Format('%.3e', [y - Expected[i]]));
  end;
end;


initialization
  RegisterTest('Special functions/Gamma function', TTestGamma);

end.

