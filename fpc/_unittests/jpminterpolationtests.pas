unit jpmInterpolationTests;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, fpcunit, testutils, testregistry;

type

  TTestInterpolations = class(TTestCase)
  private
    SinDataX: array[0..13] of Double;
    SinDataY: array[0..13] of Double;
  protected
    procedure SetUp; override;
    procedure TearDown; override;
  published
    procedure TestAkimaInterpolation;
  end;

implementation

uses
  jpmInterpolation;

procedure TTestInterpolations.SetUp;
begin
  {Input sine table
  -----------------------------------------------------------------
   Sine table values from  Handbook of mathematical functions
   by M. Abramowitz and I.A. Stegun, NBS, june 1964
  ----------------------------------------------------------------}
  SinDataX[ 0]:=0.000;   SinDataY[ 0]:=0.00000000;
  SinDataX[ 1]:=0.125;   SinDataY[ 1]:=0.12467473;
  SinDataX[ 2]:=0.217;   SinDataY[ 2]:=0.21530095;
  SinDataX[ 3]:=0.299;   SinDataY[ 3]:=0.29456472;
  SinDataX[ 4]:=0.376;   SinDataY[ 4]:=0.36720285;
  SinDataX[ 5]:=0.450;   SinDataY[ 5]:=0.43496553;
  SinDataX[ 6]:=0.520;   SinDataY[ 6]:=0.49688014;
  SinDataX[ 7]:=0.589;   SinDataY[ 7]:=0.55552980;
  SinDataX[ 8]:=0.656;   SinDataY[ 8]:=0.60995199;
  SinDataX[ 9]:=0.721;   SinDataY[ 9]:=0.66013615;
  SinDataX[10]:=0.7853981634; SinDataY[10]:=0.7071067812;
  SinDataX[11]:=0.849;   SinDataY[11]:=0.75062005;
  SinDataX[12]:=0.911;   SinDataY[12]:=0.79011709;
  SinDataX[13]:=0.972;   SinDataY[13]:=0.82601466;
end;

procedure TTestInterpolations.TearDown;
begin
  //
end;

procedure TTestInterpolations.TestAkimaInterpolation;
const
  EPS = 1E-5;
var
  x, yExpected, yActual: Double;
  err: String;
begin
  // Exactly at left edge of allowed table range
  x := 0.0;
  yExpected := sin(x);
  yActual := AkimaInterpolation(SinDataX, SinDataY, 0, -1, x, err);
  CheckEquals(yExpected, yActual, EPS, 'Akima interpolation #1 result mismatch for x = ' + FloatToStr(x));
  CheckEquals('', err, 'Akima interpolation #1 error result mismatch for x = ' + FloatToStr(x));

  // Close to left edge of allowed table range, but still ok
  x := 0.001;
  yExpected := sin(x);
  yActual := AkimaInterpolation(SinDataX, SinDataY, 0, -1, x, err);
  CheckEquals(yExpected, yActual, EPS, 'Akima interpolation #2 result mismatch for x = ' + FloatToStr(x));
  CheckEquals('', err, 'Akima interpolation #2 error result mismatch for x = ' + FloatToStr(x));

  // Some where in the center of the allowed table range
  x := 0.5;
  yExpected := sin(x);
  yActual := AkimaInterpolation(SinDataX, SinDataY, 0, -1, x, err);
  CheckEquals(yExpected, yActual, EPS, 'Akima interpolation #3 result mismatch for x = ' + FloatToStr(x));
  CheckEquals('', err, 'Akima interpolation #3 error result mismatch for x = ' + FloatToStr(x));

  // Close to right edge of allowed table range, but still ok
  x := 0.848;
  yExpected := sin(x);
  yActual := AkimaInterpolation(SinDataX, SinDataY, 0, -1, x, err);
  CheckEquals(yExpected, yActual, EPS, 'Akima interpolation #4 result mismatch for x = ' + FloatToStr(x));
  CheckEquals('', err, 'Akima interpolation #4 error result mismatch for x = ' + FloatToStr(x));

  // Exactly at right edge of allowed range
  x := SinDataX[11];
  yExpected := sin(x);
  yActual := AkimaInterpolation(SinDataX, SinDataY, 0, -1, x, err);
  CheckEquals(yExpected, yActual, EPS, 'Akima interpolation #5 result mismatch for x = ' + FloatToStr(x));
  CheckEquals('', err, 'Akima interpolation #5 error result mismatch for x = ' + FloatToStr(x));

  // Out-of-range error: Smaller than first value in table
  x := -1.0;
  yExpected := sin(x);
  yActual := AkimaInterpolation(SinDataX, SinDataY, 0, -1, x, err);
  CheckNotEquals('', err, 'Akima interpolation #6 error result mismatch for x = ' + FloatToStr(x));

  // Out-of-range error: no 3 points at the right
  x := 0.850;
  yExpected := sin(x);
  yActual := AkimaInterpolation(SinDataX, SinDataY, 0, -1, x, err);
  CheckNotEquals('', err, 'Akima interpolation #7 error result mismatch for x = ' + FloatToStr(x));

  // Out-of-range error: no 3 points at the right
  yExpected := sin(x);
  yActual := AkimaInterpolation(SinDataX, SinDataY, 0, -1, x, err);
  CheckNotEquals('', err, 'Akima interpolation #8 error result mismatch for x = ' + FloatToStr(x));

  // Out-of-range error: larger than last value in table
  yExpected := sin(x);
  yActual := AkimaInterpolation(SinDataX, SinDataY, 0, -1, x, err);
  CheckNotEquals('', err, 'Akima nterpolation error #9 result mismatch for x = ' + FloatToStr(x));
end;

initialization
  RegisterTest(TTestInterpolations);

end.

