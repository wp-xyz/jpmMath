unit jpmComplex;
{$mode objfpc}{$H+}

{-------------------------------------------------------------------------------
                               Unit jpmComplex
     Complex number arithmetic and transcendental functions.
     TComplex record stores Cartesian (re, im) form only.
-------------------------------------------------------------------------------}

interface

uses
  SysUtils, Math, jpmtypes;

type
  TComplex = record
    re: Float;  { real part }
    im: Float;  { imaginary part }
  end;

{ 16.1 Construction and conversion }
function  CSet(re, im: Float): TComplex;
function  CFromPolar(r, theta: Float): TComplex;
function  CModulus(z: TComplex): Float;
function  CArgument(z: TComplex): Float;
function  CConjugate(z: TComplex): TComplex;

{ 16.2 Arithmetic }
function  CAdd(z1, z2: TComplex): TComplex;
function  CSub(z1, z2: TComplex): TComplex;
function  CMul(z1, z2: TComplex): TComplex;
function  CDiv(z1, z2: TComplex): TComplex;
function  CNeg(z: TComplex): TComplex;
function  CScale(z: TComplex; s: Float): TComplex;

{ 16.3 Elementary functions }
function  CExp(z: TComplex): TComplex;
function  CLn(z: TComplex): TComplex;
function  CSqrt(z: TComplex): TComplex;
function  CPow(z: TComplex; n: Float): TComplex;
function  CPowC(z1, z2: TComplex): TComplex;

{ 16.4 Trigonometric }
function  CSin(z: TComplex): TComplex;
function  CCos(z: TComplex): TComplex;
function  CTan(z: TComplex): TComplex;
function  CSinh(z: TComplex): TComplex;
function  CCosh(z: TComplex): TComplex;
function  CTanh(z: TComplex): TComplex;

{ 16.5 Utility }
function  CToStr(z: TComplex; digits: integer): string;
function  CAbs2(z: TComplex): Float;

procedure self_test;

implementation

{ --------------------------------------------------------------------------- }
{ 16.1  Construction and conversion                                            }
{ --------------------------------------------------------------------------- }

function CSet(re, im: Float): TComplex;
begin
  result.re := re;
  result.im := im
end;

function CFromPolar(r, theta: Float): TComplex;
begin
  result.re := r * Cos(theta);
  result.im := r * Sin(theta)
end;

function CModulus(z: TComplex): Float;
begin
  result := Sqrt(z.re * z.re + z.im * z.im)
end;

function CArgument(z: TComplex): Float;
begin
  result := ArcTan2(z.im, z.re)
end;

function CConjugate(z: TComplex): TComplex;
begin
  result.re :=  z.re;
  result.im := -z.im
end;

{ --------------------------------------------------------------------------- }
{ 16.2  Arithmetic                                                             }
{ --------------------------------------------------------------------------- }

function CAdd(z1, z2: TComplex): TComplex;
begin
  result.re := z1.re + z2.re;
  result.im := z1.im + z2.im
end;

function CSub(z1, z2: TComplex): TComplex;
begin
  result.re := z1.re - z2.re;
  result.im := z1.im - z2.im
end;

function CMul(z1, z2: TComplex): TComplex;
begin
  result.re := z1.re * z2.re - z1.im * z2.im;
  result.im := z1.re * z2.im + z1.im * z2.re
end;

function CDiv(z1, z2: TComplex): TComplex;
var
  denom: Float;
begin
  denom := z2.re * z2.re + z2.im * z2.im;
  result.re := (z1.re * z2.re + z1.im * z2.im) / denom;
  result.im := (z1.im * z2.re - z1.re * z2.im) / denom
end;

function CNeg(z: TComplex): TComplex;
begin
  result.re := -z.re;
  result.im := -z.im
end;

function CScale(z: TComplex; s: Float): TComplex;
begin
  result.re := s * z.re;
  result.im := s * z.im
end;

{ --------------------------------------------------------------------------- }
{ 16.3  Elementary functions                                                   }
{ --------------------------------------------------------------------------- }

function CExp(z: TComplex): TComplex;
var
  er: Float;
begin
  er := Exp(z.re);
  result.re := er * Cos(z.im);
  result.im := er * Sin(z.im)
end;

function CLn(z: TComplex): TComplex;
begin
  result.re := Ln(CModulus(z));
  result.im := CArgument(z)
end;

function CSqrt(z: TComplex): TComplex;
var
  r, theta: Float;
begin
  r     := CModulus(z);
  theta := CArgument(z);
  result := CFromPolar(Sqrt(r), theta / 2.0)
end;

function CPow(z: TComplex; n: Float): TComplex;
var
  r, theta: Float;
begin
  r     := CModulus(z);
  theta := CArgument(z);
  result := CFromPolar(Power(r, n), n * theta)
end;

function CPowC(z1, z2: TComplex): TComplex;
begin
  result := CExp(CMul(z2, CLn(z1)))
end;

{ --------------------------------------------------------------------------- }
{ 16.4  Trigonometric                                                          }
{ --------------------------------------------------------------------------- }

function CSin(z: TComplex): TComplex;
begin
  result.re := Sin(z.re) * Cosh(z.im);
  result.im := Cos(z.re) * Sinh(z.im)
end;

function CCos(z: TComplex): TComplex;
begin
  result.re :=  Cos(z.re) * Cosh(z.im);
  result.im := -Sin(z.re) * Sinh(z.im)
end;

function CTan(z: TComplex): TComplex;
begin
  result := CDiv(CSin(z), CCos(z))
end;

function CSinh(z: TComplex): TComplex;
begin
  result.re := Sinh(z.re) * Cos(z.im);
  result.im := Cosh(z.re) * Sin(z.im)
end;

function CCosh(z: TComplex): TComplex;
begin
  result.re := Cosh(z.re) * Cos(z.im);
  result.im := Sinh(z.re) * Sin(z.im)
end;

function CTanh(z: TComplex): TComplex;
begin
  result := CDiv(CSinh(z), CCosh(z))
end;

{ --------------------------------------------------------------------------- }
{ 16.5  Utility                                                                }
{ --------------------------------------------------------------------------- }

function CAbs2(z: TComplex): Float;
begin
  result := z.re * z.re + z.im * z.im
end;

function CToStr(z: TComplex; digits: integer): string;
var
  fmt: string;
  imSign: string;
  absIm: Float;
begin
  fmt   := '%.' + IntToStr(digits) + 'f';
  absIm := Abs(z.im);
  if z.im < 0.0 then
    imSign := ' - '
  else
    imSign := ' + ';
  result := Format(fmt, [z.re]) + imSign + Format(fmt, [absIm]) + 'i'
end;

{ --------------------------------------------------------------------------- }
{ self_test                                                                    }
{ --------------------------------------------------------------------------- }

procedure self_test;
const
  TOL = 1e-4;

  procedure check(const name: string; computed, expected: Float);
  var
    ok: string;
  begin
    if Abs(computed - expected) < TOL then
      ok := 'OK'
    else
      ok := '*** FAIL ***';
    WriteLn(Format('  %-38s computed=%10.6f  expected=%10.6f  %s',
                   [name, computed, expected, ok]));
    if ok = '*** FAIL ***' then
      SelfTestFail(name + ': computed=' + FloatToStr(computed) + ' expected=' + FloatToStr(expected));
  end;

  procedure checkZ(const name: string; z: TComplex; ere, eim: Float);
  var
    ok: string;
  begin
    if (Abs(z.re - ere) < TOL) and (Abs(z.im - eim) < TOL) then
      ok := 'OK'
    else
      ok := '*** FAIL ***';
    WriteLn(Format('  %-38s (%9.5f, %9.5f)  expected (%9.5f, %9.5f)  %s',
                   [name, z.re, z.im, ere, eim, ok]));
    if ok = '*** FAIL ***' then
      SelfTestFail(name + ': computed=(' + FloatToStr(z.re) + ',' + FloatToStr(z.im) +
                   ') expected=(' + FloatToStr(ere) + ',' + FloatToStr(eim) + ')');
  end;

var
  z, z1, z2: TComplex;
begin
  WriteLn('=== jpmComplex self_test ===');
  WriteLn;

  { 16.1 Construction and conversion }
  WriteLn('-- 16.1  Construction and conversion --');
  z := CSet(3.0, 4.0);
  check('CModulus(3+4i)',         CModulus(z),   5.0);
  check('CArgument(3+4i)',        CArgument(z),  ArcTan2(4.0, 3.0));

  z1 := CFromPolar(5.0, ArcTan2(4.0, 3.0));
  check('CFromPolar(5,atan2(4,3)).re',  z1.re,  3.0);
  check('CFromPolar(5,atan2(4,3)).im',  z1.im,  4.0);

  checkZ('CConjugate(3-4i)',      CConjugate(CSet(3.0, -4.0)),  3.0, 4.0);
  WriteLn;

  { 16.2 Arithmetic }
  WriteLn('-- 16.2  Arithmetic --');
  z1 := CSet(1.0, 2.0);
  z2 := CSet(3.0, 4.0);
  checkZ('CAdd(1+2i, 3+4i)',      CAdd(z1, z2),   4.0,  6.0);
  checkZ('CSub(3+4i, 1+2i)',      CSub(z2, z1),   2.0,  2.0);
  checkZ('CMul(1+2i, 3+4i)',      CMul(z1, z2),  -5.0, 10.0);
  checkZ('CDiv(1+2i, 3+4i)',      CDiv(z1, z2),   0.44, 0.08);
  checkZ('CNeg(1+2i)',            CNeg(z1),       -1.0, -2.0);
  checkZ('CScale(1+2i, 3)',       CScale(z1, 3.0), 3.0,  6.0);
  WriteLn;

  { 16.3 Elementary functions }
  WriteLn('-- 16.3  Elementary functions --');
  checkZ('CExp(i*Pi)  [Euler]',   CExp(CSet(0.0, Pi)),       -1.0,  0.0);
  checkZ('CLn(1+0i)',             CLn(CSet(1.0, 0.0)),        0.0,  0.0);
  checkZ('CSqrt(-1+0i)',          CSqrt(CSet(-1.0, 0.0)),     0.0,  1.0);
  checkZ('CPow(i, 2)',            CPow(CSet(0.0, 1.0), 2.0), -1.0,  0.0);
  checkZ('CPowC(i, 2+0i)',        CPowC(CSet(0.0, 1.0), CSet(2.0, 0.0)), -1.0, 0.0);
  WriteLn;

  { 16.4 Trigonometric }
  WriteLn('-- 16.4  Trigonometric --');
  checkZ('CSin(Pi/2 + 0i)',       CSin(CSet(Pi / 2.0, 0.0)),  1.0,       0.0);
  checkZ('CCos(0+1i)',            CCos(CSet(0.0, 1.0)),        Cosh(1.0), 0.0);
  checkZ('CTan(Pi/4 + 0i)',       CTan(CSet(Pi / 4.0, 0.0)),  1.0,       0.0);
  checkZ('CSinh(0+0i)',           CSinh(CSet(0.0, 0.0)),       0.0,       0.0);
  checkZ('CCosh(0+0i)',           CCosh(CSet(0.0, 0.0)),       1.0,       0.0);
  checkZ('CTanh(0+0i)',           CTanh(CSet(0.0, 0.0)),       0.0,       0.0);
  WriteLn;

  { 16.5 Utility }
  WriteLn('-- 16.5  Utility --');
  check('CAbs2(3+4i)',            CAbs2(CSet(3.0, 4.0)), 25.0);
  z := CSet(1.5, -2.3);
  WriteLn(Format('  CToStr(1.5-2.3i, 4)  = "%s"', [CToStr(z, 4)]));
  z := CSet(-0.5, 3.7);
  WriteLn(Format('  CToStr(-0.5+3.7i, 3) = "%s"', [CToStr(z, 3)]));
  WriteLn;

  WriteLn('=== done ===')
end;

end.
