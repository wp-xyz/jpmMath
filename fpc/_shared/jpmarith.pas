unit jpmArith;
{$mode objfpc}{$H+}

interface

uses SysUtils, Math, jpmtypes;

type
  TFraction = record
    num: int64;  { numerator }
    den: int64;  { denominator, always > 0 }
  end;

  TInt64Array = array of int64;

{ 12.1 GCD and LCM }
function GCD(a, b: int64): int64;
function LCM(a, b: int64): int64;

{ 12.2 Extended GCD (Bezout coefficients) }
function ExtGCD(a, b: int64; var x, y: int64): int64;

{ 12.3 Combinatorics }
function Factorial(n: integer): Float;
function BinomCoeff(n, k: integer): Float;
function Permutation(n, k: integer): Float;

{ 12.4 Fractions }
function FracCreate(n, d: int64): TFraction;
function FracAdd(a, b: TFraction): TFraction;
function FracSub(a, b: TFraction): TFraction;
function FracMul(a, b: TFraction): TFraction;
function FracDiv(a, b: TFraction): TFraction;
function FracToFloat(a: TFraction): Float;
function FracToStr(a: TFraction): string;

{ 12.5 Base Conversion }
function IntToBase(n: int64; base: integer): string;
function BaseToInt(s: string; base: integer): int64;

{ 12.6 Prime Testing and Sieve }
function IsPrime(n: int64): boolean;
procedure PrimeSieve(limit: integer; var primes: TIntArray; var count: integer);

{ 12.7 Factorization }
procedure Factorize(n: int64; var factors: TInt64Array; var count: integer);

{ 12.8 Diophantine equation solver }
function SolveDiophantine(a, b, c: int64; var x, y: int64): boolean;

{ Self-test }
procedure self_test;

implementation

{ ------------------------------------------------------------------ }
{ 12.1  GCD — Euclidean algorithm                                     }
{ ------------------------------------------------------------------ }
function GCD(a, b: int64): int64;
var
  r: int64;
begin
  a := Abs(a);
  b := Abs(b);
  while b <> 0 do
  begin
    r := a mod b;
    a := b;
    b := r;
  end;
  result := a;
end;

{ 12.1  LCM = |a*b| / GCD(a,b) }
function LCM(a, b: int64): int64;
var
  g: int64;
begin
  g := GCD(a, b);
  if g = 0 then
    result := 0
  else
    result := (Abs(a) div g) * Abs(b);
end;

{ ------------------------------------------------------------------ }
{ 12.2  Extended GCD                                                  }
{ Returns GCD; fills x,y such that a*x + b*y = GCD(a,b)             }
{ ------------------------------------------------------------------ }
function ExtGCD(a, b: int64; var x, y: int64): int64;
var
  x1, y1, q: int64;
begin
  if b = 0 then
  begin
    x := 1;
    y := 0;
    result := a;
    exit;
  end;
  result := ExtGCD(b, a mod b, x1, y1);
  q := a div b;
  x := y1;
  y := x1 - q * y1;
end;

{ ------------------------------------------------------------------ }
{ 12.3  Combinatorics                                                 }
{ ------------------------------------------------------------------ }

{ n! as Float — exact for n<=20, Stirling/LnGamma for larger n }
function Factorial(n: integer): Float;
var
  i: integer;
  f: Float;
begin
  if n < 0 then
  begin
    result := 0;
    exit;
  end;
  if n <= 20 then
  begin
    f := 1.0;
    for i := 2 to n do
      f := f * i;
    result := f;
  end
  else
    { Stirling's approximation: ln(n!) ≈ n*ln(n) - n + 0.5*ln(2*pi*n) }
    result := Exp(n * Ln(n) - n + 0.5 * Ln(2.0 * Pi * n));
end;

{ C(n,k) = n! / (k! * (n-k)!) }
function BinomCoeff(n, k: integer): Float;
var
  i, u: integer;
  r: Float;
begin
  if (k < 0) or (k > n) then
  begin
    result := 0;
    exit;
  end;
  if k > n - k then
    u := n - k
  else
    u := k;
  r := 1.0;
  for i := 0 to u - 1 do
    r := r * (n - i) / (i + 1);
  result := r;
end;

{ P(n,k) = n! / (n-k)! }
function Permutation(n, k: integer): Float;
var
  i: integer;
  r: Float;
begin
  if (k < 0) or (k > n) then
  begin
    result := 0;
    exit;
  end;
  r := 1.0;
  for i := n - k + 1 to n do
    r := r * i;
  result := r;
end;

{ ------------------------------------------------------------------ }
{ 12.4  Fractions                                                     }
{ ------------------------------------------------------------------ }

function FracCreate(n, d: int64): TFraction;
var
  g: int64;
begin
  if d = 0 then
  begin
    result.num := 0;
    result.den := 1;
    exit;
  end;
  if d < 0 then
  begin
    n := -n;
    d := -d;
  end;
  g := GCD(Abs(n), d);
  if g = 0 then g := 1;
  result.num := n div g;
  result.den := d div g;
end;

function FracAdd(a, b: TFraction): TFraction;
begin
  result := FracCreate(a.num * b.den + b.num * a.den, a.den * b.den);
end;

function FracSub(a, b: TFraction): TFraction;
begin
  result := FracCreate(a.num * b.den - b.num * a.den, a.den * b.den);
end;

function FracMul(a, b: TFraction): TFraction;
begin
  result := FracCreate(a.num * b.num, a.den * b.den);
end;

function FracDiv(a, b: TFraction): TFraction;
begin
  result := FracCreate(a.num * b.den, a.den * b.num);
end;

function FracToFloat(a: TFraction): Float;
begin
  result := a.num / a.den;
end;

function FracToStr(a: TFraction): string;
begin
  result := IntToStr(a.num) + '/' + IntToStr(a.den);
end;

{ ------------------------------------------------------------------ }
{ 12.5  Base Conversion                                               }
{ ------------------------------------------------------------------ }

function IntToBase(n: int64; base: integer): string;
const
  Digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ';
var
  r: int64;
begin
  if (base < 2) or (base > 36) then
  begin
    result := '';
    exit;
  end;
  if n = 0 then
  begin
    result := '0';
    exit;
  end;
  result := '';
  while n > 0 do
  begin
    r := n mod base;
    result := Digits[r + 1] + result;
    n := n div base;
  end;
end;

function BaseToInt(s: string; base: integer): int64;
var
  i, j: integer;
  ch: char;
begin
  result := 0;
  s := UpperCase(s);
  for i := 1 to Length(s) do
  begin
    ch := s[i];
    if (ch >= '0') and (ch <= '9') then
      j := Ord(ch) - Ord('0')
    else if (ch >= 'A') and (ch <= 'Z') then
      j := Ord(ch) - Ord('A') + 10
    else
    begin
      result := 0;
      exit;
    end;
    result := result * base + j;
  end;
end;

{ ------------------------------------------------------------------ }
{ 12.6  Prime Testing and Sieve                                       }
{ ------------------------------------------------------------------ }

{ Trial division up to sqrt(n); correct for any n < 10^15 }
function IsPrime(n: int64): boolean;
var
  i, m: int64;
begin
  if n < 2 then
  begin
    result := false;
    exit;
  end;
  if n = 2 then
  begin
    result := true;
    exit;
  end;
  if n mod 2 = 0 then
  begin
    result := false;
    exit;
  end;
  if n = 3 then
  begin
    result := true;
    exit;
  end;
  if n mod 3 = 0 then
  begin
    result := false;
    exit;
  end;
  m := Trunc(Sqrt(n)) + 1;
  i := 5;
  result := true;
  while i <= m do
  begin
    if (n mod i = 0) or (n mod (i + 2) = 0) then
    begin
      result := false;
      exit;
    end;
    i := i + 6;
  end;
end;

{ Sieve of Eratosthenes: primes[0..count-1] <= limit }
procedure PrimeSieve(limit: integer; var primes: TIntArray; var count: integer);
var
  sieve: array of boolean;
  i, j: integer;
begin
  count := 0;
  if limit < 2 then
  begin
    SetLength(primes, 0);
    exit;
  end;
  SetLength(sieve, limit + 1);
  for i := 0 to limit do
    sieve[i] := true;
  sieve[0] := false;
  sieve[1] := false;
  i := 2;
  while i * i <= limit do
  begin
    if sieve[i] then
    begin
      j := i * i;
      while j <= limit do
      begin
        sieve[j] := false;
        j := j + i;
      end;
    end;
    Inc(i);
  end;
  for i := 2 to limit do
    if sieve[i] then
      Inc(count);
  SetLength(primes, count);
  j := 0;
  for i := 2 to limit do
    if sieve[i] then
    begin
      primes[j] := i;
      Inc(j);
    end;
end;

{ ------------------------------------------------------------------ }
{ 12.7  Factorization — trial division                                }
{ ------------------------------------------------------------------ }
procedure Factorize(n: int64; var factors: TInt64Array; var count: integer);
var
  d: int64;
begin
  count := 0;
  SetLength(factors, 0);
  if n <= 1 then
    exit;
  d := 2;
  while n > 1 do
  begin
    while n mod d = 0 do
    begin
      SetLength(factors, count + 1);
      factors[count] := d;
      Inc(count);
      n := n div d;
    end;
    if d = 2 then
      d := 3
    else
      d := d + 2;
    if d * d > n then
    begin
      if n > 1 then
      begin
        SetLength(factors, count + 1);
        factors[count] := n;
        Inc(count);
      end;
      break;
    end;
  end;
end;

{ ------------------------------------------------------------------ }
{ 12.8  Diophantine equation solver                                   }
{ Solve a*x + b*y = c. Uses ExtGCD.                                  }
{ ------------------------------------------------------------------ }
function SolveDiophantine(a, b, c: int64; var x, y: int64): boolean;
var
  g, x0, y0, scale: int64;
begin
  g := ExtGCD(a, b, x0, y0);
  if (g = 0) or (c mod g <> 0) then
  begin
    result := false;
    exit;
  end;
  scale := c div g;
  x := x0 * scale;
  y := y0 * scale;
  result := true;
end;

{ ------------------------------------------------------------------ }
{ Self-test                                                           }
{ ------------------------------------------------------------------ }
procedure self_test;
var
  x, y, g: int64;
  f1, f2, fr: TFraction;
  primes: TIntArray;
  factors: TInt64Array;
  cnt, i: integer;
begin
  WriteLn('=== jpmArith Self-Test ===');
  WriteLn;

  { GCD / LCM }
  WriteLn('GCD(48, 18)        = ', GCD(48, 18), '         [expected 6]');
  WriteLn('LCM(4, 6)          = ', LCM(4, 6),  '        [expected 12]');
  WriteLn;

  { Extended GCD }
  g := ExtGCD(35, 15, x, y);
  WriteLn('ExtGCD(35,15):  GCD=', g, '  x=', x, '  y=', y,
          '  verify 35x+15y=', 35*x + 15*y, '  [expected GCD=5, sum=5]');
  WriteLn;

  { Combinatorics }
  WriteLn('Factorial(10)      = ', Factorial(10):0:0,    '  [expected 3628800]');
  WriteLn('BinomCoeff(10,3)   = ', BinomCoeff(10,3):0:0, '   [expected 120]');
  WriteLn('Permutation(10,3)  = ', Permutation(10,3):0:0,'   [expected 720]');
  WriteLn;

  { Fractions }
  f1 := FracCreate(1, 2);
  f2 := FracCreate(1, 3);
  fr := FracAdd(f1, f2);
  WriteLn('FracAdd(1/2, 1/3)  = ', FracToStr(fr), '  [expected 5/6]');
  f1 := FracCreate(3, 4);
  f2 := FracCreate(2, 3);
  fr := FracMul(f1, f2);
  WriteLn('FracMul(3/4, 2/3)  = ', FracToStr(fr), '  [expected 1/2]');
  WriteLn;

  { Base Conversion }
  WriteLn('IntToBase(255,16)  = ', IntToBase(255, 16),    '  [expected FF]');
  WriteLn('BaseToInt("FF",16) = ', BaseToInt('FF', 16),   '  [expected 255]');
  WriteLn('IntToBase(10,2)    = ', IntToBase(10, 2),      '  [expected 1010]');
  WriteLn;

  { Primality }
  WriteLn('IsPrime(17)        = ', IsPrime(17),  '  [expected TRUE]');
  WriteLn('IsPrime(18)        = ', IsPrime(18),  '  [expected FALSE]');
  WriteLn;

  { Sieve }
  PrimeSieve(20, primes, cnt);
  Write('PrimeSieve(20)     = ');
  for i := 0 to cnt - 1 do
    Write(primes[i], ' ');
  WriteLn;
  WriteLn('                     [expected 2 3 5 7 11 13 17 19]');
  WriteLn;

  { Factorize }
  Factorize(360, factors, cnt);
  Write('Factorize(360)     = ');
  for i := 0 to cnt - 1 do
    Write(factors[i], ' ');
  WriteLn;
  WriteLn('                     [expected 2 2 2 3 3 5]');
  WriteLn;

  { Diophantine }
  if SolveDiophantine(3, 5, 1, x, y) then
    WriteLn('SolveDioph(3,5,1): x=', x, '  y=', y,
            '  verify 3x+5y=', 3*x + 5*y, '  [expected TRUE, sum=1]')
  else
    WriteLn('SolveDioph(3,5,1): no solution  [WRONG]');
  WriteLn;

  WriteLn('=== Self-Test Complete ===');
end;

end.
