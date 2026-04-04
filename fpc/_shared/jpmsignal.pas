unit jpmSignal;
{$mode objfpc}{$H+}
{-------------------------------------------------------------------------------
  Unit jpmSignal — Signal Processing
  DFT, FFT, Power Spectrum, Moving Average, Savitzky-Golay, Fourier Coefficients
  Based on Jean-Pierre Moreau's original Pascal signal-processing programs.
-------------------------------------------------------------------------------}

interface

uses SysUtils, Math, jpmtypes;

procedure DFT(var xr, xi: TFloatArray; n: integer; inverse: boolean);
procedure FFT(var xr, xi: TFloatArray; n: integer; inverse: boolean);
procedure PowerSpectrum(var xr, xi: TFloatArray; n: integer; var power: TFloatArray);
procedure MovingAverage(var data, smoothed: TFloatArray; n, window: integer);
procedure SavitzkyGolay5(var data, smoothed: TFloatArray; n: integer);
procedure FourierCoeff(var x, y: TFloatArray; ndata: integer; harmonic: integer;
  var an, bn: Float);

procedure self_test;

implementation

{ ---------------------------------------------------------------------------
  11.1  Discrete Fourier Transform (direct O(n^2) definition)
  Forward: X[k] = sum(j=0..n-1) x[j] * exp(-2*pi*i*j*k/n)
  Inverse: x[j] = (1/n) * sum(k=0..n-1) X[k] * exp(+2*pi*i*j*k/n)
  --------------------------------------------------------------------------- }
procedure DFT(var xr, xi: TFloatArray; n: integer; inverse: boolean);
var
  i, k: integer;
  sumr, sumi, theta, dir: Float;
  tr, ti: TFloatArray;
begin
  SetLength(tr, n);
  SetLength(ti, n);
  if inverse then
    dir := 1.0
  else
    dir := -1.0;
  for k := 0 to n - 1 do
  begin
    sumr := 0.0;
    sumi := 0.0;
    for i := 0 to n - 1 do
    begin
      theta := dir * 2.0 * Pi * i * k / n;
      sumr  := sumr + xr[i] * Cos(theta) - xi[i] * Sin(theta);
      sumi  := sumi + xr[i] * Sin(theta) + xi[i] * Cos(theta)
    end;
    tr[k] := sumr;
    ti[k] := sumi
  end;
  if inverse then
    for k := 0 to n - 1 do
    begin
      xr[k] := tr[k] / n;
      xi[k] := ti[k] / n
    end
  else
    for k := 0 to n - 1 do
    begin
      xr[k] := tr[k];
      xi[k] := ti[k]
    end
end;

{ ---------------------------------------------------------------------------
  Helper: test whether n is a power of 2
  --------------------------------------------------------------------------- }
function IsPowerOf2(n: integer): boolean;
begin
  result := (n > 0) and ((n and (n - 1)) = 0)
end;

{ ---------------------------------------------------------------------------
  11.2  Fast Fourier Transform — Cooley-Tukey radix-2, in-place
  Falls back to DFT when n is not a power of 2.
  Forward: exp(-2*pi*i*j*k/n); Inverse: divide result by n.
  --------------------------------------------------------------------------- }
procedure FFT(var xr, xi: TFloatArray; n: integer; inverse: boolean);
var
  i, j, k, m, halfm, bits: integer;
  tr, ti, ur, ui, wr, wi, tmp, phi: Float;
begin
  if not IsPowerOf2(n) then
  begin
    DFT(xr, xi, n, inverse);
    exit
  end;
  { bit-reversal permutation }
  j := 0;
  for i := 1 to n - 1 do
  begin
    bits := n shr 1;
    while (j and bits) <> 0 do
    begin
      j := j xor bits;
      bits := bits shr 1
    end;
    j := j xor bits;
    if i < j then
    begin
      tr := xr[i]; xr[i] := xr[j]; xr[j] := tr;
      ti := xi[i]; xi[i] := xi[j]; xi[j] := ti
    end
  end;
  { butterfly passes }
  m := 1;
  while m < n do
  begin
    halfm := m;
    m := m shl 1;
    if inverse then
      phi := Pi / halfm
    else
      phi := -Pi / halfm;
    wr := Cos(phi);
    wi := Sin(phi);
    ur := 1.0;
    ui := 0.0;
    for j := 0 to halfm - 1 do
    begin
      i := j;
      while i < n do
      begin
        k  := i + halfm;
        tr := xr[k] * ur - xi[k] * ui;
        ti := xr[k] * ui + xi[k] * ur;
        xr[k] := xr[i] - tr;
        xi[k] := xi[i] - ti;
        xr[i] := xr[i] + tr;
        xi[i] := xi[i] + ti;
        i := i + m
      end;
      tmp := ur * wr - ui * wi;
      ui  := ur * wi + ui * wr;
      ur  := tmp
    end
  end;
  if inverse then
    for i := 0 to n - 1 do
    begin
      xr[i] := xr[i] / n;
      xi[i] := xi[i] / n
    end
end;

{ ---------------------------------------------------------------------------
  11.3  Power Spectrum: power[k] = xr[k]^2 + xi[k]^2  (after FFT)
  --------------------------------------------------------------------------- }
procedure PowerSpectrum(var xr, xi: TFloatArray; n: integer; var power: TFloatArray);
var
  k: integer;
begin
  FFT(xr, xi, n, false);
  SetLength(power, n);
  for k := 0 to n - 1 do
    power[k] := xr[k] * xr[k] + xi[k] * xi[k]
end;

{ ---------------------------------------------------------------------------
  11.4  Moving Average — symmetric window of half-width `window`
  At edges, uses only available points.
  window=1 gives a 3-point average (k-1, k, k+1).
  --------------------------------------------------------------------------- }
procedure MovingAverage(var data, smoothed: TFloatArray; n, window: integer);
var
  i, j, lo, hi, cnt: integer;
  s: Float;
begin
  SetLength(smoothed, n);
  for i := 0 to n - 1 do
  begin
    lo := i - window;
    if lo < 0 then lo := 0;
    hi := i + window;
    if hi > n - 1 then hi := n - 1;
    s   := 0.0;
    cnt := 0;
    for j := lo to hi do
    begin
      s   := s + data[j];
      Inc(cnt)
    end;
    smoothed[i] := s / cnt
  end
end;

{ ---------------------------------------------------------------------------
  11.5  Savitzky-Golay smoothing — fixed 5-point, degree-2
  Coefficients: [-3, 12, 17, 12, -3] / 35
  Edge points (i < 2 or i > n-3) are copied unchanged.
  --------------------------------------------------------------------------- }
procedure SavitzkyGolay5(var data, smoothed: TFloatArray; n: integer);
var
  i: integer;
begin
  SetLength(smoothed, n);
  for i := 0 to n - 1 do
  begin
    if (i < 2) or (i > n - 3) then
      smoothed[i] := data[i]
    else
      smoothed[i] := (-3.0*data[i-2] + 12.0*data[i-1] + 17.0*data[i]
                      + 12.0*data[i+1] - 3.0*data[i+2]) / 35.0
  end
end;

{ ---------------------------------------------------------------------------
  11.6  Fourier Series Coefficients from discrete data
  Uses trapezoidal integration.
  T = x[ndata-1] - x[0],  omega = 2*pi/T
  an = (2/T) * trapezoid(y * cos(harmonic*omega*t))
  bn = (2/T) * trapezoid(y * sin(harmonic*omega*t))
  For harmonic=0: a0 = (1/T) * trapezoid(y), b0 = 0.
  --------------------------------------------------------------------------- }
procedure FourierCoeff(var x, y: TFloatArray; ndata: integer; harmonic: integer;
  var an, bn: Float);
var
  i: integer;
  T, omega, h, fi, fi1: Float;
begin
  T     := x[ndata - 1] - x[0];
  omega := 2.0 * Pi / T;
  an    := 0.0;
  bn    := 0.0;
  if harmonic = 0 then
  begin
    for i := 0 to ndata - 2 do
    begin
      h  := x[i+1] - x[i];
      an := an + 0.5 * (y[i] + y[i+1]) * h
    end;
    an := an / T
  end
  else
  begin
    for i := 0 to ndata - 2 do
    begin
      h   := x[i+1] - x[i];
      fi  := y[i]   * Cos(harmonic * omega * x[i]);
      fi1 := y[i+1] * Cos(harmonic * omega * x[i+1]);
      an  := an + 0.5 * (fi + fi1) * h
    end;
    an := an * 2.0 / T;
    for i := 0 to ndata - 2 do
    begin
      h   := x[i+1] - x[i];
      fi  := y[i]   * Sin(harmonic * omega * x[i]);
      fi1 := y[i+1] * Sin(harmonic * omega * x[i+1]);
      bn  := bn + 0.5 * (fi + fi1) * h
    end;
    bn := bn * 2.0 / T
  end
end;

{ ---------------------------------------------------------------------------
  self_test — no ReadLn, hard-coded inputs, results via WriteLn
  --------------------------------------------------------------------------- }
procedure self_test;

  procedure check(const name: string; computed, expected, tol: Float);
  var
    tag: string;
  begin
    if Abs(computed - expected) < tol then
      tag := 'OK'
    else
      tag := '*** FAIL ***';
    WriteLn(Format('  %-44s computed=%12.8f  expected=%12.8f  %s',
                   [name, computed, expected, tag]));
    if tag = '*** FAIL ***' then
      SelfTestFail(name + ': computed=' + FloatToStr(computed) + ' expected=' + FloatToStr(expected));
  end;

var
  xr, xi, xr2, xi2: TFloatArray;
  data, smoothed: TFloatArray;
  fx, fy: TFloatArray;
  an, bn, maxdiff, d: Float;
  i: integer;
begin
  WriteLn;

  { ---- Test 1: DFT of pure cosine [1, 0, -1, 0] ---- }
  WriteLn('--- Test 1: DFT of cosine [1,0,-1,0] ---');
  SetLength(xr, 4); SetLength(xi, 4);
  xr[0] :=  1.0; xr[1] := 0.0; xr[2] := -1.0; xr[3] := 0.0;
  xi[0] :=  0.0; xi[1] := 0.0; xi[2] :=  0.0; xi[3] := 0.0;
  DFT(xr, xi, 4, false);
  check('DFT xr[1] (dominant real ~ 2)',  xr[1],  2.0, 1e-10);
  check('DFT xi[1] (dominant imag ~ 0)',  xi[1],  0.0, 1e-10);
  DFT(xr, xi, 4, true);
  check('IDFT xr[0] recovers  1',  xr[0],  1.0, 1e-10);
  check('IDFT xr[1] recovers  0',  xr[1],  0.0, 1e-10);
  check('IDFT xr[2] recovers -1',  xr[2], -1.0, 1e-10);
  check('IDFT xr[3] recovers  0',  xr[3],  0.0, 1e-10);
  WriteLn;

  { ---- Test 2: FFT matches DFT for 8-point signal ---- }
  WriteLn('--- Test 2: FFT matches DFT (8 points) ---');
  SetLength(xr, 8); SetLength(xi, 8);
  SetLength(xr2, 8); SetLength(xi2, 8);
  for i := 0 to 3 do begin xr[i] := 1.0; xi[i] := 0.0 end;
  for i := 4 to 7 do begin xr[i] := 0.0; xi[i] := 0.0 end;
  for i := 0 to 7 do begin xr2[i] := xr[i]; xi2[i] := xi[i] end;
  DFT(xr,  xi,  8, false);
  FFT(xr2, xi2, 8, false);
  maxdiff := 0.0;
  for i := 0 to 7 do
  begin
    d := Abs(xr[i] - xr2[i]);
    if d > maxdiff then maxdiff := d;
    d := Abs(xi[i] - xi2[i]);
    if d > maxdiff then maxdiff := d
  end;
  check('FFT vs DFT max |diff| (< 1e-10)', maxdiff, 0.0, 1e-10);
  WriteLn;

  { ---- Test 3: Moving Average ---- }
  WriteLn('--- Test 3: Moving Average window=1 ---');
  SetLength(data, 8);
  for i := 0 to 7 do data[i] := i + 1.0;   { [1,2,3,4,5,6,7,8] }
  MovingAverage(data, smoothed, 8, 1);
  check('MovingAvg smoothed[3] = (3+4+5)/3 = 4', smoothed[3], 4.0, 1e-10);
  WriteLn;

  { ---- Test 4: Savitzky-Golay ---- }
  WriteLn('--- Test 4: Savitzky-Golay 5-point ---');
  SetLength(data, 8);
  for i := 0 to 7 do data[i] := i + 1.0;   { [1,2,3,4,5,6,7,8] linear }
  SavitzkyGolay5(data, smoothed, 8);
  check('SavitzkyGolay5 smoothed[4] = 5.0', smoothed[4], 5.0, 1e-10);
  WriteLn;

  { ---- Test 5: Fourier Coefficients of sin(x) over [0, 2*pi] ---- }
  WriteLn('--- Test 5: Fourier Coefficients of sin(x) ---');
  SetLength(fx, 5); SetLength(fy, 5);
  fx[0] := 0.0;            fy[0] := Sin(fx[0]);
  fx[1] := Pi / 2.0;       fy[1] := Sin(fx[1]);
  fx[2] := Pi;              fy[2] := Sin(fx[2]);
  fx[3] := 3.0 * Pi / 2.0; fy[3] := Sin(fx[3]);
  fx[4] := 2.0 * Pi;       fy[4] := Sin(fx[4]);
  FourierCoeff(fx, fy, 5, 1, an, bn);
  check('FourierCoeff a1 ~ 0 (cosine part)', an, 0.0, 1e-6);
  check('FourierCoeff b1 ~ 1 (sine  part)',  bn, 1.0, 1e-6);
  WriteLn;

  WriteLn('=== Done ===')
end;

end.
