unit jpmDiffeq;
{$mode objfpc}{$H+}

interface

uses SysUtils, Math, jpmtypes;

type
  TODEFunc    = procedure(t: Float; var y, dydt: TFloatArray; n: integer);

{ 7.1 Classical Runge-Kutta 4 (fixed step) }
procedure RK4Step(f: TODEFunc; t: Float; var y: TFloatArray; n: integer; h: Float;
  var yout: TFloatArray);
procedure RK4Integrate(f: TODEFunc; var y: TFloatArray; n: integer; t1, t2: Float;
  steps: integer);

{ 7.2 Runge-Kutta-Fehlberg 45 (adaptive step, Cash-Karp) }
procedure RKF45(f: TODEFunc; var y: TFloatArray; n: integer; var t: Float; tend: Float;
  hstart, hmin, hmax, tol: Float; var nOK, nBad: integer);

{ 7.3 Adams-Bashforth 4-step (explicit, fixed step) }
procedure AdamsBashforth4(f: TODEFunc; var y: TFloatArray; n: integer; var t: Float;
  tend, h: Float);

{ 7.3b Adams-Moulton 4-step (predictor-corrector, fixed step) }
procedure AdamsMoulton4(f: TODEFunc; var y: TFloatArray; n: integer; var t: Float;
  tend, h: Float);

procedure self_test;

implementation

{ ------------------------------------------------------------------ }
{  Internal helpers                                                   }
{ ------------------------------------------------------------------ }

procedure SetLen(var a: TFloatArray; n: integer);
begin
  SetLength(a, n);
end;

procedure VecCopy(var dst: TFloatArray; const src: TFloatArray; n: integer);
var
  i: integer;
begin
  for i := 0 to n - 1 do
    dst[i] := src[i];
end;

procedure VecAddScale(var dst: TFloatArray; const a: TFloatArray; scale: Float; n: integer);
var
  i: integer;
begin
  for i := 0 to n - 1 do
    dst[i] := dst[i] + scale * a[i];
end;

{ ------------------------------------------------------------------ }
{  7.1  RK4 single step                                              }
{ ------------------------------------------------------------------ }
procedure RK4Step(f: TODEFunc; t: Float; var y: TFloatArray; n: integer; h: Float;
  var yout: TFloatArray);
var
  k1, k2, k3, k4, ytmp: TFloatArray;
  i: integer;
  h2: Float;
begin
  SetLen(k1, n); SetLen(k2, n); SetLen(k3, n); SetLen(k4, n);
  SetLen(ytmp, n);
  h2 := h / 2.0;

  f(t, y, k1, n);

  for i := 0 to n - 1 do
    ytmp[i] := y[i] + h2 * k1[i];
  f(t + h2, ytmp, k2, n);

  for i := 0 to n - 1 do
    ytmp[i] := y[i] + h2 * k2[i];
  f(t + h2, ytmp, k3, n);

  for i := 0 to n - 1 do
    ytmp[i] := y[i] + h * k3[i];
  f(t + h, ytmp, k4, n);

  SetLen(yout, n);
  for i := 0 to n - 1 do
    yout[i] := y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
end;

{ ------------------------------------------------------------------ }
{  7.1  RK4 full integration                                         }
{ ------------------------------------------------------------------ }
procedure RK4Integrate(f: TODEFunc; var y: TFloatArray; n: integer; t1, t2: Float;
  steps: integer);
var
  t, h: Float;
  yout: TFloatArray;
  k: integer;
begin
  h := (t2 - t1) / steps;
  t := t1;
  SetLen(yout, n);
  for k := 1 to steps do
  begin
    RK4Step(f, t, y, n, h, yout);
    VecCopy(y, yout, n);
    t := t + h;
  end;
end;

{ ------------------------------------------------------------------ }
{  7.2  RKF45 adaptive step — Cash-Karp coefficients                 }
{ ------------------------------------------------------------------ }
procedure RKF45(f: TODEFunc; var y: TFloatArray; n: integer; var t: Float; tend: Float;
  hstart, hmin, hmax, tol: Float; var nOK, nBad: integer);
const
  { Cash-Karp Butcher tableau }
  a2  = 1.0/5.0;
  a3  = 3.0/10.0;
  a4  = 3.0/5.0;
  a5  = 1.0;
  a6  = 7.0/8.0;
  b21 = 1.0/5.0;
  b31 = 3.0/40.0;       b32 = 9.0/40.0;
  b41 = 3.0/10.0;       b42 = -9.0/10.0;    b43 = 6.0/5.0;
  b51 = -11.0/54.0;     b52 = 5.0/2.0;      b53 = -70.0/27.0;    b54 = 35.0/27.0;
  b61 = 1631.0/55296.0; b62 = 175.0/512.0;  b63 = 575.0/13824.0;
  b64 = 44275.0/110592.0; b65 = 253.0/4096.0;
  { 5th-order weights }
  c1  = 37.0/378.0;     c3  = 250.0/621.0;  c4  = 125.0/594.0;
  c6  = 512.0/1771.0;
  { error: 5th minus 4th }
  e1  = 37.0/378.0   - 2825.0/27648.0;
  e3  = 250.0/621.0  - 18575.0/48384.0;
  e4  = 125.0/594.0  - 13525.0/55296.0;
  e5  = 0.0          - 277.0/14336.0;
  e6  = 512.0/1771.0 - 1.0/4.0;
  SAFETY  = 0.9;
  ERRCON  = 1.89e-4;  { (5/SAFETY)^(1/0.2) }
  MAXSTEP = 10000;
var
  k1, k2, k3, k4, k5, k6, ytmp, yerr: TFloatArray;
  i, step: integer;
  h, hh, errmax, yscale, err, tnew: Float;
  accepted: boolean;
begin
  SetLen(k1, n); SetLen(k2, n); SetLen(k3, n);
  SetLen(k4, n); SetLen(k5, n); SetLen(k6, n);
  SetLen(ytmp, n); SetLen(yerr, n);
  h := hstart;
  nOK := 0; nBad := 0;
  step := 0;
  while (t < tend) and (step < MAXSTEP) do
  begin
    inc(step);
    if t + h > tend then
      h := tend - t;

    accepted := false;
    while not accepted do
    begin
      { k1 }
      f(t, y, k1, n);
      { k2 }
      for i := 0 to n - 1 do ytmp[i] := y[i] + h * b21 * k1[i];
      f(t + a2 * h, ytmp, k2, n);
      { k3 }
      for i := 0 to n - 1 do ytmp[i] := y[i] + h * (b31 * k1[i] + b32 * k2[i]);
      f(t + a3 * h, ytmp, k3, n);
      { k4 }
      for i := 0 to n - 1 do ytmp[i] := y[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
      f(t + a4 * h, ytmp, k4, n);
      { k5 }
      for i := 0 to n - 1 do ytmp[i] := y[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
      f(t + a5 * h, ytmp, k5, n);
      { k6 }
      for i := 0 to n - 1 do ytmp[i] := y[i] + h * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
      f(t + a6 * h, ytmp, k6, n);

      { 5th-order solution }
      for i := 0 to n - 1 do
        ytmp[i] := y[i] + h * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c6 * k6[i]);

      { error estimate }
      for i := 0 to n - 1 do
        yerr[i] := h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i]);

      { error norm }
      errmax := 0.0;
      for i := 0 to n - 1 do
      begin
        yscale := Abs(y[i]) + Abs(h * k1[i]) + 1.0e-30;
        err := Abs(yerr[i]) / (tol * yscale);
        if err > errmax then errmax := err;
      end;

      if errmax <= 1.0 then
      begin
        accepted := true;
        tnew := t + h;
        VecCopy(y, ytmp, n);
        t := tnew;
        inc(nOK);
        if errmax > ERRCON then
          h := SAFETY * h * Power(errmax, -0.2)
        else
          h := 5.0 * h;
        if h > hmax then h := hmax;
      end else
      begin
        hh := SAFETY * h * Power(errmax, -0.25);
        inc(nBad);
        if Abs(hh) < Abs(hmin) then
          h := hmin
        else
          h := hh;
        accepted := (Abs(h) <= Abs(hmin));
        if accepted then
        begin
          VecCopy(y, ytmp, n);
          t := t + h;
        end;
      end;
    end;
  end;
end;

{ ------------------------------------------------------------------ }
{  Internal: evaluate f at each of 4 past times into a history array }
{ ------------------------------------------------------------------ }

{ history[0..3, 0..n-1] — row 0 is most recent }

type
  THistory = array[0..3] of TFloatArray;

procedure ShiftHistory(var hist: THistory; var xhist: array of Float;
  xnew: Float; const ynew: TFloatArray; n: integer);
var
  i: integer;
begin
  for i := 3 downto 1 do
  begin
    VecCopy(hist[i], hist[i - 1], n);
    xhist[i] := xhist[i - 1];
  end;
  VecCopy(hist[0], ynew, n);
  xhist[0] := xnew;
end;

{ ------------------------------------------------------------------ }
{  7.3  Adams-Bashforth 4-step                                       }
{ ------------------------------------------------------------------ }
procedure AdamsBashforth4(f: TODEFunc; var y: TFloatArray; n: integer; var t: Float;
  tend, h: Float);
var
  hist: THistory;
  xhist: array[0..3] of Float;
  b: array[0..3] of TFloatArray;
  ynew, ytmp, btmp: TFloatArray;
  i, j: integer;
begin
  { allocate }
  for i := 0 to 3 do
  begin
    SetLen(hist[i], n);
    SetLen(b[i], n);
  end;
  SetLen(ynew, n);
  SetLen(ytmp, n);

  { startup: use RK4 to fill history[0..2] (we need f at t, t+h, t+2h) }
  { store f values in b[0..2] }
  f(t, y, b[0], n);
  VecCopy(hist[0], y, n);
  xhist[0] := t;

  VecCopy(ytmp, y, n);
  for j := 1 to 2 do
  begin
    RK4Step(f, t + (j - 1) * h, ytmp, n, h, ynew);
    VecCopy(ytmp, ynew, n);
    VecCopy(hist[j], ynew, n);
    xhist[j] := t + j * h;
    f(xhist[j], hist[j], b[j], n);
  end;

  { main loop: Adams-Bashforth 4-step predictor
    yn+1 = yn + h/24*(55*fn - 59*fn-1 + 37*fn-2 - 9*fn-3)
    We start at index 2 (have f at 0,1,2), first full step gives index 3 }

  { For first AB4 step we need f at index 3 but only have 3 points,
    so do one more RK4 step to get the 4th startup point }
  RK4Step(f, xhist[2], hist[2], n, h, ynew);
  VecCopy(hist[3], ynew, n);
  xhist[3] := xhist[2] + h;
  f(xhist[3], hist[3], b[3], n);

  { shift history so hist[0] = most recent }
  { rearrange: hist[0]=newest(3), hist[1]=2, hist[2]=1, hist[3]=0 }
  { b[0..3] correspond to f values at xhist[3..0] }
  { now copy into ordered arrays: index 0 = most recent }
  VecCopy(y, hist[3], n);
  t := xhist[3];

  { build working history: b[0]=f_n, b[1]=f_n-1 etc.
    We now have b[3]=f(t0), b[2]=f(t1), b[1]=f(t2), b[0]=f(t3),
    so swap to make b[0] = most recent }
  SetLen(btmp, n);
  VecCopy(btmp, b[0], n); VecCopy(b[0], b[3], n); VecCopy(b[3], btmp, n);
  VecCopy(btmp, b[1], n); VecCopy(b[1], b[2], n); VecCopy(b[2], btmp, n);

  { main integration loop }
  while t + h / 2.0 < tend do
  begin
    if t + h > tend + h * 1.0e-10 then
      h := tend - t;
    for i := 0 to n - 1 do
      ynew[i] := y[i] + h / 24.0 * (55.0 * b[0][i] - 59.0 * b[1][i]
                 + 37.0 * b[2][i] - 9.0 * b[3][i]);
    { shift b }
    VecCopy(b[3], b[2], n);
    VecCopy(b[2], b[1], n);
    VecCopy(b[1], b[0], n);
    VecCopy(y, ynew, n);
    t := t + h;
    f(t, y, b[0], n);
  end;
end;

{ ------------------------------------------------------------------ }
{  7.3b Adams-Moulton 4-step predictor-corrector                     }
{ ------------------------------------------------------------------ }
procedure AdamsMoulton4(f: TODEFunc; var y: TFloatArray; n: integer; var t: Float;
  tend, h: Float);
var
  yp, yc, ytmp, ynew: TFloatArray;
  { b[0]=f_n, b[1]=f_n-1, b[2]=f_n-2 (needed for both predictor+corrector) }
  b: array[0..2] of TFloatArray;
  bpred: TFloatArray;
  i: integer;
begin
  SetLen(yp, n); SetLen(yc, n); SetLen(ytmp, n); SetLen(ynew, n);
  SetLen(bpred, n);
  for i := 0 to 2 do SetLen(b[i], n);

  { startup: RK4 to get y at t+h and t+2h }
  f(t, y, b[2], n);   { b[2] = f(t) }

  RK4Step(f, t, y, n, h, ytmp);
  f(t + h, ytmp, b[1], n);   { b[1] = f(t+h) }

  VecCopy(ynew, ytmp, n);
  RK4Step(f, t + h, ytmp, n, h, ynew);
  f(t + 2.0 * h, ynew, b[0], n);  { b[0] = f(t+2h) }
  VecCopy(y, ynew, n);
  t := t + 2.0 * h;

  { main loop }
  while t + h / 2.0 < tend do
  begin
    if t + h > tend + h * 1.0e-10 then
      h := tend - t;

    { Adams-Bashforth predictor (3-step):
      yP_n+1 = yn + h/12*(23*fn - 16*fn-1 + 5*fn-2) }
    for i := 0 to n - 1 do
      yp[i] := y[i] + h / 12.0 * (23.0 * b[0][i] - 16.0 * b[1][i] + 5.0 * b[2][i]);

    { Adams-Moulton corrector (1 correction):
      yC_n+1 = yn + h/12*(5*f(t+h,yP) + 8*fn - fn-1) }
    f(t + h, yp, bpred, n);
    for i := 0 to n - 1 do
      yc[i] := y[i] + h / 12.0 * (5.0 * bpred[i] + 8.0 * b[0][i] - b[1][i]);

    { shift history }
    VecCopy(b[2], b[1], n);
    VecCopy(b[1], b[0], n);
    VecCopy(y, yc, n);
    t := t + h;
    f(t, y, b[0], n);
  end;
end;

{ ------------------------------------------------------------------ }
{  Self-test ODE right-hand sides                                     }
{ ------------------------------------------------------------------ }

procedure ODE_ExpDecay(t: Float; var y, dydt: TFloatArray; n: integer);
begin
  dydt[0] := -y[0];
end;

procedure ODE_Harmonic(t: Float; var y, dydt: TFloatArray; n: integer);
begin
  dydt[0] :=  y[1];
  dydt[1] := -y[0];
end;

{ ------------------------------------------------------------------ }
{  self_test                                                          }
{ ------------------------------------------------------------------ }
procedure self_test;
var
  y: TFloatArray;
  t: Float;
  nOK, nBad: integer;
  exact, err: Float;
begin
  WriteLn('=== jpmdiffeq self_test ===');
  WriteLn;

  { ---- Test 1a: Exp decay with RK4 ---- }
  SetLength(y, 1);
  y[0] := 1.0;
  RK4Integrate(@ODE_ExpDecay, y, 1, 0.0, 1.0, 100);
  exact := Exp(-1.0);
  err   := Abs(y[0] - exact);
  WriteLn(Format('Test 1a (RK4   exp decay):  y(1) = %.9f  exact = %.9f  err = %.2e',
    [y[0], exact, err]));
  if err < 1.0e-8 then WriteLn('  PASS') else begin WriteLn('  FAIL'); SelfTestFail('RK4 exp decay: err=' + FloatToStr(err)); end;

  { ---- Test 1b: Exp decay with RKF45 ---- }
  y[0] := 1.0;
  t := 0.0;
  nOK := 0; nBad := 0;
  RKF45(@ODE_ExpDecay, y, 1, t, 1.0, 0.1, 1.0e-12, 1.0, 1.0e-9, nOK, nBad);
  err := Abs(y[0] - exact);
  WriteLn(Format('Test 1b (RKF45 exp decay):  y(1) = %.9f  exact = %.9f  err = %.2e  nOK=%d nBad=%d',
    [y[0], exact, err, nOK, nBad]));
  if err < 1.0e-7 then WriteLn('  PASS')
  else begin WriteLn('  FAIL'); SelfTestFail('RKF45 exp decay: err=' + FloatToStr(err)); end;

  { ---- Test 2a: Harmonic oscillator with RK4 ---- }
  SetLength(y, 2);
  y[0] := 1.0; y[1] := 0.0;
  RK4Integrate(@ODE_Harmonic, y, 2, 0.0, Pi, 1000);
  exact := Cos(Pi);  { = -1 }
  err   := Abs(y[0] - exact);
  WriteLn(Format('Test 2a (RK4   harmonic):   y1(pi)= %.9f  exact = %.9f  err = %.2e',
    [y[0], exact, err]));
  if err < 1.0e-8 then WriteLn('  PASS')
  else begin WriteLn('  FAIL'); SelfTestFail('RK4 harmonic: err=' + FloatToStr(err)); end;

  { ---- Test 2b: Harmonic oscillator with RKF45 ---- }
  y[0] := 1.0; y[1] := 0.0;
  t := 0.0;
  nOK := 0; nBad := 0;
  RKF45(@ODE_Harmonic, y, 2, t, Pi, 0.1, 1.0e-12, 1.0, 1.0e-9, nOK, nBad);
  err := Abs(y[0] - exact);
  WriteLn(Format('Test 2b (RKF45 harmonic):   y1(pi)= %.9f  exact = %.9f  err = %.2e  nOK=%d nBad=%d',
    [y[0], exact, err, nOK, nBad]));
  if err < 1.0e-7 then WriteLn('  PASS')
  else begin WriteLn('  FAIL'); SelfTestFail('RKF45 harmonic: err=' + FloatToStr(err)); end;

  { ---- Test 3: Exp decay with Adams-Bashforth 4 ---- }
  SetLength(y, 1);
  y[0] := 1.0;
  t := 0.0;
  AdamsBashforth4(@ODE_ExpDecay, y, 1, t, 1.0, 0.01);
  exact := Exp(-1.0);
  err   := Abs(y[0] - exact);
  WriteLn(Format('Test 3  (AB4   exp decay):  y(1) = %.9f  exact = %.9f  err = %.2e',
    [y[0], exact, err]));
  if err < 1.0e-6 then WriteLn('  PASS')
  else begin WriteLn('  FAIL'); SelfTestFail('AB4 exp decay: err=' + FloatToStr(err)); end;

  { ---- Test 4: Exp decay with Adams-Moulton ---- }
  y[0] := 1.0;
  t := 0.0;
  AdamsMoulton4(@ODE_ExpDecay, y, 1, t, 1.0, 0.01);
  exact := Exp(-1.0);
  err   := Abs(y[0] - exact);
  WriteLn(Format('Test 4  (AM4   exp decay):  y(1) = %.9f  exact = %.9f  err = %.2e',
    [y[0], exact, err]));
  if err < 1.0e-7 then WriteLn('  PASS')
  else begin WriteLn('  FAIL'); SelfTestFail('AM4 exp decay: err=' + FloatToStr(err)); end;

  WriteLn;
  WriteLn('=== self_test done ===');
end;

end.
