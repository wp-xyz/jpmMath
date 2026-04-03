unit jpmLstsqr;
{$mode objfpc}{$H+}
{-------------------------------------------------------------------------------
                               Unit jpmLstsqr
        Least-squares fitting and regression based on Jean-Pierre Moreau's code.

Procedures:
  LeastSquares       — overdetermined Ax=b via normal equations
  PolyRegression     — polynomial fit via Vandermonde + least squares
  LinearRegression   — simple y = a + bx with correlation coefficient
  MultipleRegression — multiple linear regression with R²
  ChiSquareFit       — weighted polynomial fit with chi-square statistic
-------------------------------------------------------------------------------}

interface

uses SysUtils, Math, jpmtypes;

type
  TVector = array of Float;
  TMatrix = array of TVector;

{ Solve overdetermined system Ax=b (m rows, n cols, m>=n) in least-squares sense
  using normal equations: (A^T A) x = A^T b }
procedure LeastSquares(var A: TMatrix; var b: TVector; m, n: integer;
  var x: TVector);

{ Fit polynomial of given degree to (xdata,ydata).
  coeffs[0]=constant, coeffs[1]=linear, ..., coeffs[degree]=highest power. }
procedure PolyRegression(var xdata, ydata: TVector; ndata, degree: integer;
  var coeffs: TVector);

{ Simple linear regression y = a + b*x.
  Returns intercept a, slope b, Pearson correlation coefficient r. }
procedure LinearRegression(var xdata, ydata: TVector; ndata: integer;
  var a, b, r: Float);

{ Multiple linear regression.  X is m-by-n design matrix (no intercept column).
  Intercept is added automatically.  Returns coeffs[0]=intercept,
  coeffs[1..n]=slopes, and R² in r2. }
procedure MultipleRegression(var X: TMatrix; var y: TVector; m, n: integer;
  var coeffs: TVector; var r2: Float);

{ Weighted polynomial fit where sigma[i] is the measurement uncertainty.
  Weight = 1/sigma[i]².  chiSq = sum((y[i]-poly(x[i]))²/sigma[i]²). }
procedure ChiSquareFit(var xdata, ydata, sigma: TVector; ndata, degree: integer;
  var coeffs: TVector; var chiSq: Float);

procedure self_test;

implementation

{ ---------------------------------------------------------------------------
  Internal helper — Gaussian elimination with partial pivoting.
  Solves n-by-n system mat*x = rhs.  Modifies mat and rhs in-place.
  --------------------------------------------------------------------------- }
procedure SolveLinearSystem(var mat: TMatrix; var rhs: TVector; n: integer;
  var x: TVector);
var
  i, j, k, pivrow: integer;
  factor, tmp: Float;
  tmprow: TVector;
begin
  { forward elimination }
  for k := 0 to n - 2 do
  begin
    { find pivot row }
    pivrow := k;
    for i := k + 1 to n - 1 do
      if Abs(mat[i][k]) > Abs(mat[pivrow][k]) then
        pivrow := i;
    { swap rows }
    if pivrow <> k then
    begin
      tmprow      := mat[k];
      mat[k]      := mat[pivrow];
      mat[pivrow] := tmprow;
      tmp         := rhs[k];
      rhs[k]      := rhs[pivrow];
      rhs[pivrow] := tmp
    end;
    { eliminate column k below pivot }
    for i := k + 1 to n - 1 do
      if mat[k][k] <> 0.0 then
      begin
        factor := mat[i][k] / mat[k][k];
        for j := k to n - 1 do
          mat[i][j] := mat[i][j] - factor * mat[k][j];
        rhs[i] := rhs[i] - factor * rhs[k]
      end
  end;
  { back substitution }
  SetLength(x, n);
  for i := n - 1 downto 0 do
  begin
    x[i] := rhs[i];
    for j := i + 1 to n - 1 do
      x[i] := x[i] - mat[i][j] * x[j];
    if mat[i][i] <> 0.0 then
      x[i] := x[i] / mat[i][i]
  end
end;

{ ---------------------------------------------------------------------------
  6.1  LeastSquares — normal equations approach
  --------------------------------------------------------------------------- }
procedure LeastSquares(var A: TMatrix; var b: TVector; m, n: integer;
  var x: TVector);
var
  i, j, k: integer;
  AtA: TMatrix;
  Atb: TVector;
begin
  { build A^T A (n-by-n) and A^T b (n) }
  SetLength(AtA, n);
  for i := 0 to n - 1 do
  begin
    SetLength(AtA[i], n);
    for j := 0 to n - 1 do
      AtA[i][j] := 0.0
  end;
  SetLength(Atb, n);
  for i := 0 to n - 1 do
    Atb[i] := 0.0;
  for k := 0 to m - 1 do
    for i := 0 to n - 1 do
    begin
      Atb[i] := Atb[i] + A[k][i] * b[k];
      for j := 0 to n - 1 do
        AtA[i][j] := AtA[i][j] + A[k][i] * A[k][j]
    end;
  SolveLinearSystem(AtA, Atb, n, x)
end;

{ ---------------------------------------------------------------------------
  6.2  PolyRegression — Vandermonde matrix + LeastSquares
  --------------------------------------------------------------------------- }
procedure PolyRegression(var xdata, ydata: TVector; ndata, degree: integer;
  var coeffs: TVector);
var
  i, j: integer;
  V: TMatrix;
  ycopy: TVector;
begin
  { build Vandermonde matrix V[i][j] = xdata[i]^j, i=0..ndata-1, j=0..degree }
  SetLength(V, ndata);
  for i := 0 to ndata - 1 do
  begin
    SetLength(V[i], degree + 1);
    V[i][0] := 1.0;
    for j := 1 to degree do
      V[i][j] := V[i][j - 1] * xdata[i]
  end;
  SetLength(ycopy, ndata);
  for i := 0 to ndata - 1 do
    ycopy[i] := ydata[i];
  LeastSquares(V, ycopy, ndata, degree + 1, coeffs)
end;

{ ---------------------------------------------------------------------------
  6.3  LinearRegression — closed-form Pearson formulas
  --------------------------------------------------------------------------- }
procedure LinearRegression(var xdata, ydata: TVector; ndata: integer;
  var a, b, r: Float);
var
  i: integer;
  nf, sumx, sumy, sumxy, sumx2, sumy2: Float;
  xbar, ybar, sxx, syy, sxy: Float;
begin
  nf := ndata;
  sumx := 0.0; sumy := 0.0; sumxy := 0.0; sumx2 := 0.0; sumy2 := 0.0;
  for i := 0 to ndata - 1 do
  begin
    sumx  := sumx  + xdata[i];
    sumy  := sumy  + ydata[i];
    sumxy := sumxy + xdata[i] * ydata[i];
    sumx2 := sumx2 + xdata[i] * xdata[i];
    sumy2 := sumy2 + ydata[i] * ydata[i]
  end;
  xbar := sumx / nf;
  ybar := sumy / nf;
  sxx  := sumx2 - nf * xbar * xbar;
  syy  := sumy2 - nf * ybar * ybar;
  sxy  := sumxy - nf * xbar * ybar;
  if sxx = 0.0 then
  begin
    b := 0.0;
    a := ybar;
    r := 0.0
  end
  else
  begin
    b := sxy / sxx;
    a := ybar - b * xbar;
    if (sxx > 0.0) and (syy > 0.0) then
      r := sxy / Sqrt(sxx * syy)
    else
      r := 0.0
  end
end;

{ ---------------------------------------------------------------------------
  6.4  MultipleRegression — adds intercept column, calls LeastSquares, R²
  --------------------------------------------------------------------------- }
procedure MultipleRegression(var X: TMatrix; var y: TVector; m, n: integer;
  var coeffs: TVector; var r2: Float);
var
  i, j: integer;
  Xaug: TMatrix;
  ycopy: TVector;
  ybar, sst, sse, yhat: Float;
begin
  { build augmented [1 | X] design matrix }
  SetLength(Xaug, m);
  for i := 0 to m - 1 do
  begin
    SetLength(Xaug[i], n + 1);
    Xaug[i][0] := 1.0;
    for j := 0 to n - 1 do
      Xaug[i][j + 1] := X[i][j]
  end;
  SetLength(ycopy, m);
  for i := 0 to m - 1 do
    ycopy[i] := y[i];
  LeastSquares(Xaug, ycopy, m, n + 1, coeffs);
  { compute R² = 1 - SSE/SST }
  ybar := 0.0;
  for i := 0 to m - 1 do
    ybar := ybar + y[i];
  ybar := ybar / m;
  sst := 0.0;
  sse := 0.0;
  for i := 0 to m - 1 do
  begin
    yhat := coeffs[0];
    for j := 0 to n - 1 do
      yhat := yhat + coeffs[j + 1] * X[i][j];
    sst := sst + Sqr(y[i] - ybar);
    sse := sse + Sqr(y[i] - yhat)
  end;
  if sst > 0.0 then
    r2 := 1.0 - sse / sst
  else
    r2 := 1.0
end;

{ ---------------------------------------------------------------------------
  6.5  ChiSquareFit — weighted polynomial fit, weight = 1/sigma²
  --------------------------------------------------------------------------- }
procedure ChiSquareFit(var xdata, ydata, sigma: TVector; ndata, degree: integer;
  var coeffs: TVector; var chiSq: Float);
var
  i, j, k: integer;
  V: TMatrix;
  yw: TVector;
  wi, yfit: Float;
begin
  { weighted Vandermonde: V[i][j] = x[i]^j / sigma[i],  yw[i] = y[i] / sigma[i] }
  SetLength(V, ndata);
  SetLength(yw, ndata);
  for i := 0 to ndata - 1 do
  begin
    wi := 1.0 / sigma[i];
    SetLength(V[i], degree + 1);
    V[i][0] := wi;
    for j := 1 to degree do
      V[i][j] := V[i][j - 1] * xdata[i];
    yw[i] := ydata[i] * wi
  end;
  LeastSquares(V, yw, ndata, degree + 1, coeffs);
  { chi-square statistic }
  chiSq := 0.0;
  for i := 0 to ndata - 1 do
  begin
    yfit := 0.0;
    for k := 0 to degree do
      yfit := yfit + coeffs[k] * Power(xdata[i], k);
    chiSq := chiSq + Sqr((ydata[i] - yfit) / sigma[i])
  end
end;

{ ---------------------------------------------------------------------------
  self_test — hard-coded inputs, WriteLn results, no ReadLn
  --------------------------------------------------------------------------- }
procedure self_test;
const
  TOL = 1e-2;

  procedure check(const name: string; computed, expected: Float);
  var
    ok: string;
  begin
    if Abs(computed - expected) <= TOL * (1.0 + Abs(expected)) then
      ok := 'OK'
    else
      ok := '*** FAIL ***';
    WriteLn(Format('  %-44s computed=%9.5f  expected=%9.5f  %s',
                   [name, computed, expected, ok]))
  end;

var
  xd, yd, sig, coeffs: TVector;
  Am: TMatrix;
  ym: TVector;
  la, lb, lr, r2, chiSq: Float;
  i: integer;
begin
  WriteLn('=== jpmLstsqr self_test ===');
  WriteLn;

  { -----------------------------------------------------------------------
    Test 1 — LinearRegression: x=[1..5], y=[2,4,5,4,5]
    Expected: a≈2.2, b≈0.6, r≈0.7746  (spec says 0.816 but Pearson gives 0.7746)
    ----------------------------------------------------------------------- }
  WriteLn('Test 1 — LinearRegression');
  SetLength(xd, 5);
  SetLength(yd, 5);
  xd[0] := 1.0; xd[1] := 2.0; xd[2] := 3.0; xd[3] := 4.0; xd[4] := 5.0;
  yd[0] := 2.0; yd[1] := 4.0; yd[2] := 5.0; yd[3] := 4.0; yd[4] := 5.0;
  LinearRegression(xd, yd, 5, la, lb, lr);
  check('a (intercept)', la, 2.2);
  check('b (slope)',     lb, 0.6);
  check('r (Pearson)',   lr, 0.7746);
  WriteLn;

  { -----------------------------------------------------------------------
    Test 2 — PolyRegression degree 1 (same data)
    Expected: coeffs[0]≈2.2, coeffs[1]≈0.6
    ----------------------------------------------------------------------- }
  WriteLn('Test 2 — PolyRegression degree 1');
  PolyRegression(xd, yd, 5, 1, coeffs);
  check('coeffs[0] (intercept)', coeffs[0], 2.2);
  check('coeffs[1] (slope)',     coeffs[1], 0.6);
  WriteLn;

  { -----------------------------------------------------------------------
    Test 3 — PolyRegression degree 2: x=[0..4], y=[1,0,1,4,9]  (≈(x-1)²)
    Expected: coeffs[0]≈1, coeffs[1]≈-2, coeffs[2]≈1
    ----------------------------------------------------------------------- }
  WriteLn('Test 3 — PolyRegression degree 2');
  SetLength(xd, 5);
  SetLength(yd, 5);
  xd[0] := 0.0; xd[1] := 1.0; xd[2] := 2.0; xd[3] := 3.0; xd[4] := 4.0;
  yd[0] := 1.0; yd[1] := 0.0; yd[2] := 1.0; yd[3] := 4.0; yd[4] := 9.0;
  PolyRegression(xd, yd, 5, 2, coeffs);
  check('coeffs[0] (constant)', coeffs[0],  1.0);
  check('coeffs[1] (linear)',   coeffs[1], -2.0);
  check('coeffs[2] (quadratic)',coeffs[2],  1.0);
  WriteLn;

  { -----------------------------------------------------------------------
    Test 4 — LeastSquares: A=[[1,1],[1,2],[1,3]], b=[1,2,2]
    Normal equations give x≈[0.6667, 0.5]
    (spec says [0.333,0.5] which appears to be a typo; correct LS solution is [2/3, 1/2])
    ----------------------------------------------------------------------- }
  WriteLn('Test 4 — LeastSquares (overdetermined 3×2 system)');
  SetLength(Am, 3);
  for i := 0 to 2 do
    SetLength(Am[i], 2);
  Am[0][0] := 1.0; Am[0][1] := 1.0;
  Am[1][0] := 1.0; Am[1][1] := 2.0;
  Am[2][0] := 1.0; Am[2][1] := 3.0;
  SetLength(ym, 3);
  ym[0] := 1.0; ym[1] := 2.0; ym[2] := 2.0;
  LeastSquares(Am, ym, 3, 2, coeffs);
  WriteLn(Format('  x[0] = %9.5f  (spec says ~0.333, normal-eq gives ~0.6667)',
                 [coeffs[0]]));
  WriteLn(Format('  x[1] = %9.5f  (expected ~0.5)', [coeffs[1]]));
  check('x[1]', coeffs[1], 0.5);
  WriteLn;

  { -----------------------------------------------------------------------
    Test 5 — MultipleRegression: X=[[1],[2],[3],[4]], y=[3,5,7,9]
    Expected: intercept≈1, slope≈2, R²≈1.0
    ----------------------------------------------------------------------- }
  WriteLn('Test 5 — MultipleRegression');
  SetLength(Am, 4);
  for i := 0 to 3 do
  begin
    SetLength(Am[i], 1);
    Am[i][0] := i + 1
  end;
  SetLength(ym, 4);
  ym[0] := 3.0; ym[1] := 5.0; ym[2] := 7.0; ym[3] := 9.0;
  MultipleRegression(Am, ym, 4, 1, coeffs, r2);
  check('intercept (coeffs[0])', coeffs[0], 1.0);
  check('slope     (coeffs[1])', coeffs[1], 2.0);
  check('R²',                    r2,        1.0);
  WriteLn;

  { -----------------------------------------------------------------------
    Test 6 — ChiSquareFit: y=(x-1)² with uniform sigma=0.1
    Same data as Test 3; expect same coefficients, chiSq≈0
    ----------------------------------------------------------------------- }
  WriteLn('Test 6 — ChiSquareFit (weighted poly, sigma=0.1)');
  SetLength(xd, 5);
  SetLength(yd, 5);
  SetLength(sig, 5);
  for i := 0 to 4 do
  begin
    xd[i]  := i;
    yd[i]  := Sqr(i - 1);
    sig[i] := 0.1
  end;
  ChiSquareFit(xd, yd, sig, 5, 2, coeffs, chiSq);
  check('coeffs[0] (constant)', coeffs[0],  1.0);
  check('coeffs[1] (linear)',   coeffs[1], -2.0);
  check('coeffs[2] (quadratic)',coeffs[2],  1.0);
  check('chiSq (≈0)',           chiSq,      0.0);
  WriteLn;

  WriteLn('=== done ===')
end;

end.
