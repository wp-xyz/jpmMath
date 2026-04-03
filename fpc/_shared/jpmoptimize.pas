unit jpmOptimize;
{$mode objfpc}{$H+}

interface

uses
  SysUtils, Math, jpmTypes;

type
  TFloatArray  = array of Float;
  TFloatMatrix = array of TFloatArray;
  TFuncND      = function(var x: TFloatArray; n: integer): Float;

{ Golden Section Search: finds minimum of f in bracket (ax,bx,cx).
  f(bx) < f(ax) and f(bx) < f(cx). Returns f(xmin). }
function GoldenSection(f: TFunction1; ax, bx, cx: Float; tol: Float;
  var xmin: Float): Float;

{ Brent's method: parabolic interpolation + golden section fallback.
  Returns f(xmin). }
function BrentMin(f: TFunction1; ax, bx, cx: Float; tol: Float;
  var xmin: Float): Float;

{ BracketMin: starting from ax,bx, expand downhill until bracket is found.
  On exit: f(bx) < f(ax) and f(bx) < f(cx). }
procedure BracketMin(f: TFunction1; var ax, bx, cx: Float;
  var fa, fb, fc: Float);

{ Nelder-Mead downhill simplex minimization.
  p is (n+1) x n simplex; vertices p[0..n] each of length n.
  Minimizes f(x,n). Stops when simplex size < tol or maxIter reached. }
procedure NelderMead(f: TFuncND; var p: TFloatMatrix; n: integer; tol: Float;
  var nIter: integer; maxIter: integer);

{ Steepest Descent with numerical gradient (central differences, step h).
  Line search uses BrentMin. }
procedure SteepestDescent(f: TFuncND; var x: TFloatArray; n: integer;
  tol, h: Float; maxIter: integer; var nIter: integer);

procedure self_test;

implementation

const
  CGold    = 0.38196601125;   { 1 - golden ratio complement }
  CGoldR   = 0.61803398875;   { golden ratio }
  CBrentIt = 100;
  CZeps    = 1.0e-10;

{***********************************************************************
* GoldenSection
* Given bracket (ax, bx, cx) with f(bx) < f(ax), f(bx) < f(cx),
* performs golden-section search to precision tol.
***********************************************************************}
function GoldenSection(f: TFunction1; ax, bx, cx: Float; tol: Float;
  var xmin: Float): Float;
var
  c, r, x0, x1, x2, x3, f1, f2: Float;
begin
  r := CGoldR;
  c := 1.0 - r;
  x0 := ax;
  x3 := cx;
  if Abs(cx - bx) > Abs(bx - ax) then
  begin
    x1 := bx;
    x2 := bx + c * (cx - bx)
  end
  else
  begin
    x2 := bx;
    x1 := bx - c * (bx - ax)
  end;
  f1 := f(x1);
  f2 := f(x2);
  while Abs(x3 - x0) > tol * (Abs(x1) + Abs(x2)) do
  begin
    if f2 < f1 then
    begin
      x0 := x1;
      x1 := x2;
      x2 := r * x1 + c * x3;
      f1 := f2;
      f2 := f(x2)
    end
    else
    begin
      x3 := x2;
      x2 := x1;
      x1 := r * x2 + c * x0;
      f2 := f1;
      f1 := f(x1)
    end
  end;
  if f1 < f2 then
  begin
    Result := f1;
    xmin   := x1
  end
  else
  begin
    Result := f2;
    xmin   := x2
  end
end;

{***********************************************************************
* BrentMin
* Brent's method: combines golden section and parabolic interpolation.
* Bracket: (ax, bx, cx) with f(bx) < f(ax), f(bx) < f(cx).
***********************************************************************}
function BrentMin(f: TFunction1; ax, bx, cx: Float; tol: Float;
  var xmin: Float): Float;
var
  a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm: Float;
  iter: integer;
  doParabolic: boolean;
begin
  if ax < cx then
    a := ax
  else
    a := cx;
  if ax > cx then
    b := ax
  else
    b := cx;
  v  := bx;
  w  := v;
  x  := v;
  d  := 0.0;
  e  := 0.0;
  fx := f(x);
  fv := fx;
  fw := fx;
  for iter := 1 to CBrentIt do
  begin
    xm   := 0.5 * (a + b);
    tol1 := tol * Abs(x) + CZeps;
    tol2 := 2.0 * tol1;
    if Abs(x - xm) <= (tol2 - 0.5 * (b - a)) then
    begin
      xmin   := x;
      Result := fx;
      exit
    end;
    doParabolic := false;
    if Abs(e) > tol1 then
    begin
      r := (x - w) * (fx - fv);
      q := (x - v) * (fx - fw);
      p := (x - v) * q - (x - w) * r;
      q := 2.0 * (q - r);
      if q > 0.0 then
        p := -p;
      q      := Abs(q);
      etemp  := e;
      e      := d;
      if (Abs(p) < Abs(0.5 * q * etemp)) and
         (p > q * (a - x)) and
         (p < q * (b - x)) then
      begin
        d := p / q;
        u := x + d;
        if (u - a < tol2) or (b - u < tol2) then
        begin
          if xm >= x then
            d := tol1
          else
            d := -tol1
        end;
        doParabolic := true
      end
    end;
    if not doParabolic then
    begin
      if x >= xm then
        e := a - x
      else
        e := b - x;
      d := CGold * e
    end;
    if Abs(d) >= tol1 then
      u := x + d
    else
    begin
      if d >= 0.0 then
        u := x + tol1
      else
        u := x - tol1
    end;
    fu := f(u);
    if fu <= fx then
    begin
      if u >= x then
        a := x
      else
        b := x;
      v  := w;   fv := fw;
      w  := x;   fw := fx;
      x  := u;   fx := fu
    end
    else
    begin
      if u < x then
        a := u
      else
        b := u;
      if (fu <= fw) or (w = x) then
      begin
        v  := w;   fv := fw;
        w  := u;   fw := fu
      end
      else if (fu <= fv) or (v = x) or (v = w) then
      begin
        v  := u;
        fv := fu
      end
    end
  end;
  { max iterations reached }
  xmin   := x;
  Result := fx
end;

{***********************************************************************
* BracketMin
* Starting from (ax, bx), expand downhill to find bracket
* (ax, bx, cx) where f(bx) < f(ax) and f(bx) < f(cx).
***********************************************************************}
procedure BracketMin(f: TFunction1; var ax, bx, cx: Float;
  var fa, fb, fc: Float);
const
  CGold2  = 1.61803398875;
  CGlimit = 100.0;
  CTiny   = 1.0e-20;
var
  fu, r2, q2, ulim, u, tmp: Float;
begin
  fa := f(ax);
  fb := f(bx);
  if fb > fa then
  begin
    { swap so we go downhill from ax to bx }
    tmp := ax; ax := bx; bx := tmp;
    tmp := fa; fa := fb; fb := tmp
  end;
  cx := bx + CGold2 * (bx - ax);
  fc := f(cx);
  while fb >= fc do
  begin
    r2   := (bx - ax) * (fb - fc);
    q2   := (bx - cx) * (fb - fa);
    tmp  := q2 - r2;
    if Abs(tmp) < CTiny then
      tmp := CTiny;
    u    := bx - ((bx - cx) * q2 - (bx - ax) * r2) / (2.0 * tmp);
    ulim := bx + CGlimit * (cx - bx);
    if (bx - u) * (u - cx) > 0.0 then
    begin
      fu := f(u);
      if fu < fc then
      begin
        ax := bx; fa := fb;
        bx := u;  fb := fu;
        exit
      end
      else if fu > fb then
      begin
        cx := u; fc := fu;
        exit
      end;
      u  := cx + CGold2 * (cx - bx);
      fu := f(u)
    end
    else if (cx - u) * (u - ulim) > 0.0 then
    begin
      fu := f(u);
      if fu < fc then
      begin
        bx := cx; fb := fc;
        cx := u;  fc := fu;
        u  := cx + CGold2 * (cx - bx);
        fu := f(u)
      end
    end
    else if (u - ulim) * (ulim - cx) >= 0.0 then
    begin
      u  := ulim;
      fu := f(u)
    end
    else
    begin
      u  := cx + CGold2 * (cx - bx);
      fu := f(u)
    end;
    ax := bx; fa := fb;
    bx := cx; fb := fc;
    cx := u;  fc := fu
  end
end;

{***********************************************************************
* NelderMead
* Downhill simplex method (Nelder & Mead 1965).
* p[0..n] are the n+1 vertices, each a TFloatArray of length n.
* y[i] = f(p[i]) must NOT be pre-initialised by the caller —
* this routine initialises y internally.
***********************************************************************}
procedure NelderMead(f: TFuncND; var p: TFloatMatrix; n: integer; tol: Float;
  var nIter: integer; maxIter: integer);
const
  CAlpha = 1.0;
  CBeta  = 0.5;
  CGamma2 = 2.0;
var
  y, pbar, pr, prr: TFloatArray;
  i, j, ilo, ihi, inhi, mpts: integer;
  rtol, ypr, yprr: Float;
  converged: boolean;
begin
  mpts  := n + 1;
  nIter := 0;
  SetLength(y,    mpts);
  SetLength(pbar, n);
  SetLength(pr,   n);
  SetLength(prr,  n);
  for i := 0 to mpts - 1 do
    y[i] := f(p[i], n);
  converged := false;
  while (not converged) and (nIter < maxIter) do
  begin
    { find ilo, ihi, inhi }
    ilo := 0;
    if y[0] > y[1] then
    begin
      ihi  := 0;
      inhi := 1
    end
    else
    begin
      ihi  := 1;
      inhi := 0
    end;
    for i := 0 to mpts - 1 do
    begin
      if y[i] < y[ilo] then
        ilo := i;
      if y[i] > y[ihi] then
      begin
        inhi := ihi;
        ihi  := i
      end
      else if (y[i] > y[inhi]) and (i <> ihi) then
        inhi := i
    end;
    rtol := 2.0 * Abs(y[ihi] - y[ilo]) /
            (Abs(y[ihi]) + Abs(y[ilo]) + CZeps);
    if rtol < tol then
    begin
      converged := true;
      break
    end;
    Inc(nIter);
    { compute centroid pbar (excluding ihi) }
    for j := 0 to n - 1 do
      pbar[j] := 0.0;
    for i := 0 to mpts - 1 do
      if i <> ihi then
        for j := 0 to n - 1 do
          pbar[j] := pbar[j] + p[i][j];
    for j := 0 to n - 1 do
    begin
      pbar[j] := pbar[j] / n;
      pr[j]   := (1.0 + CAlpha) * pbar[j] - CAlpha * p[ihi][j]
    end;
    ypr := f(pr, n);
    if ypr <= y[ilo] then
    begin
      { try expansion }
      for j := 0 to n - 1 do
        prr[j] := CGamma2 * pr[j] + (1.0 - CGamma2) * pbar[j];
      yprr := f(prr, n);
      if yprr < y[ilo] then
      begin
        for j := 0 to n - 1 do
          p[ihi][j] := prr[j];
        y[ihi] := yprr
      end
      else
      begin
        for j := 0 to n - 1 do
          p[ihi][j] := pr[j];
        y[ihi] := ypr
      end
    end
    else if ypr >= y[inhi] then
    begin
      if ypr < y[ihi] then
      begin
        for j := 0 to n - 1 do
          p[ihi][j] := pr[j];
        y[ihi] := ypr
      end;
      { contraction }
      for j := 0 to n - 1 do
        prr[j] := CBeta * p[ihi][j] + (1.0 - CBeta) * pbar[j];
      yprr := f(prr, n);
      if yprr < y[ihi] then
      begin
        for j := 0 to n - 1 do
          p[ihi][j] := prr[j];
        y[ihi] := yprr
      end
      else
      begin
        { shrink toward best point }
        for i := 0 to mpts - 1 do
          if i <> ilo then
          begin
            for j := 0 to n - 1 do
            begin
              pr[j]     := 0.5 * (p[i][j] + p[ilo][j]);
              p[i][j]   := pr[j]
            end;
            y[i] := f(pr, n)
          end
      end
    end
    else
    begin
      for j := 0 to n - 1 do
        p[ihi][j] := pr[j];
      y[ihi] := ypr
    end
  end
end;

{***********************************************************************
* SteepestDescent
* Gradient-based minimization using numerical central-difference
* gradient and a BrentMin line search along -gradient.
***********************************************************************}

{ module-level state for the line-search closure }
var
  GSD_f:    TFuncND;
  GSD_x:    TFloatArray;
  GSD_dir:  TFloatArray;
  GSD_n:    integer;
  GSD_tmp:  TFloatArray;

function LineFunc(t: Float): Float;
var
  j: integer;
begin
  for j := 0 to GSD_n - 1 do
    GSD_tmp[j] := GSD_x[j] + t * GSD_dir[j];
  Result := GSD_f(GSD_tmp, GSD_n)
end;

procedure SteepestDescent(f: TFuncND; var x: TFloatArray; n: integer;
  tol, h: Float; maxIter: integer; var nIter: integer);
var
  grad, xp, xm: TFloatArray;
  i, j: integer;
  fx, fnew, gnorm, ax, bx, cx, fa, fb, fc, tmin: Float;
  converged: boolean;
begin
  SetLength(grad, n);
  SetLength(xp,   n);
  SetLength(xm,   n);
  { set up module-level state for LineFunc }
  GSD_f   := f;
  GSD_n   := n;
  SetLength(GSD_x,   n);
  SetLength(GSD_dir, n);
  SetLength(GSD_tmp, n);
  nIter     := 0;
  converged := false;
  while (not converged) and (nIter < maxIter) do
  begin
    fx := f(x, n);
    { numerical gradient via central differences }
    for i := 0 to n - 1 do
    begin
      for j := 0 to n - 1 do
      begin
        xp[j] := x[j];
        xm[j] := x[j]
      end;
      xp[i] := x[i] + h;
      xm[i] := x[i] - h;
      grad[i] := (f(xp, n) - f(xm, n)) / (2.0 * h)
    end;
    { compute gradient norm }
    gnorm := 0.0;
    for i := 0 to n - 1 do
      gnorm := gnorm + grad[i] * grad[i];
    gnorm := Sqrt(gnorm);
    if gnorm < tol then
    begin
      converged := true;
      break
    end;
    { direction = -grad/|grad| (unit vector along steepest descent) }
    for i := 0 to n - 1 do
      GSD_dir[i] := -grad[i] / gnorm;
    for i := 0 to n - 1 do
      GSD_x[i] := x[i];
    { bracket minimum along this direction }
    ax := 0.0;
    bx := 1.0;
    BracketMin(@LineFunc, ax, bx, cx, fa, fb, fc);
    { ensure proper bracket ordering for BrentMin }
    if fa < fb then
    begin
      tmin := GoldenSection(@LineFunc, ax, bx, cx, tol, tmin)
    end
    else
      BrentMin(@LineFunc, ax, bx, cx, tol, tmin);
    { move x in the search direction }
    for i := 0 to n - 1 do
      x[i] := x[i] + tmin * GSD_dir[i];
    fnew := f(x, n);
    Inc(nIter);
    if Abs(fnew - fx) < tol then
      converged := true
  end
end;

{***********************************************************************
* self_test
***********************************************************************}

{ 1D test functions }
function TestF(x: Float): Float;
begin
  Result := (x - 2.0) * (x - 2.0) + 1.0
end;

function TestG(x: Float): Float;
begin
  Result := x * x * x * x - 4.0 * x * x + x
end;

{ 2D test function h(x,y) = (x-1)^2 + (y-2)^2, min at (1,2) = 0 }
function TestH(var v: TFloatArray; n: integer): Float;
begin
  Result := (v[0] - 1.0) * (v[0] - 1.0) + (v[1] - 2.0) * (v[1] - 2.0)
end;

procedure self_test;
var
  xmin, fmin: Float;
  ax, bx, cx, fa, fb, fc: Float;
  p: TFloatMatrix;
  x2d: TFloatArray;
  nIter: integer;
  i: integer;
begin
  WriteLn('=== jpmOptimize Self Test ===');
  WriteLn;

  { --- GoldenSection on f(x)=(x-2)^2+1, min at x=2, f=1 --- }
  WriteLn('1) GoldenSection: f(x)=(x-2)^2+1, bracket (0,1,4)');
  xmin := 0.0;
  fmin := GoldenSection(@TestF, 0.0, 1.0, 4.0, 1e-9, xmin);
  WriteLn('   xmin = ', xmin:12:8, '  (expected 2.0)');
  WriteLn('   fmin = ', fmin:12:8, '  (expected 1.0)');
  WriteLn;

  { --- BrentMin on f(x)=(x-2)^2+1 --- }
  WriteLn('2) BrentMin: f(x)=(x-2)^2+1, bracket (0,1,4)');
  xmin := 0.0;
  fmin := BrentMin(@TestF, 0.0, 1.0, 4.0, 1e-10, xmin);
  WriteLn('   xmin = ', xmin:12:8, '  (expected 2.0)');
  WriteLn('   fmin = ', fmin:12:8, '  (expected 1.0)');
  WriteLn;

  { --- BracketMin + BrentMin on g(x)=x^4-4x^2+x, local min near x~1.3 --- }
  WriteLn('3) BracketMin+BrentMin: g(x)=x^4-4x^2+x, start (1.0,1.5)');
  ax := 1.0; bx := 1.5;
  BracketMin(@TestG, ax, bx, cx, fa, fb, fc);
  WriteLn('   bracket: ax=', ax:8:4, '  bx=', bx:8:4, '  cx=', cx:8:4);
  xmin := 0.0;
  fmin := BrentMin(@TestG, ax, bx, cx, 1e-10, xmin);
  WriteLn('   xmin = ', xmin:12:8, '  (expected ~1.3)');
  WriteLn('   fmin = ', fmin:12:8);
  WriteLn;

  { --- NelderMead on h(x,y)=(x-1)^2+(y-2)^2 --- }
  WriteLn('4) NelderMead: h(x,y)=(x-1)^2+(y-2)^2, min at (1,2), h=0');
  SetLength(p, 3);
  for i := 0 to 2 do
    SetLength(p[i], 2);
  p[0][0] :=  0.0; p[0][1] :=  0.0;
  p[1][0] :=  3.0; p[1][1] :=  0.0;
  p[2][0] :=  0.0; p[2][1] :=  4.0;
  nIter := 0;
  NelderMead(@TestH, p, 2, 1e-10, nIter, 500);
  WriteLn('   x=', p[0][0]:12:8, '  y=', p[0][1]:12:8,
          '  (expected 1.0, 2.0)');
  WriteLn('   f=', TestH(p[0], 2):12:8, '  (expected 0.0)');
  WriteLn('   iterations=', nIter);
  WriteLn;

  { --- SteepestDescent on h(x,y) --- }
  WriteLn('5) SteepestDescent: h(x,y)=(x-1)^2+(y-2)^2, start (0,0)');
  SetLength(x2d, 2);
  x2d[0] := 0.0; x2d[1] := 0.0;
  nIter := 0;
  SteepestDescent(@TestH, x2d, 2, 1e-8, 1e-5, 200, nIter);
  WriteLn('   x=', x2d[0]:12:8, '  y=', x2d[1]:12:8,
          '  (expected 1.0, 2.0)');
  WriteLn('   f=', TestH(x2d, 2):12:8, '  (expected 0.0)');
  WriteLn('   iterations=', nIter);
  WriteLn;

  WriteLn('=== End Self Test ===')
end;

end.
