
unit jpmAnneal;
{$mode objfpc}{$H+}

interface

uses SysUtils, Math, jpmtypes;

type
  TAnnealResult = record
    x: TFloatArray;      { best solution found }
    fval: Float;         { best function value }
    iterations: integer; { iterations performed }
    accepted: integer;   { accepted moves }
  end;

{ Simulated Annealing for continuous minimization
  f        : objective function (N-dimensional)
  x0       : initial point (length N)
  n        : dimension
  t0       : initial temperature
  alpha    : cooling factor (e.g. 0.99)
  maxiter  : max iterations
  step     : initial step size for random perturbation
  Returns TAnnealResult with best x and f(x) }
function SimulatedAnnealing(f: TFunctionN; var x0: TFloatArray; n: integer;
  t0, alpha, step: Float; maxiter: integer): TAnnealResult;

procedure self_test;

implementation

{ TFunctionN is function(var x: TFloatArray): Float - check jpmtypes }

function SimulatedAnnealing(f: TFunctionN; var x0: TFloatArray; n: integer;
  t0, alpha, step: Float; maxiter: integer): TAnnealResult;
var
  xcur, xbest, xnew: TFloatArray;
  fcur, fbest, fnew, temp, delta, r: Float;
  iter, acc, i: integer;
begin
  setlength(xcur,  n);
  setlength(xbest, n);
  setlength(xnew,  n);
  for i := 0 to n-1 do
  begin
    xcur[i]  := x0[i];
    xbest[i] := x0[i]
  end;
  fcur  := f(xcur);
  fbest := fcur;
  temp  := t0;
  acc   := 0;

  RandSeed := 42; { reproducible }

  for iter := 1 to maxiter do
  begin
    { generate neighbour }
    for i := 0 to n-1 do
      xnew[i] := xcur[i] + step * (2.0*Random - 1.0);

    fnew  := f(xnew);
    delta := fnew - fcur;

    { accept? }
    if delta < 0 then
      r := 1.0
    else if temp > 1e-300 then
      r := Exp(-delta / temp)
    else
      r := 0.0;

    if Random < r then
    begin
      for i := 0 to n-1 do
        xcur[i] := xnew[i];
      fcur := fnew;
      inc(acc);
      if fcur < fbest then
      begin
        fbest := fcur;
        for i := 0 to n-1 do
          xbest[i] := xnew[i]
      end
    end;

    temp := temp * alpha
  end;

  result.x          := xbest;
  result.fval       := fbest;
  result.iterations := maxiter;
  result.accepted   := acc
end;

{ Test functions }
function Sphere(var x: TFloatArray): Float;
var
  i: integer;
  s: Float;
begin
  s := 0;
  for i := 0 to length(x)-1 do
    s := s + x[i]*x[i];
  result := s
end;

function Rosenbrock(var x: TFloatArray): Float;
begin
  { f(x,y) = (1-x)^2 + 100*(y-x^2)^2, min at (1,1)=0 }
  result := Sqr(1 - x[0]) + 100*Sqr(x[1] - x[0]*x[0])
end;

procedure self_test;
var
  x0: TFloatArray;
  res: TAnnealResult;
begin
  writeln('=== jpmAnneal self_test ===');
  writeln;

  { Test 1: minimize sphere f(x,y) = x^2+y^2, min at (0,0)=0 }
  setlength(x0, 2);
  x0[0] := 3.0; x0[1] := 4.0;
  res := SimulatedAnnealing(@Sphere, x0, 2, 10.0, 0.995, 0.5, 50000);
  writeln('Sphere f(x,y)=x^2+y^2:');
  writeln('  x=', res.x[0]:8:4, ' y=', res.x[1]:8:4,
          ' fval=', res.fval:10:6, ' accepted=', res.accepted);
  if res.fval < 0.01 then
    writeln('  PASS (fval < 0.01)')
  else
  begin
    writeln('  FAIL (expected fval < 0.01)');
    SelfTestFail('Anneal Sphere: fval=' + FloatToStr(res.fval) + ' expected < 0.01');
  end;

  writeln;

  { Test 2: Rosenbrock f(x,y)=(1-x)^2+100*(y-x^2)^2, min at (1,1)=0 }
  x0[0] := -1.0; x0[1] := 1.0;
  res := SimulatedAnnealing(@Rosenbrock, x0, 2, 5.0, 0.999, 0.3, 100000);
  writeln('Rosenbrock: min at (1,1)');
  writeln('  x=', res.x[0]:8:4, ' y=', res.x[1]:8:4,
          ' fval=', res.fval:10:6, ' accepted=', res.accepted);
  if res.fval < 0.1 then
    writeln('  PASS (fval < 0.1)')
  else
  begin
    writeln('  FAIL (expected fval < 0.1, got ', res.fval:8:4, ')');
    SelfTestFail('Anneal Rosenbrock: fval=' + FloatToStr(res.fval) + ' expected < 0.1');
  end;

  writeln;
  writeln('=== done ===')
end;

end.
