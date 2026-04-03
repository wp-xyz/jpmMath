
unit jpmSimplex;
{$mode objfpc}{$H+}

interface

uses SysUtils, Math, jpmtypes;

type
  TSimplexMatrix = array of array of Float;

function SimplexMaximize(nv, nc: integer; var ts: TSimplexMatrix;
  var solution: TFloatArray): boolean;
function SimplexMinimize(nv, nc: integer; var ts: TSimplexMatrix;
  var solution: TFloatArray): boolean;
procedure self_test;

implementation

{ Simplex method: maximize row 0 objective.
  ts[0,1..nv] = objective coefficients
  ts[0,0]     = 0
  ts[i,1..nv] = constraint i coefficients (i=1..nc), <= type
  ts[i,0]     = RHS of constraint i
}
function SimplexMaximize(nv, nc: integer; var ts: TSimplexMatrix;
  var solution: TFloatArray): boolean;
var
  i, j, p1, p2, iter: integer;
  xmax, rap, v, pivot: Float;
  found: boolean;
begin
  result := false;
  setlength(solution, nv + 1);

  { negate objective row so we look for most negative }
  for j := 1 to nv do
    ts[0, j] := -ts[0, j];

  for iter := 1 to 100 do
  begin
    { find most negative in row 0 → pivot column p2 }
    xmax := 0;
    p2 := 0;
    for j := 1 to nv do
      if ts[0, j] < xmax then
      begin
        xmax := ts[0, j];
        p2 := j
      end;
    if p2 = 0 then
      break; { optimal }

    { find pivot row p1: minimum positive ratio ts[i,0]/ts[i,p2] }
    rap := 1e308;
    p1 := 0;
    for i := 1 to nc do
      if ts[i, p2] > 1e-12 then
      begin
        v := ts[i, 0] / ts[i, p2];
        if v < rap then
        begin
          rap := v;
          p1 := i
        end
      end;
    if p1 = 0 then
      exit; { unbounded }

    { pivot: divide pivot row by pivot element }
    pivot := ts[p1, p2];
    for j := 0 to nv do
      ts[p1, j] := ts[p1, j] / pivot;

    { eliminate pivot column from all other rows }
    for i := 0 to nc do
      if i <> p1 then
      begin
        v := ts[i, p2];
        for j := 0 to nv do
          ts[i, j] := ts[i, j] - v * ts[p1, j]
      end
  end;

  { extract solution: column j is basic in row i if it has exactly one 1 }
  for j := 1 to nv do
  begin
    solution[j] := 0;
    found := false;
    p1 := 0;
    for i := 1 to nc do
      if Abs(ts[i, j] - 1.0) < 1e-8 then
      begin
        if found then
        begin
          found := false;
          break
        end;
        found := true;
        p1 := i
      end
      else if Abs(ts[i, j]) > 1e-8 then
      begin
        found := false;
        break
      end;
    if found then
      solution[j] := ts[p1, 0]
  end;
  result := true
end;

function SimplexMinimize(nv, nc: integer; var ts: TSimplexMatrix;
  var solution: TFloatArray): boolean;
var
  j: integer;
begin
  for j := 1 to nv do
    ts[0, j] := -ts[0, j];
  result := SimplexMaximize(nv, nc, ts, solution)
end;

procedure self_test;
var
  ts: TSimplexMatrix;
  sol: TFloatArray;
  ok: boolean;
  obj: Float;
begin
  writeln('=== jpmSimplex self_test ===');
  writeln;

  { Maximize 15*x1 + 17*x2 + 20*x3 subject to:
      x2 - x3 <= 2
      3*x1 + 3*x2 + 5*x3 <= 15
      3*x1 + 2*x2 + x3 <= 8
    Expected: x1≈0.333, x2=3, x3=1, obj=76 }
  setlength(ts, 4);
  setlength(ts[0], 4); setlength(ts[1], 4);
  setlength(ts[2], 4); setlength(ts[3], 4);
  ts[0,0] := 0;  ts[0,1] := 15; ts[0,2] := 17; ts[0,3] := 20;
  ts[1,0] := 2;  ts[1,1] := 0;  ts[1,2] := 1;  ts[1,3] := -1;
  ts[2,0] := 15; ts[2,1] := 3;  ts[2,2] := 3;  ts[2,3] := 5;
  ts[3,0] := 8;  ts[3,1] := 3;  ts[3,2] := 2;  ts[3,3] := 1;

  ok := SimplexMaximize(3, 3, ts, sol);
  obj := 15*sol[1] + 17*sol[2] + 20*sol[3];
  writeln('Max 15x1+17x2+20x3: solved=', ok);
  writeln('  x1=', sol[1]:6:4, ' x2=', sol[2]:6:4, ' x3=', sol[3]:6:4);
  writeln('  objective=', obj:6:2, '  expected=76.00');
  if Abs(obj - 76.0) < 0.01 then
    writeln('  PASS')
  else
    writeln('  FAIL');

  writeln;
  writeln('=== done ===')
end;

end.
