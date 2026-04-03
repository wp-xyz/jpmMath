
unit jpmGeometry;
{$mode objfpc}{$H+}

interface

uses SysUtils, Math, jpmtypes;

type
  TTriangle = record
    a, b, c: Float;
    ta, tb, tc: Float;
    area: Float;
    valid: boolean;
  end;

  TArcCircle = record
    r: Float;
    s: Float;
    c: Float;
    t: Float;
    sector_area: Float;
    segment_area: Float;
    valid: boolean;
  end;

function SolveTriangle(var tri: TTriangle): boolean;
function SolveArcCircle(var arc: TArcCircle): boolean;
function PointDistance2D(x1, y1, x2, y2: Float): Float;
function PointDistance3D(x1, y1, z1, x2, y2, z2: Float): Float;
function TriangleArea(x1, y1, x2, y2, x3, y3: Float): Float;
function PolygonArea(var px, py: TFloatArray; n: integer): Float;
function CircleArea(r: Float): Float;
function CirclePerimeter(r: Float): Float;
function CircleChord(r, angle: Float): Float;
function CircleArcLen(r, angle: Float): Float;
function EllipseArea(a, b: Float): Float;
function EllipsePerimeter(a, b: Float): Float;
function HyperbolaAsymptoteAngle(a, b: Float): Float;
procedure self_test;

implementation

function SolveTriangle(var tri: TTriangle): boolean;
var
  knownSides, knownAngles: integer;
begin
  tri.valid := false;
  result := false;

  knownSides := 0;
  knownAngles := 0;
  if tri.a > 0 then inc(knownSides);
  if tri.b > 0 then inc(knownSides);
  if tri.c > 0 then inc(knownSides);
  if tri.ta > 0 then inc(knownAngles);
  if tri.tb > 0 then inc(knownAngles);
  if tri.tc > 0 then inc(knownAngles);

  { SSS }
  if knownSides = 3 then
  begin
    tri.ta := ArcCos((tri.b*tri.b + tri.c*tri.c - tri.a*tri.a) / (2*tri.b*tri.c));
    tri.tb := ArcCos((tri.a*tri.a + tri.c*tri.c - tri.b*tri.b) / (2*tri.a*tri.c));
    tri.tc := Pi - tri.ta - tri.tb;
    tri.area := 0.5 * tri.a * tri.b * Sin(tri.tc);
    tri.valid := true;
    result := true;
    exit;
  end;

  { SAS: two sides + included angle }
  if (knownSides = 2) and (knownAngles = 1) then
  begin
    if (tri.a > 0) and (tri.b > 0) and (tri.tc > 0) then
    begin
      tri.c := Sqrt(tri.a*tri.a + tri.b*tri.b - 2*tri.a*tri.b*Cos(tri.tc));
      tri.ta := ArcSin(tri.a * Sin(tri.tc) / tri.c);
      tri.tb := Pi - tri.ta - tri.tc;
      tri.area := 0.5 * tri.a * tri.b * Sin(tri.tc);
      tri.valid := true;
      result := true;
      exit;
    end;
    if (tri.a > 0) and (tri.c > 0) and (tri.tb > 0) then
    begin
      tri.b := Sqrt(tri.a*tri.a + tri.c*tri.c - 2*tri.a*tri.c*Cos(tri.tb));
      tri.ta := ArcSin(tri.a * Sin(tri.tb) / tri.b);
      tri.tc := Pi - tri.ta - tri.tb;
      tri.area := 0.5 * tri.a * tri.c * Sin(tri.tb);
      tri.valid := true;
      result := true;
      exit;
    end;
    if (tri.b > 0) and (tri.c > 0) and (tri.ta > 0) then
    begin
      tri.a := Sqrt(tri.b*tri.b + tri.c*tri.c - 2*tri.b*tri.c*Cos(tri.ta));
      tri.tb := ArcSin(tri.b * Sin(tri.ta) / tri.a);
      tri.tc := Pi - tri.ta - tri.tb;
      tri.area := 0.5 * tri.b * tri.c * Sin(tri.ta);
      tri.valid := true;
      result := true;
      exit;
    end;
  end;

  { AAS / SAA: 2 angles + 1 side }
  if (knownAngles >= 2) and (knownSides >= 1) then
  begin
    { fill missing angle }
    if tri.ta <= 0 then tri.ta := Pi - tri.tb - tri.tc;
    if tri.tb <= 0 then tri.tb := Pi - tri.ta - tri.tc;
    if tri.tc <= 0 then tri.tc := Pi - tri.ta - tri.tb;
    { use law of sines to find missing sides }
    if tri.a > 0 then
    begin
      if tri.b <= 0 then tri.b := tri.a * Sin(tri.tb) / Sin(tri.ta);
      if tri.c <= 0 then tri.c := tri.a * Sin(tri.tc) / Sin(tri.ta);
    end
    else if tri.b > 0 then
    begin
      if tri.a <= 0 then tri.a := tri.b * Sin(tri.ta) / Sin(tri.tb);
      if tri.c <= 0 then tri.c := tri.b * Sin(tri.tc) / Sin(tri.tb);
    end
    else if tri.c > 0 then
    begin
      if tri.a <= 0 then tri.a := tri.c * Sin(tri.ta) / Sin(tri.tc);
      if tri.b <= 0 then tri.b := tri.c * Sin(tri.tb) / Sin(tri.tc);
    end;
    tri.area := 0.5 * tri.a * tri.b * Sin(tri.tc);
    tri.valid := true;
    result := true;
    exit;
  end;
end;

function SolveArcCircle(var arc: TArcCircle): boolean;
var
  known: integer;
begin
  arc.valid := false;
  result := false;
  known := 0;
  if arc.r > 0 then inc(known);
  if arc.s > 0 then inc(known);
  if arc.c > 0 then inc(known);
  if arc.t > 0 then inc(known);
  if known < 2 then exit;

  { derive r and t first }
  if (arc.r > 0) and (arc.t > 0) then
  begin
    arc.s := arc.r * arc.t;
    arc.c := 2 * arc.r * Sin(arc.t / 2);
  end
  else if (arc.r > 0) and (arc.s > 0) then
  begin
    arc.t := arc.s / arc.r;
    arc.c := 2 * arc.r * Sin(arc.t / 2);
  end
  else if (arc.r > 0) and (arc.c > 0) then
  begin
    arc.t := 2 * ArcSin(arc.c / (2 * arc.r));
    arc.s := arc.r * arc.t;
  end
  else if (arc.s > 0) and (arc.t > 0) then
  begin
    arc.r := arc.s / arc.t;
    arc.c := 2 * arc.r * Sin(arc.t / 2);
  end
  else if (arc.c > 0) and (arc.t > 0) then
  begin
    arc.r := arc.c / (2 * Sin(arc.t / 2));
    arc.s := arc.r * arc.t;
  end
  else if (arc.s > 0) and (arc.c > 0) then
  begin
    { iterative: s=r*t, c=2r*sin(t/2) → t/sin(t/2) = 2s/c }
    { approximate: t ≈ 2*arcsin(c/(2r)), need Newton iteration }
    { For simplicity use: if c≈s then t is small, use t≈c/r }
    exit; { underdetermined without iteration }
  end;

  arc.sector_area  := arc.r * arc.r * arc.t / 2;
  arc.segment_area := arc.sector_area - arc.r * arc.r * Sin(arc.t) / 2;
  arc.valid := true;
  result := true;
end;

function PointDistance2D(x1, y1, x2, y2: Float): Float;
begin
  result := Sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
end;

function PointDistance3D(x1, y1, z1, x2, y2, z2: Float): Float;
begin
  result := Sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
end;

function TriangleArea(x1, y1, x2, y2, x3, y3: Float): Float;
begin
  result := Abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) / 2;
end;

function PolygonArea(var px, py: TFloatArray; n: integer): Float;
var
  i, j: integer;
  s: Float;
begin
  s := 0;
  j := n - 1;
  for i := 0 to n - 1 do
  begin
    s := s + (px[j] + px[i]) * (py[j] - py[i]);
    j := i;
  end;
  result := Abs(s) / 2;
end;

function CircleArea(r: Float): Float;
begin
  result := Pi * r * r;
end;

function CirclePerimeter(r: Float): Float;
begin
  result := 2 * Pi * r;
end;

function CircleChord(r, angle: Float): Float;
begin
  result := 2 * r * Sin(angle / 2);
end;

function CircleArcLen(r, angle: Float): Float;
begin
  result := r * angle;
end;

function EllipseArea(a, b: Float): Float;
begin
  result := Pi * a * b;
end;

function EllipsePerimeter(a, b: Float): Float;
var
  h: Float;
begin
  h := (a - b) * (a - b) / ((a + b) * (a + b));
  result := Pi * (a + b) * (1 + 3*h / (10 + Sqrt(4 - 3*h)));
end;

function HyperbolaAsymptoteAngle(a, b: Float): Float;
begin
  result := ArcTan(b / a);
end;

procedure self_test;
var
  tri: TTriangle;
  arc: TArcCircle;
  px, py: TFloatArray;
  r: boolean;
begin
  WriteLn('=== jpmGeometry self_test ===');
  WriteLn;

  { SSS: 3-4-5 right triangle }
  FillChar(tri, SizeOf(tri), 0);
  tri.a := 3; tri.b := 4; tri.c := 5;
  r := SolveTriangle(tri);
  WriteLn('SSS 3-4-5: solved=', r, ' area=', tri.area:6:4,
          ' tc=', RadToDeg(tri.tc):6:2, 'deg  expected area=6, tc=90');

  { SAS }
  FillChar(tri, SizeOf(tri), 0);
  tri.a := 5; tri.b := 7; tri.tc := Pi/3;
  r := SolveTriangle(tri);
  WriteLn('SAS a=5,b=7,tc=60: solved=', r, ' c=', tri.c:6:4,
          ' area=', tri.area:6:4, '  expected c≈6.245, area≈15.155');

  { AAS }
  FillChar(tri, SizeOf(tri), 0);
  tri.ta := Pi/6; tri.tb := Pi/4; tri.a := 4;
  r := SolveTriangle(tri);
  WriteLn('AAS ta=30,tb=45,a=4: solved=', r, ' b=', tri.b:6:4,
          ' c=', tri.c:6:4, '  expected b≈5.657, c≈7.727');

  WriteLn;

  { Arc: chord+angle }
  FillChar(arc, SizeOf(arc), 0);
  arc.c := 2; arc.t := 1;
  r := SolveArcCircle(arc);
  WriteLn('Arc c=2,t=1: solved=', r, ' r=', arc.r:6:4,
          ' s=', arc.s:6:4, '  expected r≈2.0858, s≈2.0858');

  { Arc: r+angle }
  FillChar(arc, SizeOf(arc), 0);
  arc.r := 5; arc.t := Pi/3;
  r := SolveArcCircle(arc);
  WriteLn('Arc r=5,t=Pi/3: solved=', r, ' s=', arc.s:6:4,
          ' c=', arc.c:6:4, '  expected s≈5.236, c≈5.000');

  WriteLn;

  WriteLn('PointDistance2D(0,0,3,4)    = ', PointDistance2D(0,0,3,4):6:4,
          '  expected 5.0');
  WriteLn('TriangleArea(0,0,4,0,0,3)  = ', TriangleArea(0,0,4,0,0,3):6:4,
          '  expected 6.0');

  SetLength(px, 4); SetLength(py, 4);
  px[0]:=0; px[1]:=1; px[2]:=1; px[3]:=0;
  py[0]:=0; py[1]:=0; py[2]:=1; py[3]:=1;
  WriteLn('PolygonArea unit square    = ', PolygonArea(px, py, 4):6:4,
          '  expected 1.0');

  WriteLn('CircleArea(1)              = ', CircleArea(1):8:6,
          '  expected Pi≈3.141593');
  WriteLn('EllipseArea(3,4)           = ', EllipseArea(3,4):8:4,
          '  expected 37.6991');

  WriteLn;
  WriteLn('=== done ===');
end;

end.
