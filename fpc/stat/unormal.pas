unit uNormal;

{$mode ObjFPC}{$H+}

interface

uses
  jpmTypes;
//  Classes, SysUtils, Math;

function phi(u: Float): Float;
function Normal(P: Float): Float;

implementation

{Standard Normal Probability Function}
function phi(u: Double): Float;
begin
  Result := (1.0/Sqrt(2.0 * PI)) * Exp(-0.5*u*u)
end;

{Inverse Standard Normal Function}
function Normal(P: Float): Float;
var
  chimax, chisqval, epsilon, minchisq, maxchisq: Float;
begin
  epsilon := 1e-6;
  chimax := 1e6;
  minchisq := 0.0;
  maxchisq := chimax;

  if P <= 0.0 then
  begin
    Result := maxchisq;
    exit;
  end;

  if P >= 1.0 then
  begin
    Result := 0.0;
    exit;
  end;

  chisqval := 0.5;
  while abs(maxchisq - minchisq) > epsilon do
  begin
    IF phi(chisqval) < P then
      maxchisq := chisqval
    else
      minchisq := chisqval;
    chisqval := (maxchisq + minchisq) * 0.5
  end;
  Result := chisqval;
end;

end.

