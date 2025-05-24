{-------------------------------------------------------------------------------
                                 Unit jpmStats                        
         Statistical function based on Jean-Pierre Moraus's code. 

This unit collects generally available statistical calculations from
Jean-Pierre Moraus's sample projects.
-------------------------------------------------------------------------------}

unit jpmStats;

{$mode ObjFPC}{$H+}

interface

uses
  SysUtils, jpmTypes;
//  Classes, SysUtils, Math;

function NormalDist(u: Float): Float;
function InvNormalDist(P: Float): Float;

implementation

{ Standard Normal Probability Function (Gaussian curve)
  Taken from the "normal.pas" sample project where it is named "Phi()". }
function NormalDist(u: Float): Float;
begin
  Result := 1.0/Sqrt(2.0*pi) * exp(-0.5*u*u)
end;

{ Inverse Standard Normal Function 
  Taken from the "normal.pas" sample project where it is named "Normal()". }
function InvNormalDist(P: Float): Float;
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
    IF NormalDist(chisqval) < P then
      maxchisq := chisqval
    else
      minchisq := chisqval;
    chisqval := (maxchisq + minchisq) * 0.5
  end;
  Result := chisqval;
end;

(* Alternative (direct) method of calculation
function InvNormalDist(P: Float): Float;
var
  a: Double;
begin
  if P <= 0 then
    Result := NaN
  else
  if P > 37 then  // to be more exact: replace this by InvNormalDist(MIN_FLOAT)
    Result := 0
  else
  begin
    a := -2.0 * (ln(sqrt(2.0*pi)) + ln(P));
    if a >= 0 then
      Result := sqrt(a)
    else
      Result := 0.0; //NaN;
  end;
end;               *)

end.

