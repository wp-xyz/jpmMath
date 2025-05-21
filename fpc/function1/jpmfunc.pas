unit jpmFunc;

{$mode ObjFPC}{$H+}

interface

function AkimaInterpolation(X, Y: Array of Double; xx: Double;
  var yy: Double): Boolean;

implementation

{*******************************************************
*          Akima spline fitting subroutine             *
* ---------------------------------------------------- *
* The input table is X(i), Y(i), where Y(i) is the     *
* dependant variable. The interpolation point is xx,   *
* which is assumed to be in the interval of the table  *
* with at least one table value to the left, and three *
* to the right. The interpolated returned value is yy. *
* The function result is false when an error occored.  *
* It is also assumed that the X(i) are in ascending    *
* order.                                               *
* ---------------------------------------------------- *
*   Reference: BASIC Scientific Subroutines, Vol. II   *
*   By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
*                                                      *
*         Pascal Version by J.-P. Moreau, Paris.       *
*                   (www.jpmoreau.fr)                  *
*******************************************************}
function AkimaInterpolation(X, Y: Array of Double; xx: Double;
  var yy: Double): Boolean;
var
  i, n: integer;
  a, aa, aaa, b, bb: Double;
  XM: array of Double = nil;
  Z: array of Double = nil;
begin
  Result := true;

  n := Length(X);
  if Length(Y) < n then
    n := Length(Y);

  (*              wp: why this?
  { Special case xx=0 }
  if xx = 0.0 then
  begin
    yy := 0.0;
    exit;
  end;
  *)

  { Check to see if interpolation point is correct }
  if (xx < X[0]) or (xx >= X[n-4]) then
  begin
    Result := false;
    exit;
  end;

  { Calculate Akima coefficients, a and b }
  SetLength(XM, n+3);
  SetLength(Z, n);
  for i := 0 to n-1 do
    {Shift i to i+2}
    XM[i+2] := (Y[i+1]-Y[i]) / (X[i+1]-X[i]);
  XM[n+1] := 2.0*XM[n] - XM[n-1];
  XM[n+2] := 2.0*XM[n+1] - XM[n];
  XM[1] := 2.0*XM[2] - XM[3];
  XM[0] := 2.0*XM[1] - XM[2];

  for i := 0 to n-1 do
  begin
    a := abs(XM[i+3] - XM[i+2]);
    b := abs(XM[i+1] - XM[i]);
    if a + b = 0 then
    begin
      Z[i] := (XM[i+2] + XM[i+1]) / 2.0;
      break;
    end;
    Z[i] := (a*XM[i+1] + b*XM[i+2]) / (a + b);
  end;

  {Find relevant table interval}
  i := 0;
  while (xx > X[i]) do
    inc(i);
  dec(i);
  if i = -1 then i := 0;

  {Begin interpolation}
  b := X[i+1] - X[i];
  a := xx - X[i];

  aa := a*a;
  aaa := aa*a;
  bb := b*b;
  yy := Y[i] + Z[i]*a + (3.0*XM[i+2]-2.0*Z[i]-Z[i+1])*aa/b;
  yy := yy + (Z[i]+Z[i+1]-2.0*XM[i+2])*aaa/(bb);
end;

end.

