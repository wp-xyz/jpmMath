unit jpmFunc;

{$mode ObjFPC}{$H+}

interface

uses
  jpmTypes;

function AkimaInterpolation(X, Y: Array of Float; StartIndex, EndIndex: Integer;
  xx: Float; var ErrorMsg: String): Float;

implementation

{*******************************************************
*          Akima spline fitting subroutine             *
* ---------------------------------------------------- *
* The input table is X(i), Y(i), where Y(i) is the     *
* dependent variable. The interpolation point is xx,   *
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
function AkimaInterpolation(X, Y: Array of Double; StartIndex, EndIndex: Integer;
  xx: Double; var ErrorMsg: String): Double;
var
  i, j, n: integer;
  a, aa, aaa, b, bb: Double;
  XM: array of Double = nil;
  Z: array of Double = nil;
begin
  Result := -999.0;
  ErrorMsg := '';

  // Both arrays should have the same length. If no, we use only the common parts.
  n := Length(X);
  if Length(Y) < n then
    n := Length(Y);

  // Akima interpolation requires one data point at the left of the interpolation range.
  if StartIndex < 0 then
    StartIndex := 0;

  // ... and three data points at the right
  if (EndIndex < 0) or (EndIndex >= n -3) then
    EndIndex := n - 4;

  // Not enough data points
  n := EndIndex - StartIndex + 1;
  if n < 3 then
  begin
    ErrorMsg := 'Not enough data points';
    exit;
  end;

  { Check to see if interpolation point is correct }
  if (xx < X[StartIndex]) or (xx > X[EndIndex]) then
  begin
    ErrorMsg := 'Interpolation point outside valid table range.';
    exit;
  end;

  { Calculate Akima coefficients, a and b }
  SetLength(XM, n+3);
  SetLength(Z, n);
  j := StartIndex;
  for i := 0 to n-1 do
  begin
    {Shift i to i+2}
    XM[i+2] := (Y[j+1]-Y[j]) / (X[j+1]-X[j]);
    inc(j);
  end;
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
      Z[i] := (XM[i+1] + XM[i+2]) / 2.0;
      break;
    end;
    Z[i] := (a*XM[i+1] + b*XM[i+2]) / (a + b);
  end;

  {Find relevant table interval}
  j := 0;
  while (xx > X[j]) do
    inc(j);
  dec(j);
  if j = -1 then j := 0;
  i := j - StartIndex;

  {Begin interpolation}
  b := X[j+1] - X[j];
  a := xx - X[j];
  aa := a*a;
  aaa := aa*a;
  bb := b*b;
  Result := Y[j] +
            Z[i]*a +
           (3.0*XM[i+2] - 2.0*Z[i] - Z[i+1]) * aa/b +
           (Z[i] + Z[i+1] - 2.0*XM[i+2]) * aaa/bb;
end;

end.

