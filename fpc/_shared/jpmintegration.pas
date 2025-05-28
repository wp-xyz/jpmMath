{ Integration routines extracted from Jean-Pierre Moreaus's programs }

unit jpmIntegration;

{$mode objFPC}{$H+}

interface

uses
  Classes, SysUtils,
  jpmTypes;

function RombergIntegral(AFunc: TFunction1; a, b, Prec: Float; var ObtPrec: Float;
  var Iterations: integer; MinIterations, MaxIterations: integer): Float;

implementation

{*******************************************************
* Integral of a function AFunc(X) by Romberg's method  *
* ---------------------------------------------------- *
* INPUTS:                                              *
*              a  begin value of x variable            *
*              b  end value of x variable              *
*           Prec  desired precision                    *
*  MinIterations  minimum number of iterations         *
*  MaxIterations  maximum number of iterations         *
*                                                      *
* OUTPUTS:                                             *
*        ObtPrec  obtained precision for integral      *
*     Iterations  number of iterations done            *
*                                                      *
* RETURNED VALUE  the integral of AFunc(X) from a to b *
*                                                      *
*******************************************************}
function RombergIntegral(AFunc: TFunction1; a, b, Prec: Float; var ObtPrec: Float;
  var Iterations: integer; MinIterations, MaxIterations: integer): Float;
const
  MAXITER = 15;
var
  i,j        : integer;
  pas,r,s,ta : Float;
  t          : array of array of Float;
begin
  if MaxIterations > MAXITER then
    MaxIterations := MAXITER;
  r := AFunc(a);
  ta := (r + AFunc(b)) / 2;
  Iterations := 0;
  pas := b - a;
  SetLength(t, MaxIterations+1, MaxIterations+1);
  t[0,0] := ta*pas;
  repeat
    Inc(Iterations);
    pas := pas/2;
    s := ta;
    for i := 1 to pred(1 shl Iterations) do   {2^n-1}
      s := s + AFunc(a + pas*i);
    t[0, Iterations] := s*pas;
    r := 1;
    for i := 1 to Iterations do
    begin
      r := r*4;
      j := Iterations-i;
      t[i,j] := (r*t[i-1, j+1] - t[i-1, j]) / (r-1)
    end;
    obtprec := abs(t[Iterations, 0] - t[Iterations-1, 0])
  until (Iterations >= MaxIterations) or ((ObtPrec < Prec) and (Iterations >= MinIterations));
  Result := t[Iterations, 0]
end;

end.

