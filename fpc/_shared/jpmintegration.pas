{ Integration routines extracted from Jean-Pierre Moreaus's programs }

unit jpmIntegration;

{$mode objFPC}{$H+}

interface

uses
  Classes, SysUtils,
  jpmTypes;

function GaussIntegral(AFunc: TFunction1; Order: Integer;
  x1, x2: Float): Float;
function GaussIntegral(AFunc: TFunction2; Order: Integer;
  x1, x2, y1, y2: Float): Float;
function GaussIntegral(AFunc: TFunction3; Order: Integer;
  x1, x2, y1, y2, z1, z2: Float): Float;

function RombergIntegral(AFunc: TFunction1; a, b, Prec: Float; var ObtPrec: Float;
  var Iterations: integer; MinIterations, MaxIterations: integer): Float;

implementation

{**************************************************************
* Integration by Gauss method of a real function F=F(X) or    *
* F=F(X,Y) or F=F(X,Y,Z). The integral is calculated by using *
* from 2 to 10 Gauss points.                                  *
* ----------------------------------------------------------- *
* INPUTS:                                                     *
*         AFunc  Function of 1, 2 or 3 real variables         *
*         Order  Number of Gauss points, order of polynomial  *
'        x1, x2  begin and end values of x variable           *
'        y1, y2  begin and end values of y variable           *
'        z1, z2  begin and end values of z variable           *
*                                                             *
* RETURNED VALUE: the integral of F(X) from x1 to x2 (or      *
*   F(x,y) from x1 to x2 and y1 to y2, or F(x, y, z) from     *
*   x1 to x2 and y1 to y2 and z1 to z2.                       *
* ----------------------------------------------------------- *
* Ref.: "Mécanique des vibrations linéaires By  M. Lalanne,   *
*       P. Berthier and J. Der Hagopian, Masson, Paris, 1980" *
*       [BIBLI 16].                                           *
* ----------------------------------------------------------- *
*                          TPW Version By J-P Moreau, Paris.  *
*                                          (www.jpmoreau.fr)  *
*                                                             *
* Extracted as stand-alone procedure by W. Pamler             *
**************************************************************}
procedure CalcGaussCoeffs(n: Integer; var A, H: TFloatArray);
var
  i, j: Integer;
begin
  SetLength(A, n);
  SetLength(H, n);
  case n of
    2: begin
         A[0] := -0.5773502691896257;  H[0] := 1.0;
         A[1] :=  0.5773502691896257;  H[1] := 1.0;
       end;
    3: begin
         A[0] := -0.7745966692414834;  H[0] := 0.5555555555555556;
         A[1] :=  0.0;                 H[1] := 0.8888888888888888;
         A[2] :=  0.7745966692414834;  H[2] := 0.5555555555555556;
       end;
    4: begin
         A[0] := -0.8611363115940526;  H[0] := 0.3478548451374538;
         A[1] := -0.3399810435848563;  H[1] := 0.6521451548625461;
         A[2] :=  0.3399810435848563;  H[2] := 0.6521451548625461;
         A[3] :=  0.8611363115940526;  H[3] := 0.3478548451374538;
       end;
    5: begin
         A[0] := -0.9061798459386640;  H[0] := 0.2369268850561891;
         A[1] := -0.5384693101056831;  H[1] := 0.4786286704993665;
         A[2] :=  0.0;                 H[2] := 0.5688888888888889;
         A[3] :=  0.5384693101056831;  H[3] := 0.4786286704993665;
         A[4] :=  0.9061798459386640;  H[4] := 0.2369268850561891;
       end;
    6: begin
         A[0] := -0.9324695142031521;  H[0] := 0.1713244923791704;
         A[1] := -0.6612093864662645;  H[1] := 0.3607615730481386;
         A[2] := -0.2386191860831969;  H[2] := 0.4679139345726910;
         A[3] :=  0.2386191860831969;  H[3] := 0.4679139345726910;
         A[4] :=  0.6612093864662645;  H[4] := 0.3607615730481386;
         A[5] :=  0.9324695142031521;  H[5] := 0.1713244923791704;
       end;
    7: begin
         A[0] := -0.9491079123427585;  H[0] := 0.1294849661688697;
         A[1] := -0.7415311855993945;  H[1] := 0.2797053914892766;
         A[2] := -0.4058451513773972;  H[2] := 0.3818300505051189;
         A[3] :=  0.0;                 H[3] := 0.4179591836734694;
         A[4] :=  0.4058451513773972;  H[4] := 0.3818300505051189;
         A[5] :=  0.7415311855993945;  H[5] := 0.2797053914892766;
         A[6] :=  0.9491079123427585;  H[6] := 0.1294849661688697;
       end;
    8: begin
         A[0] := -0.9602898564975363;  H[0] := 0.1012285362903763;
         A[1] := -0.7966664774136267;  H[1] := 0.2223810344533745;
         A[2] := -0.5255324099163290;  H[2] := 0.3137066458778873;
         A[3] := -0.1834346424956498;  H[3] := 0.3626837833783620;
         A[4] :=  0.1834346424956498;  H[4] := 0.3626837833783620;
         A[5] :=  0.5255324099163290;  H[5] := 0.3137066458778873;
         A[6] :=  0.7966664774136267;  H[6] := 0.2223810344533745;
         A[7] :=  0.9602898564975363;  H[7] := 0.1012285362903763;
       end;
    9: begin
         A[0] := -0.9681602395076261;  H[0] := 0.0812743883615744;
         A[1] := -0.8360311073266358;  H[1] := 0.1806481606948574;
         A[2] := -0.6133714327005904;  H[2] := 0.2606106964029354;
         A[3] := -0.3242534234038089;  H[3] := 0.3123470770400029;
         A[4] :=  0.0;                 H[4] := 0.3302393550012598;
         A[5] :=  0.3242534234038089;  H[5] := 0.3123470770400029;
         A[6] :=  0.6133714327005904;  H[6] := 0.2606106964029354;
         A[7] :=  0.8360311073266358;  H[7] := 0.1806481606948574;
         A[8] :=  0.9681602395076261;  H[8] := 0.0812743883615744;
       end;
   10: begin
         A[0] := -0.9739065285171717;  H[0] := 0.0666713443086881;
         A[1] := -0.8650633666889845;  H[1] := 0.1494513491505806;
         A[2] := -0.6794095682990244;  H[2] := 0.2190863625159820;
         A[3] := -0.4333953941292472;  H[3] := 0.2692667193099963;
         A[4] := -0.1488743389816312;  H[4] := 0.2955242247147529;
         A[5] :=  0.1488743389816312;  H[5] := 0.2955242247147529;
         A[6] :=  0.4333953941292472;  H[6] := 0.2692667193099963;
         A[7] :=  0.6794095682990244;  H[7] := 0.2190863625159820;
         A[8] :=  0.8650633666889845;  H[8] := 0.1494513491505806;
         A[9] :=  0.9739065285171717;  H[9] := 0.0666713443086881;
       end;
  end;
end;

function GaussIntegral(AFunc: TFunction1; Order:Integer;
  x1, x2: Float): Float;
var
  A: TFloatArray = nil;
  H: TFloatArray = nil;
  x: Float;
  a1, b1: Float;
  i: Integer;
begin
  Result := 0.0;
  CalcGaussCoeffs(Order, A, H);
  a1 := -(x1 - x2) / 2;
  b1 :=  (x1 + x2) / 2;
  for i := 0 to Order-1 do
  begin
    x := a1 * A[i] + b1;
    Result := Result + AFunc(x) * H[i] * a1;
  end;
end;

function GaussIntegral(AFunc: TFunction2; Order: Integer;
  x1, x2, y1, y2: Float): Float;
var
  A: Array of Float = nil;
  H: Array of Float = nil;
  x, y: Float;
  a1, a2, b1, b2: Float;
  i, j: Integer;
begin
  Result := 0.0;
  CalcGaussCoeffs(Order, A, H);
  a1 := -(x1 - x2) / 2;
  b1 :=  (x1 + x2) / 2;
  a2 := -(y1 - y2) / 2;
  b2 :=  (y1 + y2) / 2;
  for i := 0 to Order - 1 do
  begin
    x := a1 * A[i] + b1;
    for j := 0 to Order - 1 do
    begin
      y := a2 * A[j] + b2;
      Result := Result + AFunc(x,y) * H[i] * H[j] * a1 * a2;
    end;
  end;
end;

function GaussIntegral(AFunc: TFunction3; Order: Integer;
  x1, x2, y1, y2, z1, z2: Float): Float;
var
  A: Array of Float = nil;
  H: Array of Float = nil;
  x, y, z: Float;
  a1, a2, a3, b1, b2, b3: Float;
  i, j, k: Integer;
begin
  Result := 0.0;
  CalcGaussCoeffs(Order, A, H);
  a1 := -(x1 - x2) / 2;
  b1 :=  (x1 + x2) / 2;
  a2 := -(y1 - y2) / 2;
  b2 :=  (y1 + y2) / 2;
  a3 := -(z1 - z2) / 2;
  b3 :=  (z1 + z2) / 2;
  for i := 0 to Order-1 do
  begin
    x := a1 * A[i] + b1;
    for j := 0 to Order-1 do
    begin
      y := a2 * A[j] + b2;
      for k := 0 to Order-1 do
      begin
        z := a3 * A[k] + b3;
        Result := Result + AFunc(x,y,z) * H[i] * H[j] * H[k] * a1 * a2 * a3;
      end;
    end;
  end;
end;

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
  i, j: integer;
  pas, r, s, ta: Float;
  t: array of array of Float;
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

