{*************************************************************
* Program to demonstrate integration of a real function F(x) *
*                    by Romberg's method                     *
* ---------------------------------------------------------- *
* REFERENCE: "Mathematiques en Turbo-Pascal (Part 1) by      *
*             Marc Ducamp and Alain Reverchon, Eyrolles,     *
*             Paris, 1991" [BIBLI 03].                       *
* ---------------------------------------------------------- *
* SAMPLE RUN:                                                *
* (Integrate sin(x) from x=0 to x=1)                         *
*                                                            *
* Integral of a function F(X) by Romberg's method            *
*                                                            *
* Input begin and end values for x variable:                 *
*                                                            *
*         X0 = 0                                             *
*         X1 = 1                                             *
*                                                            *
* Input desired precision: 1e-10                             *
*                                                            *
*                                                            *
* Value of integral   :  4.59697694131851E-0001              *
*                                                            *
* Obtained precision  :  9.59908263986620E-0011              *
*                                                            *
* Number of iterations:  4                                   *
*                                                            *
*************************************************************}
program test_romberg_Laz;

{$mode objfpc}

uses jpmIntegration;

type
  TFunction1 = function(x: double): double;

{Given function to integrate}
function FUNC(x: double): double;
begin
  Result := sin(x)
end;

  {******************************************************
  * Integral of a function FUNC(X) by Romberg's method *
  * --------------------------------------------------- *
  * INPUTS:                                             *
  *              a  begin value of x variable           *
  *              b  end value of x variable             *
  *           Prec  desired precision                   *
  *  IterationsMin  minimum number of iterations        *
  *  IterationsMax  maximum number of iterations        *
  *                                                     *
  * OUTPUTS:                                            *
  *        ObtPrec  obtained precision for integral     *
  *     Iterations  number of iterations done           *
  *                                                     *
  * RETURNED VALUE  the integral of FUNC(X) from a to b *
  *                                                     *
  ******************************************************} 
  function RombergIntegral(AFunc: TFunction1; a, b, Prec: double; var ObtPrec: double;
    var Iterations: integer; IterationsMin, IterationsMax: integer): double;
  const
    MAXITER = 15;
  var
    i,j        : integer;
    pas,r,s,ta : double;
    t          : array of array of Double;
  begin
    if IterationsMax > MAXITER then
      IterationsMax := MAXITER;
    r := AFunc(a);
    ta := (r + AFunc(b)) / 2;
    Iterations := 0;
    pas := b - a;
    SetLength(t, IterationsMax+1, IterationsMax+1);
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
    until (Iterations >= IterationsMax) or ((ObtPrec < Prec) and (Iterations >= IterationsMin));
    Result := t[Iterations, 0]
  end;

{main program}
var
  s: String;
  res: Integer;
  x0      : double;     {begin x value}
  x1      : double;     {end   x value}
  prec    : double;     {desired precision}
  integral: double;     {result of integral}
  obtprec : double;     {obtained precision}
  niter   : integer;    {number of actual iterations}

begin
  WriteLn;
  WriteLn(' Integral of a function F(X) by Romberg''s method');
  WriteLn;
  WriteLn(' Input begin and end values for variable x:');
  WriteLn;

  Write('         X0 = [0.0] ');
  ReadLn(s);
  if s = '' then x0 := 0.0 else val(s, x0, res);
  if res <> 0 then
  begin
    WriteLn('No valid number entered. Terminating...');
    Halt;
  end;

  Write('         X1 = [1.0] ');
  ReadLn(s);
  if s = '' then x1 := 1.0 else val(s, x1, res);
  if res <> 0 then
  begin
    WriteLn('No valid number entered. Terminating...');
    Halt;
  end;

  Write(' Input desired precision: [1e-10] ');
  ReadLn(s);
  if s = '' then prec := 1e-10 else val(s, prec, res);
  if (res <> 0) or (prec <= 0) then
  begin
    WriteLn('Input not valid. Terminating...');
    Halt;
  end;

  integral := RombergIntegral(@func, x0, x1, prec, obtprec, niter, 1, 50);

  WriteLn;
  WriteLn(' Value of integral   : ', integral);
  WriteLn(' Obtained precision  : ', obtprec);
  WriteLn(' Number of iterations:  ', niter);
  WriteLn;

  Write('Press ENTER to quit...');
  ReadLn;
end.

