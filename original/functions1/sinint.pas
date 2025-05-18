{*************************************************************
* Function sinintegral(x) and cosintegral(x) by Simpson's    *
* method                                                     *
* ---------------------------------------------------------- *
*   REFERENCE:  (for Simpson's method)                       * 
*              "Mathematiques en Turbo-Pascal (Part 1) by    *
*               Marc Ducamp and Alain Reverchon, Eyrolles,   *
*               Paris, 1987" [BIBLI 03].                     *
* ---------------------------------------------------------- *
*                                                            *
*                                TPW Version By J-P Moreau.  *
*                                    (www.jpmoreau.fr)       *
*************************************************************}
Unit Sinint;

Interface

  Const gamma = 0.57721566490153286;   {Euler's constant}

  Function sinintegral(x:double): double;
  Function cosintegral(x:double): double;

Implementation

{Given function to integrate (2 kinds) }
Function FUNC(kind:integer; t: double): Double;
Begin
  Case kind of
    1: if abs(t)<1e-10 then
         FUNC := 1.0
       else
         FUNC := sin(t)/t;      {for sinintegral}
    2: if abs(t)<1e-10 then
         FUNC := 0.0
       else
         FUNC := (cos(t)-1.0)/t {for cosintegral}
  end
End;

{******************************************************
* Integral of a function FUNC(X) by Simpson's method  *
* --------------------------------------------------- *
* INPUTS:                                             *
*          kind   =1 for sinintegral                  *
*                 =2 for cosintegral                  * 
*          a      begin value of x variable           *
*          b      end value of x variable             *
*          n      number of integration steps         *
*                                                     *
* OUTPUT:                                             *
*          res    the integral of FUNC(X) from a to b *
*                                                     *
******************************************************}
Procedure Integral_Simpson(kind:integer; a:double; b:double; n:longint; var res:double);
Var
  i: longint; step,r: double;
Begin
  step:=(b-a)/2/n;
  r:=FUNC(kind,a);
  res:=(r+FUNC(kind,b))/2;
  for i:=1 to 2*n-1 do
  begin
    r:=FUNC(kind,a+i*step);
    if (i Mod 2) <> 0 then res := res + r+r
    else res := res + r
  end;
  res := res * step*2/3
End;

{******************************************************
*     Si(x) = Integral of sin(u)/u from u=0 to u=x    *
******************************************************}          
Function sinintegral(x: double): double;
Var
   x0:double;      {begin x value
   x                end x value}
   nstep:longint;  {number of integration steps}
   result:double;  {result of integral}
Begin

  x0:=0.0;
  nstep := Round(200*abs(x));   {this should ensure about 14 exact digits}
  
  Integral_Simpson(1,x0,x,nstep,result);     {kind=1}

  sinintegral := result

End;

{******************************************************
* Calcualte Ci(x)= gamma + Ln(x) + I(x) with I(x) =   *
* Integral of (cos(u)-1)/u from u=0 to u=x (x > 0).   *
* Gamma = 0.57721566... (Euler's constant).           *
******************************************************}          
Function cosintegral(x: double): double;
Var
   x0:double;      {begin x value
   x                end x value}
   nstep:longint;  {number of integration steps}
   result:double;  {result of integral}
Begin

  x0:=0.0;
  nstep := Round(200*abs(x));   {this should ensure about 14 exact digits}
  
  Integral_Simpson(2,x0,x,nstep,result);     {kind=2}

  cosintegral := gamma + Ln(x) + result

End;

End.

{End of file sinint.pas}