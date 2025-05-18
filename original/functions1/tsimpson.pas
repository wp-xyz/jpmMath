{*************************************************************
* Program to demonstrate integration of a real function F(x) *
*                    by Simpson's method                     *
* ---------------------------------------------------------- *
* REFERENCE: "Mathematiques en Turbo-Pascal (Part 1)         *
*             By Marc Ducamp and Alain Reverchon, Eyrolles,  *
*             Paris, 1991" [BIBLI 03].                       *
* ---------------------------------------------------------- *
* SAMPLE RUN:                                                *
* (Integrate sin(x) from x=0 to x=1)                         *
*                                                            *
* Integral of a function F(X) by Simpson's method            *
*                                                            *
* Input begin and end values for x variable:                 *
*                                                            *
*         X0 = 0                                             *
*         X1 = 1                                             *
*                                                            *
* Number of integration steps: 100                           *
*                                                            *
*                                                            *
* Value of integral:  4.59697694133457E-0001                 *
*                                                            *
*************************************************************}
Program Test_Romberg;
Uses WinCrt;

Var
        x0      : double;     {begin x value}
        x1      : double;     {end   x value}
        nstep   : integer;    {number of integration steps}
        result  : double;     {result of integral}


  {Given function to integrate}
  Function FUNC(x:double): double;
  Begin
    FUNC := SIN(x)
  End;


  {******************************************************
  * Integral of a function FUNC(X) by Simpson's method  *
  * --------------------------------------------------- *
  * INPUTS:                                             *
  *          a      begin value of x variable           *
  *          b      end value of x variable             *
  *          n      number of integration steps         *
  *                                                     *
  * OUTPUT:                                             *
  *          res    the integral of FUNC(X) from a to b *
  *                                                     *
  ******************************************************} 
  Procedure Integral_Simpson(a,b:double; n:Integer; VAR res:double);
  Var
       i:integer; step,r: double;
  Begin
    step:=(b-a)/2/n;
    r:=FUNC(a);
    res:=(r+FUNC(b))/2;
    for i:=1 to 2*n-1 do
    begin
      r:=FUNC(a+i*step);
      if Odd(i) then res:=res+r+r
      else res:=res+r
    end;
    res:=res*step*2/3
  End;

          
{main program}
BEGIN
  writeln;
  writeln(' Integral of a function F(X) by Simpson''s method');
  writeln;
  writeln(' Input begin and end values for x variable:');
  writeln;
  write('         X0 = '); read(x0);
  write('         X1 = '); read(x1);
  writeln;
  write(' Number of integration steps: '); read(nstep);

  Integral_Simpson(x0,x1,nstep,result);

  writeln;
  writeln;
  writeln(' Value of integral: ', result);
  writeln;
  ReadKey; DoneWinCrt

END.

{End of file tsimpson.pas}