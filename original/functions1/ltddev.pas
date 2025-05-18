{********************************************************
*  Calculate a limited development of a real function   *
*  f(x) at point x0 with step h (to order 5).           *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
*                                                       *
* Limited development of f(x) at x0:                    *
*                                                       *
* Function to develop: 1/(1-x)                          *
*                                                       *
* At point x0= 0                                        *
* Value of step h= 0.1                                  *
*                                                       *
* a1 =  1.0000416250                                    *
* a2 =  1.0000416250                                    *
* a3 =  1.0011207928                                    *
* a4 =  1.0011207928                                    *
* a5 =  1.0136468470                                    *
*                                                       *
* (The correct answer is a1=a2=a3=a4=a5=1)              *         
*                                                       *
* Function to develop: Ln(x)                            *
*                                                       *
* At point x0= 1                                        *
* Value of step h= 0.1                                  *
*                                                       *
* a1 =  1.0000057552                                    *
* a2 = -0.5000050520                                    *
* a3 =  0.3334508142                                    *
* a4 = -0.2501062334                                    *
* a5 =  0.2011301593                                    *
*                                                       *
* (the correct answer is a1=1, a2=-0.5, a3=1/3,         *
*  a4=-1/4 and a5=1/5)                                  *   
* ----------------------------------------------------- *
* Ref.: "Mathematiques en Turbo Pascal By Alain         *
*        Reverchon and Marc Ducamp, Editions Eyrolles,  *
*        Paris, 1991" [BIBLI 03].                       *
*                                                       *
*                            TPW version by J-P Moreau. *
*                                (www.jpmoreau.fr)      *
*********************************************************
Note: as you can see this direct method is not suitable to
      calculate a limited development in an accurate way.
      The next programs (ltddev1 to ltddev3) will allow us
      to calculate limited developments at x=0 up to a high
      order with a good accuracy for functions f*g, f/g or
      gof, knowing the limited developments of f(x) and g(x)
      at x=0.
-----------------------------------------------------------}
Program Ltddev;
Uses WinCrt;

Var
    d,h,r,s,x0: Double;
    n: Integer;


{sample function f(x)=1/(1-x) }
Function f(x:Double):Double;
begin
  if abs(x-1) >1e-12 then
    f := 1.0/(1.0-x)
  else
    f := 1e12
end;

{sample function f(x)=Ln(x)
Function f(x:Double):Double;
begin
  if x >1e-12 then
    f := Ln(x)
  else
    f := -1e12
end;                      }

{returns the value of nth derivative of function
 f(x) at point x0 with increment step h (n=1 to 5).
 Note that for polynomial functions h must not be
 too small. Here h=1 gives best results than h=0.01}
Function ar_nderiv(x:double;n:integer;h:double):Double;
Const table:Array[1..5,0..6] of Integer =
    ((3, 45, -9, 1, 0, 0, 60),
     (3, 270, -27, 2, 0, 0, 180),
     (4, -488, 338, -72, 7, 0, 240),
     (4, -1952, 676, -96, 7, 0, 240),
     (5, 1938, -1872, 783, -152, 13, 288));
Var
  i:Integer;
  fa,fb,fc,t:Double;
Begin
  if(n>5) then
  begin
    ar_nderiv:=0;
    exit
  end;
  fa := f(x);
  t:=0;
  for i:=1 to table[n,0] do
  begin
    fb := f(x-i*h);
    fc := f(x+i*h);
    if Odd(n) then
      t := t + table[n,i]*(fc-fb)
    else
      t := t + table[n,i]*((fc-fa) + (fb-fa))
  end;
  t := t / table[n,6];
  for i:=1 to n do  t := t / h;
  ar_nderiv:= t
End;

{main program}
BEGIN
  writeln;
  writeln(' Limited development of a real function f(x):');
  writeln;
  writeln(' Function to develop: 1/(1-x)');
  writeln;
  write(' At point x0= '); readln(x0);
  write(' Value of step h= '); readln(h);
  writeln;

  s:=1;

  for n:=1 to 5 do
  begin
    s := s * n;
    r := ar_nderiv(x0,n,h) / s;
    writeln(' a',n,' = ',r:13:10)
  end;
  writeln;

  Readkey; DoneWinCrt
END.

{end of file ltddev.pas}