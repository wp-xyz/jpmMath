{********************************************************
*  Estimate the nth derivative of a real function f(x)  *
*  at point x with step h (n from 1 to 5).              *
* ----------------------------------------------------- *
* SAMPLE RUN:                                           *
*                                                       *
* Nth derivative of a real function f(x):               *
*                                                       *
* Function to derivate: x^3                             *
*                                                       *
* At point x0= 10                                       *
* Order of derivative (1 to 5): 1                       *
* Value of step h= 0.1                                  *
*                                                       *
* f'(x0)=300.000000                                     *  
*                                                       *
* At point x0= 10                                       *
* Order of derivative (1 to 5): 2                       *
* Value of step h= 0.1                                  *
*                                                       *
* f"(x0)= 60.000000                                     *
*                                                       *
* ----------------------------------------------------- *
* Ref.: "Analyse numérique en C By Alain Reverchon and  *
*        Marc Ducamp, Armand Colin Editeur, Paris,      *
*        1990" [BIBLI 14].                              *
*                                                       *
*                         Pascal version By J-P Moreau. *
*                              (www.jpmoreau.fr)        *
********************************************************}
Program Nderiv;
Uses WinCrt;

Var
    d,h,x: Double;
    n: Integer;


{sample function f(x)=x^3}
Function f(x:Double):Double;
begin
  f := x*x*x
end;

{returns the value of nth derivative of function
 f(x) at point x0 with increment step h
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
  writeln(' Nth derivative of a real function f(x):');
  writeln;
  writeln(' Function to derivate: x^3');
  writeln;
  write(' At point x0= '); readln(x);
  write(' Order of derivative (1 to 5): '); readln(n);
  write(' Value of step h= '); readln(h);
  writeln;

  d := ar_nderiv(x,n,h);

  Case n of
    1: writeln(' f''(x0)=',d:10:6);
    2: writeln(' f"(x0)=',d:10:6);
    3: writeln(' f"''(x0)=',d:10:6);
    4: writeln(' f""(x0)=',d:10:6);
    5: writeln(' f""''(x0)=',d:10:6)
  End;
  if (n<1) or (n>5) then
    writeln(' n must be between 1 and 5 !');
  Readkey; DoneWinCrt
END.

{end of file nderiv.pas}