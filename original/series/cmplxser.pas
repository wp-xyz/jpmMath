{***************************************************
* Program to demonstrate complex series evaluation *
* It is assumed that the coefficients are obtained *
* from a subroutine.                               *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].*
*                                                  *
*          Pascal Version by J.-P. Moreau, Paris.  *
*                    (www.jpmoreau.fr)             *
* ------------------------------------------------ *
* Reference: BASIC Scientific Subroutines, Vol. II *
* by F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981.    *
* ------------------------------------------------ *
* SAMPLE RUN:                                      *
*                                                  *
* Input the complex number as prompted:            *
*                                                  *
*   Real part    = 1                               *
*   Complex part = 0                               *
*                                                  *
* Results are:                                     *
*                                                  *
*   Z1 =  32.0000                                  *
*   Z2 =   0.0000                                  *
*                                                  *
***************************************************}
PROGRAM CMPLXSER;
Uses WinCrt;

VAR
        A : Array[0..5] of DOUBLE;
        u,u1,v,v1,x,y,z1,z2 : DOUBLE;
        m,n : INTEGER;

{Coefficients subroutine}
PROCEDURE Coeff;
Begin
  m:=5;
  A[0] := 1;
  A[1] := 5;
  A[2] := 10;
  A[3] := 10;
  A[4] := 5;
  A[5] := 1
End;

Procedure ATAN(Numerateur,denominateur:DOUBLE;
Var Phase:DOUBLE);
{Returns a phase value between 0 and 2*PI}
BEGIN
if (Abs(Denominateur)<1.0E-6) then
   begin
   if (numerateur < 0) then Phase:=-PI/2
      else Phase:=PI/2;
   exit;
   end
else begin
     Phase:=Arctan(numerateur/denominateur);
     if(denominateur > 0 ) then
         begin
         if (numerateur >=0 ) then exit
             else begin
                  Phase:=phase+2*PI;
                  exit;
                  end;
         end
         else Phase:=phase+PI;
     end;       {else}
END;            {ATAN}

Function Power(x:DOUBLE; n:integer): DOUBLE;
{ Calculates x power n}
var i,m : integer;
    result :DOUBLE;
begin
result := 1.0;
if n=0 then
   begin
   Power:=result;
   exit;
   end;
m:=  n;
if n<0 then m:=-n;
for i:=1 to m do result :=x*result;
Power :=result;
if n<0 then Power:=1.0/result;
end;

{Rectangular to polar conversion}
PROCEDURE RectPol;
Begin
  u := sqrt(x*x+y*y);
  {Guard against ambiguous vector}
  if y=0 then y := 1e-30;
  {Guard against divide by zero}
  if x=0 then x := 1e-30;
  ATAN(y,x,v)  {see math.pas}
End;

{Polar to rectangular conversion}
PROCEDURE PolRect;
Begin
  x := u*cos(v);
  y := u*sin(v)
End;

{Polar power}
PROCEDURE PolPower;
Begin
  u1 := Power(u,n);  
  v1 := n * v;
  v1 := v1 - 2*pi * int(v1/2.0/pi)
End;

{Rectangular complex number power}
PROCEDURE RectPower;
Begin
  {Rectangular to polar conversion}
  RectPol;
  {Polar power}
  PolPower;
  {Change variable for conversion}
  u:=u1 ; v:=v1;
  {Polar to rectangular conversion}
  PolRect
End;

{*******************************************************
*        Complex series evaluation subroutine          *
* ---------------------------------------------------- *
* The series coefficients are A(i), assumed real, the  *
* order of the polynomial is m. The subroutine uses    *
* repeated calls to the nth power (Z^N) complex number *
* subroutine. Inputs to the subroutine are x, y, m and *
* the A(i). Outputs are z1 (real) and z2 (imaginary).  *
*******************************************************}
PROCEDURE Complex_Series;
Var a1, a2 : DOUBLE;
Begin
  z1 := A[0] ; z2 := 0.0;
  {Store x and y}
  a1 := x ; a2 := y;
  for n := 1 to m do
  begin
    {Recall original x and y}
    x := a1 ; y := a2;
    {Call Z^N subroutine}
    RectPower;
    {Form partial sum}
    z1 := z1 + A[n] * x;
    z2 := z2 + A[n] * y
  end;
  {Restore x and y}
  x := a1 ; y := a2
End;


{main^program}
BEGIN
  {Get coefficients}
  Coeff;
  {input complex number}
  clrscr;
  writeln;
  writeln(' Input the complex number as prompted:');
  writeln;
  write('   Real part    = '); read(x);
  write('   Complex part = '); read(y);

  Complex_Series;

  writeln;
  writeln(' Results are:');
  writeln;
  writeln('   Z1 = ', z1:8:4);
  writeln('   Z2 = ', z2:8:4);
  ReadKey; DoneWincrt
END.

{End of file cmplxser.pas}