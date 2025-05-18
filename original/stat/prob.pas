{************************************************************
* Evaluate probabilities by parsing a string expression     *
*        and performing the implied calculations            *
* --------------------------------------------------------- *
* SAMPLE RUN:                                               *
*                                                           *
* Type expression to be evaluated (no syntax checking).     *
* Example: 1-c(10,0,0.2)-c(10,1,0.2)                        *
*                                                           *
* 'c(10,0,.1)+c(10,1,0.1)+c(10,2,0.1)'                      *
*                                                           *
* C(10,0) = 1                                               *
* C(10,1) = 10                                              *
* C(10,2) = 45                                              *
* probability =    0.9298                                   *
*                                                           *
* --------------------------------------------------------- *
* Ref.: "Problem Solving with Fortran 90 By David R.Brooks, *
*        Springer-Verlag New York, 1997".                   *
*                                                           *
*                        TPW Release By J-P Moreau, Paris.  *
*                               (www.jpmoreau.fr)           *
*************************************************************
!Explanations:
!------------
!   A manufacturer's experience has shown that 10% of all integrated
! circuits (ICs) will be defective. Because of this high failure rate,
! a quality control engineer monitors the manufacturing process by
! testing random samples every day. What is the probability that:
!   (a) exactly two ICs in a sample of 10 will be defective?
!   (b) at least two will be defective?
!   (c) no more than two will be defective?
! Answers:
! (a) the probability that a particular sample of 10 will contain
!                                     2     8
! exactly two defective ICs is:  (0.1) (0.9) = 0.004305. However,
! there are C(10,2)=45 possible combinations of two defective and eight
! good ICs. From probability theory, the number of combinations of n
! things taken k at a time is:  C(n,k) = n!/[k!(n-k)!]  where
! ! indicates the factorial function (n! = 1 x 2 x3 x ....n). Therefore,
! the probability that a sample of 10 will contain exactly two defects
! is:                          2     8
!          P(=2) = C(10,2)(0.1) (0.9)  = 45 x 0.004305 = 0.1937 
!
! (b) the probability of finding at least two defective ICs is equal to
!  1 minus the probability of 0 defective IC minus the probability of
!  1 defective IC:              0     10              1     9
!      P(>=2) = 1 - C(10,0)(0.1) (0.9)  - C(10,0)(0.1) (0.9) = 0.2639
!
! (Reemeber that 0! = 1 by definition.)
!
! (c) the probability of finding no more than two defective ICs is:
!                           0     10              1     9
!      P(<=2) = C(10,0)(0.1) (0.9)  + C(10,1)(0.1) (0.9)
!                             2     8
!               + C(10,2)(0.1) (0.9) = 0.9298   
!
! For example, for part (b) of the problem, the user will type the
! string:  '1-c(10,0,.1)-c(10,1,.1)'.
!----------------------------------------------------------------------}
Program Prob;

Uses WinCrt;

Var
    a: String;     {string to be calculated}
    s,s1: String[20];
    error, i, ierror, j, k, len, n, sign, left: Integer;
    probability, prob_a: Real;

Function Fact(x:real):Real;
{ Calculate x! }
Var
    prod: Real;
    i,ix: Integer;
Begin
  prod := 1.0;
  ix:=Round(x);
  For i:=2 to ix do prod := prod * i;
  Fact := prod
End;

Function C(n,k:integer): Integer;
{ Calculate combinations of n things taken k at a time }
Var temp, denom: real;
Begin
  denom := Fact(1.0*k)*Fact(1.0*(n-k));
  temp := Fact(n) / denom;
  C:=Round(temp)
End;

Function Power(x:real; n:integer): Real;
{calculates x power n}
var result : real;
    i: integer;
begin
  result := 1.0;
  if n=0 then
  begin
    Power:=result;
    exit
  end
  else
    for i:=1 to n do
      result := x * result;
  Power :=result
end;

{main program}
Begin
  Writeln;
  Writeln(' Type expression to be evaluated (no syntax checking).');
  Writeln(' Example: 1-c(10,0,0.2)-c(10,1,0.2)');
  writeln;
  write(' '); readln(a); 

  writeln;

  len := Length(a);
  probability := 0.0;
  if a[1]='1' then probability := 1.0;
  sign:=1;    {a leading + sign is optional}
  For i:=1 to len do
  begin
    if a[i]='+' then sign:=1;
    if a[i]='-' then sign:=-1;
    if a[i]='(' then left:=i;
    if a[i]=')' then
    begin
      {Read(a(left+1:i-1),*) n, k, prob_a}
      s:='';
      For j:=left+1 to i do s:=s+a[j];
      s1:=''; j:=1;
      while (s[j]<>',') and (s[j]<>')') do
      begin s1:=s1+s[j]; Inc(j) end;
      Val(s1,n,ierror); error:=error+ierror;
      s1:=''; Inc(j);
      while (s[j]<>',') and (s[j]<>')') do
      begin s1:=s1+s[j]; Inc(j) end;
      Val(s1,k,ierror); error:=error+ierror;
      s1:=''; Inc(j);
      while (s[j]<>',') and (s[j]<>')') do
      begin s1:=s1+s[j]; Inc(j) end;
      Val(s1,prob_a,ierror); error:=error+ierror;
      if ierror<>0 then writeln(' Error in read n,k,prob.');
      probability:=probability + sign*C(n,k)*Power(prob_a,k)*Power(1.0-prob_a,n-k);
      writeln(' C(',n,',',k,') = ',C(n,k)) 
    end
  end;
  writeln(' probability = ', probability:10:4);	  

  ReadKey;
  DoneWinCrt

End.

{end of file prob.pas}