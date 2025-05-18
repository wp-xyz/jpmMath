{*****************************************
*  Unit for basic statistical Functions  *
*****************************************}
Unit Stats;

Interface

Uses WinCrt1;

Const SIZE = 500;

Type
     pVec = ^VEC;
     VEC = Array[1..SIZE] of Real;


     Procedure NormalStats(a:pVec; n:Integer; flag:Char; Var avg,
                           std_dev: Real);

     Procedure NormalArray(a:pVec; n:Integer);

     Procedure LinearReg(x,y:pVec; n:Integer; flag:Char;
                         Var a,b,s_yx,r: Real);


Implementation

  Function Dot_Product(n:Integer; x,y:pVec): Real;
  Var i:Integer;
      temp:Real;
  Begin
    temp:=0.0;
    For i:=1 to n do temp:=temp + x^[i]*y^[i];
    Dot_Product := temp
  End;

  Function Sum(n:Integer; x:pVec): Real;
  Var i:Integer;
      temp:Real;
  Begin
    temp:=0.0;
    For i:=1 to n do temp:=temp + x^[i];
    Sum := temp
  End;

{---------------------------------------------------------
! Basic description statistics for normally distributed
! data. Calculates either sample or population statistics.
! Sets std_dev to -1 if error condition is detected.
!--------------------------------------------------------}
  Procedure NormalStats(a:pVec; n:Integer; flag:Char; Var avg,
                        std_dev: Real);
  Var
      sum_, sum_sq, variance: Real;
      i: Integer;
  Begin
    sum_ := Sum(n,a);
    sum_sq := DOT_PRODUCT(n,a,a);

    if (flag='p') or (flag='P') then
      variance := (sum_sq-Sqr(sum_)/(1.0*n))/(1.0*n)
    else if (flag='s') or (flag='S') then
      variance := (sum_sq-Sqr(sum_)/(1.0*(n-1)))/(1.0*(n-1))
    else
    begin
      writeln(' From NormalStats: Flag Error, <P> assumed.');
      variance := (sum_sq-Sqr(sum_)/(1.0*n))/(1.0*n)
    end;
	
    If variance < 0.0 Then  {an error exists}
    begin
      writeln(' From NormalStats: negative variance ', variance);
      std_dev := -1.0
    end
    Else
      std_dev := SQRT(variance);

    avg := sum_/n

  End; {NormalStats}
	  	 
{----------------------------------------------------------
! For data to be represented by y=ax+b, calculates linear 
! regression coefficients, sample standard error of y on x,
! and sample correlation coefficients. Sets r=0 if an error
! exists. If the intercept coefficient a is set to 0 on
! input, the regression is forced through (0,0).
!---------------------------------------------------------}
  Procedure LinearReg(x,y:pVec; n:Integer; flag:Char;
                     Var a,b,s_yx,r: Real);
  Var
      avg, std_dev: Real;
      sum_x,sum_y,sum_xy,sum_xx,sum_yy,temp: Real;
  Begin
    sum_x := SUM(n,x);
    sum_y := SUM(n,y);
    sum_xy := DOT_PRODUCT(n,x,y);
    sum_xx := DOT_PRODUCT(n,x,x);
    sum_yy := DOT_PRODUCT(n,y,y);

    If a <> 0.0 Then  {calculate full expression}
    begin
      temp := n*sum_xx - Sqr(sum_x);
      a := (sum_y*sum_xx - sum_x*sum_xy)/temp;
      b := (n*sum_xy - sum_x*sum_y)/temp;
      s_yx := SQRT((sum_yy - a*sum_y - b*sum_xy)/n)
    end
    Else   {just calculate slope}
    begin
      b := sum_y/sum_x;
      s_yx := SQRT((sum_yy - 2.0*b*sum_xy + b*b*sum_xx)/n)
    end;

    if (flag='s') or (flag='S') then
      s_yx := s_yx * SQRT((1.0*n)/(1.0*(n-2)));
	
{ Use NormalStats to get standard deviation of y }
    NormalStats(y,n,flag,avg,std_dev);
	
    If std_dev > 0.0 Then
    begin
      temp := 1.0 - Sqr(s_yx/std_dev);
      If temp > 0.0 Then
        r := SQRT(temp)
      Else  {an error exists}
      begin
        r := 0.0;
        writeln(' From LinearReg: error in temp ', temp)
      end
    end
    Else   {an error exists}
      r := 0.0
	
  End; {LinearReg} 	        		 

{-----------------------------------------------------------
! Generates an array of normal random numbers from pairs of
! uniform random numbers in range [0,1].
!----------------------------------------------------------}
  Procedure NormalArray(a:pVec; n:Integer);
  Var
      i: Integer;
      u1,u2: Real;
  Begin

    {fills array with uniform random}
    For i:=1 to n do a^[i]:=Random;

    i:=1;
    Repeat
      u1 := a^[i];
      u2 := a^[i+1];
      If u1=0.0 Then u1 := 1e-12;  {u must not be zero}
      If u2=0.0 Then u2 := 1e-12;
      a^[i] := SQRT(-2.0*Ln(u1))*COS(2.0*pi*u2);
      a^[i+1] := SQRT(-2.0*Ln(u2))*SIN(2.0*pi*u2);
      Inc(i,2);
    Until i>=n;
	
    If (n MOD 2) <> 0 Then   {there is one extra element}
    begin
      If a^[n] = 0.0 Then a^[n] := 1e-12;  	   
      a^[n] := SQRT(-2.0*Ln(a^[n]))*SIN(2.0*pi*a^[n])
    end

  End; {NormalArray}


End. {Unit Stats

end of file Stats.pas}