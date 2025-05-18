{*********************************************
* Solving a diophantian equation ax+by = c   *
*      (a,b,c,x,y are integer numbers)       *
* ------------------------------------------ *
* Ref.: "Mathématiques en Turbo-Pascal       *
*        By M. Ducamp and A. Reverchon (2),  *
*        Eyrolles, Paris, 1988" [BIBLI 05].  *
* ------------------------------------------ *
* Sample run:                                *
*                                            *
* SOLVING IN Z EQUATION AX + BY = C          *
*                                            *
*    A = 3                                   *
*    B = -2                                  *
*    C = 7                                   *
*                                            *
* Solutions are:                             *
*                                            *
*   X =    1 +     2*K                       *
*   Y =   -2 +     3*K                       *
*                                            *       
*********************************************}
Program Diophantian_Equation;
Uses WinCrt;

VAR     a,b,c,x0,y0,p,q : REAL;


{return greatest common divisor of two integer numbers}
Function GCD(a,b:REAL):REAL;
Var  x,r: REAL;
Begin
  a:=INT(ABS(a));
  b:=INT(ABS(b));
  gcd:=1;
  if (a>1E10) or (b>1E10) then exit;
  if (a=0) or (b=0) then exit;
  if a<b then
  begin
    x:=a; a:=b; b:=x
  end;
  Repeat
    r:=a-b*INT(a/b);
    a:=b; b:=r
  Until abs(r)<1E-10;
  gcd:=a
End;


{**********************************************************
* Solving equation ax+by=c, a,b,c,x,y are integer numbers *
* ------------------------------------------------------- *
* INPUT:   a,b,c       coefficients of equation           *
* OUTPUT:  solutions are x0+kp and y0-kq, with k=0,1,2... *
*          or k=-1,-2,-3...                               *
* The function returns TRUE if solutions exist (that is,  *
* if the GCD of a,b is also a divisor of c).              *
**********************************************************}                                                         
Function Diophantian(a,b,c:REAL; VAR x0,y0,p,q:REAL):BOOLEAN;
Var  pg,x1,x2,y1,y2 : REAL;
     found : BOOLEAN;
Begin
  Diophantian:=FALSE;
  if a=0 then exit; if b=0 then exit;
  pg:=GCD(a,b);
  a:=a/pg; b:=b/pg; c:=c/pg;
  if c<>INT(c) then exit; {pg must be also a divisor of c}
  x1:=0; y2:=0;
  found:=FALSE;
  Repeat
    y1:=(c-a*x1)/b;
    if y1=INT(y1) then
    begin
      x0:=x1; y0:=y1;
      found:=TRUE;
    end
    else
    begin      
      x1:=-x1; if x1>=0 then x1:=x1+1;
      x2:=(c-b*y2)/a;
      if x2=INT(x2) then
      begin
        x0:=x2; y0:=y2; found:=TRUE;
      end
      else
      begin
        y2:=-y2; if y2>=0 then y2:=y2+1
      end
    end
  Until found;
  p:=a; q:=b;
  Diophantian:=TRUE
End;


{main program}
BEGIN

  Writeln;
  Writeln(' SOLVING IN Z EQUATION AX + BY = C');
  Writeln;
  Write('    A = '); Readln(a);
  Write('    B = '); Readln(b);
  Write('    C = '); Readln(c);
  Writeln;

  if Diophantian(a,b,c,x0,y0,p,q) then
  begin
    writeln(' Solutions are:');
    writeln;
    writeln('   X=',x0:6:0,' +',abs(q):6:0,' *K');
    write('   Y=',y0:6:0);
    if p*q>0 then write(' -') else write(' +');
    writeln(abs(p):6:0,' *K');
  end
  else
    writeln(' No solutions.');
  writeln;
  ReadKey; DoneWinCrt      

END.

{end of file diophan.pas}