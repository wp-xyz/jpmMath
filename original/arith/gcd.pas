{*********************************************
* Geatest common divisor and Smallest common *
*     multiple of several integer numbers    *
* ------------------------------------------ *
* Ref.: "Mathématiques en Turbo-Pascal       *
*        By M. Ducamp and A. Reverchon (2),  *
*        Eyrolles, Paris, 1988" [BIBLI 05].  *
* ------------------------------------------ *
* Sample run:                                *
*                                            *
* GCD and SCM of integer numbers:            *
*                                            *
* First number: 9936                         *
* Next number : 414                          *
* Next number : 3174                         *
* Next number : 0                            *
*                                            *
* GCD =      138                             *
* SCM =   228528                             *
*                                            *
*********************************************}
PROGRAM GCD_SCM;
Uses WinCrt;

Var  x,pg,pp : REAL;


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


BEGIN

  clrscr;
  writeln;
  writeln(' GCD and SCM of integer numbers:');
  writeln;
  write(' First number: '); readln(pg); pp:=pg;
  Repeat
    write(' Next number : '); readln(x);
    if x>0 then
    begin
      pg:=gcd(pg,x); pp:=pp*x/gcd(pp,x)
    end
  Until x=0;
  writeln;
  writeln(' GCD = ',pg:12:0);
  writeln(' SCM = ',pp:12:0);
  writeln;
  Readkey; Donewincrt

END.

{end of file gcd.pas}