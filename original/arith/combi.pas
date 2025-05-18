{*********************************************
* This program allows user performing combi- *
* natory analysis, such as Factorial N,      *
* Combination C(n,p) and Permutation A(n,p). *
* ------------------------------------------ *
* Ref.: "Mathématiques en Turbo-Pascal       *
*        By M. Ducamp and A. Reverchon (2),  *
*        Eyrolles, Paris, 1988" [BIBLI 05].  *
* ------------------------------------------ *
* Sample runs:                               *
*                                            *
* COMBINATORY ANALYSIS                       *
*                                            *
*   1:  Factorial n!                         *
*   2:  Combination C n,p                    *
*   3:  Permutation A n,p                    *
*   0:  Quit                                 *
*                                            *
*   Your choice (0 to 3): 1                  *
*                                            *
*   N = 100                                  *
*   N! = 9.3248476268  10^ 157               *
*                                            *
*   Your choice (0 to 3): 2                  *
*                                            *
*   N = 7                                    *
*   P = 3                                    *
*   Cnp =           35                       *
*                                            *
*   Your choice (0 to 3): 3                  *
*                                            *
*   N = 10                                   *
*   P = 6                                    *
*   Anp =       151200                       *
*                                            *
* ------------------------------------------ *
* Function Cnp: bug corrected 03/13/2004.    * 
*********************************************}
PROGRAM COMBI_ANALYSIS;
Uses WinCrt;

Const MAXREAL = 1E4000;

Var   n,p,choice: INTEGER; r1,r2,s: REAL;


  {Factorial n!}
  Function Factorial(n:INTEGER;VAR mantissa,exponent:REAL):REAL;
  Var  i:INTEGER; fa:REAL;
  Begin
    fa:=1; Factorial:=0;
    if (n<25) then
    begin
      for i:=1 to n do fa:=fa*i;
      Factorial:=fa
    end
    else
    begin
      fa:=(Ln(2*PI*n)/2 + n*Ln(n)-n)/Ln(10);
      exponent:=INT(fa); mantissa:=Exp(Frac(fa)*Ln(10))
    end
  End;

  {Combination Cn,p}
  Function Cnp(n,p:INTEGER):REAL;
  Var  i,u:INTEGER; r,x:REAL;
  Begin
    Cnp:=0;
    if p>n then exit;
    r:=1;
    if p>n-p then u:=n-p else u:=p;
    for i:=0 to u-1 do
    begin
      x:=(n-i)/(u-i);
      if r>MAXREAL/x then exit;
      r:=r*x
    end;
    Cnp:=r
  End;

  {Permutation An,p}
  Function Anp(n,p:INTEGER):REAL;
  Var  i:INTEGER; r:REAL;
  Begin
    Anp:=0;
    if p>n then exit;
    r:=1;
    for i:=n-p+1 to n do
    begin
      if r>MAXREAL/i then exit;
      r:=r*i
    end;
    Anp:=r
  End;


BEGIN     {main program}
  Repeat
    Clrscr;
    writeln;
    writeln(' COMBINATORY ANALYSIS');
    writeln;
    writeln('  1:  Factorial n!');
    writeln('  2:  Combination Cn,p');
    writeln('  3:  Permutation An,p');
    writeln('  0:  Quit');
    writeln;
    write('  Your choice (0 to 3): '); readln(choice);
    writeln;
    if choice<0 then choice:=0;
    if choice>3 then choice:=3;
    Case choice of
      0:  DoneWinCrt;     {quit}
      1:  begin           {Factorial n}
            write('  N = '); readln(n);
            s:=Factorial(n,r1,r2);
            write('  N! = ');
            if s>0 then writeln(s:12:0)
                   else writeln(r1:12:10,'   10^',r2:4:0)
          end;
      2:  begin           {Combination n,p}
            write('  N = '); readln(n);
            write('  P = '); readln(p);
            s:=Cnp(n,p);
            write('  Cnp = ');
            if s>0 then writeln(s:12:0)
                   else writeln('  Quantity impossible to evaluate')
          end;
      3:  begin           {Permutation n,p}
            write('  N = '); readln(n);
            write('  P = '); readln(p);
            s:=Anp(n,p);
            write('  Anp = ');
            if s>0 then writeln(s:12:0)
                   else writeln('  Quantity impossible to evaluate')
          end;
    End;
    Readkey  {Pause}
  Until choice=0
END.

{end of file combi.pas}