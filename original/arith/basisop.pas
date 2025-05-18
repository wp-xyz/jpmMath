{*********************************************
*  This program allows making arithmetic     *
*  operations in any base between 2 and 36.  *
* ------------------------------------------ *
* Ref.: "Mathématiques en Turbo-Pascal       *
*        By M. Ducamp and A. Reverchon (2),  *
*        Eyrolles, Paris, 1988" [BIBLI 05].  *
* ------------------------------------------ *
* Sample runs:                               *
*                                            *
* ARITHMETIC OPERATIONS IN ANY BASE          *
* BETWEEN 2 AND 36                           *
*                                            *
* Enter base (2 to 36): 25                   *
* Enter 1st number : AM3G                    *
*                                            *
* Enter next number: 28                      *
* Operator (+-*/^) : +                       *
* Result           : AM50                    *
* Enter next number: J                       *
* Operator (+-*/^) : *                       *
* Result           : 86MD6                   *
* Enter next number: 2B                      *
* Operator (+-*/^) : /                       *
* Result           : 39JM                    *
*                                            *
* ARITHMETIC OPERATIONS IN ANY BASE          *
* BETWEEN 2 AND 36                           *
*                                            *
* Enter base (2 to 36): 2                    *
* Enter 1st number : 11001010                *
*                                            *
* Enter next number: 1100                    *
* Operator (+-*/^) : -                       *
* Result           : 10111110                *
* Enter next number: 1110                    *
* Operator (+-*/^) : -                       *
* Result           : 10110000                *
* Enter next number: 10101                   *
* Operator (+-*/^) : +                       *
* Result           : 11000101                *
* Enter next number: 10                      *
* Operator (+-*/^) : /                       *
* Result           : 1100010                 *
* Enter next number: 11                      *
* Operator (+-*/^) : *                       *
* Result           : 100100110               *
*                                            *
* ARITHMETIC OPERATIONS IN ANY BASE          *
* BETWEEN 2 AND 36                           *
*                                            *
* Enter base (2 to 36): 16                   *
* Enter 1st number : FEAA                    *
*                                            *
* Enter next number: FF                      *
* Operator (+-*/^) : +                       *
* Result           : FFA9                    *
* Enter next number: 1799                    *
* Operator (+-*/^) : -                       *
* Result           : E810                    *
* Enter next number: 2                       *
* Operator (+-*/^) : /                       *
* Result           : 7408                    *
*                                            *
* Note: Enter null string to exit.           *
*                                            *
*                  TPW Version By J-P Moreau *
*                      (www.jpmoreau.fr)     *
**********************************************
Explanations:
------------
The letters A,B,...,Z  represent 10,11,...,36
in the base > 10.

The numbers are first converted from base b to
base 10 by Function Decodebase, then  operation
is made in base 10 and result is converted from
base 10 to base b by Function Codebase.

Note: assuming that a real number has a mantissa of
      40 digits in base 2, the number of useful
      digits udg in base b is given by formula:
	     
	   udg = integer(40*log(2)/log(b))

           b=10: udg = 12
	   b=16: udg = 10
	   b=36: udg = 7

--------------------------------------------------}
PROGRAM base;
Uses WinCrt;

CONST  MAXREAL = 1E40;

VAR  b    : INTEGER;
     x,y  : STRING;
     op   : CHAR;
     r1,r2: REAL;

{Convert a number from base b to base 10. The function returns
 FALSE if b not in [2..36] or if string x contains invalid
 characters in base b or if result y is too big}
Function Decodebase(x:STRING;b:INTEGER;VAR y:REAL): BOOLEAN;
Var  i,j,long: INTEGER;
     mult: REAL; ch:CHAR;
Begin
  Decodebase:=FALSE;
  if (b<2) or (b>36) then
  begin
    writeln(' base must be between 2 and 36 !');
    exit
  end;
  y:=0; mult:=1;
  long:=length(x);
  for i:=1 to long do
  begin
    ch:=UpCase(x[long+1-i]);
    if (ch<'0') or (ch>'Z') or ((ch>'9') and (ch<'A')) then  exit;
    if ch<='9' then j:=Ord(ch)-Ord('0')
               else j:=Ord(ch)-Ord('A')+10;
    if j>=b then exit;
    y:=y+mult*j;
    if mult>MAXREAL/b then exit;
    mult:=mult*b
  end;
  Decodebase:=TRUE
End;

{Convert a number from base 10 to base b. The function returns
 FALSE if b not in [2..36] or if string x contains invalid
 characters in base 10 or if number x is too big}
Function Codebase(x:REAL;b:INTEGER;VAR y:STRING): BOOLEAN;
Var n: INTEGER;
Begin
  Codebase:=FALSE;
  if (b<2) or (b>36) then
  begin
    writeln(' base must be between 2 and 36 !');
    exit
  end;
  y:='';
  While x>0 do
  begin
    n:=TRUNC(x-b*INT(x/b));
    if n<10 then y:=Chr(Ord('0')+n)+y
            else y:=Chr(Ord('A')+n-10) +y;
    x:=INT(x/b)           
  end;
  Codebase:=TRUE
End;

{main program}
BEGIN
  writeln;
  writeln(' ARITHMETIC OPERATIONS IN ANY BASE BETWEEN 2 AND 36');
  writeln;
  Repeat
    Write(' Enter base (2 to 36): '); readln(b);
    write(' Enter 1st number : '); readln(x)
  Until Decodebase(x,b,r1);
  Repeat
    writeln;
    write(' Enter next number: '); readln(x);
    if x<>'' then
    if Decodebase(x,b,r2) then
    begin
      write(' Operator (+-*/^) : '); readln(op);
      Case op of
        '+': r1:=r1+r2;
        '-': r1:=ABS(r1-r2);
        '*': r1:=r1*r2;
        '/': r1:=INT(r1/r2);
        '^': r1:=INT(EXP(r2*LN(r1)))
      End;
      if Codebase(r1,b,y) then writeln(' Result           : ',y)
                           else writeln(' Error in coding result.')
    end
    else
      writeln(' Number not valid or decoding error.')                           
  Until x='';
  DoneWinCrt
END.

{end of file basisop.pas}