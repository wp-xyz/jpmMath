{*********************************************
* This program converts a number in base a   *
* into a number in base b (a and b must be   * 
* between 2 and 36).                         * 
* ------------------------------------------ *
* Ref.: "Mathématiques en Turbo-Pascal       *
*        By M. Ducamp and A. Reverchon (2),  *
*        Eyrolles, Paris, 1988" [BIBLI 05].  *
* ------------------------------------------ *
* Sample runs:                               *
*                                            *
* BASE CONVERSION                            *
*                                            *
* Start  base (2 to 36): 5                   *
* Arival base (2 to 36): 3                   *
*                                            *
* Enter number in start base: 3421           *
* In base 3: 200000                          *
* Enter number in start base: 3420           *
* In base 3: 122222                          *
*                                            *
*                                            *
* BASE CONVERSION                            *
*                                            *
* Start  base (2 to 36): 10                  *
* Arival base (2 to 36): 16                  *
*                                            *
* Enter number in start base: 65535          *
* In base 3: FFFF                            *
* Enter number in start base: 100            *
* In base 3: 64                              *
*                                            *
* BASE CONVERSION                            *
*                                            *
* Start  base (2 to 36): 16                  *
* Arival base (2 to 36): 2                   *
*                                            *
* Enter number in start base: FF             *
* In base 3: 11111111                        *
* Enter number in start base: 3A             *
* In base 3: 111010                          *
* Enter number in start base: E2             *
* In base 3: 11100010                        *
*                                            *
* Note: Enter null string to exit.           *
*                                            *
*                  TPW version by J-P Moreau *
*                      (www.jpmoreau.fr)     *
**********************************************
Explanations:
------------
The letters A,B,...,Z  represent 10,11,...,36
in the base > 10.

The number is first converted from base a to
base 10 by Function Decodebase, then  converted
from base 10 to base b by Function Codebase.
----------------------------------------------}
PROGRAM base;
Uses WinCrt;

CONST  MAXREAL = 1E40;

VAR  ba,bd: INTEGER;
     x,y: STRING;
     r: REAL;

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
  writeln(' BASE CONVERSION');
  writeln;
  Write(' Start  base (2 to 36): '); readln(bd);
  Write(' Arival base (2 to 36): '); readln(ba);
  Repeat
    writeln;
    write(' Enter number in start base: '); readln(x);
    if x<>'' then
      if Decodebase(x,bd,r) then
        if Codebase(r,ba,y) then
          writeln(' In base ',ba:2,': ',y)
        else
          writeln(' Error in coding number.')
      else
        writeln(' Errorin decoding number.')
  Until  x='';
  DoneWinCrt
END.

{end of file base.pas}