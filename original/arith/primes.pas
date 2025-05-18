{*********************************************
*           TABLE OF PRIME NUMBERS           *
* ------------------------------------------ *
* This small program writes a table of prime *
* numbers from 1 to N (here N=2000).         *
* ------------------------------------------ *
*                                            *
*             Pascal version by J-P Moreau.  *
*                   (www.jpmoreau.fr)        *
*********************************************}
PROGRAM Table_of_primes;
Uses WinCrt;

LABEL   10;

VAR
        i, j, m, N, prem : INTEGER;
        Premier : Array[1..400] of INTEGER;
        Raye    : Array[1..2000] of INTEGER;

BEGIN

  N := 2000;

  for i:=1 to N do  Raye[i]:=0;
  Raye[1] := 1;

  prem := 1;
  while prem*prem <= N do
  begin
    while 1<2 do  {permanent condition}
    begin 
      prem:=prem+1;
      if (prem > N) or (Raye[prem]=0) then goto 10;
    end;
10: i := 2*prem;
    while i<=N do
    begin
      Raye[i] := 1;
      i := i + prem;
    end
  end;
  
  j := 0;
  for i:=1 to N do
    if (Raye[i]=0) then
    begin
      j:=j+1; 
      Premier[j] := i
    end;

  m:=0;
  writeln(' Between 1 and ',N,', there are ',j,' prime numbers:');
  writeln;
  for i:=1 to j do
  begin
    m:=m+1;
    write(Premier[i]:5);
    if (m=15) then begin writeln; m:=0 end
  end;
  writeln;
  Readkey; DoneWinCrt

END.

{end of file primes.pas}