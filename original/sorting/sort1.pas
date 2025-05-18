{**************************************************
* Sorting an array with straight insertion method *
* ----------------------------------------------- *
* REFERENCE:                                      *
*                                                 *
* "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,*
*  S.A. Teukolsky and W.T. Vetterling, Cambridge  *
*  University Press, 1986" [BIBLI 08].            *
*                                                 *
*                   Pascal Version By J-P Moreau  *
*                        (www.jpmoreau.fr)        *
* ----------------------------------------------- *
* SAMPLE RUN:                                     *
*                                                 *
* Table to be sorted:                             *
*                                                 *
*   241   338   602   630   211                   *
*   608    13   674   317   361                   *
*   255   873   159   541   865                   *
*   222   125   862   505   637                   *
*    94   120   652   190   338                   *
*                                                 *
* Sorted table (straight insertion):              *
*                                                 *
*    13    94   120   125   159                   *
*   190   211   222   241   255                   *
*   317   338   338   361   505                   *
*   541   602   608   630   637                   *
*   652   674   862   865   873                   *
*                                                 *
**************************************************}
PROGRAM SORT1;
Uses WinCrt;

CONST   SIZE = 100;                   {maximum size of table}

TYPE
        Table = Array[1..SIZE] of real;

VAR
        A         : Table;               {Table to be sorted}
        MAX_VALUE : real;            {Maximum value of table}
        i,N       : integer;


{****************************************************
* Sorts an array ARR of length N in ascending order *
* by straight insertion.                            *
* ------------------------------------------------- *
* INPUTS:                                           *
*	   N	  size of table ARR                 *
*          ARR	  table to be sorted                *
* OUTPUT:                                           *
*	   ARR   table sorted in ascending order    *
*                                                   *
* NOTE: Straight insertion is a N² routine and      *
*       should only be used for relatively small    *
*       arrays (N<100).                             *
****************************************************}         
Procedure PIKSRT(N:INTEGER;VAR ARR:Table);
Label 10;
Var i,j:INTEGER; a:REAL;
Begin
  for j:=2 to N do
  begin
    a:=ARR[j];
    for i:=j-1 downto 1 do
    begin
      if ARR[i]<=a then goto 10;
      ARR[i+1]:=ARR[i]
    end;
    i:=0;
10:  ARR[i+1]:=a
  end
End;

{write table of size N to standard output}
Procedure TWRIT(N:INTEGER;ARR:Table);
Var i:integer;
begin
  writeln;
  for i:=1 to N do
  begin
    write(ARR[i]:6:0);
    if (i MOD 5)=0 then writeln
  end
end;


{main program}
BEGIN
  {initialize random number generator}
  Randomize;

  N:=25;    {initialize size of table}
  MAX_VALUE := 1000.0;

  {generate random table of numbers (between 0 and MAX_VALUE) }
  for i:=1 to N do A[i]:=MAX_VALUE*Random;

  writeln;
  writeln(' Table to be sorted:');
  TWRIT(N,A);
 
  {call sorting subroutine}
  PIKSRT(N,A);

  writeln;
  writeln(' Sorted table (straight insertion):');
  TWRIT(N,A);

  writeln;
  ReadKey; DoneWinCrt

END.

{end of file sort1.pas}