{**************************************************************
*          Sorting an array with the Heapsort method          *
* ----------------------------------------------------------- *
* REFERENCE:                                                  *
*                                                             *
*   "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,          *
*    S.A. Teukolsky and W.T. Vetterling, Cambridge            *
*    University Press, 1986" [BIBLI 08].                      *
*                                                             *
*                              Pascal Version By J-P Moreau   *
*                                    (www.jpmoreau.fr)        *
* ----------------------------------------------------------- *
* SAMPLE RUN:                                                 *
*                                                             *
* Table to be sorted:                                         *
*                                                             *
*   798   805   701   112   378   561   508   525   956   633 *
*   915    66     5   582   881   942    25   357   847   643 *
*   549   284   815   348   572   798   330   139   433   543 *
*   594   757   338   581   492   720   731   788   369   248 *
*   700   822   946   457   188   425   141   701   500   282 *
*   672   110   574   170   497   350   921   447   393   986 *
*   782   929   205   398   559   441   808   520   757   375 *
*   305   503   590   595   516   177   430   551   911   645 *
*                                                             *
* Sorted table (Heapsort method):                             *
*                                                             *
*     5    25    66   110   112   139   141   170   177   188 *
*   205   248   282   284   305   330   338   348   350   357 *
*   369   375   378   393   398   425   430   433   441   447 *
*   457   492   497   500   503   508   516   520   525   543 *
*   549   551   559   561   572   574   581   582   590   594 *
*   595   633   643   645   672   700   701   701   720   731 *
*   757   757   782   788   798   798   805   808   815   822 *
*   847   881   911   915   921                               *
*                                                             *
**************************************************************}
PROGRAM SORT3;
Uses WinCrt;

CONST   SIZE = 100;                   {maximum size of table}

TYPE
        Table = Array[1..SIZE] of real;

VAR
        A         : Table;               {Table to be sorted}
        MAX_VALUE : real;            {Maximum value of table}
        i,N       : integer;

{****************************************************
*  Sorts an array RA of length N in ascending order *
*                by the Heapsort method             *
* ------------------------------------------------- *
* INPUTS:                                           *
*	   N	  size of table RA                  *
*          RA	  table to be sorted                *
* OUTPUT:                                           *
*	   RA    table sorted in ascending order    *
*                                                   *
* NOTE: The Heapsort method is a N Log2 N routine,  *
*       and can be used for very large arrays.      *
****************************************************}         
PROCEDURE HPSORT(N:INTEGER;VAR RA:Table);
Label 10, 20;
Var
  i,ir,j,l:INTEGER;
  rra:REAL;
Begin
  L:=(N DIV 2)+1;
  IR:=N;
  {The index L will be decremented from its initial value during the
   "hiring" (heap creation) phase. Once it reaches 1, the index IR 
   will be decremented from its initial value down to 1 during the
   "retirement-and-promotion" (heap selection) phase.                 }
10: if L > 1 then
  begin 
    L:=L-1;
    RRA:=RA[L]
  end
  else
  begin
    RRA:=RA[IR];
    RA[IR]:=RA[1];
    IR:=IR-1;
    if IR=1 then
    begin
      RA[1]:=RRA;
      exit
    end;
  end;
  I:=L;
  J:=L+L;
20: if J <= IR then
  begin
    if J < IR then
      if RA[J] < RA[J+1] then J:=J+1;
    if RRA < RA[J] then
    begin
      RA[I]:=RA[J];
      I:=J; J:=J+J
    end
    else
      J:=IR+1;
    goto 20
  end;
  RA[I]:=RRA;
  goto 10

END;


{write table of size N to standard output}
Procedure TWRIT(N:INTEGER;ARR:Table);
Var i:integer;
begin
  writeln;
  for i:=1 to N do
  begin
    write(ARR[i]:6:0);
    if (i MOD 10)=0 then writeln
  end
end;


{main program}
BEGIN
  {initialize random number generator}
  Randomize;

  N:=80;    {initialize size of table}
  MAX_VALUE := 1000.0;

  {generate random table of numbers (between 0 and MAX_VALUE) }
  for i:=1 to N do A[i]:=MAX_VALUE*Random;

  writeln;
  writeln(' Table to be sorted:');
  TWRIT(N,A);
 
  {call sorting subroutine}
  HPSORT(N,A);

  writeln;
  writeln(' Sorted table (Heapsort method):');
  TWRIT(N,A);

  ReadKey; DoneWinCrt

END.

{end of file sort3.pas}