{**************************************************************
*           Sorting an array with the Shell method            *
* ----------------------------------------------------------- *
* REFERENCE:                                                  *
*                                                             *
*   "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,          *
*    S.A. Teukolsky and W.T. Vetterling, Cambridge            *
*    University Press, 1986" [BIBLI 08].                      *
*                                                             *
*                              Pascal Version By J-P Moreau   *
*                                   (www.jpmoreau.fr)         *
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
* Sorted table (Shell method):                                *
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
PROGRAM SORT2;
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
*            by the Shell-Mezgar method             *
* ------------------------------------------------- *
* INPUTS:                                           *
*	   N	  size of table ARR                 *
*          ARR	  table to be sorted                *
* OUTPUT:                                           *
*	   ARR   table sorted in ascending order    *
*                                                   *
* NOTE: The Shell method is a N^3/2 routine and can *
*       be used for relatively large arrays.        *
****************************************************}         
PROCEDURE SHELL(N:integer;VAR ARR:Table);
Label 10;
Const
      ALN2I=1.0/0.69314718;
      TINY =1E-5;
var
      LOGNB2, t : REAL;
      i,j,k,l,m,nn : INTEGER;
Begin
  LOGNB2:=INT(LN(N)*ALN2I+TINY);
  m:=n;
  for nn:=1 to Round(LOGNB2) do
  begin
    m:=m DIV 2; k:=n-m;
    for j:=1 to k do
    begin
      i:=j;
10:   l:=i+m;
      if ARR[l] < ARR[i] then
      begin
        t:=ARR[i];
	ARR[i]:=ARR[l];
	ARR[l]:=t;
	i:=i-m;
	if i >= 1 then GOTO 10
      end
    end
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
  SHELL(N,A);

  writeln;
  writeln(' Sorted table (Shell method):');
  TWRIT(N,A);

  ReadKey; DoneWinCrt

END.

{end of file sort2.pas}