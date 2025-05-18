{*************************************************************************
*    Calculate the median value of an array with the Heapsort method     *
* ---------------------------------------------------------------------- *
* REFERENCE:                                                             *
*              "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,          *
*               S.A. Teukolsky and W.T. Vetterling, Cambridge            *
*               University Press, 1986" [BIBLI 08].                      *
* ---------------------------------------------------------------------- *
* SAMPLE RUN:                                                            *
*                                                                        *
* Given table:                                                           *
*                                                                        *
*    1.0   63.8  278.9  406.2  546.8  657.7  638.4  324.6  745.5  852.3  *
*  165.0  950.6  142.1  319.3  120.4  587.6  166.4  736.8  451.7  656.9  *
*  605.7  312.7  565.0  614.3  326.3  660.0  933.0  494.3  349.6  559.1  *
*  964.5  299.4  252.3  575.6  455.5   48.1  986.1  225.1  346.4   41.6  *
*  283.1  288.0  999.4   44.4  815.1   20.3  452.0  699.7  460.0  584.8  *
*  886.0  413.1  638.8  815.3   90.1  712.9    4.3  488.6  648.7  592.4  *
*  172.8  456.4  991.6  237.6  981.1  854.4   94.6  631.3  678.4   16.2  *
*  590.3   23.1  452.9  501.9  396.4  288.3  173.8  432.6  822.4  271.6  *
*                                                                        *
* Median value: 455.9250                                                 *
*                                                                        *
* Sorted table (Heapsort method):                                        *
*                                                                        *
*    1.0    4.3   16.2   20.3   23.1   41.6   44.4   48.1   63.8   90.1  *
*   94.6  120.4  142.1  165.0  166.4  172.8  173.8  225.1  237.6  252.3  *
*  271.6  278.9  283.1  288.0  288.3  299.4  312.7  319.3  324.6  326.3  *
*  346.4  349.6  396.4  406.2  413.1  432.6  451.7  452.0  452.9  455.5  *
*  456.4  460.0  488.6  494.3  501.9  546.8  559.1  565.0  575.6  584.8  *
*  587.6  590.3  592.4  605.7  614.3  631.3  638.4  638.8  648.7  656.9  *
*  657.7  660.0  678.4  699.7  712.9  736.8  745.5  815.1  815.3  822.4  *
*  852.3  854.4  886.0  933.0  950.6  964.5  981.1  986.1  991.6  999.4  *
*                                                                        *
*                                                                        *
*                                 Pascal Release By J-P Moreau, Paris.   *
*                                          (www.jpmoreau.fr)             *
*************************************************************************}
PROGRAM TMDIAN;
Uses WinCrt;

Const NMAX = 1024;          {max. number of items}
      MAX_VALUE:REAL = 1000;

Type  pVEC = ^VEC;
       VEC = Array[1..NMAX] of REAL;

Var   A:pVEC;              {pointer to input random Table}
      I,N:Integer;
      AMED:REAL;


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
PROCEDURE HPSORT(N:INTEGER; RA:pVEC);
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
    RRA:=RA^[L]
  end
  else
  begin
    RRA:=RA^[IR];
    RA^[IR]:=RA^[1];
    IR:=IR-1;
    if IR=1 then
    begin
      RA^[1]:=RRA;
      exit
    end;
  end;
  I:=L;
  J:=L+L;
20: if J <= IR then
  begin
    if J < IR then
      if RA^[J] < RA^[J+1] then J:=J+1;
    if RRA < RA^[J] then
    begin
      RA^[I]:=RA^[J];
      I:=J; J:=J+J
    end
    else
      J:=IR+1;
    goto 20
  end;
  RA^[I]:=RRA;
  goto 10

END;


{******************************************************
* Given an array X of N numbers, returns their median *
* value XMED. The array X is modified and returned    *
* sorted in ascending order.                          *
******************************************************}
Procedure MDIAN(X:pVEC; N:Integer; Var XMED: REAL);
Var N2:Integer;
Begin
  hpsort(N,X);
  N2:=N Div 2;
  if 2*N2 = N then
    XMED := 0.5*(X^[N2]+X^[N2+1])
  else
    XMED := X^[N2+1]
End;


{write table of size N to standard output}
Procedure TWRIT(N:Integer; ARR:pVEC);
Var i: Integer;
Begin
  writeln;
  For i:=1 to N do
  begin
    write(' ',ARR^[i]:6:1);
    if (i MOD 10) = 0 then writeln
  end;
End;


{main program}
BEGIN

  New(A);

  N:=80;    {initialize size of table}

  {generate random table of numbers (from 0 to 1000) }
  For i:=1 to N do A^[i]:=1.0 + MAX_VALUE*Random;

  writeln(' Given table:');
  TWRIT(N,A);
 
  {call MDIAN procedure}
  MDIAN(A,N,AMED);

  writeln;
  writeln(' Median value: ', AMED:8:4);

  writeln;
  writeln(' Sorted table (Heapsort method):');
  TWRIT(N,A);

  writeln;
  ReadKey;
  Dispose(A);
  DoneWincrt

END.

{end of file tmdian.pas}