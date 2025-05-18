{*********************************************************
*   Demonstration Program of In-place Merge Sorting      * 
*  (about n*n comparisons used with no extra storage).   *
* ------------------------------------------------------ *
* Reference: After a java algorithm By Jason Harrison,   * 
*            1995 University of British Columbia.        *
*                                                        *
*                    TPW Version By J-P Moreau, Paris.   *
*                            (www.jpmoreau.fr)           *
* ------------------------------------------------------ *
* SAMPLE RUN:                                            *
*                                                        *
* Initial table A:                                       *
*  4.00   3.00   1.00  67.00  55.00   8.00   0.00   4.00 *
* -5.00  37.00   7.00   4.00   2.00   9.00   1.00  -1.00 *
*                                                        *
* Sorted table A:                                        *
* -5.00  -1.00   0.00   1.00   1.00   2.00   3.00   4.00 *
*  4.00   4.00   7.00   8.00   9.00  37.00  55.00  67.00 *
*                                                        *
*********************************************************}
Program MSORT;

Uses WinCrt;

Const NMAX = 1024;

Type  VEC = Array[0..NMAX-1] of Double;

Var   a: VEC;
      i, n: Integer;
       

  {*****************************************************
  * In-place sorting of a table a[] in ascending order *
  * -------------------------------------------------- *
  * Inputs:                                            *
  *         a :  Table of n real numbers               *
  *         lo:  Starting index (>=0)                  *
  *         hi:  Ending index (<=n-1), with lo<hi.     *
  * Output:                                            *
  *         a :  Table sorted in ascending order.      *
  *****************************************************}
  Procedure sort(Var a: VEC; lo, hi: Integer);
  Label return;
  Var end_lo,k,mid,start_hi: Integer;
      T: Double;
  Begin

    if lo>=hi then goto return;  {no action}  
    mid := (lo + hi) Div 2;

    {Partition the list into two lists and sort them recursively}
    sort(a, lo, mid);
    sort(a, mid + 1, hi);

    {Merge the two sorted lists}
    start_hi := mid + 1;
    end_lo := mid;
    while (lo <= end_lo) and (start_hi <= hi) do
    begin
      if a[lo] < a[start_hi] then
        Inc(lo)
      else
      begin
      { a[lo] >= a[start_hi]
        The next element comes from the second list, 
        move the a[start_hi] element into the next 
        position and shuffle all the other elements up }
        T := a[start_hi];
        for k := start_hi - 1 Downto lo do a[k+1] := a[k];
        a[lo] := T;
        Inc(lo);
        Inc(end_lo);
        Inc(start_hi)
      end
    end;

  return: End;


  {main program}
  BEGIN

    n:=16;
    
    a[0] := 4; a[1] := 3; a[2] :=1; a[3] :=67;
    a[4] :=55; a[5] := 8; a[6] :=0; a[7] := 4;
    a[8] :=-5; a[9] :=37; a[10]:=7; a[11]:= 4;
    a[12]:= 2; a[13]:= 9; a[14]:=1; a[15]:=-1;

    writeln;
    writeln(' Initial table A:');
    for i:=0 to n-1 do
    begin
      write(a[i]:6:2,' ');
      if i=7 then writeln
    end;

    sort(a,0, n-1);

    writeln;
    writeln;
    writeln(' Sorted table A:');
    for i:=0 to n-1 do
    begin
      write(a[i]:6:2,' ');
      if i=7 then writeln
    end;
    writeln;
    readKey;
    DoneWinCrt

  END.

{end of file msort.pas}