{*********************************************
*   Demonstration program of Merge sorting   * 
*      (about n*log(n) comparisons used).    *
* ------------------------------------------ *
* Reference: "A book on C By Al Kelley and   *
* Ira Pohl, The Benjamin/Cummings Publishing *
* Company, Inc, 1984" [BIBLI 09].            *
*                                            *
*             Pascal version by J-P Moreau.  *
*                   (www.jpmoreau.fr)        *
* ------------------------------------------ *
* SAMPLE RUN:                                *
*                                            *
* Initial table A:                           *
* 4 3 1 67 55 8 0 4 -5 37 7 4 2 9 1 -1       *
*                                            *
* Sorted table A:                            *
* -5 -1 0 1 1 2 3 4 4 4 7 8 9 37 55 67       *
*                                            *
*********************************************}
Program Merge_Sorting;
Uses WinCrt;

Const MAXSIZE=1023;

Type
     ptab = ^TAB;
     TAB  = Array[0..MAXSIZE] of Integer;

Var  i,error,n: Integer;
     a: ptab; 

{ Merge Sorting of a[] of size m and b[] of size n
 into c[] of size m+n }
Procedure Merge(a,b,c:ptab; m,n:Integer);
Var i,j,k:Integer;
Begin
  i:=0; j:=0; k:=0;
  while (i<m) and (j<n) do
  begin
    if a^[i]<b^[j] then
    begin
      c^[k]:=a^[i];
      Inc(i)
    end
    else
    begin
      c^[k]:=b^[j];
      Inc(j)
    end;
    Inc(k)
  end;
  while i<m do  {pick up any remainder}
  begin
    c^[k]:=a^[i];
    Inc(k); Inc(i)
  end;
  while j<n do
  begin
    c^[k]:=b^[j];
    Inc(k); Inc(j)
  end
End;

{use function Merge() to sort an array key[] of size n
 In this version, n is a power of two }
Procedure MergeSort(key:ptab; n:Integer; VAR error:Integer);
Var j,k,m,m1:Integer; w:ptab;  {w is a work space}
Begin
  New(w);
  for j:=0 to MAXSIZE do w^[j]:=0;
  m:=1;
  for m1:=1 to n-1 do
    while m<n do m:=m*2;
  if (n=m) and (n<=MAXSIZE) then   {n is a power of two <= MAXSIZE}
  begin
    k:=1;
    while k<n do
    begin
      j:=0;  
      while j<n-k do
      begin
        Merge(@key^[j],@key^[j+k],@w^[j],k,k);
        j:=j+2*k
      end;
      k:=k*2;
      for j:=0 to n-1 do    {write w back into key}
        key^[j] := w^[j];
    end
  end
  else
  begin
    if n<>m then
    begin
      writeln(' Error #1: size n of array not a power of 2!');
      error:=1
    end;
    if n>MAXSIZE then
    begin
      writeln(' Error #2: size of array is two bog!');
      error:=2
    end
  end;
  Dispose(w);
End; {mergesort}


BEGIN

  n:=16; New(a);

  for i:=0 to MAXSIZE do a^[i]:=0;

  a^[0]:=4;  a^[1]:=3;   a^[2]:=1;  a^[3]:=67;
  a^[4]:=55; a^[5]:=8;   a^[6]:=0;  a^[7]:=4;
  a^[8]:=-5; a^[9]:=37; a^[10]:=7; a^[11]:=4;
  a^[12]:=2; a^[13]:=9;  a^[14]:=1; a^[15]:=-1;


  writeln;
  writeln(' Initial table A:');
  for i:=0 to n-1 do write(' ',a^[i]);
  writeln;

  MergeSort(a,n,error);

  if error=0 then
  begin
    writeln(' Sorted table A:');
    for i:=0 to n-1 do write(' ',a^[i]);
  end;
  writeln;

  Dispose(a);
  Readkey; Donewincrt

END.

{end of file merge.pas}

