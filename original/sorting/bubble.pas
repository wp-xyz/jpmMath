{*********************************************
*   Demonstration program of Bubble sorting  * 
*        (about n*n comparisons used).       *
* ------------------------------------------ *
* Reference: "A book on C  By Al Kelley and  *
* Ira Pohl, The Benjamin/Cummings Publishing *
* Company, Inc, 1984" [BIBLI 09].            *
*                                            *
*             Pascal version by J-P Moreau.  *
*                  (www.jpmoreau.fr)         *
* ------------------------------------------ *
* SAMPLE RUN:                                *
*                                            *
* Initial table A:                           *
* 7  3  66  3  -5  22  -77  2  36  -12       *
*                                            *
* Sorted table A:                            *
* -77  -12  -5  2  3  3  7  22  36  66       *
*                                            *
*********************************************}
Program Bubble_sort;
Uses WinCrt;

Const MAX = 1024;

Type TAB = Array[1..MAX] of Integer;

Var i,n:Integer; A:TAB;

{return p,q in ascending order}
Procedure Order(VAR p,q:Integer);
Var temp:Integer;
Begin
  if p>q then
  begin
    temp:=p;
    p:=q;
    q:=temp
  end
End;

{Buuble sorting of integer array A}
Procedure Bubble(VAR A:TAB; n:Integer);
Var i,j:Integer;
Begin
  for i:=1 to n do
    for j:=n downto i+1 do
      Order(A[j-1], A[j])
End;

BEGIN
  n:=10;
  A[1]:=7;  A[2]:=3;   A[3]:=66; A[4]:=3;  A[5]:=-5;
  A[6]:=22; A[7]:=-77; A[8]:=2;  A[9]:=36; A[10]:=-12;

  writeln;
  writeln(' Initial table A:');
  for i:=1 to n do
    write(' ',A[i]);
  writeln;

  Bubble(A,n);

  writeln;
  writeln(' Sorted table A:');
  for i:=1 to n do
    write(' ',A[i]);
  writeln;

  writeln;
  ReadKey; DoneWinCrt
END.

{end of file bubble.pas}

