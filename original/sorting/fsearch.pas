{******************************************
!* Searching procedures in list of names  *
!* to be used by program Search.pas.      *
!*****************************************} 
UNIT FSearch;

INTERFACE

TYPE
     Str20 = String[20];

VAR
     a: Array[1..100] of Str20;
     target: Str20;

     Procedure FindFirst(size:Integer; Var where:Integer);
     Procedure FindAll(size:Integer; Var how_many:Integer);
     Procedure Binary(low,high:Integer; Var where:Integer);


IMPLEMENTATION

Procedure FindFirst(size:Integer; Var where:Integer);
{ Linear search of a list of names for first occurence
  of target value.  }
Begin  
  where:=1;
  WHILE (a[where] <> target) AND (where < size) DO Inc(where);
  if a[where] <> target then where:=0
End;

Procedure FindAll(size:Integer; Var how_many:Integer);
{ Linear search of a list of names for all occurences
  of target value.  }
Var i: INTEGER;
Begin
  how_many:=0;
  For i:=1 to size do
    IF a[i] = target THEN
    begin
      WRITELN(' ',target,' at position ', i);
      Inc(how_many)
    end
End;

Procedure Binary(low,high:Integer; Var where:Integer);
{ Binary search of an ordered list for one occurence of
  specified target value. Assumes low < high, i.e. there
  is something in the list to look for.   }
Var
    mid, lo, hi: Integer;
Begin
  lo:=low; hi:=high; where:=0;
  WHILE (lo <= hi) AND (where = 0) DO
  begin
    mid:=(lo+hi) Div 2;
    IF a[mid] = target THEN
      where:=mid
    ELSE IF a[mid] > target THEN
      hi:=mid-1
    ELSE
      lo:=mid+1
  end
End;

END.

{end of file fsearch.pas}