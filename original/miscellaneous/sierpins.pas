{*************************************************************
*     One-dimensional cellular automata with rule that       *
*           generates a Sierpinski triangle                  *
* ---------------------------------------------------------- *
* Ref.: "Problem Solving with Fortran 90 By David R. Brooks, *
*        Springer-Verlag New York, 1997".                    *
*                                                            *
*                         TPW Release By J-P Moreau, Paris.  *
*                                (www.jpmoreau.fr)           *
*************************************************************} 
Program Sierpinski;

Uses WinCrt1;

Const Size = 60;

Type  Vec = Array[1..Size] of Boolean;

Var
    a, old_a: Vec;
    a_1, a0, a1: Boolean;

    b: Array[1..Size] of char;

    i, j, cycle, n_cycles: Integer;

Begin

{ Data n_cycles,a,b /15,Size*.FALSE.,Size*' '/ }

  n_cycles := 15;

  For i := 1 to Size do
  begin
    a[i]:=False;
    b[i]:=' '
  end;

  cycle := 0;
  a[Size Div 2] := TRUE;   {start with a single live cell}
  b[Size Div 2] := '*';
  
  writeln;
  writeln('                                Sierpinski  Triangle');
  writeln;
  write(' Generation ',cycle:2);
  For i:=1 to Size do write(b[i]);
  writeln;

{ Generate more cycles }

  For i := 1 to n_cycles do
  begin
    old_a := a;
    For j := 2 to Size-1 do
    begin
      a_1 := old_a[j-1];
      a0 := old_a[j];
      a1 := old_a[j+1];
      a[j] := (a_1 and (not a0) and (not a1))
        or ((not a_1) and a1)
    end;
    For j := 1 to Size do
      If (a[j]) Then  b[j] := '*'
                Else  b[j] := ' ';
    write(' Generation ',i:2);
    For j:=1 to Size do write(b[j]);
    writeln;
  End;

  writeln;
  ReadKey;
  DoneWinCrt

End.

{end of file sierpins.pas}