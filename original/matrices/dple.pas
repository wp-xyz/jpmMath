{**********************************************************************
* Solve AX = B using a partial pivoting algorithm and reduced storage *
* ------------------------------------------------------------------- *
* SAMPLE RUN:                                                         *
*                                                                     *
* System to solve:                                                    *
*   2.0000  -1.0000   1.0000   7.0000 -12.5400    5.0000              *
*   1.0000   5.0000  -2.0000  -8.0000 100.0000    1.0000              *
*   3.0000  -2.0000   3.0000  45.0000  27.3333    3.0000              *
*  11.0000   0.5500  -2.0000  -4.0000   1.0000    4.0000              *
*  33.0000   2.0000  -3.0000   5.0000   7.3333  -10.0000              *
*                                                                     *
* Solution is:                                                        *
*   2.11149597961869E+0000                                            *
*  -2.58290267820056E+0001                                            *
*   8.17423194407132E+0000                                            *
*  -2.52146730210577E+0000                                            *
*   1.24210363401706E+0000                                            *
*                                                                     *
* ------------------------------------------------------------------- *
* Ref.: "Wassyng, A. - Solving Ax = b: A method with reduced storage  *
*        requirements, SIAM J. Numerical Analysis, vol.19 (1982),     *
*        pp. 197-204".                                                *
*                                                                     *
*                               Pascal Release By J-P Moreau, Paris.  *
*                                        (www.jpmoreau.fr)            *
**********************************************************************}
PROGRAM Test_DPLE;                                                    

Uses WinCrt;

Const SIZE = 25;
      zero = 0.0;

Type  pVec = ^Vec;
      Vec = Array[1..SIZE] of Double;

Var   a, b, soln: pVec;
      error, temp: Double;
      i, j, ierr, row: Integer;


Procedure rowk(n, k: integer; r:pVec); Forward;

Procedure Swap(Var a,b: Double);
Var temp:Double;
Begin
  temp:=a;
  a:=b;
  b:=temp
End;

Procedure dple(n:Integer; b, c:pVec; Var ierr:Integer);
 
{ ******************************************************************
          SOLUTION OF LINEAR EQUATIONS WITH REDUCED STORAGE
  ******************************************************************

{ Uses the Henderson-Wassyng partial pivot algorithm.
  Wassyng, A. 'Solving Ax = b: A method with reduced storage requirements',
  SIAM J. Numerical Analysis, vol.19 (1982), pp. 197-204.

  The user must provide a routine ROWK to return the requested row of the
  matrix A.     }
Label 130, Return;
Var  bk, cj, dkj: Double;
     wk: pVec;
     i, iflag, ij, ijold, ik, j, k, kjold, km1, kp1,
     last, lastm1, lcol, lcolp1, m, maxwk, mjold, nm1, np1: Integer;
     iwk:Array[1..SIZE] of Integer;
Begin

  New(wk);

{ SET THE NECESSARY CONSTANTS }

  ierr := 0;
  maxwk := (n * n Div 4) + n + 3;
  np1 := n + 1;
  k := 1;
  iflag := -1;

{ GET THE FIRST COLUMN OF THE TRANSPOSED SYSTEM }

  rowk(n, 1, c);
  bk := b^[1];

  IF n <= 1 THEN
  begin
    IF c^[1] = zero then GOTO 130;
    c^[1] := bk / c^[1];
    goto Return
  end;

{ FIND THE PIVOT FOR COLUMN 1 }

  m := 1;
  For  i := 2 to n do
    IF ABS(c^[m]) < ABS(c^[i]) then m := i;

  iwk[1] := m;
  Swap(c^[m],c^[1]);

  IF c^[1] <> zero THEN
  begin
  
{ FIND THE FIRST ELEMENTARY MATRIX AND STORE IT IN D }
  
    For i := 2 to n do
      wk^[i-1] := -c^[i] / c^[1];
    wk^[n] := bk / c^[1];
  
{ K LOOP - EACH K FOR A NEW COLUMN OF THE TRANSPOSED SYSTEM }
  
    For k := 2 to n do
    begin
      kp1 := k + 1;
      km1 := k - 1;
    
{ GET COLUMN K }
    
      rowk(n, k, c);
      For j := 1 to km1 do
      begin
        m := iwk[j];
        Swap(c^[j],c^[m])
      end;
      bk := b^[k];
    
      iflag := -iflag;
      lcol := np1 - k;
      lcolp1 := lcol + 1;
      lastm1 := 1;
      last := maxwk - n + k;
      IF k <> 2 THEN
      begin
        lastm1 := maxwk - n + km1;
        IF iflag < 0 then last := last - n + k - 2;
        IF iflag > 0 then lastm1 := lastm1 - n + k - 3
      end;
    
{ J LOOP - EFFECT OF COLUMNS 1 TO K-1 OF L-INVERSE }
    
      For j := 1 to km1 do
      begin
        cj := c^[j];
        ij := (j-1) * lcolp1;
        IF j = km1 then ij := lastm1 - 1;
      
{ I LOOP - EFFECT OF L-INVERSE ON ROWS K TO N+1 }
      
        For i := k to n do
        begin
          ij := ij + 1;
          c^[i] := c^[i] + wk^[ij] * cj
        end;
        bk := bk - wk^[ij+1] * cj
      end;
    
{ K=N CASE }
    
      m := k;
      IF k >= n THEN
      begin
        IF c^[k] = zero then GOTO 130;
        wk^[last] := bk / c^[k]
      end
      ELSE
      begin

{ FIND THE PIVOT }
      
        For i := kp1 to n do
          IF ABS(c^[m]) < ABS(c^[i]) then m := i;
      
        iwk[k] := m;
        Swap(c^[m],c^[k]);
        IF c^[k] = zero then GOTO 130;
      
{ FIND THE K-TH ELEMENTARY MATRIX }
      
        ik := last;
        For i := kp1 to n do
        begin
          wk^[ik] := -c^[i] / c^[k];
          ik := ik + 1
        end;
        wk^[ik] := bk / c^[k]
      end;
    
{ FORM THE PRODUCT OF THE ELEMENTARY MATRICES }
    
      For j := 1 to km1 do
      begin
        kjold := j * lcolp1 + k - np1;
        mjold := kjold + m - k;
        ij := (j-1) * lcol;
        ijold := ij + j;
        IF j = km1 THEN
        begin
          kjold := lastm1;
          mjold := lastm1 + m - k;
          ijold := lastm1
        end;
      
        ik := last - 1;
        dkj := wk^[mjold];
        wk^[mjold] := wk^[kjold];
        For i := kp1 to np1 do
        begin
          ij := ij + 1;
          ijold := ijold + 1;
          ik := ik + 1;
          wk^[ij] := wk^[ijold] + wk^[ik] * dkj
        end
      end
    end; {k loop}
  
    last := maxwk;
    IF iflag < 0 then last := maxwk - 2;
    wk^[n] := wk^[last];
  
{ INSERT THE SOLUTION IN C }
  
    For i:=1 to n do c^[i] := wk^[i];
  
    nm1 := n - 1;
    For i := 1 to nm1 do
    begin
      k := n - i;
      m := iwk[k];
      Swap(c^[k],c^[m])
    end;
    goto RETURN
  end; {if C^[1]<>zero}

{ THE SYSTEM IS SINGULAR }

130: ierr := k;
Return: Dispose(wk)
END; {dple}


Procedure rowk(n, k: integer; r:pVec);
Var i:integer; temp:double;
Begin
  if k=1 then
  begin
    r^[1]:=2.0; r^[2]:=-1.0; r^[3]:=1.0; r^[4]:=7.0; r^[5]:=-12.54
  end
  else if k=2 then
  begin
    r^[1]:=1.0; r^[2]:=5.0; r^[3]:=-2.0; r^[4]:=-8.0; r^[5]:=100.0
  end
  else if k=3 then
  begin
    r^[1]:=3.0; r^[2]:=-2.0; r^[3]:=3.0; r^[4]:=45.0; r^[5]:=27.3333
  end
  else if k=4 then
  begin
    r^[1]:=11.0; r^[2]:=0.55; r^[3]:=-2.0; r^[4]:=-4.0; r^[5]:=1.0
  end
  else if k=5 then
  begin
    r^[1]:=33.0; r^[2]:=2.0; r^[3]:=-3.0; r^[4]:=5.0; r^[5]:=7.3333
  end
End; {rowk}


{main program}
BEGIN

  New(a); New(b); New(soln);

{ Define the right-hand side (n=5)  }
  b^[1]:=5.0; b^[2]:=1.0; b^[3]:=3.0; b^[4]:=4.0; b^[5]:=-10.0;

  writeln;
  writeln(' System to solve:');
  For i:=1 to 5 do
  begin
    rowk(5,i,soln);
    For j:=1 to 5 do write(' ',soln^[j]:8:4);
    writeln('  ',b^[i]:8:4)
  end;
  writeln;

  dple(5, b, soln, ierr);

  IF ierr = 0 THEN
  begin
    writeln(' Solution is:');
    For i:=1 to 5 do writeln('  ',soln^[i])
  end
  ELSE
    WRITE(' Error = ', ierr);

  ReadKey;
  Dispose(a); Dispose(b); Dispose(soln);
  DoneWinCrt

END.

{end of file dple.pas}