unit jpmSort;
{$mode objfpc}{$H+}

{-------------------------------------------------------------------------------
                                 Unit jpmSort
         Sorting and searching routines based on Jean-Pierre Moreau's code.

  Provides: BubbleSort, InsertionSort, SelectionSort, QuickSort, MergeSort,
            HeapSort, ShellSort, BinarySearch, IntSort, RankSort.
-------------------------------------------------------------------------------}

interface

uses SysUtils, Math, jpmtypes;

procedure BubbleSort(var a: TFloatArray; n: integer);
procedure InsertionSort(var a: TFloatArray; n: integer);
procedure SelectionSort(var a: TFloatArray; n: integer);
procedure QuickSort(var a: TFloatArray; lo, hi: integer);
procedure MergeSort(var a: TFloatArray; n: integer);
procedure HeapSort(var a: TFloatArray; n: integer);
procedure ShellSort(var a: TFloatArray; n: integer);
function  BinarySearch(var a: TFloatArray; n: integer; key: Float): integer;
procedure IntSort(var a: TIntArray; n: integer);
procedure RankSort(var a: TFloatArray; n: integer; var rank: TIntArray);

procedure self_test;

implementation

{ ---- BubbleSort ----------------------------------------------------------- }
{ Standard O(n^2) bubble sort with early-exit optimisation. }
procedure BubbleSort(var a: TFloatArray; n: integer);
var
  i, j : integer;
  tmp  : Float;
  swapped : boolean;
begin
  for i := 0 to n - 2 do
  begin
    swapped := false;
    for j := 0 to n - 2 - i do
    begin
      if a[j] > a[j + 1] then
      begin
        tmp      := a[j];
        a[j]     := a[j + 1];
        a[j + 1] := tmp;
        swapped  := true
      end
    end;
    if not swapped then
      break
  end
end;

{ ---- InsertionSort --------------------------------------------------------- }
procedure InsertionSort(var a: TFloatArray; n: integer);
var
  i, j : integer;
  key  : Float;
begin
  for i := 1 to n - 1 do
  begin
    key := a[i];
    j   := i - 1;
    while (j >= 0) and (a[j] > key) do
    begin
      a[j + 1] := a[j];
      dec(j)
    end;
    a[j + 1] := key
  end
end;

{ ---- SelectionSort --------------------------------------------------------- }
procedure SelectionSort(var a: TFloatArray; n: integer);
var
  i, j, minIdx : integer;
  tmp          : Float;
begin
  for i := 0 to n - 2 do
  begin
    minIdx := i;
    for j := i + 1 to n - 1 do
      if a[j] < a[minIdx] then
        minIdx := j;
    if minIdx <> i then
    begin
      tmp        := a[i];
      a[i]       := a[minIdx];
      a[minIdx]  := tmp
    end
  end
end;

{ ---- QuickSort helpers ----------------------------------------------------- }
{ Median-of-three pivot: sort a[lo], a[mid], a[hi] and return a[mid]. }
function MedianOfThree(var a: TFloatArray; lo, hi: integer): Float;
var
  mid : integer;
  tmp : Float;
begin
  mid := lo + (hi - lo) div 2;
  if a[lo] > a[mid] then
  begin
    tmp   := a[lo];
    a[lo] := a[mid];
    a[mid]:= tmp
  end;
  if a[lo] > a[hi] then
  begin
    tmp   := a[lo];
    a[lo] := a[hi];
    a[hi] := tmp
  end;
  if a[mid] > a[hi] then
  begin
    tmp    := a[mid];
    a[mid] := a[hi];
    a[hi]  := tmp
  end;
  result := a[mid]
end;

{ ---- QuickSort ------------------------------------------------------------ }
procedure QuickSort(var a: TFloatArray; lo, hi: integer);
var
  i, j : integer;
  pivot, tmp : Float;
begin
  if lo >= hi then
    exit;
  pivot := MedianOfThree(a, lo, hi);
  i := lo;
  j := hi;
  repeat
    while a[i] < pivot do inc(i);
    while a[j] > pivot do dec(j);
    if i <= j then
    begin
      tmp  := a[i];
      a[i] := a[j];
      a[j] := tmp;
      inc(i);
      dec(j)
    end
  until i > j;
  if lo < j then QuickSort(a, lo, j);
  if i < hi then QuickSort(a, i, hi)
end;

{ ---- MergeSort helpers ----------------------------------------------------- }
procedure MergeSortRec(var a: TFloatArray; var tmp: TFloatArray; lo, hi: integer);
var
  mid, i, j, k : integer;
begin
  if lo >= hi then
    exit;
  mid := lo + (hi - lo) div 2;
  MergeSortRec(a, tmp, lo, mid);
  MergeSortRec(a, tmp, mid + 1, hi);
  { merge a[lo..mid] and a[mid+1..hi] into tmp, then copy back }
  i := lo;
  j := mid + 1;
  k := lo;
  while (i <= mid) and (j <= hi) do
  begin
    if a[i] <= a[j] then
    begin
      tmp[k] := a[i];
      inc(i)
    end
    else
    begin
      tmp[k] := a[j];
      inc(j)
    end;
    inc(k)
  end;
  while i <= mid do
  begin
    tmp[k] := a[i];
    inc(i); inc(k)
  end;
  while j <= hi do
  begin
    tmp[k] := a[j];
    inc(j); inc(k)
  end;
  for i := lo to hi do
    a[i] := tmp[i]
end;

{ ---- MergeSort ------------------------------------------------------------ }
procedure MergeSort(var a: TFloatArray; n: integer);
var
  tmp : TFloatArray;
begin
  setlength(tmp, n);
  MergeSortRec(a, tmp, 0, n - 1);
  setlength(tmp, 0)
end;

{ ---- HeapSort helpers ------------------------------------------------------ }
{ Sift down element at position i in a max-heap of size n. }
procedure SiftDown(var a: TFloatArray; i, n: integer);
var
  largest, left, right : integer;
  tmp                  : Float;
begin
  largest := i;
  left    := 2 * i + 1;
  right   := 2 * i + 2;
  if (left < n) and (a[left] > a[largest]) then
    largest := left;
  if (right < n) and (a[right] > a[largest]) then
    largest := right;
  if largest <> i then
  begin
    tmp        := a[i];
    a[i]       := a[largest];
    a[largest] := tmp;
    SiftDown(a, largest, n)
  end
end;

{ ---- HeapSort ------------------------------------------------------------- }
procedure HeapSort(var a: TFloatArray; n: integer);
var
  i   : integer;
  tmp : Float;
begin
  { build max-heap }
  for i := n div 2 - 1 downto 0 do
    SiftDown(a, i, n);
  { extract elements one by one }
  for i := n - 1 downto 1 do
  begin
    tmp  := a[0];
    a[0] := a[i];
    a[i] := tmp;
    SiftDown(a, 0, i)
  end
end;

{ ---- ShellSort ------------------------------------------------------------ }
{ Uses Knuth gap sequence: 1, 4, 13, 40, 121, ... }
procedure ShellSort(var a: TFloatArray; n: integer);
var
  gap, i, j : integer;
  tmp       : Float;
begin
  gap := 1;
  while gap < n div 3 do
    gap := gap * 3 + 1;
  while gap >= 1 do
  begin
    for i := gap to n - 1 do
    begin
      tmp := a[i];
      j   := i;
      while (j >= gap) and (a[j - gap] > tmp) do
      begin
        a[j] := a[j - gap];
        dec(j, gap)
      end;
      a[j] := tmp
    end;
    gap := gap div 3
  end
end;

{ ---- BinarySearch --------------------------------------------------------- }
function BinarySearch(var a: TFloatArray; n: integer; key: Float): integer;
var
  lo, hi, mid : integer;
begin
  lo := 0;
  hi := n - 1;
  result := -1;
  while lo <= hi do
  begin
    mid := lo + (hi - lo) div 2;
    if a[mid] = key then
    begin
      result := mid;
      exit
    end
    else if a[mid] < key then
      lo := mid + 1
    else
      hi := mid - 1
  end
end;

{ ---- IntSort helpers ------------------------------------------------------- }
procedure IntQuickSort(var a: TIntArray; lo, hi: integer);
var
  i, j, mid : integer;
  pivot, tmp : integer;
begin
  if lo >= hi then
    exit;
  mid   := lo + (hi - lo) div 2;
  pivot := a[mid];
  i     := lo;
  j     := hi;
  repeat
    while a[i] < pivot do inc(i);
    while a[j] > pivot do dec(j);
    if i <= j then
    begin
      tmp  := a[i];
      a[i] := a[j];
      a[j] := tmp;
      inc(i);
      dec(j)
    end
  until i > j;
  if lo < j then IntQuickSort(a, lo, j);
  if i < hi then IntQuickSort(a, i, hi)
end;

{ ---- IntSort -------------------------------------------------------------- }
procedure IntSort(var a: TIntArray; n: integer);
begin
  if n > 1 then
    IntQuickSort(a, 0, n - 1)
end;

{ ---- RankSort ------------------------------------------------------------- }
{ Indirect sort: fill rank[] so a[rank[0]] <= a[rank[1]] <= ...
  Uses insertion sort on the rank array (stable, simple). }
procedure RankSort(var a: TFloatArray; n: integer; var rank: TIntArray);
var
  i, j, r : integer;
begin
  setlength(rank, n);
  for i := 0 to n - 1 do
    rank[i] := i;
  for i := 1 to n - 1 do
  begin
    r := rank[i];
    j := i - 1;
    while (j >= 0) and (a[rank[j]] > a[r]) do
    begin
      rank[j + 1] := rank[j];
      dec(j)
    end;
    rank[j + 1] := r
  end
end;

{ ---- Helpers for self_test ------------------------------------------------ }
function ArraysEqual(var a, b: TFloatArray; n: integer): boolean;
var i : integer;
begin
  result := true;
  for i := 0 to n - 1 do
    if a[i] <> b[i] then
    begin
      result := false;
      exit
    end
end;

procedure PrintFloatArray(var a: TFloatArray; n: integer);
var i : integer;
begin
  for i := 0 to n - 1 do
    write(a[i]:4:1, ' ');
  writeln
end;

procedure PrintIntArray(var a: TIntArray; n: integer);
var i : integer;
begin
  for i := 0 to n - 1 do
    write(a[i]:3, ' ');
  writeln
end;

procedure ResetArray(var a: TFloatArray; const src: array of Float);
var i : integer;
begin
  for i := 0 to high(src) do
    a[i] := src[i]
end;

{ ---- self_test ------------------------------------------------------------- }
procedure self_test;
const
  INPUT : array[0..9] of Float = (5, 3, 8, 1, 9, 2, 7, 4, 6, 0);
  N     = 10;
var
  a, expected : TFloatArray;
  ia          : TIntArray;
  rank        : TIntArray;
  pass        : boolean;
  i           : integer;
begin
  writeln('=== jpmsort self_test ===');
  writeln;

  setlength(a, N);
  setlength(expected, N);
  for i := 0 to N - 1 do
    expected[i] := i;

  { --- Test 1: All sort algorithms ---------------------------------------- }
  writeln('Test 1 — Sort algorithms on [5,3,8,1,9,2,7,4,6,0]:');

  ResetArray(a, INPUT);
  BubbleSort(a, N);
  pass := ArraysEqual(a, expected, N);
  write('  BubbleSort:    '); PrintFloatArray(a, N);
  writeln('  Result: ', BoolToStr(pass, 'PASS', 'FAIL'));
  SelfTestCheck(pass, 'BubbleSort');

  ResetArray(a, INPUT);
  InsertionSort(a, N);
  pass := ArraysEqual(a, expected, N);
  write('  InsertionSort: '); PrintFloatArray(a, N);
  writeln('  Result: ', BoolToStr(pass, 'PASS', 'FAIL'));
  SelfTestCheck(pass, 'InsertionSort');

  ResetArray(a, INPUT);
  SelectionSort(a, N);
  pass := ArraysEqual(a, expected, N);
  write('  SelectionSort: '); PrintFloatArray(a, N);
  writeln('  Result: ', BoolToStr(pass, 'PASS', 'FAIL'));
  SelfTestCheck(pass, 'SelectionSort');

  ResetArray(a, INPUT);
  QuickSort(a, 0, N - 1);
  pass := ArraysEqual(a, expected, N);
  write('  QuickSort:     '); PrintFloatArray(a, N);
  writeln('  Result: ', BoolToStr(pass, 'PASS', 'FAIL'));
  SelfTestCheck(pass, 'QuickSort');

  ResetArray(a, INPUT);
  MergeSort(a, N);
  pass := ArraysEqual(a, expected, N);
  write('  MergeSort:     '); PrintFloatArray(a, N);
  writeln('  Result: ', BoolToStr(pass, 'PASS', 'FAIL'));
  SelfTestCheck(pass, 'MergeSort');

  ResetArray(a, INPUT);
  HeapSort(a, N);
  pass := ArraysEqual(a, expected, N);
  write('  HeapSort:      '); PrintFloatArray(a, N);
  writeln('  Result: ', BoolToStr(pass, 'PASS', 'FAIL'));
  SelfTestCheck(pass, 'HeapSort');

  ResetArray(a, INPUT);
  ShellSort(a, N);
  pass := ArraysEqual(a, expected, N);
  write('  ShellSort:     '); PrintFloatArray(a, N);
  writeln('  Result: ', BoolToStr(pass, 'PASS', 'FAIL'));
  SelfTestCheck(pass, 'ShellSort');

  writeln;

  { --- Test 2: BinarySearch ----------------------------------------------- }
  writeln('Test 2 — BinarySearch on sorted [0..9]:');
  for i := 0 to N - 1 do
    a[i] := i;
  i := BinarySearch(a, N, 7);
  writeln('  BinarySearch(7) = ', i, '  Result: ', BoolToStr(i = 7, 'PASS', 'FAIL'));
  SelfTestCheck(i = 7, 'BinarySearch(7)');
  i := BinarySearch(a, N, 11);
  writeln('  BinarySearch(11) = ', i, '  Result: ', BoolToStr(i = -1, 'PASS', 'FAIL'));
  SelfTestCheck(i = -1, 'BinarySearch(11) not found');
  writeln;

  { --- Test 3: IntSort ----------------------------------------------------- }
  writeln('Test 3 — IntSort on [5,3,8,1,9,2,7,4,6,0]:');
  setlength(ia, N);
  ia[0] := 5; ia[1] := 3; ia[2] := 8; ia[3] := 1; ia[4] := 9;
  ia[5] := 2; ia[6] := 7; ia[7] := 4; ia[8] := 6; ia[9] := 0;
  IntSort(ia, N);
  write('  IntSort result: '); PrintIntArray(ia, N);
  pass := true;
  for i := 0 to N - 1 do
    if ia[i] <> i then pass := false;
  writeln('  Result: ', BoolToStr(pass, 'PASS', 'FAIL'));
  SelfTestCheck(pass, 'IntSort');
  writeln;

  { --- Test 4: RankSort ---------------------------------------------------- }
  writeln('Test 4 — RankSort on [3.0, 1.0, 4.0, 1.0, 5.0]:');
  setlength(a, 5);
  a[0] := 3.0; a[1] := 1.0; a[2] := 4.0; a[3] := 1.0; a[4] := 5.0;
  RankSort(a, 5, rank);
  write('  rank indices: '); PrintIntArray(rank, 5);
  write('  sorted values via rank: ');
  for i := 0 to 4 do
    write(a[rank[i]]:4:1, ' ');
  writeln;
  { verify non-decreasing order via rank }
  pass := true;
  for i := 0 to 3 do
    if a[rank[i]] > a[rank[i + 1]] then pass := false;
  writeln('  Result: ', BoolToStr(pass, 'PASS', 'FAIL'));
  SelfTestCheck(pass, 'RankSort');
  writeln;

  writeln('=== self_test complete ===')
end;

end.
