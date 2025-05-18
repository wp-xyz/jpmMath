{************************************************
*    Dynamic table of size greater than 64 Ko   *
*        (used by program rocketv2.pas)         *
* --------------------------------------------- *
* The elements can have any type: REAL, DOUBLE, *
* INTEGER etc.                                  *
************************************************}
Unit LargArra;

(**) INTERFACE (**)
Uses WinCrtMy,Strings;

TYPE
  ByteArrayPtr = ^ByteArray;
  ByteArray = ARRAY[0..65520] OF Byte;

  MultArrayPtr = ^MultArray;
  MultArray = ARRAY[0..127] OF ByteArrayPtr;

  LgArrayPtr = ^LgArray;
  LgArray = OBJECT
    data          : MultArrayPtr;
    NumItems      : LongInt;
    ItemSize,
    NumPerArray,
    Shift,
    NumArrays,
    OneArraySize,
    LastArraySize : Word;
    CONSTRUCTOR Init(iNum : LongInt; iSize : Word);
    DESTRUCTOR Done;
    PROCEDURE ZeroOut;
    FUNCTION Index(N : LongInt) : Pointer;
  END;

(**) IMPLEMENTATION (**)
CONST
  Masks : ARRAY[0..15] OF Word =
   ((1 SHL  0) - 1, (1 SHL  1) - 1,
    (1 SHL  2) - 1, (1 SHL  3) - 1,
    (1 SHL  4) - 1, (1 SHL  5) - 1,
    (1 SHL  6) - 1, (1 SHL  7) - 1,
    (1 SHL  8) - 1, (1 SHL  9) - 1,
    (1 SHL 10) - 1, (1 SHL 11) - 1,
    (1 SHL 12) - 1, (1 SHL 13) - 1,
    (1 SHL 14) - 1, (1 SHL 15) - 1);

  CONSTRUCTOR LgArray.Init(iNum : LongInt; iSize : Word);
  VAR
    N        : Word;
    ShortMem : Boolean;
  BEGIN
    NumItems := iNum;
    ItemSize := iSize;
    NumPerArray := (65520 DIV ItemSize);
    N := 32768; Shift := 15;
    WHILE N > NumPerArray DO
      BEGIN
        N := N DIV 2;
        Dec(Shift);
      END;
    NumPerArray := N;
    NumArrays := Succ(NumItems SHR Shift);
    GetMem(data, NumArrays * SizeOf(Pointer));
    OneArraySize := NumPerArray * ItemSize;
    LastArraySize := (NumItems AND Masks[shift]) *
                      ItemSize;
    N := 0; ShortMem := FALSE;
    REPEAT
      IF N < Pred(NumArrays) THEN
        IF MaxAvail < OneArraySize THEN ShortMem := TRUE
        ELSE GetMem(data^[N], OneArraySize)
      ELSE
        IF MaxAvail < LastArraySize THEN ShortMem := TRUE
        ELSE GetMem(data^[N], LastArraySize);
      Inc(N);
    UNTIL (N >= NumArrays) OR ShortMem;
    IF ShortMem THEN
      BEGIN
        Dec(N);
        IF N > 0 THEN
          REPEAT
            Dec(N);
            FreeMem(Data^[N], OneArraySize);
          UNTIL N = 0;
        FreeMem(Data, NumArrays * SizeOf(Pointer));
        Fail;
      END;

  END;

  DESTRUCTOR LgArray.Done;
  VAR N : Byte;
  BEGIN
    N := Pred(NumArrays);
    FreeMem(Data^[N], LastArraySize);
    IF N > 0 THEN
      REPEAT
        Dec(N);
        FreeMem(Data^[N], OneArraySize);
      UNTIL N = 0;
    FreeMem(Data, NumArrays * SizeOf(Pointer));
  END;

  PROCEDURE LgArray.ZeroOut;
  VAR N : Word;
  BEGIN
    N := Pred(NumArrays);
    FillChar(Data^[N]^, LastArraySize, 0);
    IF N > 0 THEN
      REPEAT
        Dec(N);
        FillChar(Data^[N]^, OneArraySize, 0);
      UNTIL N = 0;
  END;

  FUNCTION LgArray.Index(N : LongInt) : Pointer;
  BEGIN
    Index := @Data^[N SHR Shift]^
             [(N AND masks[shift])*ItemSize];
  END;

END.

{end of file Largarra.pas}