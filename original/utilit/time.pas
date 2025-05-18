{*********************************************
*                UNIT  TIME                  *
* ------------------------------------------ *
* To display computing time on PC.           *
*                                            *
* INSTRUCTIONS:                              *
*                                            *
* After a call to procedure StartTiming and  *
* a call to StopTiming, the function Elapsed *
* returns the time between the two calls in  *
* seconds.                                   *
*                                            *
*********************************************}
Unit Time;
Interface

TYPE
  TTimeString = String[20];
VAR
  TickCount : LongInt ABSOLUTE $0040:$006C;
  Tstart, Ttime : LongInt;

  PROCEDURE StartTiming;
  PROCEDURE StopTiming;
  FUNCTION  Elapsed : TTimeString;

Implementation

  PROCEDURE StartTiming;
  BEGIN
    TStart := TickCount;
    {start at the beginning of a tick!}
    REPEAT UNTIL TStart <> TickCount;
    TStart := TickCount;
  END;

  PROCEDURE StopTiming;
  BEGIN TTime := TickCount - TStart; END;

  FUNCTION Elapsed : TTimeString;
  VAR Temp : TTimeString;
    Sec10 : LongInt;
  BEGIN
    Sec10 := TTime * 2470 DIV 4497;
    Str(Sec10:2,Temp);
    IF Temp[1] = ' ' THEN Temp[1] := '0';
    Inc(Temp[0]);
    Temp[length(Temp)] := Temp[pred(length(Temp))];
    Temp[pred(length(Temp))] := '.';
    Elapsed := Temp;
  END;

END.

{end of file time.pas}

