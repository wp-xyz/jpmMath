{***************************************************
* Procedures to Read/Save a black and white screen *
* stored on disk.                                  *            
*                                                  *
*                   J-P Moreau   november 1994.    *
***************************************************}
UNIT SaveCrt;

INTERFACE

  USES WinProcs, WinTypes, Strings, Type_def;

  PROCEDURE WCrttoFile( P : HDC; name : string);
  PROCEDURE WLoadCrt  ( P : HDC; name : string);

IMPLEMENTATION

  {Save black and white screen to disk}
  PROCEDURE WCrttoFile(P:HDC;name:string);
  TYPE      PTab  = ^Table;
            Table = Array[0..300,0..200] of Byte;    {screen part < 64 Kb}
  VAR
    T         : FILE of Table;
    I,J,L,M   : Word;
    finx,finy : integer;
    Blanc     : LongInt;
    Pt        : Array[0..3] of PTab;
  BEGIN
    for i:=0 to 3 do
    begin
      New(Pt[i]);
      for j:=0 to 300 do
        for l:=0 to 200 do
          Pt[i]^[j][l]:=0
    end;
    Blanc:=RGB(255,255,255);
    FOR I:=0 TO 300 DO
      FOR J:=0 TO 200 DO
	IF GetPixel(P,I,J)<>Blanc THEN Pt[0]^[i][j]:=1;
    l:=301;
    IF MaxX < 600 THEN finx:=MaxX ELSE finx:=600;
    FOR I:=l TO finx DO
      FOR J:=0 TO 200 DO
	IF GetPixel(P,I,J)<>Blanc THEN Pt[1]^[i-l][j]:=1;
    m:=201;
    IF MaxY < 445 THEN finy:=MaxY-45 ELSE finy:=400;
    FOR I:=0 TO 300 DO
      FOR J:=m TO finy DO
	IF GetPixel(P,I,J)<>Blanc THEN Pt[2]^[i][j-m]:=1;
    FOR I:=l TO finx DO
      FOR J:=m TO finy DO
	IF GetPixel(P,I,J)<>Blanc THEN Pt[3]^[i-l][j-m]:=1;
    {IO-}
    Assign(T, name);
    ReWrite(T);
    {IO+}
    Write(T,Pt[0]^,Pt[1]^,Pt[2]^,Pt[3]^);
    Close(T);
    for i:=0 to 3 do Dispose(Pt[i]);
  END;

  {Load from disk black and white screen}
  PROCEDURE WLoadCrt(P:HDC;name:string);
  TYPE      PTab  = ^Table;
            Table = Array[0..300,0..200] of Byte;
  VAR
    T         : FILE of Table;
    I,J,L,M   : Word;
    finx,finy : integer;
    Bleu      : LongInt;
    Pt        : Array[0..3] of PTab;
  BEGIN
    for i:=0 to 3 do New(Pt[i]);
    {IO-}
    Assign(T, name);
    Reset(T);
    {IO+}
    if IOResult<>0 then
    begin
      MessageBeep(0);
      MessageBox(P,'File not found !',
			'ERROR',mb_Ok);
      exit
    end;
    Bleu:=RGB(0,0,255);
    Read(T,Pt[0]^,Pt[1]^,Pt[2]^,Pt[3]^);
    Close(T);
    FOR I:=0 TO 300 DO
      FOR J:=0 TO 200 DO
        IF (Pt[0]^[I][J]<>0) THEN SetPixel(P,I,J,Bleu);
    l:=301;
    IF MaxX < 600 THEN finx:=MaxX ELSE finx:=600;
    FOR I:=l TO finx DO
      FOR J:=0 TO 200 DO
        IF (Pt[1]^[I-l][J]<>0) THEN SetPixel(P,I,J,Bleu);
    m:=201;
    IF MaxY < 445 THEN finy:=MaxY-45 ELSE finy:=400;
    FOR I:=0 TO 300 DO
      FOR J:=m TO finy DO
        IF (Pt[2]^[I][J-m]<>0) THEN SetPixel(P,I,J,Bleu);
    FOR I:=l TO finx DO
      FOR J:=m TO finy DO
        IF (Pt[3]^[I-l][J-m]<>0) THEN SetPixel(P,I,J,Bleu);
    for i:=0 to 3 do Dispose(Pt[i])
  END;

END. {of file savecrt.pas}

