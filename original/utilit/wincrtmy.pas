{*******************************************************}
{                                                       }
{       Turbo Pascal for Windows Runtime Library        }
{       Windows CRT Interface Unit                      }
{                                                       }
{       Copyright (c) 1992 Borland International        }
{                                                       }
{*******************************************************}
unit WinCrtMy;

{$S-}

interface

uses WinTypes, WinProcs, WinDos;

const
  WindowOrg: TPoint =                       { CRT window origin }
    (X: cw_UseDefault; Y: cw_UseDefault);
  WindowSize: TPoint =                      { CRT window size }
    (X: cw_UseDefault; Y: cw_UseDefault);
  ScreenSize: TPoint = (X: 80; Y: 30);      { Screen buffer dimensions }
  Cursor: TPoint = (X: 0; Y: 0);            { Cursor location }
  Origin: TPoint = (X: 0; Y: 0);            { Client area origin }
  InactiveTitle: PChar = '(Inactive %s)';   { Inactive window title }
  AutoTracking: Boolean = True;             { Track cursor on Write? }
  CheckEOF: Boolean = False;                { Allow Ctrl-Z for EOF? }
  CheckBreak: Boolean = True;               { Allow Ctrl-C for break? }

var
  WindowTitle: array[0..79] of Char;        { CRT window title }
  CrtWindow: HWnd;                          { CRT window handle }

procedure InitWinCrt;
procedure DoneWinCrt;

procedure WriteBuf(Buffer: PChar; Count: Word);
procedure WriteChar(Ch: Char);

function KeyPressed: Boolean;
function ReadKey: Char;
function ReadBuf(Buffer: PChar; Count: Word): Word;

procedure GotoXY(X, Y: Integer);
function WhereX: Integer;
function WhereY: Integer;
procedure ClrScr;
procedure ClrEol;

procedure CursorTo(X, Y: Integer);
procedure ScrollTo(X, Y: Integer);
procedure TrackCursor;

procedure AssignCrt(var F: Text);

implementation

{ Double word record }

type
  LongRec = record
    Lo, Hi: Integer;
  end;

{ MinMaxInfo array }

type
  PMinMaxInfo = ^TMinMaxInfo;
  TMinMaxInfo = array[0..4] of TPoint;

{ Scroll key definition record }

type
  TScrollKey = record
    Key: Byte;
    Ctrl: Boolean;
    SBar: Byte;
    Action: Byte;
  end;

{ CRT window procedure }

function CrtWinProc(Window: HWnd; Message, WParam: Word;
  LParam: Longint): Longint; export; forward;

{ CRT window class }

const
  CrtClass: TWndClass = (
    style: cs_HRedraw + cs_VRedraw;
    lpfnWndProc: @CrtWinProc;
    cbClsExtra: 0;
    cbWndExtra: 0;
    hInstance: 0;
    hIcon: 0;
    hCursor: 0;
    hbrBackground: 0;
    lpszMenuName: nil;
    lpszClassName: 'TPWinCrt');

const
  FirstLine: Integer = 0;               { First line in circular buffer }
  KeyCount: Integer = 0;                { Count of keys in KeyBuffer }
  Created: Boolean = False;       	{ CRT window created? }
  Focused: Boolean = False;             { CRT window focused? }
  Reading: Boolean = False;             { Reading from CRT window? }
  Painting: Boolean = False;            { Handling wm_Paint? }

var
  SaveExit: Pointer;                    { Saved exit procedure pointer }
  ScreenBuffer: PChar;                  { Screen buffer pointer }
  ClientSize: TPoint;                   { Client area dimensions }
  Range: TPoint;                        { Scroll bar ranges }
  CharSize: TPoint;                     { Character cell size }
  CharAscent: Integer;                  { Character ascent }
  DC: HDC;                              { Global device context }
  PS: TPaintStruct;                     { Global paint structure }
  SaveFont: HFont;                      { Saved device context font }
  KeyBuffer: array[0..63] of Char;      { Keyboard type-ahead buffer }

{ Scroll keys table }

const
  ScrollKeyCount = 12;
  ScrollKeys: array[1..ScrollKeyCount] of TScrollKey = (
    (Key: vk_Left;  Ctrl: False; SBar: sb_Horz; Action: sb_LineUp),
    (Key: vk_Right; Ctrl: False; SBar: sb_Horz; Action: sb_LineDown),
    (Key: vk_Left;  Ctrl: True;  SBar: sb_Horz; Action: sb_PageUp),
    (Key: vk_Right; Ctrl: True;  SBar: sb_Horz; Action: sb_PageDown),
    (Key: vk_Home;  Ctrl: False; SBar: sb_Horz; Action: sb_Top),
    (Key: vk_End;   Ctrl: False; SBar: sb_Horz; Action: sb_Bottom),
    (Key: vk_Up;    Ctrl: False; SBar: sb_Vert; Action: sb_LineUp),
    (Key: vk_Down;  Ctrl: False; SBar: sb_Vert; Action: sb_LineDown),
    (Key: vk_Prior; Ctrl: False; SBar: sb_Vert; Action: sb_PageUp),
    (Key: vk_Next;  Ctrl: False; SBar: sb_Vert; Action: sb_PageDown),
    (Key: vk_Home;  Ctrl: True;  SBar: sb_Vert; Action: sb_Top),
    (Key: vk_End;   Ctrl: True;  SBar: sb_Vert; Action: sb_Bottom));

{ Return the smaller of two integer values }

function Min(X, Y: Integer): Integer;
begin
  if X < Y then Min := X else Min := Y;
end;

{ Return the larger of two integer values }

function Max(X, Y: Integer): Integer;
begin
  if X > Y then Max := X else Max := Y;
end;

{ Allocate device context }

procedure InitDeviceContext;
begin
  if Painting then
    DC := BeginPaint(CrtWindow, PS) else
    DC := GetDC(CrtWindow);
  SaveFont := SelectObject(DC, GetStockObject(System_Fixed_Font));
end;

{ Release device context }

procedure DoneDeviceContext;
begin
  SelectObject(DC, SaveFont);
  if Painting then
    EndPaint(CrtWindow, PS) else
    ReleaseDC(CrtWindow, DC);
end;

{ Show caret }

procedure ShowCursor;
begin
  CreateCaret(CrtWindow, 0, CharSize.X, 2);
  SetCaretPos((Cursor.X - Origin.X) * CharSize.X,
    (Cursor.Y - Origin.Y) * CharSize.Y + CharAscent);
  ShowCaret(CrtWindow);
end;

{ Hide caret }

procedure HideCursor;
begin
  DestroyCaret;
end;

{ Update scroll bars }

procedure SetScrollBars;
begin
  SetScrollRange(CrtWindow, sb_Horz, 0, Max(1, Range.X), False);
  SetScrollPos(CrtWindow, sb_Horz, Origin.X, True);
  SetScrollRange(CrtWindow, sb_Vert, 0, Max(1, Range.Y), False);
  SetScrollPos(CrtWindow, sb_Vert, Origin.Y, True);
end;

{ Terminate CRT window }

procedure Terminate;
begin
  if Focused and Reading then HideCursor;
  Halt(255);
end;

{ Set cursor position }

procedure CursorTo(X, Y: Integer);
begin
  Cursor.X := Max(0, Min(X, ScreenSize.X - 1));
  Cursor.Y := Max(0, Min(Y, ScreenSize.Y - 1));
end;

{ Scroll window to given origin }

procedure ScrollTo(X, Y: Integer);
begin
  if Created then
  begin
    X := Max(0, Min(X, Range.X));
    Y := Max(0, Min(Y, Range.Y));
    if (X <> Origin.X) or (Y <> Origin.Y) then
    begin
      if X <> Origin.X then SetScrollPos(CrtWindow, sb_Horz, X, True);
      if Y <> Origin.Y then SetScrollPos(CrtWindow, sb_Vert, Y, True);
      ScrollWindow(CrtWindow,
	(Origin.X - X) * CharSize.X,
	(Origin.Y - Y) * CharSize.Y, nil, nil);
      Origin.X := X;
      Origin.Y := Y;
      UpdateWindow(CrtWindow);
    end;
  end;
end;

{ Scroll to make cursor visible }

procedure TrackCursor;
begin
  ScrollTo(Max(Cursor.X - ClientSize.X + 1, Min(Origin.X, Cursor.X)),
    Max(Cursor.Y - ClientSize.Y + 1, Min(Origin.Y, Cursor.Y)));
end;

{ Return pointer to location in screen buffer }

function ScreenPtr(X, Y: Integer): PChar;
begin
  Inc(Y, FirstLine);
  if Y >= ScreenSize.Y then Dec(Y, ScreenSize.Y);
  ScreenPtr := @ScreenBuffer[Y * ScreenSize.X + X];
end;

{ Update text on cursor line }

procedure ShowText(L, R: Integer);
begin
  if L < R then
  begin
    InitDeviceContext;
    TextOut(DC, (L - Origin.X) * CharSize.X,
      (Cursor.Y - Origin.Y) * CharSize.Y,
      ScreenPtr(L, Cursor.Y), R - L);
    DoneDeviceContext;
  end;
end;

{ Write text buffer to CRT window }

procedure WriteBuf(Buffer: PChar; Count: Word);
var
  L, R: Integer;

procedure NewLine;
begin
  ShowText(L, R);
  L := 0;
  R := 0;
  Cursor.X := 0;
  Inc(Cursor.Y);
  if Cursor.Y = ScreenSize.Y then
  begin
    Dec(Cursor.Y);
    Inc(FirstLine);
    if FirstLine = ScreenSize.Y then FirstLine := 0;
    FillChar(ScreenPtr(0, Cursor.Y)^, ScreenSize.X, ' ');
    ScrollWindow(CrtWindow, 0, -CharSize.Y, nil, nil);
    UpdateWindow(CrtWindow);
  end;
end;

begin
  InitWinCrt;
  L := Cursor.X;
  R := Cursor.X;
  while Count > 0 do
  begin
    case Buffer^ of
      #32..#255:
	begin
	  ScreenPtr(Cursor.X, Cursor.Y)^ := Buffer^;
	  Inc(Cursor.X);
	  if Cursor.X > R then R := Cursor.X;
	  if Cursor.X = ScreenSize.X then NewLine;
	end;
      #13:
	NewLine;
      #8:
	if Cursor.X > 0 then
	begin
	  Dec(Cursor.X);
	  ScreenPtr(Cursor.X, Cursor.Y)^ := ' ';
	  if Cursor.X < L then L := Cursor.X;
	end;
      #7:
        MessageBeep(0);
    end;
    Inc(Buffer);
    Dec(Count);
  end;
  ShowText(L, R);
  if AutoTracking then TrackCursor;
end;

{ Write character to CRT window }

procedure WriteChar(Ch: Char);
begin
  WriteBuf(@Ch, 1);
end;

{ Return keyboard status }

function KeyPressed: Boolean;
var
  M: TMsg;
begin
  InitWinCrt;
  while PeekMessage(M, 0, 0, 0, pm_Remove) do
  begin
    if M.Message = wm_Quit then Terminate;
    TranslateMessage(M);
    DispatchMessage(M);
  end;
  KeyPressed := KeyCount > 0;
end;

{ Read key from CRT window }

function ReadKey: Char;
begin
  TrackCursor;
  if not KeyPressed then
  begin
    Reading := True;
    if Focused then ShowCursor;
    repeat until KeyPressed;
    if Focused then HideCursor;
    Reading := False;
  end;
  ReadKey := KeyBuffer[0];
  Dec(KeyCount);
  Move(KeyBuffer[1], KeyBuffer[0], KeyCount);
end;

{ Read text buffer from CRT window }

function ReadBuf(Buffer: PChar; Count: Word): Word;
var
  Ch: Char;
  I: Word;
begin
  I := 0;
  repeat
    Ch := ReadKey;
    case Ch of
      #8:
	if I > 0 then
	begin
	  Dec(I);
	  WriteChar(#8);
	end;
      #32..#255:
	if I < Count - 2 then
	begin
	  Buffer[I] := Ch;
	  Inc(I);
	  WriteChar(Ch);
	end;
    end;
  until (Ch = #13) or (CheckEOF and (Ch = #26));
  Buffer[I] := Ch;
  Inc(I);
  if Ch = #13 then
  begin
    Buffer[I] := #10;
    Inc(I);
    WriteChar(#13);
  end;
  TrackCursor;
  ReadBuf := I;
end;

{ Set cursor position }

procedure GotoXY(X, Y: Integer);
begin
  CursorTo(X - 1, Y - 1);
end;

{ Return cursor X position }

function WhereX: Integer;
begin
  WhereX := Cursor.X + 1;
end;

{ Return cursor Y position }

function WhereY: Integer;
begin
  WhereY := Cursor.Y + 1;
end;

{ Clear screen }

procedure ClrScr;
begin
  InitWinCrt;
  FillChar(ScreenBuffer^, ScreenSize.X * ScreenSize.Y, ' ');
  Longint(Cursor) := 0;
  Longint(Origin) := 0;
  SetScrollBars;
  InvalidateRect(CrtWindow, nil, True);
  UpdateWindow(CrtWindow);
end;

{ Clear to end of line }

procedure ClrEol;
begin
  InitWinCrt;
  FillChar(ScreenPtr(Cursor.X, Cursor.Y)^, ScreenSize.X - Cursor.X, ' ');
  ShowText(Cursor.X, ScreenSize.X);
end;

{ wm_Create message handler }

procedure WindowCreate;
begin
  Created := True;
  GetMem(ScreenBuffer, ScreenSize.X * ScreenSize.Y);
  FillChar(ScreenBuffer^, ScreenSize.X * ScreenSize.Y, ' ');
  if not CheckBreak then
    EnableMenuItem(GetSystemMenu(CrtWindow, False), sc_Close,
      mf_Disabled + mf_Grayed);
end;

{ wm_Paint message handler }

procedure WindowPaint;
var
  X1, X2, Y1, Y2: Integer;
begin
  Painting := True;
  InitDeviceContext;
  X1 := Max(0, PS.rcPaint.left div CharSize.X + Origin.X);
  X2 := Min(ScreenSize.X,
    (PS.rcPaint.right + CharSize.X - 1) div CharSize.X + Origin.X);
  Y1 := Max(0, PS.rcPaint.top div CharSize.Y + Origin.Y);
  Y2 := Min(ScreenSize.Y,
    (PS.rcPaint.bottom + CharSize.Y - 1) div CharSize.Y + Origin.Y);
  while Y1 < Y2 do
  begin
    TextOut(DC, (X1 - Origin.X) * CharSize.X, (Y1 - Origin.Y) * CharSize.Y,
      ScreenPtr(X1, Y1), X2 - X1);
    Inc(Y1);
  end;
  DoneDeviceContext;
  Painting := False;
end;

{ wm_VScroll and wm_HScroll message handler }

procedure WindowScroll(Which, Action, Thumb: Integer);
var
  X, Y: Integer;

function GetNewPos(Pos, Page, Range: Integer): Integer;
begin
  case Action of
    sb_LineUp: GetNewPos := Pos - 1;
    sb_LineDown: GetNewPos := Pos + 1;
    sb_PageUp: GetNewPos := Pos - Page;
    sb_PageDown: GetNewPos := Pos + Page;
    sb_Top: GetNewPos := 0;
    sb_Bottom: GetNewPos := Range;
    sb_ThumbPosition: GetNewPos := Thumb;
  else
    GetNewPos := Pos;
  end;
end;

begin
  X := Origin.X;
  Y := Origin.Y;
  case Which of
    sb_Horz: X := GetNewPos(X, ClientSize.X div 2, Range.X);
    sb_Vert: Y := GetNewPos(Y, ClientSize.Y, Range.Y);
  end;
  ScrollTo(X, Y);
end;

{ wm_Size message handler }

procedure WindowResize(X, Y: Integer);
begin
  if Focused and Reading then HideCursor;
  ClientSize.X := X div CharSize.X;
  ClientSize.Y := Y div CharSize.Y;
  Range.X := Max(0, ScreenSize.X - ClientSize.X);
  Range.Y := Max(0, ScreenSize.Y - ClientSize.Y);
  Origin.X := Min(Origin.X, Range.X);
  Origin.Y := Min(Origin.Y, Range.Y);
  SetScrollBars;
  if Focused and Reading then ShowCursor;
end;

{ wm_GetMinMaxInfo message handler }

procedure WindowMinMaxInfo(MinMaxInfo: PMinMaxInfo);
var
  X, Y: Integer;
  Metrics: TTextMetric;
begin
  InitDeviceContext;
  GetTextMetrics(DC, Metrics);
  CharSize.X := Metrics.tmMaxCharWidth;
  CharSize.Y := Metrics.tmHeight + Metrics.tmExternalLeading;
  CharAscent := Metrics.tmAscent;
  X := Min(ScreenSize.X * CharSize.X + GetSystemMetrics(sm_CXVScroll),
    GetSystemMetrics(sm_CXScreen)) + GetSystemMetrics(sm_CXFrame) * 2;
  Y := Min(ScreenSize.Y * CharSize.Y + GetSystemMetrics(sm_CYHScroll) +
    GetSystemMetrics(sm_CYCaption), GetSystemMetrics(sm_CYScreen)) +
    GetSystemMetrics(sm_CYFrame) * 2;
  MinMaxInfo^[1].x := X;
  MinMaxInfo^[1].y := Y;
  MinMaxInfo^[3].x := CharSize.X * 16 + GetSystemMetrics(sm_CXVScroll) +
    GetSystemMetrics(sm_CXFrame) * 2;
  MinMaxInfo^[3].y := CharSize.Y * 4 + GetSystemMetrics(sm_CYHScroll) +
    GetSystemMetrics(sm_CYFrame) * 2 + GetSystemMetrics(sm_CYCaption);
  MinMaxInfo^[4].x := X;
  MinMaxInfo^[4].y := Y;
  DoneDeviceContext;
end;

{ wm_Char message handler }

procedure WindowChar(Ch: Char);
begin
  if CheckBreak and (Ch = #3) then Terminate;
  if KeyCount < SizeOf(KeyBuffer) then
  begin
    KeyBuffer[KeyCount] := Ch;
    Inc(KeyCount);
  end;
end;

{ wm_KeyDown message handler }

procedure WindowKeyDown(KeyDown: Byte);
var
  CtrlDown: Boolean;
  I: Integer;
begin
  if CheckBreak and (KeyDown = vk_Cancel) then Terminate;
  CtrlDown := GetKeyState(vk_Control) < 0;
  for I := 1 to ScrollKeyCount do
    with ScrollKeys[I] do
      if (Key = KeyDown) and (Ctrl = CtrlDown) then
      begin
	WindowScroll(SBar, Action, 0);
	Exit;
      end;
end;

{ wm_SetFocus message handler }

procedure WindowSetFocus;
begin
  Focused := True;
  if Reading then ShowCursor;
end;

{ wm_KillFocus message handler }

procedure WindowKillFocus;
begin
  if Reading then HideCursor;
  Focused := False;
end;

{ wm_Destroy message handler }

procedure WindowDestroy;
begin
  FreeMem(ScreenBuffer, ScreenSize.X * ScreenSize.Y);
  Longint(Cursor) := 0;
  Longint(Origin) := 0;
  PostQuitMessage(0);
  Created := False;
end;

{ CRT window procedure }

function CrtWinProc(Window: HWnd; Message, WParam: Word;
  LParam: Longint): Longint;
begin
  CrtWinProc := 0;
  CrtWindow := Window;
  case Message of
    wm_Create: WindowCreate;
    wm_Paint: WindowPaint;
    wm_VScroll: WindowScroll(sb_Vert, WParam, LongRec(LParam).Lo);
    wm_HScroll: WindowScroll(sb_Horz, WParam, LongRec(LParam).Lo);
    wm_Size: WindowResize(LongRec(LParam).Lo, LongRec(LParam).Hi);
    wm_GetMinMaxInfo: WindowMinMaxInfo(PMinMaxInfo(LParam));
    wm_Char: WindowChar(Char(WParam));
    wm_KeyDown: WindowKeyDown(Byte(WParam));
    wm_SetFocus: WindowSetFocus;
    wm_KillFocus: WindowKillFocus;
    wm_Destroy: WindowDestroy;
  else
    CrtWinProc := DefWindowProc(Window, Message, WParam, LParam);
  end;
end;

{ Text file device driver output function }

function CrtOutput(var F: TTextRec): Integer; far;
begin
  if F.BufPos <> 0 then
  begin
    WriteBuf(PChar(F.BufPtr), F.BufPos);
    F.BufPos := 0;
    KeyPressed;
  end;
  CrtOutput := 0;
end;

{ Text file device driver input function }

function CrtInput(var F: TTextRec): Integer; far;
begin
  F.BufEnd := ReadBuf(PChar(F.BufPtr), F.BufSize);
  F.BufPos := 0;
  CrtInput := 0;
end;

{ Text file device driver close function }

function CrtClose(var F: TTextRec): Integer; far;
begin
  CrtClose := 0;
end;

{ Text file device driver open function }

function CrtOpen(var F: TTextRec): Integer; far;
begin
  if F.Mode = fmInput then
  begin
    F.InOutFunc := @CrtInput;
    F.FlushFunc := nil;
  end else
  begin
    F.Mode := fmOutput;
    F.InOutFunc := @CrtOutput;
    F.FlushFunc := @CrtOutput;
  end;
  F.CloseFunc := @CrtClose;
  CrtOpen := 0;
end;

{ Assign text file to CRT device }

procedure AssignCrt(var F: Text);
begin
  with TTextRec(F) do
  begin
    Handle := $FFFF;
    Mode := fmClosed;
    BufSize := SizeOf(Buffer);
    BufPtr := @Buffer;
    OpenFunc := @CrtOpen;
    Name[0] := #0;
  end;
end;

{ Create CRT window if required }

procedure InitWinCrt;
begin
  if not Created then
  begin
    CrtWindow := CreateWindow(
      CrtClass.lpszClassName,
      WindowTitle,
      ws_OverlappedWindow + ws_HScroll + ws_VScroll,
      WindowOrg.X, WindowOrg.Y,
      WindowSize.X, WindowSize.Y,
      0,
      0,
      HInstance,
      nil);
    ShowWindow(CrtWindow, CmdShow);
    UpdateWindow(CrtWindow);
  end;
end;

{ Destroy CRT window if required }

procedure DoneWinCrt;
begin
  if Created then DestroyWindow(CrtWindow);
  Halt(0);
end;

{ WinCrt unit exit procedure }

procedure ExitWinCrt; far;
var
  P: PChar;
  Message: TMsg;
  Title: array[0..127] of Char;
begin
  ExitProc := SaveExit;
  if Created and (ErrorAddr = nil) then
  begin
    P := WindowTitle;
    WVSPrintF(Title, InactiveTitle, P);
    SetWindowText(CrtWindow, Title);
    EnableMenuItem(GetSystemMenu(CrtWindow, True), sc_Close, mf_Enabled);
    CheckBreak := False;
    while GetMessage(Message, 0, 0, 0) do
    begin
      TranslateMessage(Message);
      DispatchMessage(Message);
    end;
  end;
end;

begin
  if HPrevInst = 0 then
  begin
    CrtClass.hInstance := HInstance;
    CrtClass.hIcon := LoadIcon(0, idi_Application);
    CrtClass.hCursor := LoadCursor(0, idc_Arrow);
    CrtClass.hbrBackground := GetStockObject(White_Brush);
    RegisterClass(CrtClass);
  end;
  AssignCrt(Input);
  Reset(Input);
  AssignCrt(Output);
  Rewrite(Output);
  GetModuleFileName(HInstance, WindowTitle, SizeOf(WindowTitle));
  SaveExit := ExitProc;
  ExitProc := @ExitWinCrt;
end.