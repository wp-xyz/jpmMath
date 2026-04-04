program lstsqr_demo;
{$mode objfpc}{$H+}

uses
 {$IFDEF MSWINDOWS}Windows,{$ENDIF}
  jpmlstsqr;

begin
 {$IFDEF MSWINDOWS}
  SetConsoleOutputCP(CP_UTF8);
  SetMultiByteConversionCodePage(CP_UTF8);
  SetTextCodePage(Output, CP_UTF8);
 {$ENDIF}

  jpmlstsqr.self_test;

 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
