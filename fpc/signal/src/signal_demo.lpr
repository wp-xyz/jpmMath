program signal_demo;
{$mode objfpc}{$H+}
uses
  jpmsignal;
begin
  jpmsignal.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
