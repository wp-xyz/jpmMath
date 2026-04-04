program optimize_demo;
{$mode objfpc}{$H+}
uses
  jpmoptimize;
begin
  jpmoptimize.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
