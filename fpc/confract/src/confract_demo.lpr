
program confract_demo;
{$mode objfpc}{$H+}
uses jpmcontinued;
begin
  jpmcontinued.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
