program complex_demo;
{$mode objfpc}{$H+}
uses
  jpmcomplex;
begin
  jpmcomplex.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
