program arith_demo;
{$mode objfpc}{$H+}
uses jpmarith;
begin
  jpmarith.self_test;
  {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
  {$ENDIF}
end.
