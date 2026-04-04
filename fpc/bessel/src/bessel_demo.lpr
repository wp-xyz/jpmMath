program bessel_demo;
{$mode objfpc}{$H+}
uses
  jpmbessel;
begin
  jpmbessel.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
