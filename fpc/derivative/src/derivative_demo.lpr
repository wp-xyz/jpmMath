program derivative_demo;

{$mode objfpc}{$H+}

uses jpmderivative;

begin
  jpmderivative.self_test;

 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
