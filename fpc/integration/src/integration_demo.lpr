program integration_demo;
{$mode objfpc}{$H+}

uses
  jpmintegration;

begin
  jpmintegration.self_test;

 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
