program special_demo;
{$mode objfpc}{$H+}

uses
  jpmspecial;

begin
  jpmspecial.self_test;

 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
