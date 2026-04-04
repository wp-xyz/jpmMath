program roots_demo;

{$mode objfpc}{$H+}

uses
  jpmRoots;

begin
  jpmRoots.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
