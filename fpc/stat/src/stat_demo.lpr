program stat_demo;
{$mode objfpc}{$H+}

uses
  jpmstats;

begin
  jpmstats.self_test;

 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
