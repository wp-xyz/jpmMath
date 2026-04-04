program anneal_demo;
{$mode objfpc}{$H+}

uses
  jpmanneal;
begin
  jpmanneal.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close... ');
  ReadLn;
 {$ENDIF}
end.
