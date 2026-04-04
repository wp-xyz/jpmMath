program diffeq_demo;
{$mode objfpc}{$H+}

uses
  jpmdiffeq;

begin
  jpmdiffeq.self_test;

 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
