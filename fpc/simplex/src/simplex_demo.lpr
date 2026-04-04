program simplex_demo;
{$mode objfpc}{$H+}

uses
  jpmsimplex;
begin
  jpmsimplex.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
