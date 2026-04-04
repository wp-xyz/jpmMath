program matrices_demo;
{$mode objfpc}{$H+}
uses
  jpmmatrices;
begin
  jpmmatrices.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
