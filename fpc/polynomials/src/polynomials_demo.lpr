program polynomials_demo;
{$mode objfpc}{$H+}
uses
  jpmpolynomials;
begin
  jpmpolynomials.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
