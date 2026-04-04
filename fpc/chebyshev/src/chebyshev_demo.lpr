program chebyshev_demo;
{$mode objfpc}{$H+}
uses jpmchebyshev;
begin
  jpmchebyshev.self_test;
 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
