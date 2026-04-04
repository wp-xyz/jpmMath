program geometry_demo;
{$mode objfpc}{$H+}
uses
  {$IFDEF MSWINDOWS}Windows,{$ENDIF}
  jpmgeometry;
begin
  {$IFDEF MSWINDOWS}
   SetConsoleOutputCP(CP_UTF8);
  {$ENDIF}

  jpmgeometry.self_test;

 {$IFDEF MSWINDOWS}
  Write('Press ENTER to close...');
  ReadLn;
 {$ENDIF}
end.
