program jpmTestsGUI;

{$mode objfpc}{$H+}

uses
  Interfaces, Forms, GuiTestRunner, jpmInterpolationTests;

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TGuiTestRunner, TestRunner);
  Application.Run;
end.

