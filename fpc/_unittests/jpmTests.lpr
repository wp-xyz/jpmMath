program jpmTests;

{$mode objfpc}{$H+}

uses
  Interfaces, Forms, GuiTestRunner, InterpolationTests;

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TGuiTestRunner, TestRunner);
  Application.Run;
end.

