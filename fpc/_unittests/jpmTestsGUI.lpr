program jpmTestsGUI;

{$mode objfpc}{$H+}

uses
  Interfaces, Forms, GuiTestRunner,
  jpmInterpolationTests, jpmIntegrationTests, jpmSpecialFuncsTests, 
  jpmAppointmentTests, jpmSelfTests;

{$R *.res}

begin
  Application.Initialize;
  Application.CreateForm(TGuiTestRunner, TestRunner);
  Application.Run;
end.

