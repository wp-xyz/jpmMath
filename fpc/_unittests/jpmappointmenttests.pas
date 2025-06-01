unit jpmAppointmentTests;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, fpcunit, testutils, testregistry,
  jpmAppointment;

type

  TTestAppointment = class(TTestCase)
  published
    procedure TestAppointment;
  end;

implementation

procedure TTestAppointment.TestAppointment;
const
  EPS = 1E-9;
var
  Solver: TAppointmentSolver;
  RegretMatrix: TMatrix = nil;
  SatisfactionMatrix: TMatrix = nil;
  Results: tappointmentresults;
  i, j: integer;
  TotalSatisfaction: double;
  TestProblemSize: integer;
  ExpectedApplicantIndex: array of Integer = nil;
  ExpectedJobIndex: array of Integer = nil;
  ExpectedSatisfaction: array of Double = nil;
  ExpectedTotalSatisfaction: Integer;
  ExpectedResultsCount: Integer;
begin
  Solver := TAppointmentSolver.Create; // create an instance of the solver
  try
    TestProblemSize := 4; // define the problem size for the self-test

    // Initialize the regret matrix
    SetLength(RegretMatrix, TestProblemSize + 1, TestProblemSize + 1);
    RegretMatrix[0, 0] := 10; RegretMatrix[0, 1] := 20; RegretMatrix[0, 2] := 80; RegretMatrix[0, 3] := 60;
    RegretMatrix[1, 0] := 10; RegretMatrix[1, 1] := 30; RegretMatrix[1, 2] := 70; RegretMatrix[1, 3] := 20;
    RegretMatrix[2, 0] := 60; RegretMatrix[2, 1] := 30; RegretMatrix[2, 2] := 80; RegretMatrix[2, 3] := 20;
    RegretMatrix[3, 0] := 50; regretMatrix[3, 1] := 60; RegretMatrix[3, 2] := 80; RegretMatrix[3, 3] := 40;

    // Initialize the original satisfaction matrix
    Setlength(SatisfactionMatrix, TestProblemSize + 1, TestProblemSize + 1);
    SatisfactionMatrix[0, 0] := 90; SatisfactionMatrix[0, 1] := 80; SatisfactionMatrix[0, 2] := 20; SatisfactionMatrix[0, 3] := 40;
    SatisfactionMatrix[1, 0] := 90; SatisfactionMatrix[1, 1] := 70; SatisfactionMatrix[1, 2] := 30; SatisfactionMatrix[1, 3] := 80;
    SatisfactionMatrix[2, 0] := 40; SatisfactionMatrix[2, 1] := 70; SatisfactionMatrix[2, 2] := 20; SatisfactionMatrix[2, 3] := 80;
    SatisfactionMatrix[3, 0] := 50; SatisfactionMatrix[3, 1] := 40; SatisfactionMatrix[3, 2] := 20; SatisfactionMatrix[3, 3] := 60;

    // Tun the solver with the defined matrices for np = 4
    Results := Solver.Solve(TestProblemSize, RegretMatrix, SatisfactionMatrix);

    // expected output from appoint.txt for verification:
    // trainee #0 ==> job #1  (s=80)         // adapted for 0-based matrices (1-based in original demo)
    // trainee #1 ==> job #0  (s=90)
    // trainee #2 ==> job #3  (s=80)
    // trainee #3 ==> job #2  (s=20)
    // total satisfaction: 270

    // Prepare expected results
    ExpectedResultsCount := TestProblemSize;
    Setlength(ExpectedApplicantIndex, TestProblemSize);
    Setlength(ExpectedJobIndex, TestProblemSize);
    Setlength(ExpectedSatisfaction, TestProblemSize);
    ExpectedApplicantIndex[0] := 0;  ExpectedJobIndex[0] := 1;  ExpectedSatisfaction[0] := 80;
    ExpectedApplicantIndex[1] := 1;  ExpectedJobIndex[1] := 0;  ExpectedSatisfaction[1] := 90;
    ExpectedApplicantIndex[2] := 2;  ExpectedJobIndex[2] := 3;  ExpectedSatisfaction[2] := 80;
    ExpectedApplicantIndex[3] := 3;  ExpectedJobIndex[3] := 2;  ExpectedSatisfaction[3] := 20;
    ExpectedTotalSatisfaction := 270;

    CheckEquals(ExpectedResultsCount, Length(Results), 'Results count mismatch');
    TotalSatisfaction := 0;
    for i := Low(Results) to High(Results) do
    begin
      CheckEquals(ExpectedApplicantIndex[i], results[i].ApplicantIndex, '#' + IntToStr(i) + ': ApplicantIndex mismatch');
      CheckEquals(ExpectedJobIndex[i], results[i].JobIndex, '#' + IntToStr(i) + ': JobIndex mismatch');
      CheckEquals(ExpectedSatisfaction[i], results[i].Satisfaction, EPS, '#' + IntToStr(i) + ': Satisfaction mismatch');
      TotalSatisfaction := TotalSatisfaction + Results[i].Satisfaction;
    end;
    CheckEquals(ExpectedTotalSatisfaction, ExpectedTotalSatisfaction, 'Total satisfaction mismatch');

  finally
    Solver.Free; // Ensure the solver instance is freed
  end;
end;


initialization

  RegisterTest(TTestAppointment);
end.

