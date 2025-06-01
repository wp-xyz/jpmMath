{*******************************************************************************
                           Unit jpmAppointment

  Implements the Appointment method which aims to find the optimal assignment
  configuration among a set of possible assignments (Hungarian Algorithm).

  The code, basically, was extracted from Jean-Pierre Morau's "appoint" demo
  by AI.
*******************************************************************************}

unit jpmAppointment;

interface

uses
  SysUtils;

type
  TRealArray = array of double;
  TMatrix = array of TRealarray; // Dynamic array of dynamic arrays for 2D matrix (used for cost/satisfaction)

  // Enumerated type for clarity of job assignment status in fjobapplicantmatrix
  TJobAssignmentStatus = (
    jsUnassigned = 0,    // Cell is not marked or appointed
    jsAppointed = 1,     // Cell represents an optimal appointment
    jsCoveredZero = -1   // Cell contains a zero covered by an appointment in the same row/column
  );
  // Define a specific matrix type for job assignment status, using the enumerated type
  TJobAssignmentStatusMatrix = array of array of TJobAssignmentStatus;

  TAppointmentResult = record
    ApplicantIndex: integer; // 0-based index of the applicant
    JobIndex: integer;       // 0-based index of the job
    Satisfaction: double;    // Original satisfaction value for this appointment
  end;

  TAppointmentResults = array of TAppointmentResult; // Dynamic array of result records

  TAppointmentSolver = class
  private
    FProblemSize: integer; // Stores the number of jobs/applicants for the current problem instance
    FCostMatrix: TMatrix; // The regret/cost matrix
    FJobApplicantMatrix: TJobAssignmentStatusMatrix; // Stores appointment status (as per TJobAssignmentStatus)
    FRowMarked: array of boolean; // Boolean array to mark rows during matrix reduction steps
    FColMarked: array of boolean; // Boolean array to mark columns during matrix reduction steps

    // Private helper methods, refactored from original procedures to operate on class fields
    procedure Zeroes;
    procedure Appoint;
    procedure Mark;
    procedure SubAdd;
    function CalculateCurrentappointments: integer; // Counts current valid appointments (FJobApplicantMatrix[I,J]=jsAppointed)
    procedure ClearAppointmentData; // Initializes FJobApplicantMatrix, FRowMarked, FColMarked before an appointment pass

  public
    constructor Create;
    destructor Destroy; override;

    // Main solver method: takes problem size, regret matrix, and original satisfaction matrix
    // Performs input validation for robustness.
    function Solve(AProblemSize: integer; const ACostMatrix: TMatrix; const ASatisfactionMatrix: TMatrix): TAppointmentResults;

    { -- wp: moved to unit tests
    // Self-test method: demonstrates usage with static example data
    class procedure SelfTest;
    }
  end;

implementation

constructor TAppointmentSolver.Create;
begin
  inherited;
  // Dynamic arrays are automatically initialized to empty (nil).
  // Memory will be allocated dynamically when SetLength is called within solve.
end;

destructor TAppointmentSolver.Destroy;
begin
  // Explicitly set lengths to 0 to deallocate memory held by dynamic arrays
  SetLength(FCostMatrix, 0); // Deallocates outer array, and then inner arrays are also deallocated.
  SetLength(FJobApplicantMatrix, 0); // Deallocates outer array, and then inner arrays are also deallocated.
  SetLength(FRowMarked, 0);
  SetLength(FColMarked, 0);

  inherited Destroy;
end;

procedure TAppointmentSolver.ClearAppointmentData;
var
  i, j: integer;
begin
  // Allocate or resize dynamic arrays to accommodate 0-based indexing up to FProblemSize
  SetLength(FJobApplicantMatrix, FProblemSize, FProblemSize);
  SetLength(FRowMarked, FProblemSize);
  SetLength(FColMarked, FProblemSize);

  // Initialize all elements.
  // FRowMarked and FColMarked arrays are used as flags for rows/columns.
  // FJobApplicantmatrix elements are initialized to jsUnassigned (0).
  for i := 0 to FProblemSize-1 do // Loop through all allocated indices
  begin
    FRowMarked[i] := false; // all rows are initially unmarked
    FColMarked[i] := false; // all columns are initially unmarked
    for j := 0 to FProblemSize-1 do
      FJobApplicantMatrix[i, j] := jsUnassigned; // jsunassigned (0) indicates no mark/appointment
  end;
end;

procedure TAppointmentSolver.Zeroes;
var
  i, j: integer;
  xmin: double;
  hasZero: boolean;
begin
  // step 1: row reduction
  // for each row, find the minimum element. If no zero exists, subtract this minimum from all elements in the row.
  for i := 0 to FProblemSize-1 do
  begin
    xmin := fcostmatrix[i, 0]; // initialize xmin with the first element of the row
    hasZero := false;          // flag to check if a zero already exists in the row

    // find min element and check for existing zeroes in the current row
    for j := 0 to FProblemSize-1 do
    begin
      if FCostMatrix[i, j] = 0 then
        hasZero := true; // a zero is found
      if FCostMatrix[i, j] < xmin then
        xmin := FCostMatrix[i, j]; // update minimum
    end;

    // If no zero was found, subtract xmin from all elements in the row to create at least one zero
    if not haszero then
    begin
      for j := 0 to FProblemSize-1 do
        FCostMatrix[i, j] := FCostMatrix[i, j] - xmin;
    end;
  end;

  // step 2: column reduction
  // similar to row reduction, but applied to columns.
  for j := 0 to FProblemSize-1 do
  begin
    xmin := FCostMatrix[1, j]; // initialize xmin with the first element of the column
    hasZero := false;          // flag to check if a zero already exists in the column

    // find min element and check for existing zeroes in the current column
    for i := 0 to FProblemSize-1 do
    begin
      if FCostMatrix[i, j] = 0 then
        hasZero := true; // a zero is found
      if FCostMatrix[i, j] < xmin then
        xmin := FCostMatrix[i, j]; // update minimum
    end;

    // if no zero was found, subtract xmin from all elements in the column
    if not hasZero then
    begin
      for i := 0 to FProblemSize-1 do
        FCostMatrix[i, j] := FCostMatrix[i, j] - xmin;
    end;
  end;
end;

procedure TAppointmentSolver.Appoint;
var
  i, j, k: integer;
  xcase: double; // stores the minimum number of zeroes found in a column for selection
  nz: integer;   // current count of zeroes in a column
  zi, zj: integer; // coordinates (row, col) of the zero selected for appointment
begin
  ClearAppointmentData; // reset appointment status (FJobApplicantMatrix) and marking arrays

  // Attempt to make initial appointments by selecting zeroes
  // This process iterates through rows, trying to find an optimal zero for appointment.
  for i := 0 to FProblemSize-1 do // iterate through each applicant (row)
  begin
    xcase := 1E304; // Initialize with a very large value to find the minimum nz
    zi := -1; // Initialize selected zero coordinates to invalid values
    zj := -1;

    for j := 0 to FProblemSize-1 do // iterate through each job (column)
    begin
      // consider only current zeroes that are not yet "covered" or appointed
      if (FCostMatrix[i, j] = 0) and (FJobApplicantMatrix[i, j] = jsUnassigned) then
      begin
        nz := 0; // reset zero count for the current column

        // count how many zeroes are in the current column (including the current zero)
        // this is a heuristic to prioritize columns with fewer zeroes, aiming for unique assignments.
        for k := 0 to FProblemSize-1 do
          if FCostMatrix[k, j] = 0 then
            inc(nz);

        // select the zero that minimizes 'nz' (heuristic: prioritize zeroes in columns with fewer zeroes)
        if nz < xcase then
        begin
          xcase := nz;
          zi := i; // store row of selected zero
          zj := j; // store column of selected zero
        end;
      end;
    end;

    // if a suitable zero was found in the current row 'i' to make an appointment
    if (zi <> -1) and (zj <> -1) then
    begin
      FJobApplicantMatrix[zi, zj] := jsAppointed; // mark this zero as an official appointment

      // After making an appointment, mark other zeroes in the same row/column as 'covered'
      // these 'covered' zeroes (marked with jsCoveredZero) cannot be selected for other appointments in this pass.
      for k := 0 to FProblemSize-1 do
      begin
        // Mark other zeroes in the same column 'zj' as covered (-1)
        if (FCostMatrix[k, zj] = 0) and (FJobApplicantMatrix[k, zj] = jsUnassigned) then
          FJobApplicantMatrix[k, zj] := jsCoveredZero;
        // Mark other zeroes in the same row 'zi' as covered (-1)
        if (FCostMatrix[zi, k] = 0) and (FJobApplicantMatrix[zi, k] = jsUnassigned) then
          FJobApplicantMatrix[zi, k] := jsCoveredZero;
      end;
    end;
  end;
end;

function TAppointmentSolver.CalculateCurrentAppointments: integer;
var
  i, j: integer;
  count: integer;
begin
  count := 0;
  // Iterate through the FJobApplicantMatrix to count all entries marked as 'jsAppointed' (an appointment)
  for i := 0 to FProblemSize-1 do
    for j := 0 to FProblemSize-1 do
      if FJobApplicantMatrix[i, j] = jsAppointed then
        inc(count);
  Result := count;
end;

procedure TAppointmentSolver.Mark;
var
  i, j: integer;
  m_flag: boolean; // flag to indicate if any new mark was made in an iteration pass
  n_flag: boolean; // helper flag used to check if a row has an appointed zero
begin
  // loop until no new marks are made in a complete pass
  repeat
    m_flag := false; // reset the flag for changes for the current pass

    // rule 1: mark all rows which do not contain any assigned zero (FJobApplicantMatrix[I,J]=jsAppointed).
    for i := 0 to FProblemSize-1 do
    begin
      n_flag := false; // assume no appointed zero in this row initially
      for j := 0 to FProblemSize-1 do
      begin
        if FJobApplicantMatrix[i, j] = jsAppointed then
        begin
          n_flag := true; // an appointed zero is found in this row
          break; // no need to check further in this row
        end;
      end;
      // If no appointed zero was found in the row and the row is not already marked
      if (not n_flag) and (not FRowMarked[i]) then
      begin
        FRowMarked[i] := true; // mark the row
        m_flag := true;       // indicate that a new mark was made
      end;
    end;

    // rule 2: mark all columns that contain a 'covered' zero (fjobapplicantmatrix[I,J]=jsCoveredZero)
    //         if that zero is in a row that is already marked (frowmarked[I]).
    for j := 0 to FProblemSize-1 do
      for i := 0 to FProblemSize-1 do
        if (FJobApplicantMatrix[i, j] = jsCoveredZero) and FRowMarked[i] and (not FColMarked[j]) then
        begin
          FColMarked[j] := true; // mark the column
          m_flag := true;        // indicate that a new mark was made
        end;

    // rule 3: mark all rows that contain an appointed zero (fjobapplicantmatrix[I,J]=jsAppointed)
    //         if that zero is in a column that is already marked (fcolmarked[J]).
    for i := 0 to FProblemSize-1 do
      for j := 0 to FProblemSize-1 do
        if (FJobApplicantMatrix[i, j] = jsAppointed) and FColMarked[j] and (not FRowMarked[i]) then
        begin
          FRowMarked[i] := true; // mark the row
          m_flag := true;        // indicate that a new mark was made
        end;
  until not m_flag; // continue looping until no new marks are added in an entire pass
end;

procedure TAppointmentSolver.SubAdd;
var
  i, j: integer;
  xmin: double;
begin
  xmin := 1E304; // Initialize with a very large value to find the minimum

  // Step 1: find xmin
  // xmin is the minimum value among all elements that are in a marked row and an unmarked column.
  // this corresponds to elements covered by a single line (marked row).
  for i := 0 to FProblemSize-1 do
    for j := 0 to FProblemSize-1 do
    begin
      // only consider elements in marked rows and unmarked columns
      if FRowMarked[i] and (not FColMarked[j]) then
      begin
        if FCostMatrix[i, j] < xmin then
          xmin := FCostMatrix[i, j];
      end;
    end;

  // Step 2: apply transformations to the cost matrix based on xmin
  for i := 0 to FProblemSize-1 do
    for j := 0 to FProblemSize-1 do
    begin
      // rule a: subtract xmin from elements in marked rows and unmarked columns.
      // these elements are covered by only the row line.
      if FRowMarked[i] and (not FColMarked[j]) then
        FCostMatrix[i, j] := FCostMatrix[i, j] - xmin;

      // rule b: add 2*xmin to elements in unmarked rows and marked columns.
      // these elements are covered by only the column line.
      // note: this is a specific variant of the hungarian method used in the original source,
      // where 2*xmin is added instead of just xmin for elements at the intersection of unmarked rows and marked columns.
      if (not FRowMarked[i]) and FColMarked[j] then
        FCostMatrix[i, j] := FCostMatrix[i, j] + 2 * xmin;

      // note: elements in unmarked rows and unmarked columns (uncovered) are not changed by this rule.
      // elements in marked rows and marked columns (intersections, covered by two lines) are also not changed by this rule set.
    end;
end;

function TAppointmentSolver.Solve(AProblemSize: integer; const ACostMatrix: TMatrix;
  const ASatisfactionMatrix: TMatrix): TAppointmentResults;
var
  i, j: integer;
  CurrentAppointments: integer;
  CurrentResultIndex: integer;
begin
  // --- Input validation ---
  // 1. Check if AProblemSize is a positive integer
  if AProblemSize <= 0 then
    raise EArgumentException.Create('Problem Size must be a positive integer.');

  // 2. Check if input matrices are nil (empty)
  if (ACostMatrix = nil) or (ASatisfactionMatrix = nil) then
    raise EArgumentException.Create('Input matrices (ACostMatrix and ASatisfactionMatrix) cannot be empty (nil).');

  // 3. Check dimensions of the input matrices.
  //    Length(matrix) should be equal to AProblemSize, and Length(matrix[i]) as well.
  //    Check the outer dimension (rows)
  if (Length(ACostmatrix) < AProblemSize) then
    raise EArgumentException.CreateFmt('Cost matrix (ACostMatrix) row count mismatch. Expected at least %d rows, got %d.', [AProblemSize, Length(ACostMatrix)]);

  if (Length(ASatisfactionMatrix) < AProblemSize) then
    raise EArgumentException.CreateFmt('Satisfaction matrix (ASatisfactionMatrix) row count mismatch. Expected at least %d rows, got %d.', [AProblemSize, Length(ASatisfactionMatrix)]);

  for i := 0 to AProblemSize-1 do // check each inner array (column count)
  begin
    // It's possible for an inner array to be nil if not properly initialized, handle this.
    if (ACostMatrix[i] = nil) then
      raise EArgumentException.CreateFmt('Cost matrix (ACostmatrix) row %d is nil.', [i]);
    if (Length(ACostMatrix[i]) < AProblemSize) then
      raise EArgumentException.CreateFmt('Cost matrix (ACostMatrix) column count mismatch for row %d. Expected at least %d columns, got %d.', [i, AProblemSize, Length(ACostMatrix[i])]);

    if (ASatisfactionMatrix[i] = nil) then
      raise EArgumentException.CreateFmt('Satisfaction matrix (ASatisfactionMatrix) row %d is nil.', [i]);
    if (Length(ASatisfactionMatrix[i]) < AProblemSize) then
      raise EArgumentException.CreateFmt('Satisfaction matrix (ASatisfactionMatrix) column count mismatch for row %d. Expected at least %d columns (for 1-based index up to %d), got %d.', [i, AProblemSize, Length(ASatisfactionMatrix[i])]);
  end;
  // --- End input validation ---

  // Store the number of jobs/applicants for the current problem instance
  FProblemSize := AProblemSize;

  // Initialize and copy the input regret/cost matrix.
  SetLength(FCostmatrix, FProblemSize, FProblemSize);
  for i := 0 to FProblemSize-1 do
    for j := 0 to FProblemSize-1 do
      FCostMatrix[i, j] := ACostMatrix[i, j];

  // Main loop of the hungarian algorithm implementation (adapted from original 'main' procedure)
  repeat
    Zeroes; // step 1: reduce rows and columns to create zeroes
    Appoint; // step 2: attempt initial assignments and mark covered zeroes

    CurrentAppointments := CalculateCurrentappointments; // get the number of assignments found

    // If the number of current appointments is less than fproblem_size,
    // it means an optimal assignment hasn't been found yet.
    // Apply marking and matrix adjustment steps, then re-attempt appointments.
    if CurrentAppointments <> FProblemSize then
    begin
      Mark;    // step 3: mark rows and columns
      SubAdd ; // step 4: adjust the matrix using xmin based on marked lines
      Appoint; // step 5: re-attempt appointments after matrix adjustment
    end;
    CurrentAppointments := CalculateCurrentAppointments; // recalculate after the second 'appoint' or if no 'mark'/'sub_add' was needed
  until CurrentAppointments = FProblemSize; // loop continues until FProblemSize optimal appointments are found

  // Collect the final appointment results
  SetLength(Result, FProblemSize);  // allocate space for FProblemSize results
  CurrentResultIndex := 0;
  for i := 0 to FProblemSize-1 do
    for j := 0 to FProblemSize-1 do
      if FJobApplicantMatrix[i, j] = jsAppointed then // check if this (i,j) pair is an appointed job
      begin
        // Store the result: applicant index, job index, and original satisfaction
        Result[CurrentResultIndex].ApplicantIndex := i;
        Result[CurrentResultIndex].Jobindex := j;
        Result[CurrentResultIndex].Satisfaction := ASatisfactionMatrix[i, j];
        inc(CurrentResultIndex); // move to the next slot in the result array
      end;
end;

(*
class procedure TAppointmentSolver.SelfTest;
var
  solver: TAppointmentSolver;
  regret_matrix: TMatrix = nil;
  satisfaction_matrix: TMatrix = nil;
  results: TAppointmentResults;
  i, j: integer;
  total_satisfaction: double;
  test_problem_size: integer;
begin
  writeln('--- starting tappointmentsolver self-test ---');

  solver := TAppointmentSolver.Create; // create an instance of the solver
  try
    test_problem_size := 4; // define the problem size for the self-test

    // initialize the regret matrix
    // the problem example uses 4x4 matrix
    SetLength(regret_matrix, test_problem_size, test_problem_size);
    regret_matrix[0, 0] := 10; regret_matrix[0, 1] := 20; regret_matrix[0, 2] := 80; regret_matrix[0, 3] := 60;
    regret_matrix[1, 0] := 10; regret_matrix[1, 1] := 30; regret_matrix[1, 2] := 70; regret_matrix[1, 3] := 20;
    regret_matrix[2, 0] := 60; regret_matrix[2, 1] := 30; regret_matrix[2, 2] := 80; regret_matrix[2, 3] := 20;
    regret_matrix[3, 0] := 50; regret_matrix[3, 1] := 60; regret_matrix[3, 2] := 80; regret_matrix[3, 3] := 40;

    // initialize the original satisfaction matrix (1-indexed)
    SetLength(satisfaction_matrix, test_problem_size, test_problem_size);
    satisfaction_matrix[0, 0] := 90; satisfaction_matrix[0, 1] := 80; satisfaction_matrix[0, 2] := 20; satisfaction_matrix[0, 3] := 40;
    satisfaction_matrix[1, 0] := 90; satisfaction_matrix[1, 1] := 70; satisfaction_matrix[1, 2] := 30; satisfaction_matrix[1, 3] := 80;
    satisfaction_matrix[2, 0] := 40; satisfaction_matrix[2, 1] := 70; satisfaction_matrix[2, 2] := 20; satisfaction_matrix[2, 3] := 80;
    satisfaction_matrix[3, 0] := 50; satisfaction_matrix[3, 1] := 40; satisfaction_matrix[3, 2] := 20; satisfaction_matrix[3, 3] := 60;

    // run the solver with the defined matrices for np = 4
    results := solver.solve(test_problem_size, regret_matrix, satisfaction_matrix);

    // print the results in the format specified by the problem
    writeln('appointments:');
    total_satisfaction := 0;
    for i := low(results) to high(results) do
    begin
      writeln(format('      applicant #%d => job #%d  (s=%.0f)', [results[i].applicantindex, results[i].jobindex, results[i].satisfaction]));
      total_satisfaction := total_satisfaction + results[i].satisfaction;
    end;
    writeln;
    writeln(format('      total satisfaction: %.0f', [total_satisfaction]));

    // expected output from appoint.txt for verification:
    // trainee #1 ==> job #2  (s=80)
    // trainee #2 ==> job #1  (s=90)
    // trainee #3 ==> job #4  (s=80)
    // trainee #4 ==> job #3  (s=20)
    // total satisfaction: 270

  finally
    solver.free; // ensure the solver instance is freed
  end;
  writeln('--- self-test finished ---');
end;
*)

end.
