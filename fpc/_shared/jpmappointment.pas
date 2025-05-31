
unit jpmappointment;

interface

uses
  SysUtils, Math;

type
  trealarray = array of double;
  tmatrix = array of trealarray; // Dynamic array of dynamic arrays for 2D matrix (used for cost/satisfaction)

  // Enumerated type for clarity of job assignment status in fjobapplicantmatrix
  TJobAssignmentStatus = (
    jsUnassigned = 0,    // Cell is not marked or appointed
    jsAppointed = 1,     // Cell represents an optimal appointment
    jsCoveredZero = -1   // Cell contains a zero covered by an appointment in the same row/column
  );
  // Define a specific matrix type for job assignment status, using the enumerated type
  TJobAssignmentStatusMatrix = array of array of TJobAssignmentStatus;

  tappointmentresult = record
    applicantindex: integer; // 1-based index of the applicant
    jobindex: integer;     // 1-based index of the job
    satisfaction: double;  // Original satisfaction value for this appointment
  end;

  tappointmentresults = array of tappointmentresult; // Dynamic array of result records

  TAppointmentSolver = class
  private
    fproblem_size: integer; // Stores the number of jobs/applicants for the current problem instance
    fcostmatrix: tmatrix; // The regret/cost matrix, internally uses 1-based indexing (size fproblem_size x fproblem_size)
    fjobapplicantmatrix: TJobAssignmentStatusMatrix; // Stores appointment status (as per TJobAssignmentStatus), uses 1-based indexing
    frowmarked: array of boolean; // Boolean array to mark rows during matrix reduction steps
    fcolmarked: array of boolean; // Boolean array to mark columns during matrix reduction steps

    // Private helper methods, refactored from original procedures to operate on class fields
    procedure zeroes;
    procedure appoint;
    procedure mark;
    procedure sub_add;
    function calculatecurrentappointments: integer; // Counts current valid appointments (fjobapplicantmatrix[I,J]=jsAppointed)
    procedure clearappointmentdata; // Initializes fjobapplicantmatrix, frowmarked, fcolmarked before an appointment pass

  public
    constructor create;
    destructor destroy; override;

    // Main solver method: takes problem size, regret matrix, and original satisfaction matrix
    // Performs input validation for robustness.
    function solve(aproblem_size: integer; const acostmatrix: tmatrix; const asatisfactionmatrix: tmatrix): tappointmentresults;

    // Self-test method: demonstrates usage with static example data
    class procedure self_test;
  end;

implementation

constructor TAppointmentSolver.create;
begin
  inherited;
  // Dynamic arrays are automatically initialized to empty (nil).
  // Memory will be allocated dynamically when setlength is called within solve.
end;

destructor TAppointmentSolver.destroy;
begin
  // Explicitly set lengths to 0 to deallocate memory held by dynamic arrays
  setlength(fcostmatrix, 0); // Deallocates outer array, and then inner arrays are also deallocated.
  setlength(fjobapplicantmatrix, 0); // Deallocates outer array, and then inner arrays are also deallocated.
  setlength(frowmarked, 0);
  setlength(fcolmarked, 0);

  inherited destroy;
end;

procedure TAppointmentSolver.clearappointmentdata;
var
  i, j: integer;
begin
  // Allocate or resize dynamic arrays to accommodate 1-based indexing up to fproblem_size
  // +1 is for convenience to use 1..fproblem_size and ignore index 0, mimicking original logic.
  setlength(fjobapplicantmatrix, fproblem_size + 1, fproblem_size + 1);
  setlength(frowmarked, fproblem_size + 1);
  setlength(fcolmarked, fproblem_size + 1);

  // Initialize all elements.
  // frowmarked and fcolmarked arrays are used as flags for rows/columns.
  // fjobapplicantmatrix elements are initialized to jsUnassigned (0).
  for i := 0 to fproblem_size do // Loop through all allocated indices, including 0
  begin
    frowmarked[i] := false; // all rows are initially unmarked
    fcolmarked[i] := false; // all columns are initially unmarked
    for j := 0 to fproblem_size do
      fjobapplicantmatrix[i, j] := jsUnassigned; // jsunassigned (0) indicates no mark/appointment
  end;
end;

procedure TAppointmentSolver.zeroes;
var
  i, j: integer;
  xmin: double;
  haszero: boolean;
begin
  // step 1: row reduction
  // for each row, find the minimum element. if no zero exists, subtract this minimum from all elements in the row.
  for i := 1 to fproblem_size do
  begin
    xmin := fcostmatrix[i, 1]; // initialize xmin with the first element of the row
    haszero := false;          // flag to check if a zero already exists in the row

    // find min element and check for existing zeroes in the current row
    for j := 1 to fproblem_size do
    begin
      if fcostmatrix[i, j] = 0 then
        haszero := true; // a zero is found
      if fcostmatrix[i, j] < xmin then
        xmin := fcostmatrix[i, j]; // update minimum
    end;

    // if no zero was found, subtract xmin from all elements in the row to create at least one zero
    if not haszero then
    begin
      for j := 1 to fproblem_size do
        fcostmatrix[i, j] := fcostmatrix[i, j] - xmin;
    end;
  end;

  // step 2: column reduction
  // similar to row reduction, but applied to columns.
  for j := 1 to fproblem_size do
  begin
    xmin := fcostmatrix[1, j]; // initialize xmin with the first element of the column
    haszero := false;          // flag to check if a zero already exists in the column

    // find min element and check for existing zeroes in the current column
    for i := 1 to fproblem_size do
    begin
      if fcostmatrix[i, j] = 0 then
        haszero := true; // a zero is found
      if fcostmatrix[i, j] < xmin then
        xmin := fcostmatrix[i, j]; // update minimum
    end;

    // if no zero was found, subtract xmin from all elements in the column
    if not haszero then
    begin
      for i := 1 to fproblem_size do
        fcostmatrix[i, j] := fcostmatrix[i, j] - xmin;
    end;
  end;
end;

procedure TAppointmentSolver.appoint;
var
  i, j, k: integer;
  xcase: double; // stores the minimum number of zeroes found in a column for selection
  nz: integer;   // current count of zeroes in a column
  zi, zj: integer; // coordinates (row, col) of the zero selected for appointment
begin
  clearappointmentdata; // reset appointment status (fjobapplicantmatrix) and marking arrays

  // attempt to make initial appointments by selecting zeroes
  // this process iterates through rows, trying to find an optimal zero for appointment.
  for i := 1 to fproblem_size do // iterate through each applicant (row)
  begin
    xcase := maxdouble; // initialize with a very large value to find the minimum nz
    zi := -1; // initialize selected zero coordinates to invalid values
    zj := -1;

    for j := 1 to fproblem_size do // iterate through each job (column)
    begin
      // consider only current zeroes that are not yet "covered" or appointed
      if (fcostmatrix[i, j] = 0) and (fjobapplicantmatrix[i, j] = jsUnassigned) then
      begin
        nz := 0; // reset zero count for the current column

        // count how many zeroes are in the current column (including the current zero)
        // this is a heuristic to prioritize columns with fewer zeroes, aiming for unique assignments.
        for k := 1 to fproblem_size do
          if fcostmatrix[k, j] = 0 then
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
      fjobapplicantmatrix[zi, zj] := jsAppointed; // mark this zero as an official appointment

      // after making an appointment, mark other zeroes in the same row/column as 'covered'
      // these 'covered' zeroes (marked with jsCoveredZero) cannot be selected for other appointments in this pass.
      for k := 1 to fproblem_size do
      begin
        // mark other zeroes in the same column 'zj' as covered (-1)
        if (fcostmatrix[k, zj] = 0) and (fjobapplicantmatrix[k, zj] = jsUnassigned) then
          fjobapplicantmatrix[k, zj] := jsCoveredZero;
        // mark other zeroes in the same row 'zi' as covered (-1)
        if (fcostmatrix[zi, k] = 0) and (fjobapplicantmatrix[zi, k] = jsUnassigned) then
          fjobapplicantmatrix[zi, k] := jsCoveredZero;
      end;
    end;
  end;
end;

function TAppointmentSolver.calculatecurrentappointments: integer;
var
  i, j: integer;
  count: integer;
begin
  count := 0;
  // iterate through the fjobapplicantmatrix to count all entries marked as 'jsAppointed' (an appointment)
  for i := 1 to fproblem_size do
    for j := 1 to fproblem_size do
      if fjobapplicantmatrix[i, j] = jsAppointed then
        inc(count);
  result := count;
end;

procedure TAppointmentSolver.mark;
var
  i, j: integer;
  m_flag: boolean; // flag to indicate if any new mark was made in an iteration pass
  n_flag: boolean; // helper flag used to check if a row has an appointed zero
begin
  // loop until no new marks are made in a complete pass
  repeat
    m_flag := false; // reset the flag for changes for the current pass

    // rule 1: mark all rows which do not contain any assigned zero (fjobapplicantmatrix[I,J]=jsAppointed).
    for i := 1 to fproblem_size do
    begin
      n_flag := false; // assume no appointed zero in this row initially
      for j := 1 to fproblem_size do
      begin
        if fjobapplicantmatrix[i, j] = jsAppointed then
        begin
          n_flag := true; // an appointed zero is found in this row
          break; // no need to check further in this row
        end;
      end;
      // if no appointed zero was found in the row and the row is not already marked
      if (not n_flag) and (not frowmarked[i]) then
      begin
        frowmarked[i] := true; // mark the row
        m_flag := true;       // indicate that a new mark was made
      end;
    end;

    // rule 2: mark all columns that contain a 'covered' zero (fjobapplicantmatrix[I,J]=jsCoveredZero)
    //         if that zero is in a row that is already marked (frowmarked[I]).
    for j := 1 to fproblem_size do
      for i := 1 to fproblem_size do
        if (fjobapplicantmatrix[i, j] = jsCoveredZero) and frowmarked[i] and (not fcolmarked[j]) then
        begin
          fcolmarked[j] := true; // mark the column
          m_flag := true;        // indicate that a new mark was made
        end;

    // rule 3: mark all rows that contain an appointed zero (fjobapplicantmatrix[I,J]=jsAppointed)
    //         if that zero is in a column that is already marked (fcolmarked[J]).
    for i := 1 to fproblem_size do
      for j := 1 to fproblem_size do
        if (fjobapplicantmatrix[i, j] = jsAppointed) and fcolmarked[j] and (not frowmarked[i]) then
        begin
          frowmarked[i] := true; // mark the row
          m_flag := true;        // indicate that a new mark was made
        end;
  until not m_flag; // continue looping until no new marks are added in an entire pass
end;

procedure TAppointmentSolver.sub_add;
var
  i, j: integer;
  xmin: double;
begin
  xmin := maxdouble; // initialize with a very large value to find the minimum

  // step 1: find xmin
  // xmin is the minimum value among all elements that are in a marked row and an unmarked column.
  // this corresponds to elements covered by a single line (marked row).
  for i := 1 to fproblem_size do
    for j := 1 to fproblem_size do
    begin
      // only consider elements in marked rows and unmarked columns
      if frowmarked[i] and (not fcolmarked[j]) then
      begin
        if fcostmatrix[i, j] < xmin then
          xmin := fcostmatrix[i, j];
      end;
    end;

  // step 2: apply transformations to the cost matrix based on xmin
  for i := 1 to fproblem_size do
    for j := 1 to fproblem_size do
    begin
      // rule a: subtract xmin from elements in marked rows and unmarked columns.
      // these elements are covered by only the row line.
      if frowmarked[i] and (not fcolmarked[j]) then
        fcostmatrix[i, j] := fcostmatrix[i, j] - xmin;

      // rule b: add 2*xmin to elements in unmarked rows and marked columns.
      // these elements are covered by only the column line.
      // note: this is a specific variant of the hungarian method used in the original source,
      // where 2*xmin is added instead of just xmin for elements at the intersection of unmarked rows and marked columns.
      if (not frowmarked[i]) and fcolmarked[j] then
        fcostmatrix[i, j] := fcostmatrix[i, j] + 2 * xmin;

      // note: elements in unmarked rows and unmarked columns (uncovered) are not changed by this rule.
      // elements in marked rows and marked columns (intersections, covered by two lines) are also not changed by this rule set.
    end;
end;

function TAppointmentSolver.solve(aproblem_size: integer; const acostmatrix: tmatrix; const asatisfactionmatrix: tmatrix): tappointmentresults;
var
  i, j: integer;
  currentappointments: integer;
  result_array: tappointmentresults;
  current_result_index: integer;
begin
  // --- input validation ---
  // 1. check if aproblem_size is a positive integer
  if aproblem_size <= 0 then
    raise eargumentexception.create('problem size (aproblem_size) must be a positive integer.');

  // 2. check if input matrices are nil (empty)
  if (acostmatrix = nil) or (asatisfactionmatrix = nil) then
    raise eargumentexception.create('input matrices (acostmatrix and asatisfactionmatrix) cannot be empty (nil).');

  // 3. check dimensions of the input matrices.
  //    assuming 1-based indexing for input matrices as per the self-test and original code.
  //    so, length(matrix) should be aproblem_size + 1 (to include index aproblem_size),
  //    and length(matrix[i]) should be aproblem_size + 1.
  //    check the outer dimension (rows)
  if (length(acostmatrix) < aproblem_size + 1) then
    raise eargumentexception.createfmt('cost matrix (acostmatrix) row count mismatch. expected at least %d rows (for 1-based index up to %d), got %d.', [aproblem_size + 1, aproblem_size, length(acostmatrix)]);

  if (length(asatisfactionmatrix) < aproblem_size + 1) then
    raise eargumentexception.createfmt('satisfaction matrix (asatisfactionmatrix) row count mismatch. expected at least %d rows (for 1-based index up to %d), got %d.', [aproblem_size + 1, aproblem_size, length(asatisfactionmatrix)]);

  for i := 0 to aproblem_size do // check each inner array (column count)
  begin
    // it's possible for an inner array to be nil if not properly initialized, handle this.
    if (acostmatrix[i] = nil) then
      raise eargumentexception.createfmt('cost matrix (acostmatrix) row %d is nil.', [i]);
    if (length(acostmatrix[i]) < aproblem_size + 1) then
      raise eargumentexception.createfmt('cost matrix (acostmatrix) column count mismatch for row %d. expected at least %d columns (for 1-based index up to %d), got %d.', [i, aproblem_size + 1, aproblem_size, length(acostmatrix[i])]);

    if (asatisfactionmatrix[i] = nil) then
      raise eargumentexception.createfmt('satisfaction matrix (asatisfactionmatrix) row %d is nil.', [i]);
    if (length(asatisfactionmatrix[i]) < aproblem_size + 1) then
      raise eargumentexception.createfmt('satisfaction matrix (asatisfactionmatrix) column count mismatch for row %d. expected at least %d columns (for 1-based index up to %d), got %d.', [i, aproblem_size + 1, aproblem_size, length(asatisfactionmatrix[i])]);
  end;
  // --- end input validation ---

  // store the number of jobs/applicants for the current problem instance
  fproblem_size := aproblem_size;

  // initialize and copy the input regret/cost matrix. use 1-based indexing internally.
  setlength(fcostmatrix, fproblem_size + 1, fproblem_size + 1);
  for i := 1 to fproblem_size do
    for j := 1 to fproblem_size do
      fcostmatrix[i, j] := acostmatrix[i, j];

  // main loop of the hungarian algorithm implementation (adapted from original 'main' procedure)
  repeat
    zeroes; // step 1: reduce rows and columns to create zeroes
    appoint; // step 2: attempt initial assignments and mark covered zeroes

    currentappointments := calculatecurrentappointments; // get the number of assignments found

    // if the number of current appointments is less than fproblem_size,
    // it means an optimal assignment hasn't been found yet.
    // apply marking and matrix adjustment steps, then re-attempt appointments.
    if currentappointments <> fproblem_size then
    begin
      mark;    // step 3: mark rows and columns
      sub_add; // step 4: adjust the matrix using xmin based on marked lines
      appoint; // step 5: re-attempt appointments after matrix adjustment
    end;
    currentappointments := calculatecurrentappointments; // recalculate after the second 'appoint' or if no 'mark'/'sub_add' was needed
  until currentappointments = fproblem_size; // loop continues until fproblem_size optimal appointments are found

  // collect the final appointment results
  setlength(result_array, fproblem_size); // allocate space for fproblem_size results
  current_result_index := 0;
  for i := 1 to fproblem_size do
    for j := 1 to fproblem_size do
      if fjobapplicantmatrix[i, j] = jsAppointed then // check if this (i,j) pair is an appointed job
      begin
        // store the result: applicant index, job index, and original satisfaction
        result_array[current_result_index].applicantindex := i;
        result_array[current_result_index].jobindex := j;
        result_array[current_result_index].satisfaction := asatisfactionmatrix[i, j];
        inc(current_result_index); // move to the next slot in the result array
      end;

  result := result_array; // return the dynamic array of appointment results
end;

class procedure TAppointmentSolver.self_test;
var
  solver: TAppointmentSolver;
  regret_matrix: tmatrix;
  satisfaction_matrix: tmatrix;
  results: tappointmentresults;
  i, j: integer;
  total_satisfaction: double;
  test_problem_size: integer;
begin
  writeln('--- starting tappointmentsolver self-test ---');

  solver := tappointmentsolver.create; // create an instance of the solver
  try
    test_problem_size := 4; // define the problem size for the self-test

    // initialize the regret matrix (1-indexed to match problem description)
    // the problem example uses 4x4 matrix, so allocate 5x5 to allow 1..4 indexing
    setlength(regret_matrix, test_problem_size + 1, test_problem_size + 1);
    regret_matrix[1, 1] := 10; regret_matrix[1, 2] := 20; regret_matrix[1, 3] := 80; regret_matrix[1, 4] := 60;
    regret_matrix[2, 1] := 10; regret_matrix[2, 2] := 30; regret_matrix[2, 3] := 70; regret_matrix[2, 4] := 20;
    regret_matrix[3, 1] := 60; regret_matrix[3, 2] := 30; regret_matrix[3, 3] := 80; regret_matrix[3, 4] := 20;
    regret_matrix[4, 1] := 50; regret_matrix[4, 2] := 60; regret_matrix[4, 3] := 80; regret_matrix[4, 4] := 40;

    // initialize the original satisfaction matrix (1-indexed)
    setlength(satisfaction_matrix, test_problem_size + 1, test_problem_size + 1);
    satisfaction_matrix[1, 1] := 90; satisfaction_matrix[1, 2] := 80; satisfaction_matrix[1, 3] := 20; satisfaction_matrix[1, 4] := 40;
    satisfaction_matrix[2, 1] := 90; satisfaction_matrix[2, 2] := 70; satisfaction_matrix[2, 3] := 30; satisfaction_matrix[2, 4] := 80;
    satisfaction_matrix[3, 1] := 40; satisfaction_matrix[3, 2] := 70; satisfaction_matrix[3, 3] := 20; satisfaction_matrix[3, 4] := 80;
    satisfaction_matrix[4, 1] := 50; satisfaction_matrix[4, 2] := 40; satisfaction_matrix[4, 3] := 20; satisfaction_matrix[4, 4] := 60;

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

end.
