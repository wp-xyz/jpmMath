program run_self_tests_console;
{$mode objfpc}{$H+}
uses
  SysUtils,
  jpmtypes,
  jpmanneal, jpmarith, jpmbessel, jpmchebyshev, jpmcomplex,
  jpmcontinued, jpmderivative, jpmdiffeq, jpmgeometry, jpmintegration,
  jpmlstsqr, jpmmatrices, jpmoptimize, jpmpolynomials, jpmroots,
  jpmsignal, jpmsimplex, jpmsort, jpmspecial, jpmstats;

var
  passed, failed: integer;

procedure RunTest(const name: string; proc: TProcedure);
begin
  try
    proc();
    Inc(passed);
    WriteLn('[PASS] ', name);
  except
    on E: Exception do
    begin
      Inc(failed);
      WriteLn('[FAIL] ', name, ': ', E.Message);
    end;
  end;
end;

begin
  passed := 0; failed := 0;
  WriteLn('=== jpmath self_test suite ===');
  WriteLn;
  RunTest('jpmanneal',      @jpmanneal.self_test);
  RunTest('jpmarith',       @jpmarith.self_test);
  RunTest('jpmbessel',      @jpmbessel.self_test);
  RunTest('jpmchebyshev',   @jpmchebyshev.self_test);
  RunTest('jpmcomplex',     @jpmcomplex.self_test);
  RunTest('jpmcontinued',   @jpmcontinued.self_test);
  RunTest('jpmderivative',  @jpmderivative.self_test);
  RunTest('jpmdiffeq',      @jpmdiffeq.self_test);
  RunTest('jpmgeometry',    @jpmgeometry.self_test);
  RunTest('jpmintegration', @jpmintegration.self_test);
  RunTest('jpmlstsqr',      @jpmlstsqr.self_test);
  RunTest('jpmmatrices',    @jpmmatrices.self_test);
  RunTest('jpmoptimize',    @jpmoptimize.self_test);
  RunTest('jpmpolynomials', @jpmpolynomials.self_test);
  RunTest('jpmroots',       @jpmroots.self_test);
  RunTest('jpmsignal',      @jpmsignal.self_test);
  RunTest('jpmsimplex',     @jpmsimplex.self_test);
  RunTest('jpmsort',        @jpmsort.self_test);
  RunTest('jpmspecial',     @jpmspecial.self_test);
  RunTest('jpmstats',       @jpmstats.self_test);
  WriteLn;
  WriteLn('Results: ', passed, ' passed, ', failed, ' failed');
  if failed > 0 then
    Halt(1);
end.
