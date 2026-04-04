unit jpmSelfTests;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, fpcunit, testregistry;

type
  TTestSelfTests = class(TTestCase)
  published
    procedure TestAnneal;
    procedure TestArith;
    procedure TestBessel;
    procedure TestChebyshev;
    procedure TestComplex;
    procedure TestContinued;
    procedure TestDerivative;
    procedure TestDiffEq;
    procedure TestGeometry;
    procedure TestIntegration;
    procedure TestLstSqr;
    procedure TestMatrices;
    procedure TestOptimize;
    procedure TestPolynomials;
    procedure TestRoots;
    procedure TestSignal;
    procedure TestSimplex;
    procedure TestSort;
    procedure TestSpecial;
    procedure TestStats;
  end;

implementation

uses
  jpmTypes,
  jpmanneal, jpmarith, jpmbessel, jpmchebyshev, jpmcomplex,
  jpmcontinued, jpmderivative, jpmdiffeq, jpmgeometry, jpmintegration,
  jpmlstsqr, jpmmatrices, jpmoptimize, jpmpolynomials, jpmroots,
  jpmsignal, jpmsimplex, jpmsort, jpmspecial, jpmstats;

procedure TTestSelfTests.TestAnneal;      begin jpmanneal.self_test;     end;
procedure TTestSelfTests.TestArith;       begin jpmarith.self_test;      end;
procedure TTestSelfTests.TestBessel;      begin jpmbessel.self_test;     end;
procedure TTestSelfTests.TestChebyshev;   begin jpmchebyshev.self_test;  end;
procedure TTestSelfTests.TestComplex;     begin jpmcomplex.self_test;    end;
procedure TTestSelfTests.TestContinued;   begin jpmcontinued.self_test;  end;
procedure TTestSelfTests.TestDerivative;  begin jpmderivative.self_test; end;
procedure TTestSelfTests.TestDiffEq;      begin jpmdiffeq.self_test;     end;
procedure TTestSelfTests.TestGeometry;    begin jpmgeometry.self_test;   end;
procedure TTestSelfTests.TestIntegration; begin jpmintegration.self_test;end;
procedure TTestSelfTests.TestLstSqr;      begin jpmlstsqr.self_test;     end;
procedure TTestSelfTests.TestMatrices;    begin jpmmatrices.self_test;   end;
procedure TTestSelfTests.TestOptimize;    begin jpmoptimize.self_test;   end;
procedure TTestSelfTests.TestPolynomials; begin jpmpolynomials.self_test;end;
procedure TTestSelfTests.TestRoots;       begin jpmroots.self_test;      end;
procedure TTestSelfTests.TestSignal;      begin jpmsignal.self_test;     end;
procedure TTestSelfTests.TestSimplex;     begin jpmsimplex.self_test;    end;
procedure TTestSelfTests.TestSort;        begin jpmsort.self_test;       end;
procedure TTestSelfTests.TestSpecial;     begin jpmspecial.self_test;    end;
procedure TTestSelfTests.TestStats;       begin jpmstats.self_test;      end;

initialization
  RegisterTest('Self Tests', TTestSelfTests);

end.
