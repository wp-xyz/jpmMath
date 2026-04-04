{ Linear algebra unit: LU, Gauss-Seidel, Cholesky, Tridiagonal, Jacobi }
{$mode objfpc}{$H+}

unit jpmMatrices;

interface

uses
  SysUtils, Math, jpmtypes;

{ --- 5.1 LU Decomposition --- }
procedure LUDecompose(var A: TMatrix; n: integer; var indx: TIntArray; var d: Float);
procedure LUSolve(var A: TMatrix; n: integer; var indx: TIntArray; var b: TFloatArray);
procedure LUInverse(var A: TMatrix; n: integer; var Ainv: TMatrix);
function  LUDeterminant(var A: TMatrix; n: integer): Float;

{ --- 5.2 Gauss-Seidel --- }
function GaussSeidel(var A: TMatrix; var b, x: TFloatArray; n: integer;
  maxIter: integer; tol: Float; var iter: integer): boolean;

{ --- 5.3 Cholesky --- }
function Cholesky(var A: TMatrix; n: integer; var L: TMatrix): boolean;
procedure CholeskySolve(var L: TMatrix; n: integer; var b, x: TFloatArray);

{ --- 5.4 Tridiagonal --- }
procedure TriDiagSolve(var a, b, c, r: TFloatArray; n: integer; var u: TFloatArray);

{ --- 5.5 Jacobi eigenvalues --- }
procedure JacobiEigen(var A: TMatrix; n: integer; var eigenvalues: TFloatArray;
  var eigenvectors: TMatrix; var nRot: integer);

{ --- 5.6 Helpers --- }
procedure MatVecMul(var A: TMatrix; var x: TFloatArray; n: integer; var res: TFloatArray);
procedure MatMatMul(var A, B: TMatrix; n: integer; var C: TMatrix);
procedure PrintMatrix(var A: TMatrix; n: integer; name: string);

procedure self_test;

implementation

{ ======================================================================
  5.1  LU DECOMPOSITION  (Crout's method with partial pivoting)
  0-based indexing.  A is replaced in-place by its LU factors.
  d = +1 or -1 (parity of row swaps).
  ====================================================================== }
procedure LUDecompose(var A: TMatrix; n: integer; var indx: TIntArray; var d: Float);
const
  TINY = 1.5e-16;
var
  i, imax, j, k : integer;
  amax, dum, sum : Float;
  vv : TFloatArray;
begin
  SetLength(vv, n);
  d := 1.0;
  for i := 0 to n - 1 do
  begin
    amax := 0.0;
    for j := 0 to n - 1 do
      if Abs(A[i][j]) > amax then
        amax := Abs(A[i][j]);
    if amax < TINY then
      amax := TINY;
    vv[i] := 1.0 / amax;
  end;
  for j := 0 to n - 1 do
  begin
    for i := 0 to j - 1 do
    begin
      sum := A[i][j];
      for k := 0 to i - 1 do
        sum := sum - A[i][k] * A[k][j];
      A[i][j] := sum;
    end;
    amax := 0.0;
    imax := j;
    for i := j to n - 1 do
    begin
      sum := A[i][j];
      for k := 0 to j - 1 do
        sum := sum - A[i][k] * A[k][j];
      A[i][j] := sum;
      dum := vv[i] * Abs(sum);
      if dum >= amax then
      begin
        amax := dum;
        imax := i;
      end;
    end;
    if j <> imax then
    begin
      for k := 0 to n - 1 do
      begin
        dum := A[imax][k];
        A[imax][k] := A[j][k];
        A[j][k] := dum;
      end;
      d := -d;
      vv[imax] := vv[j];
    end;
    indx[j] := imax;
    if Abs(A[j][j]) < TINY then
      A[j][j] := TINY;
    if j <> n - 1 then
    begin
      dum := 1.0 / A[j][j];
      for i := j + 1 to n - 1 do
        A[i][j] := A[i][j] * dum;
    end;
  end;
end;

{ Forward- then back-substitution using in-place LU }
procedure LUSolve(var A: TMatrix; n: integer; var indx: TIntArray; var b: TFloatArray);
var
  i, ii, j, ip : integer;
  sum : Float;
begin
  ii := -1;
  for i := 0 to n - 1 do
  begin
    ip  := indx[i];
    sum := b[ip];
    b[ip] := b[i];
    if ii >= 0 then
      for j := ii to i - 1 do
        sum := sum - A[i][j] * b[j]
    else if sum <> 0.0 then
      ii := i;
    b[i] := sum;
  end;
  for i := n - 1 downto 0 do
  begin
    sum := b[i];
    for j := i + 1 to n - 1 do
      sum := sum - A[i][j] * b[j];
    b[i] := sum / A[i][i];
  end;
end;

{ Matrix inverse via LU: Ainv is set to A^(-1) }
procedure LUInverse(var A: TMatrix; n: integer; var Ainv: TMatrix);
var
  i, j    : integer;
  d       : Float;
  col     : TFloatArray;
  indx    : TIntArray;
  LU      : TMatrix;
begin
  { copy A into LU }
  SetLength(LU, n);
  for i := 0 to n - 1 do
  begin
    SetLength(LU[i], n);
    for j := 0 to n - 1 do
      LU[i][j] := A[i][j];
  end;
  SetLength(indx, n);
  SetLength(col,  n);
  SetLength(Ainv, n);
  for i := 0 to n - 1 do
    SetLength(Ainv[i], n);
  LUDecompose(LU, n, indx, d);
  for j := 0 to n - 1 do
  begin
    for i := 0 to n - 1 do
      col[i] := 0.0;
    col[j] := 1.0;
    LUSolve(LU, n, indx, col);
    for i := 0 to n - 1 do
      Ainv[i][j] := col[i];
  end;
end;

{ Determinant via LU }
function LUDeterminant(var A: TMatrix; n: integer): Float;
var
  i       : integer;
  d       : Float;
  indx    : TIntArray;
  LU      : TMatrix;
  j       : integer;
begin
  SetLength(LU, n);
  for i := 0 to n - 1 do
  begin
    SetLength(LU[i], n);
    for j := 0 to n - 1 do
      LU[i][j] := A[i][j];
  end;
  SetLength(indx, n);
  LUDecompose(LU, n, indx, d);
  for i := 0 to n - 1 do
    d := d * LU[i][i];
  result := d;
end;

{ ======================================================================
  5.2  GAUSS-SEIDEL ITERATION
  Solves A*x = b.  Returns true if converged within maxIter iterations.
  x must be pre-initialised (zero is fine).
  ====================================================================== }
function GaussSeidel(var A: TMatrix; var b, x: TFloatArray; n: integer;
  maxIter: integer; tol: Float; var iter: integer): boolean;
var
  i, j   : integer;
  xnew, sigma, diff, maxdiff : Float;
begin
  iter := 0;
  result := false;
  repeat
    maxdiff := 0.0;
    for i := 0 to n - 1 do
    begin
      sigma := 0.0;
      for j := 0 to n - 1 do
        if j <> i then
          sigma := sigma + A[i][j] * x[j];
      if Abs(A[i][i]) < 1e-20 then
        exit;
      xnew := (b[i] - sigma) / A[i][i];
      diff := Abs(xnew - x[i]);
      if diff > maxdiff then
        maxdiff := diff;
      x[i] := xnew;
    end;
    Inc(iter);
    if maxdiff < tol then
    begin
      result := true;
      exit;
    end;
  until iter >= maxIter;
end;

{ ======================================================================
  5.3  CHOLESKY DECOMPOSITION
  Factors symmetric positive-definite A = L * L^T.
  L is lower-triangular and stored separately.
  Returns false if A is not positive definite.
  ====================================================================== }
function Cholesky(var A: TMatrix; n: integer; var L: TMatrix): boolean;
const
  EPSQUAD = 1e-30;
var
  i, j, k : integer;
  sum      : Float;
begin
  SetLength(L, n);
  for i := 0 to n - 1 do
  begin
    SetLength(L[i], n);
    for j := 0 to n - 1 do
      L[i][j] := 0.0;
  end;
  result := true;
  for i := 0 to n - 1 do
  begin
    for j := 0 to i do
    begin
      sum := A[i][j];
      for k := 0 to j - 1 do
        sum := sum - L[i][k] * L[j][k];
      if i = j then
      begin
        if sum < EPSQUAD then
        begin
          result := false;
          exit;
        end;
        L[i][j] := Sqrt(sum);
      end
      else
        L[i][j] := sum / L[j][j];
    end;
  end;
end;

{ Solve A*x = b using Cholesky factor L (A = L*L^T).
  Forward substitution: L*y = b, then back substitution: L^T*x = y. }
procedure CholeskySolve(var L: TMatrix; n: integer; var b, x: TFloatArray);
var
  i, k : integer;
  sum  : Float;
  y    : TFloatArray;
begin
  SetLength(y, n);
  { Forward: L * y = b }
  for i := 0 to n - 1 do
  begin
    sum := b[i];
    for k := 0 to i - 1 do
      sum := sum - L[i][k] * y[k];
    y[i] := sum / L[i][i];
  end;
  { Back: L^T * x = y }
  SetLength(x, n);
  for i := n - 1 downto 0 do
  begin
    sum := y[i];
    for k := i + 1 to n - 1 do
      sum := sum - L[k][i] * x[k];
    x[i] := sum / L[i][i];
  end;
end;

{ ======================================================================
  5.4  TRIDIAGONAL SOLVER  (Thomas algorithm)
  a = sub-diagonal  (a[0] unused)
  b = main diagonal
  c = super-diagonal (c[n-1] unused)
  r = right-hand side
  u = solution  (0-based, length n)
  ====================================================================== }
procedure TriDiagSolve(var a, b, c, r: TFloatArray; n: integer; var u: TFloatArray);
var
  j    : integer;
  bet  : Float;
  gam  : TFloatArray;
begin
  SetLength(gam, n);
  SetLength(u, n);
  if Abs(b[0]) < 1e-30 then
    exit;
  bet    := b[0];
  u[0]   := r[0] / bet;
  for j := 1 to n - 1 do
  begin
    gam[j] := c[j - 1] / bet;
    bet    := b[j] - a[j] * gam[j];
    if Abs(bet) < 1e-30 then
      exit;
    u[j] := (r[j] - a[j] * u[j - 1]) / bet;
  end;
  for j := n - 2 downto 0 do
    u[j] := u[j] - gam[j + 1] * u[j + 1];
end;

{ ======================================================================
  5.5  JACOBI EIGENVALUE METHOD  (cyclic Jacobi for symmetric matrices)
  A is modified in place (becomes diagonal).
  eigenvectors[i] = i-th eigenvector (row i of the accumulation matrix).
  nRot = total number of rotations performed.
  0-based indexing.
  ====================================================================== }
procedure JacobiEigen(var A: TMatrix; n: integer; var eigenvalues: TFloatArray;
  var eigenvectors: TMatrix; var nRot: integer);
const
  MAXITER = 50;
  TOL     = 1.0e-10;
var
  iter, p, q, i : integer;
  theta, t, c, s, tau : Float;
  g, h, aOffDiag : Float;
  done : boolean;
  b, z : TFloatArray;
begin
  { initialise eigenvectors to identity }
  SetLength(eigenvectors, n);
  for i := 0 to n - 1 do
  begin
    SetLength(eigenvectors[i], n);
    for p := 0 to n - 1 do
      eigenvectors[i][p] := 0.0;
    eigenvectors[i][i] := 1.0;
  end;
  SetLength(eigenvalues, n);
  SetLength(b, n);
  SetLength(z, n);
  for i := 0 to n - 1 do
  begin
    eigenvalues[i] := A[i][i];
    b[i]           := A[i][i];
    z[i]           := 0.0;
  end;
  nRot := 0;
  for iter := 1 to MAXITER do
  begin
    { sum of off-diagonal elements }
    aOffDiag := 0.0;
    for p := 0 to n - 2 do
      for q := p + 1 to n - 1 do
        aOffDiag := aOffDiag + Abs(A[p][q]);
    if aOffDiag < TOL then
      exit;
    done := true;
    for p := 0 to n - 2 do
    begin
      for q := p + 1 to n - 1 do
      begin
        g := 100.0 * Abs(A[p][q]);
        if (iter > 4) and
           (Abs(eigenvalues[p]) + g = Abs(eigenvalues[p])) and
           (Abs(eigenvalues[q]) + g = Abs(eigenvalues[q])) then
          A[p][q] := 0.0
        else if Abs(A[p][q]) > TOL then
        begin
          done := false;
          h := eigenvalues[q] - eigenvalues[p];
          if Abs(h) + g = Abs(h) then
            t := A[p][q] / h
          else
          begin
            theta := 0.5 * h / A[p][q];
            t := 1.0 / (Abs(theta) + Sqrt(1.0 + theta * theta));
            if theta < 0.0 then
              t := -t;
          end;
          c   := 1.0 / Sqrt(1.0 + t * t);
          s   := t * c;
          tau := s / (1.0 + c);
          h   := t * A[p][q];
          z[p] := z[p] - h;
          z[q] := z[q] + h;
          eigenvalues[p] := eigenvalues[p] - h;
          eigenvalues[q] := eigenvalues[q] + h;
          A[p][q] := 0.0;
          { rotate rows 0..p-1 }
          for i := 0 to p - 1 do
          begin
            g := A[i][p];
            h := A[i][q];
            A[i][p] := g - s * (h + g * tau);
            A[i][q] := h + s * (g - h * tau);
          end;
          { rotate rows p+1..q-1 }
          for i := p + 1 to q - 1 do
          begin
            g := A[p][i];
            h := A[i][q];
            A[p][i] := g - s * (h + g * tau);
            A[i][q] := h + s * (g - h * tau);
          end;
          { rotate rows q+1..n-1 }
          for i := q + 1 to n - 1 do
          begin
            g := A[p][i];
            h := A[q][i];
            A[p][i] := g - s * (h + g * tau);
            A[q][i] := h + s * (g - h * tau);
          end;
          { accumulate eigenvectors }
          for i := 0 to n - 1 do
          begin
            g := eigenvectors[i][p];
            h := eigenvectors[i][q];
            eigenvectors[i][p] := g - s * (h + g * tau);
            eigenvectors[i][q] := h + s * (g - h * tau);
          end;
          Inc(nRot);
        end;
      end;
    end;
    for i := 0 to n - 1 do
    begin
      b[i] := b[i] + z[i];
      eigenvalues[i] := b[i];
      z[i] := 0.0;
    end;
    if done then
      exit;
  end;
end;

{ ======================================================================
  5.6  HELPER UTILITIES
  ====================================================================== }
procedure MatVecMul(var A: TMatrix; var x: TFloatArray; n: integer; var res: TFloatArray);
var
  i, j : integer;
  sum  : Float;
begin
  SetLength(res, n);
  for i := 0 to n - 1 do
  begin
    sum := 0.0;
    for j := 0 to n - 1 do
      sum := sum + A[i][j] * x[j];
    res[i] := sum;
  end;
end;

procedure MatMatMul(var A, B: TMatrix; n: integer; var C: TMatrix);
var
  i, j, k : integer;
  sum     : Float;
begin
  SetLength(C, n);
  for i := 0 to n - 1 do
  begin
    SetLength(C[i], n);
    for j := 0 to n - 1 do
    begin
      sum := 0.0;
      for k := 0 to n - 1 do
        sum := sum + A[i][k] * B[k][j];
      C[i][j] := sum;
    end;
  end;
end;

procedure PrintMatrix(var A: TMatrix; n: integer; name: string);
var
  i, j : integer;
begin
  WriteLn('Matrix ', name, ':');
  for i := 0 to n - 1 do
  begin
    for j := 0 to n - 1 do
      Write(A[i][j]:10:5, ' ');
    WriteLn;
  end;
end;

{ ======================================================================
  SELF TEST
  ====================================================================== }
procedure self_test;
var
  i, j, n   : integer;
  A, Ainv, L, evecs, C : TMatrix;
  b, x, u, evals  : TFloatArray;
  a_tri, b_tri, c_tri, r_tri : TFloatArray;
  indx  : TIntArray;
  d     : Float;
  ok    : boolean;
  iter  : integer;
  nRot  : integer;
  det   : Float;
  err   : Float;
  evMin, evMid, evMax : Float;

  { allocate n x n TMatrix }
  procedure AllocMat(var M: TMatrix; sz: integer);
  var
    ii, jj : integer;
  begin
    SetLength(M, sz);
    for ii := 0 to sz - 1 do
    begin
      SetLength(M[ii], sz);
      for jj := 0 to sz - 1 do
        M[ii][jj] := 0.0;
    end;
  end;

begin
  WriteLn('=== jpmMatrices Self Test ===');
  WriteLn;

  { ------------------------------------------------------------------ }
  { TEST 1 : LU solve 3x3 system                                       }
  { A = [[2,1,-1],[−3,−1,2],[−2,1,2]]  b=[8,−11,−3]  x=[2,3,-1]       }
  { ------------------------------------------------------------------ }
  WriteLn('--- Test 1: LU solve 3x3 ---');
  n := 3;
  AllocMat(A, n);
  A[0][0] :=  2; A[0][1] :=  1; A[0][2] := -1;
  A[1][0] := -3; A[1][1] := -1; A[1][2] :=  2;
  A[2][0] := -2; A[2][1] :=  1; A[2][2] :=  2;
  SetLength(b, n);
  b[0] :=  8; b[1] := -11; b[2] := -3;
  SetLength(indx, n);
  LUDecompose(A, n, indx, d);
  LUSolve(A, n, indx, b);
  WriteLn('  Expected : x = [2, 3, -1]');
  WriteLn('  Computed : x = [', b[0]:6:3, ', ', b[1]:6:3, ', ', b[2]:6:3, ']');
  SelfTestCheck(Abs(b[0] - 2.0) < 1e-8, 'LUSolve x[0]=2');
  SelfTestCheck(Abs(b[1] - 3.0) < 1e-8, 'LUSolve x[1]=3');
  SelfTestCheck(Abs(b[2] - (-1.0)) < 1e-8, 'LUSolve x[2]=-1');
  WriteLn;

  { ------------------------------------------------------------------ }
  { TEST 2 : LU inverse — verify A * A^-1 ≈ I                         }
  { ------------------------------------------------------------------ }
  WriteLn('--- Test 2: LU inverse 3x3 ---');
  AllocMat(A, n);
  A[0][0] :=  2; A[0][1] :=  1; A[0][2] := -1;
  A[1][0] := -3; A[1][1] := -1; A[1][2] :=  2;
  A[2][0] := -2; A[2][1] :=  1; A[2][2] :=  2;
  LUInverse(A, n, Ainv);
  MatMatMul(A, Ainv, n, C);
  err := 0.0;
  for i := 0 to n - 1 do
    for j := 0 to n - 1 do
    begin
      d := C[i][j];
      if i = j then d := d - 1.0;
      if Abs(d) > err then err := Abs(d);
    end;
  WriteLn('  Max |A*A^-1 - I| = ', err:20:15, '  (expect < 1e-13)');
  SelfTestCheck(err < 1e-10, 'LUInverse A*A^-1≈I');
  WriteLn;

  { ------------------------------------------------------------------ }
  { TEST 3 : LU determinant                                            }
  { det([[2,1,-1],[-3,-1,2],[-2,1,2]]) = -2*2 + ... = ?              }
  { computed manually = -4 (use cofactor expansion)                   }
  { Let us just print it — sign matters, abs value = 4                 }
  { ------------------------------------------------------------------ }
  WriteLn('--- Test 3: LU determinant ---');
  AllocMat(A, n);
  A[0][0] :=  2; A[0][1] :=  1; A[0][2] := -1;
  A[1][0] := -3; A[1][1] := -1; A[1][2] :=  2;
  A[2][0] := -2; A[2][1] :=  1; A[2][2] :=  2;
  det := LUDeterminant(A, n);
  WriteLn('  det(A) = ', det:10:4, '  (expected -1.0000)');
  SelfTestCheck(Abs(det - (-1.0)) < 1e-6, 'LUDeterminant=-1');
  WriteLn;

  { ------------------------------------------------------------------ }
  { TEST 4 : Gauss-Seidel  4x4 diagonally dominant                    }
  { A*x = b  exact solution x=[1,2,3,4]                               }
  { ------------------------------------------------------------------ }
  WriteLn('--- Test 4: Gauss-Seidel 4x4 ---');
  n := 4;
  AllocMat(A, n);
  A[0][0] := 10; A[0][1] :=  1; A[0][2] :=  1; A[0][3] :=  1;
  A[1][0] :=  1; A[1][1] := 10; A[1][2] :=  1; A[1][3] :=  1;
  A[2][0] :=  1; A[2][1] :=  1; A[2][2] := 10; A[2][3] :=  1;
  A[3][0] :=  1; A[3][1] :=  1; A[3][2] :=  1; A[3][3] := 10;
  { b = A * [1,2,3,4]^T }
  SetLength(b, n);
  b[0] := 10*1 +  1*2 +  1*3 +  1*4;
  b[1] :=  1*1 + 10*2 +  1*3 +  1*4;
  b[2] :=  1*1 +  1*2 + 10*3 +  1*4;
  b[3] :=  1*1 +  1*2 +  1*3 + 10*4;
  SetLength(x, n);
  for i := 0 to n - 1 do x[i] := 0.0;
  ok := GaussSeidel(A, b, x, n, 200, 1e-10, iter);
  WriteLn('  Expected : x = [1, 2, 3, 4]');
  WriteLn('  Computed : x = [', x[0]:6:3, ', ', x[1]:6:3, ', ',
    x[2]:6:3, ', ', x[3]:6:3, ']');
  WriteLn('  Converged: ', ok, '  iterations: ', iter);
  SelfTestCheck(ok, 'GaussSeidel converged');
  SelfTestCheck(Abs(x[0] - 1.0) < 1e-6, 'GaussSeidel x[0]=1');
  SelfTestCheck(Abs(x[1] - 2.0) < 1e-6, 'GaussSeidel x[1]=2');
  SelfTestCheck(Abs(x[2] - 3.0) < 1e-6, 'GaussSeidel x[2]=3');
  SelfTestCheck(Abs(x[3] - 4.0) < 1e-6, 'GaussSeidel x[3]=4');
  WriteLn;

  { ------------------------------------------------------------------ }
  { TEST 5 : Cholesky  3x3 SPD                                        }
  { A = [[4,2,2],[2,5,3],[2,3,6]]  b=[8,13,13]  x=[1,1,1]            }
  { ------------------------------------------------------------------ }
  WriteLn('--- Test 5: Cholesky 3x3 ---');
  n := 3;
  AllocMat(A, n);
  A[0][0] := 4; A[0][1] := 2; A[0][2] := 2;
  A[1][0] := 2; A[1][1] := 5; A[1][2] := 3;
  A[2][0] := 2; A[2][1] := 3; A[2][2] := 6;
  SetLength(b, n);
  b[0] :=  8; b[1] := 10; b[2] := 11;  { = A*[1,1,1] }
  ok := Cholesky(A, n, L);
  WriteLn('  Cholesky factored: ', ok);
  SetLength(x, n);
  CholeskySolve(L, n, b, x);
  WriteLn('  Expected : x = [1, 1, 1]');
  WriteLn('  Computed : x = [', x[0]:6:3, ', ', x[1]:6:3, ', ', x[2]:6:3, ']');
  SelfTestCheck(ok, 'Cholesky factored');
  SelfTestCheck(Abs(x[0] - 1.0) < 1e-8, 'CholeskySolve x[0]=1');
  SelfTestCheck(Abs(x[1] - 1.0) < 1e-8, 'CholeskySolve x[1]=1');
  SelfTestCheck(Abs(x[2] - 1.0) < 1e-8, 'CholeskySolve x[2]=1');
  WriteLn;

  { ------------------------------------------------------------------ }
  { TEST 6 : Tridiagonal  4x4                                         }
  {  [2 -1  0  0][x0]   [1]                                           }
  {  [-1 2 -1  0][x1] = [0]  =>  x=[3/5, 4/5, 4/5, 3/5]? use simple  }
  {  [0 -1  2 -1][x2]   [0]  b=[1,0,0,1]  => x=[1,1,1,1]? (below)    }
  {  [0  0 -1  2][x3]   [1]                                           }
  { A*[1,1,1,1] = [2-1, -1+2-1, -1+2-1, -1+2] = [1,0,0,1] check!    }
  { ------------------------------------------------------------------ }
  WriteLn('--- Test 6: Tridiagonal 4x4 ---');
  n := 4;
  SetLength(a_tri, n);
  SetLength(b_tri, n);
  SetLength(c_tri, n);
  SetLength(r_tri, n);
  { main diagonal b = 2 }
  b_tri[0] := 2; b_tri[1] := 2; b_tri[2] := 2; b_tri[3] := 2;
  { sub-diagonal a (a[0] unused) }
  a_tri[0] := 0; a_tri[1] := -1; a_tri[2] := -1; a_tri[3] := -1;
  { super-diagonal c (c[n-1] unused) }
  c_tri[0] := -1; c_tri[1] := -1; c_tri[2] := -1; c_tri[3] := 0;
  { rhs = [1,0,0,1] }
  r_tri[0] := 1; r_tri[1] := 0; r_tri[2] := 0; r_tri[3] := 1;
  SetLength(u, n);
  TriDiagSolve(a_tri, b_tri, c_tri, r_tri, n, u);
  WriteLn('  Expected : u = [1, 1, 1, 1]');
  WriteLn('  Computed : u = [', u[0]:6:3, ', ', u[1]:6:3, ', ',
    u[2]:6:3, ', ', u[3]:6:3, ']');
  SelfTestCheck(Abs(u[0] - 1.0) < 1e-8, 'TriDiagSolve u[0]=1');
  SelfTestCheck(Abs(u[1] - 1.0) < 1e-8, 'TriDiagSolve u[1]=1');
  SelfTestCheck(Abs(u[2] - 1.0) < 1e-8, 'TriDiagSolve u[2]=1');
  SelfTestCheck(Abs(u[3] - 1.0) < 1e-8, 'TriDiagSolve u[3]=1');
  WriteLn;

  { ------------------------------------------------------------------ }
  { TEST 7 : Jacobi eigenvalues  3x3 symmetric                        }
  { A = [[4,1,0],[1,3,1],[0,1,2]]                                      }
  { eigenvalues approx 1.268, 3.000, 4.732  (sorted ascending)        }
  { ------------------------------------------------------------------ }
  WriteLn('--- Test 7: Jacobi eigenvalues 3x3 ---');
  n := 3;
  AllocMat(A, n);
  A[0][0] := 4; A[0][1] := 1; A[0][2] := 0;
  A[1][0] := 1; A[1][1] := 3; A[1][2] := 1;
  A[2][0] := 0; A[2][1] := 1; A[2][2] := 2;
  JacobiEigen(A, n, evals, evecs, nRot);
  WriteLn('  nRot = ', nRot);
  WriteLn('  Eigenvalues (unsorted):');
  for i := 0 to n - 1 do
    WriteLn('    lambda[', i, '] = ', evals[i]:10:6);
  WriteLn('  Expected ~1.267949, 3.000000, 4.732051');
  { verify eigenvalues are present (order may vary) }
  evMin := evals[0]; evMax := evals[0];
  for j := 1 to n-1 do
  begin
    if evals[j] < evMin then evMin := evals[j];
    if evals[j] > evMax then evMax := evals[j];
  end;
  evMid := (evals[0] + evals[1] + evals[2]) - evMin - evMax;
  SelfTestCheck(Abs(evMin - 1.267949) < 1e-4, 'JacobiEigen lambda_min');
  SelfTestCheck(Abs(evMid - 3.0) < 1e-4, 'JacobiEigen lambda_mid');
  SelfTestCheck(Abs(evMax - 4.732051) < 1e-4, 'JacobiEigen lambda_max');
  WriteLn;

  WriteLn('=== End jpmMatrices Self Test ===');
end;

end.
