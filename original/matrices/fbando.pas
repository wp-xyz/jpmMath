{ -------------------------- UNIT fbando.pas ------------------------ }
{       Translated from C language By J.-P. Moreau, March, 1999.      }
{   Reference: Numerical Algotithms with C, Springer, 1996 [BIBLI 11] }
{                                                                     }
{                           (www.jpmoreau.fr)                         }
{ ------------------------------------------------------------------- }
UNIT fbando;
 
INTERFACE
Uses WinCrt, Basis_r;

PROCEDURE bando                { Linear banded system without using pivots. }
          (
           mode : INTEGER;         { Modus: 0, 1, 2 ..................}
           n    : INTEGER;         { size of system ..................}
           ld   : INTEGER;         { # of lower co-diagonals .........}
           ud   : INTEGER;         { # of upper co-diagonals .........}
           VAR pmat : pVECT;       { condensed input matrix ..........}
           VAR b    : pVECT;       { right hand side .................}
           VAR code : INTEGER      { return code ( 0 = OK ) ..........}
          );
PROCEDURE banodec             { Decompose a banded matrix .................}
          (
           n  : INTEGER;           { size of system ..................}
           ld : INTEGER;           { # of lower co-diagonals .........}
           ud : INTEGER;           { # of upper co-diagonals .........}
           VAR pmat : pVECT;       { condensed input matrix ..........}
           VAR code : INTEGER      { return code ( 0 = OK ) ..........}
          );
PROCEDURE banosol             { Solve a banded system .....................}
          (
           n  : INTEGER;           { size of system ..................}
           ld : INTEGER;           { # of lower co-diagonals .........}
           ud : INTEGER;           { # of upper co-diagonals .........}
           VAR pmat : pVECT;       { condensed input matrix ..........}
           VAR b    : pVECT;       { right hand side .................}
           VAR code : INTEGER      { return code ( 0 = OK ) ..........}
          );

IMPLEMENTATION

PROCEDURE bando;
 {=====================================================================
 *                                                                    *
 *  The procedure bando solves a linear banded system:  pmat * x = b. *
 *  Here pmat is a nonsingular n x n matrix in condensed form, i.e.   *
 *  represented in an ld+ud+1 x n matrix for its ld lower and ud upper*
 *  co-diagonals. b denotes the right hand side of the system, and x  *
 *  is the solution.                                                  *
 *                                                                    *
 *  bando uses the Gauss algorithm without column pivot search.       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Applications:                                                    *
 *   =============                                                    *
 *      Solve linear systems with nonsingular banded system matrices. *
 *      Particularly useful for large sparse and banded and diagonally*
 *      dominant matrices.                                            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Control parameter:                                               *
 *   ==================                                               *
 *      mode     integer mode                                         *
 *               calling modus for band:                              *
 *       = 0     factor matrix and solve linear system                *
 *       = 1     factor only, store factors in pmat                   *
 *       = 2     solve linear system only; for this call, the factors *
 *               must already be available in pmat, such as when      *
 *               many systems are solved for differing right hand     *
 *               sides and the same system matrix.                    *
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer n;  ( n > 2 )                                *
 *               Dimension of pmat, size of b                         *
 *      ld       integer ld; ( ld >= 0 )                              *
 *               number of lower co-diagonals                         *
 *      ud       integer ud; ( ud >= 0 )                              *
 *               number of upper co-diagonals                         *
 *      pmat     REAL pmat(n,n) stored in a vector of type VECT.      *
 *               mode = 0, 1:                                         *
 *               Matrix of the system in comndensed form.             *
 *               Each row has length at least ld + 1 + ud + min(ld,ud)*
 *               where the columns 0, .., ld-1 denote the lower       *
 *               co-diagonals, column ld stores the diagonal and the  *
 *               columns ld+1, .., ld+ud contain the upper            *
 *               co-diagonals.                                        *
 *               If A is the original uncondensed band matrix, then : *
 *               A[i][k] = pmat[i][ld+k-i],                           *
 *                           for k,i inside the band                  *
 *               mode = 2:                                            *
 *               LU factors in condensed form.                        *
 *      b        REAL   b[n];          ( for mode = 0, 2 )            *
 *               Right hand side                                      *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      pmat     REAL pmat(n,n);       ( for mode = 0, 1 )            *
 *               LU factorization in condensed form                   *
 *      b        REAL   b[n];          ( for mode = 0, 2 )            *
 *               solution vector for the system                       *
 *                                                                    *
 *   Return value rc:                                                 *
 *   ================                                                 *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *      = 2      lack of space                                        *
 *      = 3      Matrix is numerically singular                       *
 *      = 4      wrong calling modus                                  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Procedures called :                                              *
 *   =================                                                *
 *      banodec : Factor matrix                                       *
 *      banosol : solve linear system                                 *
 *                                                                    *
 *====================================================================}
VAR rc : INTEGER;
Begin
  if((n < 3) or (ld < 0) or (ud < 0) or (n < ld + 1 + ud)) then
  begin
    rc := 1;
    exit
  end;

  CASE mode OF
    0: { factor and solve system ...............................}
       begin
            banodec (n, ld, ud, pmat, rc);
            if rc = 0 then
              banosol (n, ld, ud, pmat, b, code)
            else
              code:=rc
       end;
    1: { factor only ...........................................}
            banodec (n, ld, ud, pmat, code);

    2: { solve only ............................................}
            banosol (n, ld, ud, pmat, b, code)
    else
            code := 3
  END;
End;


PROCEDURE banodec;
{======================================================================
 *                                                                    *
 *  The procedure banodec factors a condensed banded matrix pmat.     *
 *  banddec uses the Gauss algorithm without column pivot search.     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer n;  ( n > 2 )                                *
 *               Dimension of pmat, size of b                         *
 *      ld       integer ld; ( ld >= 0 )                              *
 *               number of lower co-diagonals                         *
 *      ud       integer ud; ( ud >= 0 )                              *
 *               number of upper co-diagonals                         *
 *      pmat     REAL pmat(n,n) stored in a vector of type VECT.      *
 *               Matrix of the system in comndensed form.             *
 *               Each row has length at least ld + 1 + ud + min(ld,ud)*
 *               where the columns 0, .., ld-1 denote the lower       *
 *               co-diagonals, column ld stores the diagonal and the  *
 *               columns ld+1, .., ld+ud contain the upper            *
 *               co-diagonals.                                        *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      pmat     REAL pmat(n,n) stored in a vector of type VECT.      *
 *               LU factorization in condensed form                   *
 *                                                                    *
 *   Return value code:                                               *
 *   =================                                                *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *      = 2      Matrix is numerically singular                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :     MACH_EPS                                    *
 *   ================                                                 *
 *                                                                    *
 *   Functions:  min (user defined)                                   *
 *   =========                                                        *
 *                                                                    *
 *====================================================================}
VAR
  kend, kjend, jm, jk, k, j, i, m : INTEGER;
Begin
  m := ud + ld + 1;                  { dimension #2 of pmat          }

  if ld = 0 then                     { Matrix already in upper       }
  begin                              { triangular form               }
    code := 0;
    exit
  end;

  for i := 0 to n - 2 do                    { loop over all rws      }
  begin
    kend  := min (ld + 1, n - i);
    kjend := min (ud + 1, n - i);

    if ABS(pmat^[i*m+ld]) < MACH_EPS then    { LU decompsition does   }
    begin                                   { not exist              }
      code := 2;
      exit
    end;

    for k:=1 to kend-1 do                   { loop over all rows     }
    begin                                   { below row i            }
      pmat^[(k+i)*m+ld-k] := pmat^[(k+i)*m+ld-k] / pmat^[i*m+ld];

      for j := 1 to kjend-1 do
      begin
        jk := j + ld - k;
        jm := j + ld;
        pmat^[(k+i)*m+jk] := pmat^[(k+i)*m+jk] - pmat^[(k+i)*m+ld-k] * pmat^[i*m+jm]
      end  
    end   { for k }
  end;  { for i }

  code := 0
End;


PROCEDURE banosol;
{======================================================================
 *                                                                    *
 *  The procedure banosol solves a factored linear banded system in   *
 *  condensed form using banodec.                                     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer n;  ( n > 2 )                                *
 *               Dimension of pmat, size of b                         *
 *      ld       integer ld; ( ld >= 0 )                              *
 *               number of lower co-diagonals                         *
 *      ud       integer ud; ( ud >= 0 )                              *
 *               number of upper co-diagonals                         *
 *      pmat     REAL pmat(n,n) stored in a vector of type VECT.      *
 *               Matrices for the factored system in comndensed form. *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      b        REAL b[n];                                           *
 *               solution vector of linear system                     *
 *                                                                    *
 *   Return value code:                                               *
 *   =================                                                *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used :     None                                        *
 *   ================                                                 *
 *                                                                    *
 *   Functions:  min (user defined)                                   *
 *   =========                                                        *
 *                                                                    *
 *====================================================================}
VAR
  kend, i, k, m : INTEGER;
Begin
  m := ld + ud + 1;
                             {  Invalid input parameter   }
  if ((n < 3) or (ld < 0) or (ud < 0) or (m > n)) then
  begin
    code := 1;
    exit
  end;

  for i := 0 to n - 2 do
  begin
    kend := min (ld + 1, n - i);
    for k := 1 to kend-1 do
      b^[k+i] := b^[k+i] - pmat^[(k+i)*m+ld-k] * b^[i]
  end;

  for i := n - 1 downto 0 do     { back substitution        }
  begin
    kend := min (ud + 1, n - i);
    for k := 1 to kend-1 do
      b^[i] := b^[i] - pmat^[i*m+k+ld] * b^[i+k];
    b^[i] := b^[i] / pmat^[i*m+ld]
  end;

  code := 0
End;

END.
{ ------------------------ END fbando.pas ------------------------- }