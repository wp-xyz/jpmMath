{ ------------------------- fband.pas UNIT ------------------------- }
{       Translated from C language By J.-P. Moreau, March, 1999.     }
{  Reference: Numerical Algotithms with C, Springer, 1996 [BIBLI 11] }
{ ------------------------------------------------------------------ }
UNIT fband;           {Pascal for Windows version in simple precision}
                                   {www.jpmoreau.fr}

INTERFACE
Uses WinCrt, Basis_r;

PROCEDURE band           { Linear systems with banded matrices .......}
         (
          mode : INTEGER;         { Modus: 0, 1, 2 .................. }
          n    : INTEGER;         { size of system .................. }
          ld   : INTEGER;         { # of lower co-diagonals ......... }
          ud   : INTEGER;         { # of upper co-diagonals ......... }
          VAR pmat : pVECT;       { condensed input matrix .......... }
          VAR b    : pVECT;       { right hand side ................. }
          VAR perm : IVECT;       { row permutation vector .......... }
          VAR signd: INTEGER;     { sign of perm .................... }
          VAR code : INTEGER      { return code ..................... }
         );

PROCEDURE banddec       { Factor a banded matrix .................... }
            (
             n  : INTEGER;        { size of system .................. }
             ld : INTEGER;        { # of lower co-diagonals ......... }
             ud : INTEGER;        { # of upper co-diagonals ......... }
             VAR pmat : pVECT;    { condensed input matrix .......... }
             VAR perm : IVECT;    { row permutation vector .......... }
             VAR signd: INTEGER;  { sign of perm .................... }
             VAR code : INTEGER   { return code ..................... }
            );

PROCEDURE bandsol       { Solve a banded system ..................... }
            (
             n  : INTEGER;        { size of system .................. }
             ld : INTEGER;        { # of lower co-diagonals ......... }
             ud : INTEGER;        { # of upper co-diagonals ......... }
             VAR pmat : pVECT;    { condensed input matrix .......... }
             VAR b    : pVECT;    { right hand side ................. }
             VAR perm : IVECT;    { row permutation vector .......... }
             VAR code : INTEGER   { return code ..................... }
            );


IMPLEMENTATION
 
 PROCEDURE band;
 {=====================================================================
 *                                                                    *
 *  The procedure band solves a linear banded system:  pmat * x = b.  *
 *  Here pmat is a nonsingular n x n matrix in condensed form, i.e.   *
 *  represented in an ld+ud+1 x n matrix for its ld lower and ud upper*
 *  co-diagonals. b denotes the right hand side of the system, and x  *
 *  is the solution.                                                  *
 *                                                                    *
 *  band uses the Gauss algorithm with column pivot search.           *
 *  The result of pivoting are min( ud, ld) additional columns,       *
 *  so that pmat needs all in all a n x (ld+1+ud+min(ld,ud)) matrix.  *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Applications:                                                    *
 *   -------------                                                    *
 *      Solve linear systems with nonsingular banded system matrices. *
 *      Particularly useful for large sparse and banded matrices with *
 *      n >> ld+1+ud.                                                 *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Control parameter:                                               *
 *   ------------------                                               *
 *      mode     mode:INTEGER;                                        *
 *               calling modus for band:                              *
 *       = 0     factor matrix and solve linear system                *
 *       = 1     factor only, store factors in pmat                   *
 *       = 2     solve linear system only; for this call, the factors *
 *               must already be available in pmat, such as when      *
 *               many systems are solved for differing right hand     *
 *               sides and the same system matrix.                    *
 *                                                                    *
 *   Input parameters:                                                *
 *   -----------------                                                *
 *      n        integer n;  ( n > 2 )                                *
 *               Dimension of pmat, size of b                         *
 *      ld       integer ld; ( ld >= 0 )                              *
 *               number of lower co-diagonals                         *
 *      ud       integer ud; ( ud >= 0 )                              *
 *               number of upper co-diagonals                         *
 *      pmat     REAL  pmat(n,n) stored in a vector of VECT type.     *
 *               mode = 0, 1:                                         *
 *               Matrix of the system in condensed form.              *
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
 *      perm     integer perm[n];      ( for mode = 2 )               *
 *               row permutation vector                               *
 *      signd    integer signd;        ( for mode = 2 )               *
 *               sign of perm. The determinant of A can be computed   *
 *               as the product of the diagonal entries times signd.  *
 *                                                                    *
 *   Output parameters:                                               *
 *   ------------------                                               *
 *      pmat     REAL pmat[n];         ( for mode = 0, 1 )            *
 *               LU factorization in condensed form                   *
 *      perm     integer perm[n];      ( for mode = 0, 1 )            *
 *               row permutation vector                               *
 *      b        REAL   b[n];          ( for mode = 0, 2 )            *
 *               solution vector for the system                       *
 *      signd    integer signd;        ( for mode = 0, 1 )            *
 *               sign of perm                                         *
 *                                                                    *
 *   Return value code:                                               *
 *   -----------------                                                *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *      = 2      lack of space ( not used here )                      *
 *      = 3      Matrix is numerically singular                       *
 *      = 4      wrong calling modus                                  *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Procedures called :                                              *
 *   -----------------                                                *
 *      banddec : Factor matrix                                       *
 *      bandsol : solve linear system                                 *
 *                                                                    *
 =====================================================================}
VAR
  k,rc : INTEGER;
BEGIN
  if ((n < 1) or (ld < 0) or (ud < 0)) then {wrong parameters}
  begin
    code:=1;
    exit
  end;

  CASE mode OF
    0: { factor and solve system ...............................}
       begin
            banddec (n, ld, ud, pmat, perm, signd,rc);
            for k:=0 to n-1 do write(b^[k]:9:6); writeln;
            for k:=0 to n-1 do write(perm[k]:2); writeln;
            writeln(' signd=',signd:2);
            if rc = 0 then
            begin
              bandsol (n, ld, ud, pmat, b, perm, code);
              for k:=0 to n-1 do write(b^[k]:9:6); writeln;
              for k:=0 to n-1 do write(perm[k]:2); writeln; readkey
            end
            else
              code:=rc
       end;
    1: { factor only ...........................................}
            banddec (n, ld, ud, pmat, perm, signd, code);

    2: { solve only ............................................}
            bandsol (n, ld, ud, pmat, b, perm, code)
    else
            code := 4
  END;

END; {band}

PROCEDURE banddec;
 {=====================================================================
 *                                                                    *
 *  The procedure banddec factors a condensed banded matrix pmat.     *
 *  banddec uses the Gauss algorithm with column pivot search.        *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Input parameters:                                                *
 *   ----------------                                                 *
 *      n        integer n;  ( n > 2 )                                *
 *               Dimension of pmat, size of b                         *
 *      ld       integer ld; ( ld >= 0 )                              *
 *               number of lower co-diagonals                         *
 *      ud       integer ud; ( ud >= 0 )                              *
 *               number of upper co-diagonals                         *
 *      pmat     REAL pmat(n,n) stored in a vector of VECT type.      *
 *               Matrix of the system in comndensed form.             *
 *               Each row has length at least ld + 1 + ud + min(ld,ud)*
 *               where the columns 0, .., ld-1 denote the lower       *
 *               co-diagonals, column ld stores the diagonal and the  *
 *               columns ld+1, .., ld+ud contain the upper            *
 *               co-diagonals.                                        *
 *      perm     integer perm[n] stored in vector of IVECT type.      *
 *               row permutation vector                               *
 *      signd    integer signd;                                       *
 *               sign of perm. The determinant of A can be computed   *
 *               as the product of the diagonal entries times signd.  *
 *                                                                    *
 *   Output parameters:                                               *
 *   ------------------                                               *
 *      pmat     REAL pmat(n,n) stored in a vector of VECT type.      *
 *               LU factorization in condensed form                   *
 *      perm     int perm[n] stored in a vector of IVECT type.        *
 *               row permutation vector                               *
 *      signd    integer signd;                                       *
 *               sign of perm                                         *
 *                                                                    *
 *   Return value code:                                               *
 *   -----------------                                                *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *      = 2      not used                                             *
 *      = 3      Matrix is numerically singular                       *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Constants used :     MACH_EPS                                    *
 *   ----------------                                                 *
 *                                                                    *
 *   Functions used:   min, max, swap (user defined)                  *
 *   --------------                                                   *
 *                                                                    *
 =====================================================================}
VAR
  j0, mm, up, istart, iend, step, kstart,kend, kjend, km, jm, jk:INTEGER;
  k, j, i, ii, m : INTEGER;
  piv : REAL;
  tmp : pVECT;
  fp  : TEXT;
BEGIN
  New(tmp);

  m := ld + ud + 1;                       { second dimension of original pmat }

  if ((ld < 0) or (ud < 0) or (n < 1)) then
  begin
    code:=1;                             { Invalid input parameter    }
    exit
  end;

  mm := ld + 1 + ud + min (ld, ud);     { second dimension of extended pmat }

  if ld <= ud then up:=1                { up = 0 ==> transform into   }
  else up:=0;                           { lower triangular matrix     }

  for i := 0 to n-1 do
  { initialize needed extra columns using a temporary matrix}
    for k:=0 to mm-1 do
      if k<m then tmp^[i*mm+k]:=pmat^[i*m+k] else tmp^[i*mm+k]:=0.0;;

  Copy_vector(pmat,tmp,n*mm);

  signd := 1;                                { initialize signd       }
  if up=1 then 
  begin
    istart := 0; iend := n-1; step := 1;     { Find start, end and    }
    kstart := 1;                             { directions of loops    }
  end                                        { depending of up        }
  else
  begin
    istart := n-1; iend := 0; step := -1;
    kstart := -1;
  end;

  i:=istart;
  while i<>iend do
  begin
    {kend := (up ? min (ld+1, n-i) : max (-i-1, -ud-1)); }
    if up = 1 then kend := min (ld+1, n-i)
              else kend := max (-i-1, -ud-1);
    j0 := 0;
    piv := ABS (pmat^[i*mm+ld]);                { choose pivot          }

    k:=kstart;
    while k <> kend do
    begin
      if ABS(pmat^[(k+i)*mm+ld-k]) > piv  then
      begin
        piv := ABS(pmat^[(k+i)*mm+ld-k]);
        j0 := k
      end;
      k := k + step
    end;

    if piv < MACH_EPS then
    begin
      code := 3;                          { If piv = 0, matrix is  }
      writeln(fp,' BANDDEC - code=',code);
      close(fp);
      exit                                { singular               }
    end;

    perm[i] := j0;
    {kjend :=  (up ? min (j0+ud+1, n-i) : max ( -i-1, j0-ld-1)); }
    if up=1 then kjend := min (j0+ud+1, n-i)
            else kjend := max ( -i-1, j0-ld-1);

    if j0 <> 0 then
    begin
      signd := - signd;                   { swap rows               }
      k := 0;
      while k <> kjend do
      begin
        km := k + ld;
        if km < 0 then km := km + mm;

        SWAP(pmat^[i*mm+km], pmat^[(i+j0)*mm+k+ld-j0]);

        k := k + step
      end
    end;

    k:=kstart;
    while k<>kend do
    begin                                        { below row i         }
      pmat^[(k+i)*mm+ld-k] := pmat^[(k+i)*mm+ld-k]/ pmat^[i*mm+ld];
      j:=kstart;
      while j <> kjend do
      begin                                      { loop over all       }
        jk := j + ld - k;                        { columns right of i  }
        jm := j + ld;
                                       { additional columns from pivoting }
        if jk < 0 then jk := jk + mm;  { are stored to right of column    }
        if jm < 0 then jm := jm + mm;  {  ud+ld+1                         }

        pmat^[(k+i)*mm+jk] := pmat^[(k+i)*mm+jk] - pmat^[(k+i)*mm+ld-k] * pmat^[i*mm+jm];

        j := j + step
      end;    { while j }
      k := k + step
    end;    { while k }
    i := i + step
  end;    { while i }

  if up=1 then ii := n-1 else ii:=0;
  piv := ABS (pmat^[ii*mm+ld]);                { choose pivot            }
                                               { If piv = 0, matrix is   }
  if piv < MACH_EPS then                       { singular                }
  begin
    code:=3;
    writeln(fp,' BANDDEC - code=',code);
    close(fp);
    exit
  end;
  perm[iend] := 0;
  code := 0
END;  {banddec}


PROCEDURE bandsol;
 {=====================================================================
 *                                                                    *
 *  The procedure bandsol solves a factored linear banded system in   *
 *  condensed form using banddec.                                     *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Input parameters:                                                *
 *   -----------------                                                *
 *      n        integer n;  ( n > 2 )                                *
 *               Dimension of pmat, size of b                         *
 *      ld       integer ld; ( ld >= 0 )                              *
 *               number of lower co-diagonals                         *
 *      ud       integer ud; ( ud >= 0 )                              *
 *               number of upper co-diagonals                         *
 *      pmat     REAL pmat(n,n) stored in a vector of type VECT.      *
 *               Matrices for the factored system in comndensed form. *
 *      perm     integer perm[n] stored in a vector of type IVECT.    *
 *               row permutation vector                               *
 *                                                                    *
 *   Output parameters:                                               *
 *   ------------------                                               *
 *      b        REAL b[n];                                           *
 *               solution vector of linear system                     *
 *                                                                    *
 *   Return value code:                                               *
 *   -----------------                                                *
 *      = 0      all ok                                               *
 *      = 1      n < 3 or other incorrect input parameter             *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Functions used:  min, max, SWAP (user defined)                   *
 *   --------------                                                   *
 *                                                                    *
 =====================================================================}
Var
  i, k, s, mm, up, istart, iend, m, step, kstart, kend, km : INTEGER;
Begin

  m := ld + ud + 1;                        { second dimension of pmat }
                                       
  if ((ld < 0) or (ud < 0) or (n < 1)) then { invalid input           }
  begin
    code := 1;
    exit
  end;

  mm := ld + ud + 1 + min (ld, ud);     { mm := max. column number    }

  if ld <= ud then up:=1                { up = 0 ==> transform into   }
  else up:=0;                           { lower triangular matrix     }

  if up=1 then
  begin
    istart := 0; iend := n-1; step := 1;    { determine bounds and      }
    kstart := 1; s := -1;                   { direction of loop         }
  end                                       { depending on up           }
  else
  begin
    istart := n-1; iend := 0; step := -1;
    kstart := -1; s := 1;
  end;

  i:=istart;

  while i <> iend do
  begin
    if perm[i] <> 0 then
      SWAP( b^[i], b^[i+perm[i]]);

    {kend = (up ? min (ld+1, n-i) : max (-i-1, -ud-1)); }
    if up=1 then kend := min(ld+1, n-i)
            else kend := max(-i-1, -ud-1);

    k:=kstart;
    while k<>kend do
    begin
      b^[k+i] := b^[k+i] - pmat^[(k+i)*mm+ld-k] * b^[i];
      k:=k+step
    end;
    i:=i+step
  end;

  i:=iend;
  while i <> istart+s do
  begin  
    {kend :=  (up ? min (ld+ud+1, n-i) : max (-i-1, -ud-ld)); }
    if up=1 then kend := min(ld+ud+1, n-i)
            else kend := max(-i-1, -ud-ld);

    k:=kstart;
    while k <> kend do
    begin
      km := k + ld;                     { update and                  }
      if km < 0 then km := km + mm;     { back substitute             }
      b^[i] := b^[i] - pmat^[i*mm+km] * b^[i+k];
      k:=k+step
    end;
    if ABS(pmat^[i*mm+ld]) > MACH_EPS then
      b^[i] := b^[i] / pmat^[i*mm+ld]
    else
      b^[i]:=0.0;
    i:=i-step
  end;
  code := 0
End; {bandsol}

END.
{ --------------------------- END fband.pas ------------------------- }