{****************************************************************
*         Gauss algorithm for solving linear equations          *
* ------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].    *
*                                                               *
*                           TPW Release By J-P Moreau, Paris.   *
*                                  (www.jpmoreau.fr)            *
****************************************************************} 
UNIT FGAUSS;

INTERFACE

Const
      MACH_EPS = 1.2e-16;
      MAXROOT  = 1e16;
      NMAX = 10;
      ONE  = 1.0;
      ZERO = 0.0;

Type
      pMAT = ^MATRIX;
      MATRIX = Array[0..NMAX,0..NMAX] of Double;
      pVEC = ^VECTOR;
      VECTOR = Array[0..NMAX] of Double;
      pIVEC  = ^IVEC;
      IVEC   = Array[0..NMAX] of Integer;

      Procedure gauss(mode,n: integer; mat: pMAT; lumat: pMAT;
                      perm: pIVEC; b:pVEC; x: pVEC; var signd: integer;
                      var code: integer);

      Procedure mgauss            {Gauss for multiple right hand sides}
           (
            n,                      {Dimension of system .............}
            k: Integer;             {number of right hand sides ......}
            mat:  pMAT;             {original matrix .................}
            rmat: pMAT;             {Right hand sides/solutions ......}
        var code: Integer           {Error code ......................}
           );

IMPLEMENTATION

{headers of functions used by gauss}

Procedure  gaudec         {Gauss decomposition .......................}
           (
            n:     integer;        { size of matrix ..................}
            mat:   pMAT;           { Input matrix ....................}
            lumat: pMAT;           { matrix decomposition ............}
            perm:  pIVEC;          { row interchanges ................}
        var signd: integer;        { sign of perm ....................}
        var rc:    integer         { error code ......................}
           ); Forward;

Procedure  gausol        { Gauss solution ............................}
           (
            n:     integer;        { size of matrix ..................}
            lumat: pMAT;           { decomposed matrix (LU) ..........}
            perm:  pIVEC;          { row permutation vector ..........}
            b:     pVEC;           { Right hand side .................}
            x:     pVEC;           { solution ........................}
        var rc:    integer         { error code ......................}
           ); Forward;

        Procedure ISWAP(var a,b: integer);
        Var tmp:integer;
        Begin
          tmp:=b; b:=a; a:=tmp
        End;

        Procedure SWAP(var a,b: double);
        Var tmp:double;
        Begin
          tmp:=b; b:=a; a:=tmp
        End;

Procedure gauss          {Gauss algorithm for solving linear equations}
          (
           mode,                   { Modus: 0, 1, 2, 3 ...............}
           n:      integer;        { Dimension of matrix .............}
           mat:    pMAT;           { Input matrix ....................}
           lumat:  pMAT;           { LU decomposition ................}
           perm:   pIVEC;          { row remutation vector ...........}
           b:      pVEC;           { right hand side .................}
           x:      pVEC;           { solution of the system ..........}
       var signd:  integer;        { sign of the permutation .........}
       var code:   integer
          );
 {====================================================================*
 *                                                                    *
 *  The procedure gauss solves a linear system :  mat * x = b.        *
 *  Here mat is the nonsingular system matrix, b the right hand side  *
 *  of the system and x the solution vector.                          *
 *                                                                    *
 *  gauss uses the Gauss algorithm and computes a triangular factori- *
 *  zation of mat and scaled column pivot search.  (Crout method with *
 *  row swaps).                                                       *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Application:                                                     *
 *   ============                                                     *
 *      Solve general linear system with a nonsingular coefficient    *
 *      matrix.                                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Control parameter:                                               *
 *   ==================                                               *
 *      mode     integer;                                             *
 *               calling modus for gauss:                             *
 *       = 0     Find factorization and solve linear system           *
 *       = 1     Find factorization only.                             *
 *       = 2     Solve linear system only; the factorization is       *
 *               already available in lumat. This saves work when     *
 *               solving a linear system repeatedly for several right *
 *               hand sides and the same system matrix such as when   *
 *               inverting the matrix.                                *
 *       = 3     as under 2, additionally we improve the solution     *
 *               via iterative refinement (not available here).       *
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer;  (n > 0)                                    *
 *               Dimension of mat and lumat,                          *
 *               size of the vector b, the right hand side, the       *
 *               solution x and the permutation vector perm.          *
 *      mat      type pMAT; pointer to matrix of the linear system.   *
 *      lumat    type pMAT;        (for mode = 2, 3)                  *
 *               LU factors of mat                                    *
 *               lumat can be stored in the space of mat.             *
 *      perm     type pIVEC;       (for mode = 2, 3)                  *
 *               Permutation vector, of the row interchangfes in lumat*
 *      b        type pVEC;        (for mode = 0, 2, 3)               *
 *               Right hand side of the system.                       *
 *      signd    type Integer;     (for mode = 2, 3)                  *
 *               sign of the permutation in perm; the determinant of  *
 *               mat can be computed as the product of the diagonal   *
 *               entries of lumat times signd.                        *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      lumat    type pMAT;        (for mode = 0, 1)                  *
 *               LU factorization of mat.                             *
 *      perm     type pIVEC;       (for mode = 0, 1)                  *
 *               row ermutation vektor                                *
 *      x        type pVEC;        (for mode = 0, 2, 3)               *
 *               pointer to solution vector.                          *
 *      signd    type Integer;     (for mode = 0, 1)                  *
 *               sign of perm.                                        *
 *                                                                    *
 *   Return value (code):                                             *
 *   ===================                                              *
 *      =-1      Max. number (MAXITER) of iterative refinements       *
 *               reached (MAXITER) while mod = 3                      *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or other invalid input                         *
 *      = 2      lack of memory                                       *
 *      = 3      Matrix singular                                      *
 *      = 4      Matrix numerically singular                          *
 *      = 5      incorrect call                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Procedures used:                                                 *
 *   ================                                                 *
 *                                                                    *
 *      gaudec: determines  LU decomposition                          *
 *      gausol: solves the linear system                              *
 *                                                                    *
 *====================================================================}
Label fin;  {to exit procedure}
Var rc: integer;
Begin
  if n < 1 then
  begin
    code := 1;
    goto fin
  end;

  Case mode of

    0: begin  {Find factorization and solve system ...................}
         gaudec (n, mat, lumat, perm, signd,rc);
         if rc = 0 then
         begin
           gausol(n, lumat, perm, b, x,rc);
           code:=rc
         end
         else
           code:=rc;
         goto fin
       end;

    1: begin  {Find factorization only ...............................}
         gaudec(n, mat, lumat, perm, signd,rc);
         code:=rc;
         goto fin
       end;

    2: begin  {Solve only ............................................}
         gausol(n, lumat, perm, b, x,rc);
         code:=rc;
         goto fin
       end;

    3: begin  {solve and then use iterative refinement ...............}
         writeln(' fgauss: gausoli not implemented.');
	 code:=5;
         goto fin
       end
  End;

  code := 5;                                               {Wrong call}
fin:End;


Procedure  gaudec         {Gauss decomposition .......................}
           (
            n:     integer;        { size of matrix ..................}
            mat:   pMAT;           { Input matrix ....................}
            lumat: pMAT;           { matrix decomposition ............}
            perm:  pIVEC;          { row interchanges ................}
        var signd: integer;        { sign of perm ....................}
        var rc:    integer         { error code ......................}
           );
 {====================================================================*
 *                                                                    *
 *  gaudec decomposes a nonsingular n x n matrix into a product of a  *
 *  unit lower and an upper triangular matrix. Both triangular factors*
 *  are stored in lumat (minus the unit diagonal, of course).         *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Input parameter:                                                 *
 *   ================                                                 *
 *      n        integer;  (n > 0)                                    *
 *               Dimension of  mat and lumat,                         *
 *               size of  b , x and perm.                             *
 *      mat      type pMAT;                                           *
 *               pointer to original system matrix in vector form.    *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      lumat    type pMAT;                                           *
 *               pointer to LU factorization                          *
 *      perm     type pIVEC;                                          *
 *               pointer to row permutation vector for lumat          *
 *      signd    TYPE integer;                                        *
 *               sign of perm. The determinant of mat can be computed *
 *               as the product of the diagonal entreis of lumat times*
 *               signd.                                               *
 *                                                                    *
 *   Return value (rc):                                               *
 *   =================                                                *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or invalid input                               *
 *      = 2      lack of memory                                       *
 *      = 3      Matrix is singular                                   *
 *      = 4      Matrix numerically singular                          *
 *                                                                    *
 *====================================================================}
Label fin; {to exit procedure}
Var
  m, j, i, j0: Integer;
  piv, tmp, zmax: Double;
  d: VECTOR;  {scaling vector for pivoting}

Begin

  if n < 1 then
  begin
    rc:=1;                              {Invalid parameters}
    goto fin
  end;

  {copy mat to lumat}
  for i:=0 to n-1 do
    for j:=0 to n-1 do
      lumat^[i,j] := mat^[i,j];

  for i := 0 to n-1 do
  begin
    perm^[i] := i;                       { Initialize perm}
    zmax := ZERO;
    for j := 0 to n-1 do                 {find row maxima}
    begin
      tmp := ABS(lumat^[i,j]);
      if tmp > zmax then zmax := tmp
    end;

    if zmax = ZERO then                  {mat is singular}
    begin
      rc := 3;
      goto fin
    end;
    d[i] := ONE / zmax
  end;

  signd := 1;                            {initialize sign of perm}

  for i := 0 to n-1 do
  begin
    piv := ABS(lumat^[i,i]) * d[i];
    j0 := i;                             {Search for pivot element}
    for j := i + 1 to n-1 do
    begin
      tmp := ABS(lumat^[j,i]) * d[j];
      if piv < tmp then
      begin
        piv := tmp;                      {Mark pivot element and}
        j0 := j                          {its location}
      end
    end;

    if piv < MACH_EPS then               {If piv is small, mat is}
    begin                                {nearly singular}
      signd := 0;
      rc :=4;
      goto fin
    end;

    if j0 <> i then
    begin
      signd := - signd;                  {update signd}
      ISWAP(perm^[j0], perm^[i]);        {swap pivot entries}
      SWAP(d[j0], d[i]);                 {swap scaling vector}
      for j:=0 to n-1 do                 
        SWAP(lumat^[j0,j], lumat^[i,j])  {swap j0-th and i-th rows of lumat}
    end;

    for j := i + 1 to n-1 do             {Gauss elimination}
    begin
      if lumat^[j,i] <> ZERO then
      begin
        lumat^[j,i] := lumat^[j,i] / lumat^[i,i];
        tmp := lumat^[j,i];
        for m := i + 1 to  n-1 do
          lumat^[j,m] := lumat^[j,m] - tmp * lumat^[i,m]
      end
    end
  end; {i loop}

  rc:=0;  {all ok}
fin:End;

Procedure  gausol        { Gauss solution ............................}
           (
            n:     integer;        { size of matrix ..................}
            lumat: pMAT;           { decomposed matrix (LU) ..........}
            perm:  pIVEC;          { row permutation vector ..........}
            b:     pVEC;           { Right hand side .................}
            x:     pVEC;           { solution ........................}
        var rc:    integer         { error code ......................}
           );
 {====================================================================*
 *                                                                    *
 *  gausol  finds the solution x of the linear system  lumat * x = b  *
 *  for the product matrix lumat, that describes an LU decomposition, *
 *  as produced by gaudec.                                            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer;  (n > 0)                                    *
 *               Dimension of lumat,                                  *
 *      lumat    type pMAT;                                           *
 *               pointer to LU factorization, as produced from gaudec *
 *      perm     type pIVEC;                                          *
 *               pointer to row permutation vector for lumat          *
 *      b        type pVEC;                                           *
 *               pointer to Right hand side of the system.            *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      x        type pVEC;                                           *
 *               pointer to solution vector                           *
 *                                                                    *
 *   Return value (rc):                                               *
 *   =================                                                *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or other invalid input parameter               *
 *      = 3      improper LU decomposition ( zero diagonal entry)     *
 *                                                                    *
 *====================================================================}
Label fin; {to exit procedure}
Var
  j,k: integer;
  sum: double;
Begin
  if n < 1 then
  begin
     rc:=1;                                   {Invalid input parameter}
     goto fin
  end;

  for k := 0 to n-1 do                                       {update b}
  begin
    x^[k] := b^[perm^[k]];
    for j := 0 to k-1 do
      x^[k] := x^[k] - lumat^[k,j] * x^[j]
  end;

  for k := n - 1 Downto 0 do                          {back substitute}
  begin
    sum := ZERO;
    for j := k + 1 to n-1 do
      sum := sum + lumat^[k,j] * x^[j];

    if ABS(lumat^[k,k]) < MACH_EPS then
    begin
      rc:=3;
      goto fin
    end;
    x^[k] := (x^[k] - sum) / lumat^[k,k]
  end;

  rc := 0; {all ok}
fin:End;

Function det              {Determinant  ..............................}
           (
            n: integer;             {Dimension of the matrix .........}
            mat:pMAT                {matrix ..........................}
           ): Double;
 {====================================================================*
 *                                                                    *
 *  det computes the determinant of an n x n real matrix mat          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameter:                                                 *
 *   ================                                                 *
 *      n        integer;  (n > 0)                                    *
 *               Dimension of mat                                     *
 *      mat      pointer to (type pMAT)                               *
 *               n x n matrix                                         *
 *                                                                    *
 *   Return value:                                                    *
 *   =============                                                    *
 *      REAL     Determinant of mat.                                  *
 *               If the return value = 0, then the matrix is singular *
 *               or the storage is insufficient                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use  :                                              *
 *   ===================                                              *
 *                                                                    *
 *      int gaudec ():    LU decomposition of mat                     *
 *                                                                    *
 *====================================================================}
Label fin;  {to exit function}
Var
  i, rc, signd: Integer;
  perm: pIVEC;
  lu:   pMAT;
  tmpdet: Double;
Begin

  New(perm); New(lu);

  if n < 1 then
  begin
    det:=ZERO;
    goto fin
  end;
 
  gaudec(n, mat, lu, perm, signd,rc);                       {decompose}

  if (rc <> 0) or (signd = 0) then
  begin
    det:=ZERO;
    goto fin
  end;

  tmpdet := 1.0*signd;

  for i := 0 to n-1 do
    if ABS(tmpdet) < MACH_EPS then
    begin
      det:=ZERO;
      goto fin
    end
    else if (ABS(tmpdet) > MAXROOT) or (ABS(lu^[i,i]) > MAXROOT) then
    begin
      det:=MAXROOT;
      goto fin
    end
    else
      tmpdet := tmpdet * lu^[i,i];                        {compute det}

  det:=tmpdet;

fin: Dispose(perm); Dispose(lu)

End;

Procedure mgauss                  {Gauss for multiple right hand sides}
           (
            n,                      {Dimension of system .............}
            k: Integer;             {number of right hand sides ......}
            mat:  pMAT;             {original matrix .................}
            rmat: pMAT;             {Right hand sides/solutions ......}
        var code: Integer           {Error code ......................}
           );
 {====================================================================*
 *                                                                    *
 *  mgauss  finds the solution matrix x for the linear system         *
 *  mat * x = rmat with an  n x n coefficient matrix mat and a        *
 *  n x k matrix rmat of right hand sides. Here mat must be           *
 *  nonsingular.                                                      *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer;  (n > 0)                                    *
 *               Dimension of mat.                                    *
 *      k        integer k;  (k > 0)                                  *
 *               number of right hand sides                           *
 *      mat      type pMAT; pointer to a                              *
 *               n x n original system matrix                         *
 *      rmat     type pMAT;                                           *
 *               pointer to matrix of right hand sides                *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      rmat     type pMAT;                                           *
 *               solution matrix for the system.                      *
 *               The input right hand sides are lost.                 *
 *                                                                    *
 *   Return value (code):                                             *
 *   ===================                                              *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or k < 1 or invalid input parameter            *
 *      = 2      lack of memory (not used here)                       *
 *      = 3      mat is numerically singular.                         *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   procedures used:                                                 *
 *   ================                                                 *
 *                                                                    *
 *      gaudec:  LU decomposition of mat.                             *
 *                                                                    *
 *====================================================================}
Label fin; {to exit procedure}
Var
  i, j: integer;
  m, signd, rc: integer;
  perm: pIVEC;
  lu: pMAT;
  x:  pVEC;
  sum: Double;
Begin

  New(perm); New(lu); New(x);

  if (n < 1) or (k < 1) then                        {Invalid parameter}
  begin
    code:=1;                   
    goto fin
  end;

  gaudec (n, mat, lu, perm, signd, rc);         {compute factorization}
                                                {in matrix lu         }
  if (rc <> 0) or (signd = 0) then              {if not possible      }
  begin
    code:=3;                                    {exit with code=3     }
    goto fin
  end;

  for m := 0 to k-1 do                 {Loop over the right hand sides}
  begin
    for i := 0 to n-1 do                             {Updating the b's}
    begin
      x^[i] := rmat^[perm^[i],m];
      for j := 0 to i-1 do
        x^[i] := x^[i] - lu^[i,j] * x^[j]
    end;

    for i := n - 1 Downto 0 do                      {back substitution}
    begin
      sum := ZERO;
      for j := i + 1 to n-1 do
        sum := sum + lu^[i,j] * x^[j];

      if ABS(lu^[i,i]) < MACH_EPS then       {invalid LU decomposition}
      begin
        code := 2;
        goto fin
      end;
      x^[i] := (x^[i] - sum) / lu^[i,i];
    end;

    for j := 0 to n-1 do                                 {Save result}
      rmat^[j,m] := x^[j]
  end;

  code := 0;

fin: Dispose(perm); Dispose(lu); Dispose(x)

End;


END. {of unit}

{ ------------------------ END fgauss.pas ------------