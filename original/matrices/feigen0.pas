{ ------------------------ UNIT feigen.pas ------------------------- *
 * Reference: "Numerical Algorithms with C  By G. Engeln-Mueller and  *
 *             F. Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
 *                                                                    *
 *                              Pascal version By J-P Moreau, Paris.  *
 *                                       (www.jpmoreau.fr)            *
 * ------------------------------------------------------------------ }
UNIT FEIGEN0;          {version with index starting from zero}
                       {and static allocations}
INTERFACE

Uses WinCrt, Type_def; {for REAL_AR}

CONST
      Maxc=6; Nsol=6; Ndiag=6;
      ONE=1.0; TWO=2.0; ZERO=0.0;
      MACH_EPS = 2.22e-16; MAXIT = 50;

TYPE
  Square_Matrix   = ARRAY[0..Maxc, 0..Maxc]  OF REAL_AR;  {declared double in Type_def}
  Rect_Matrix     = ARRAY[0..Nsol, 0..Maxc]  OF REAL_AR;
  Matrix_n_mPlus1 = ARRAY[0..Maxc, 0..Ndiag] OF REAL_AR;
  Matrix_3x3      = ARRAY[0..2,0..2]         OF REAL_AR;
  Real_Vector     = ARRAY[0..Maxc]           OF REAL_AR;
  Integer_Vector  = ARRAY[0..Maxc]           OF INTEGER;

{visible from calling program}
PROCEDURE eigen (         {Compute all evalues/evectors of a matrix ..}
           vec:    Boolean;         {switch for computing evectors ...}
           ev_norm:Boolean;         {normalize Eigenvectors? .........}
           n:      integer;         {size of matrix ..................}
           mat:    Square_Matrix;   {input matrix ....................}
       VAR eivec:  Square_Matrix;   {Eigenvectors ....................}
       VAR valre:  Real_Vector;     {real parts of eigenvalues .......}
       VAR valim:  Real_Vector;     {imaginary parts of eigenvalues ..}
       VAR cnt:    Integer_Vector;  {Iteration counter ...............}
       VAR rc:     Integer          {return code .....................}
          );

IMPLEMENTATION

Procedure RSWAP(VAR a,b:REAL_AR);
Var t:REAL_AR;
Begin
  t:=a; a:=b; b:=t
End;

 {--------------------------------------------------------------------*
 * Auxiliary procedures for procedure eigen                           *
 *--------------------------------------------------------------------}

PROCEDURE balance          {balance a matrix .........................}
                 (n:   integer;         {size of matrix ..............}
              VAR mat:  Square_Matrix;  {input matrix ................}
              VAR scal: Real_Vector;    {Scaling data ................}
              VAR low:  Integer;        {first relevant row index ....}
              VAR high: Integer         {last relevant row index .....}
                 );
 {====================================================================*
 *                                                                    *
 *  balance balances the matrix mat so that the rows with zero entries*
 *  off the diagonal are isolated and the remaining columns and rows  *
 *  are resized to have one norm close to 1.                          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer;  ( n > 0 )                                  *
 *               Dimension of mat                                     *
 *      mat      Square_Matrix;                                       *
 *               n x n input matrix                                   *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      mat      Square_Matrix;                                       *
 *               scaled matrix                                        *
 *      low      integer;                                             *
 *      high     integer;                                             *
 *               the rows 0 to low-1 and those from high to n-1       *
 *               contain isolated eigenvalues (only nonzero entry on  *
 *               the diagonal)                                        *
 *      scal     Real_Vector;                                         *
 *               the vector scal contains the isolated eigenvalues in *
 *               the positions 0 to low-1 and high to n-1, its other  *
 *               components contain the scaling factors for           *
 *               transforming mat.                                    *
 *                                                                    *
 *====================================================================}
CONST
  basis = 2;     {Floating point basis of used CPU}
VAR
  i, j: Integer;
  iter, k, m: Integer;
  b2, r, c, f, g, s: REAL_AR;
Begin
  b2 := basis * basis;
  m := 0;
  k := n - 1;

  Repeat
    iter := 0;
    for j := k downto 0 do
    begin
      r := ZERO;
      for i := 0 to k do
        if i <> j then  r := r + ABS (mat[j,i]);

      if r = ZERO then
      begin
        scal[k] := j;
        if j <> k then
        begin
          for i := 0 to k do RSWAP(mat[i,j], mat[i,k]);
          for i := m to n-1 do RSWAP(mat[j,i], mat[k,i])
        end;
        Dec(k);
        iter := 1
      end
    end {of j loop}
  Until iter=0;

  Repeat
    iter := 0;
    for j := m to k do
    begin
      c := ZERO;
      for i := m to k do
        if i <> j then c := c + ABS(mat[i,j]);
      if c = ZERO then
      begin
        scal[m] := j;
        if j <> m then
        begin
          for i := 0 to k do RSWAP(mat[i,j], mat[i,m]);
          for i := m to n-1 do RSWAP(mat[j,i], mat[m,i])
        end;
        Inc(m);
        iter := 1
      end
    end {of j loop}
  Until iter=0;

  low := m;
  high := k;
  for i := m to k do  scal[i] := ONE;

  Repeat
    iter := 0;
    for i := m to k do
    begin
      c:=ZERO; r:=ZERO;
      for j := m to k do
        if j <> i then
        begin
          c := c + ABS(mat[j,i]);
          r := r + ABS(mat[i,j])
        end;
      g := r / basis;
      f := ONE;
      s := c + r;

      while c < g do
      begin
        f := f * basis;
        c := c * b2
      end;

      g := r * basis;
      while c >= g do
      begin
        f := f / basis;
        c := c / b2
      end;

      if (c + r) / f < 0.95 * s then
      begin
        g := ONE / f;
        scal[i] := scal[i] * f;
        iter := 1;
        for j := m to n-1 do  mat[i,j] := mat[i,j] * g;
        for j := 0 to k do  mat[j,i] := mat[j,i] * f
      end
    end {of i loop}
  Until iter=0
End;  {balance}


Procedure balback            {reverse balancing ........................}
                  (n:      integer;      { Dimension of matrix .........}
                   low:    integer;      { first nonzero row ...........}
                   high:   integer;      { last nonzero row ............}
                   scal:   Real_Vector;  { Scaling data ................}
               VAR eivec:  Square_Matrix { Eigenvectors ................}
                  );
 {====================================================================*
 *                                                                    *
 *  balback reverses the balancing of balance for the eigenvactors.   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ----------------                                                 *
 *      n        integer;  ( n > 0 )                                  *
 *               Dimension of mat                                     *
 *      low      integer;                                             *
 *      high     integer;   see balance                               *
 *      eivec    Square_Matrix;                                       *
 *               Matrix of eigenvectors, as computed in  qr2          *
 *      scal     Real_Vector;                                         *
 *               Scaling data from  balance                           *
 *                                                                    *
 *   Output parameter:                                                *
 *   ----------------                                                 *
 *      eivec    Square_Matrix;                                       *
 *               Non-normalized eigenvectors of the original matrix   *
 *                                                                    *
 *   Functions:  RSWAP                                                *
 *   ---------                                                        *
 *                                                                    *
 *====================================================================}
Var
  i, j, k: Integer;
  s: REAL_AR;
Begin
  for i := low to high do
  begin
    s := scal[i];
    for j := 0 to n-1 do eivec[i,j] := eivec[i,j] * s
  end;

  for i := low-1 downto 0 do
  begin
    k := Round(scal[i]);
    if k <> i then
      for j := 0 to n-1 do RSWAP(eivec[i,j], eivec[k,j])
  end;

  for i := high + 1 to n-1 do
  begin
    k := Round(scal[i]);
    if k <> i then
      for j := 0 to n-1 do RSWAP(eivec[i,j], eivec[k,j])
  end
End;


PROCEDURE elmhes        {reduce matrix to upper Hessenberg form ....}
                (n:    integer;       {Dimension of matrix .........}
                 low:  integer;       {first nonzero row ...........}
                 high: integer;       {last nonzero row ............}
             VAR mat:  Square_Matrix; {input/output matrix .........}
             VAR perm: Integer_Vector {Permutation vector ..........}
                );
 {====================================================================*
 *                                                                    *
 *  elmhes transforms the matrix mat to upper Hessenberg form.        *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer;  ( n > 0 )                                  *
 *               Dimension of mat                                     *
 *      low      integer;                                             *
 *      high     integer; see  balance                                *
 *      mat      Square_Matrix;                                       *
 *               n x n matrix                                         *
 *                                                                    *
 *   Output parameter:                                                *
 *   =================                                                *
 *      mat      Square_Matrix;                                       *
 *               upper Hessenberg matrix; additional information on   *
 *               the transformation is stored in the lower triangle   *
 *      perm     Integer_Vector;                                      *
 *               Permutation vector for elmtrans                      *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions:  RSWAP                                                *
 *   =========                                                        *
 *                                                                    *
 *====================================================================}
Var
  i, j, m: Integer;
  x, y: REAL_AR;
Begin
  for m := low + 1 to high-1 do
  begin
    i := m;
    x := ZERO;
    for j := m to high do
      if ABS(mat[j,m-1]) > ABS (x) then
      begin
        x := mat[j,m-1];
        i := j
      end;

    perm[m] := i;
    if i <> m then
    begin
      for j := m - 1 to n-1 do  RSWAP(mat[i,j], mat[m,j]);
      for j := 0 to high do  RSWAP( mat[j,i], mat[j,m])
    end;

    if x <> ZERO then
    begin
      for i := m + 1 to high do
      begin
        y := mat[i,m-1];
        if y <> ZERO then
        begin
          y := y / x;
          mat[i,m-1] := y;
          for j := m to n-1 do  mat[i,j] := mat[i,j] - y * mat[m,j];
          for j := 0 to high do  mat[j,m] := mat[j,m] + y * mat[j,i]
        end
      end {of i loop}
    end {of x <> ZERO}
  end {of m loop}

End;

PROCEDURE elmtrans             {copy to Hessenberg form .................}
                  (n:      integer;        {Dimension of matrix .........}
                   low:    integer;        {first nonzero row ...........}
                   high:   integer;        {last nonzero row ............}
                   mat:    Square_Matrix;  {input matrix ................}
                   perm:   Integer_Vector; {row permutations ............}
               VAR h:      Square_Matrix   {Hessenberg matrix ...........}
                   );
 {====================================================================*
 *                                                                    *
 *  elmtrans copies the Hessenberg matrix stored in mat to h.         *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        integer;  ( n > 0 )                                  *
 *               Dimension of  mat and eivec                          *
 *      low      integer;                                             *
 *      high     integer; see  balance                                *
 *      mat      Square_Matrix;                                       *
 *               n x n input matrix                                   *
 *      perm     Integer_Vector;                                      *
 *               Permutation data from  elmhes                        *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      h        Square_Matrix;                                       *
 *               Hessenberg matrix                                    *
 *                                                                    *
 *====================================================================}
Var
  k, i, j: Integer;
Begin
  for i := 0 to n-1 do
  begin
    for k := 0  to n-1 do  h[i,k] := ZERO;
    h[i,i] := ONE
  end;

  for i := high - 1 downto low+1 do
  begin
    j := perm[i];
    for k := i + 1 to high do h[k,i] := mat[k,i-1];
    if i <> j then
    begin
      for k := i to high do
      begin
        h[i,k] := h[j,k];
        h[j,k] := ZERO
      end;
      h[j,i] := ONE
    end
  end
End;

Procedure Comdiv           { Complex division ..........................}
           (
            ar:      REAL_AR;        { Real part of numerator ..........}
            ai:      REAL_AR;        { Imaginary part of numerator .....}
            br:      REAL_AR;        { Real part of denominator ........}
            bi:      REAL_AR;        { Imaginary part of denominator ...}
        VAR cr:      REAL_AR;        { Real part of quotient ...........}
        VAR ci:      REAL_AR;        { Imaginary part of quotient ......}
        VAR rc:      Integer         { return code .....................}
           );
 {====================================================================*
 *                                                                    *
 *  Complex division  c = a / b                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      ar,ai    double   ar, ai;                                     *
 *               real, imaginary parts of numerator                   *
 *      br,bi    double   br, bi;                                     *
 *               real, imaginary parts of denominator                 *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      cr,ci    double   *cr, *ci;                                   *
 *               real , imaginary parts of the quotient               *
 *                                                                    *
 *   Return code value:                                               *
 *   =================                                                *
 *      = 0      ok                                                   *
 *      = 1      division by 0                                        *
 *                                                                    *
 *====================================================================}
Label 10;
Var tmp: REAL_AR;
Begin
  if (br = ZERO) AND (bi = ZERO) then
  begin
    rc := 1;
    goto 10
  end;
  if abs(br) > abs(bi) then
  begin
    tmp  := bi / br;
    br   := tmp * bi + br;
    cr  := (ar + tmp * ai) / br;
    ci  := (ai - tmp * ar) / br
  end
  else
  begin
    tmp  := br / bi;
    bi   := tmp * br + bi;
    cr  := (tmp * ar + ai) / bi;
    ci  := (tmp * ai - ar) / bi
 end;
 rc := 0;
10: End; {Comdiv}

FUNCTION Comabs           {Complex absolute value ....................}
              (
               ar: REAL_AR;         {Real part .......................}
               ai: REAL_AR          {Imaginary part ..................}
              ): REAL_AR;
 {====================================================================*
 *                                                                    *
 *  Complex absolute value of   a                                     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      ar,ai    REAL_AR;                                             *
 *               Real, imaginary parts of  a                          *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      Absolute value of a (REAL_AR)                                 *
 *                                                                    *
 *   Procedure used:    RSWAP                                         *
 *   ==============                                                   *
 *                                                                    *
 *====================================================================}
Label 10;
Begin
  if (ar = ZERO) and (ai = ZERO) then
  begin
    Comabs := ZERO;
    goto 10  {return}
  end;

  ar := ABS(ar);
  ai := ABS(ai);

  if ai > ar then                                   {Switch  ai and ar}
    RSWAP(ai, ar);

  if ai = ZERO then
    Comabs := ar
  else
    Comabs := ar * SQRT(ONE + ai / ar * ai / ar);
10:End;


Procedure hqrvec          {compute eigenvectors ......................}
                  (n:      integer;       {Dimension of matrix .......}
                   low:    integer;       {first nonzero row .........}
                   high:   integer;       {last nonzero row ..........}
                   h:      square_matrix; {upper Hessenberg matrix ...}
                   wr:     real_vector;   {Real parts of evalues .....}
                   wi:     real_vector;   {Imaginary parts of evalues }
               VAR eivec:  square_matrix; {Eigenvectors ..............}
               VAR rc:     integer        {return code ...............}
                  );
 {====================================================================*
 *                                                                    *
 *  hqrvec computes the eigenvectors for the eigenvalues found in hqr2*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of  mat and eivec, number of eigenvalues.  *
 *      low      int low;                                             *
 *      high     int high; see  balance                               *
 *      h        REAL   *h[n];                                        *
 *               upper Hessenberg matrix                              *
 *      wr       REAL   wr[n];                                        *
 *               Real parts of the n eigenvalues.                     *
 *      wi       REAL   wi[n];                                        *
 *               Imaginary parts of the n eigenvalues.                *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      eivec    REAL   *eivec[n];                                    *
 *               Matrix, whose columns are the eigenvectors           *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      =  0     all ok                                               *
 *      =  1     h is the zero matrix.                                *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   function in use  :                                               *
 *   ==================                                               *
 *                                                                    *
 *      int   comdiv(): complex division                              *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used:    MACH_EPS                                      *
 *   ==============                                                   *
 *                                                                    *
 *====================================================================}
Label 10;
Var
  code, i, j, k: integer;
  l, m, en, na: integer;
  p, q, r, s, t, w, x, y, z, ra, sa, vr, vi, norm, temp: REAL_AR;

Begin

  r:=ZERO; s:=ZERO; z:=ZERO; norm:=ZERO;

  for i := 0 to n-1 do                        {find norm of h}
    for j := i to n-1 do
      norm := norm + ABS(h[i,j]);

  if norm = ZERO then
  begin
    rc := 1;                                  {zero matrix}
    goto 10                                   {return}
  end;

  for en := n-1 downto 0 do                   {transform back}
  begin
    p := wr[en];
    q := wi[en];
    na := en - 1;
    if q = ZERO then
    begin
      m := en;
      h[en,en] := ONE;
      for i := na downto 0 do
      begin
        w := h[i,i] - p;
        r := h[i,en];
        for j := m to na do  r := r + h[i,j] * h[j,en];
        if wi[i] < ZERO then
        begin
          z := w;
          s := r
        end
        else
        begin
          m := i;
          if wi[i] = ZERO then
          begin
            if w <> ZERO then temp:=w else temp:=MACH_EPS * norm;
            h[i,en] := -r/temp            
          end
          else
          begin  {Solve the linear system:
                | w   x |  | h[i][en]   |   | -r |
                |       |  |            | = |    |
                | y   z |  | h[i+1][en] |   | -s |  }
             x := h[i,i+1];
             y := h[i+1,i];
             q := SQR(wr[i] - p) + SQR(wi[i]);
             h[i,en] := (x * s - z * r) / q;
             t := h[i,en];
            if ABS(x) > ABS(z) then temp:=(-r -w * t) / x else temp:=(-s -y * t) / z;
            h[i+1,en] := temp
          end
        end {end case wi[i] >= 0}
      end {end i loop}
    end  {end if q = 0}
    else if q < ZERO then
    begin
      m := na;
      if ABS(h[en,na]) > ABS(h[na,en]) then
      begin
        h[na,na] := - (h[en,en] - p) / h[en,na];
        h[na,en] := - q / h[en,na];
      end
      else
        Comdiv(-h[na,en], ZERO, h[na,na]-p, q, h[na,na], h[na,en],code);

      h[en,na] := ONE;
      h[en,en] := ZERO;
      for i := na - 1 downto 0 do
      begin
        w := h[i,i] - p;
        ra := h[i,en];
        sa := ZERO;
        for j := m to na do
        begin
          ra := ra + h[i,j] * h[j,na];
          sa := sa + h[i,j] * h[j,en]
        end;

        if wi[i] < ZERO then
        begin
          z := w;
          r := ra;
          s := sa
        end
        else
        begin
          m := i;
          if wi[i] = ZERO then
            Comdiv (-ra, -sa, w, q, h[i,na], h[i,en],code)
          else
          begin
            {solve complex linear system:
            | w+i*q     x | | h[i][na] + i*h[i][en]  |   | -ra+i*sa |
            |             | |                        | = |          |
            |   y    z+i*q| | h[i+1][na]+i*h[i+1][en]|   | -r+i*s   |  }
            x := h[i,i+1];
            y := h[i+1,i];
            vr := SQR (wr[i] - p) + SQR (wi[i]) - SQR (q);
            vi := TWO * q * (wr[i] - p);
            if (vr = ZERO) AND (vi = ZERO) then
              vr := MACH_EPS * norm * (ABS (w) + ABS (q) + ABS (x) + ABS (y) + ABS (z));

            Comdiv (x * r - z * ra + q * sa, x * s - z * sa -q * ra, vr, vi, h[i,na], h[i,en], code);
            if ABS (x) > ABS (z) + ABS (q) then
            begin
              h[i+1,na] := (-ra - w * h[i,na] + q * h[i,en]) / x;
              h[i+1,en] := (-sa - w * h[i,en] - q * h[i,na]) / x
            end
            else
              Comdiv (-r - y * h[i,na], -s - y * h[i,en], z, q, h[i+1,na], h[i+1,en],code)
          end {wi[i] > 0}
        end {wi[i] >= 0}
      end {of i loop}
    end {of else if q < 0}
  end; {of en loop}

  for i := 0 to n-1 do                {Eigenvectors for the evalues for}
    if (i < low) or (i > high) then   {rows < low  and rows  > high    }
      for k := i + 1 to n-1 do
        eivec[i,k] := h[i,k];

  for j := n-1 downto low do
  begin
    if j<=high then m:=j else j:=high;
    if wi[j] < ZERO then
    begin
      l:=j-1;
      for i := low to high do
      begin
        y:=ZERO; z:=ZERO;
        for k := low to m do
        begin
          y := y + eivec[i,k] * h[k,l];
          z := z + eivec[i,k] * h[k,j]
        end;
        eivec[i,l] := y;
        eivec[i,j] := z
      end
    end
    else
      if wi[j] = ZERO then
      begin
        for i := low to high do
        begin
          z:=ZERO;
          for k := low to m do
            z := z + eivec[i,k] * h[k,j];
          eivec[i,j] := z
        end
      end
  end; {of j loop}

  rc := 0;
10: End;

{utility procedure used only for debugging file}
PROCEDURE WriteMat(VAR fp:TEXT;caption:string;A:Square_Matrix;n:Integer);
Var i,j: Integer;
Begin
  Writeln(fp,' ',caption);
  For i:=0 to n-1 do
  begin
    For j:=0 to n-1 do write(fp,A[i,j]:10:6);
    Writeln(fp)
  end
End;

PROCEDURE hqr2             {compute eigenvalues .......................}
                (vec:    Boolean;        {switch for computing evectors}
                 n:      integer;        {Dimension of matrix .........}
                 low:    integer;        {first nonzero row ...........}
                 high:   integer;        {last nonzero row ............}
                 h:      Square_Matrix;  {Hessenberg matrix ...........}
             VAR wr:     Real_Vector;    {Real parts of eigenvalues ...}
             VAR wi:     Real_Vector;    {Imaginary parts of evalues ..}
             VAR eivec:  Square_Matrix;  {Matrix of eigenvectors ......}
             VAR cnt:    Integer_Vector; {Iteration counter ...........}
             VAR rc:     Integer         {return code .................}
                );
 {*********************************************************************
 *                                                                    *
 * hqr2 computes the eigenvalues and (if vec = True) the eigenvectors *
 * of an  n * n upper Hessenberg matrix.                              *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Control parameter:                                               *
 *   -----------------                                                *
 *      vec      Boolean;                                             *
 *       = False compute eigenvalues only                             *
 *       = True  compute all eigenvalues and eigenvectors             *
 *                                                                    *
 *   Input parameters:                                                *
 *   ----------------                                                 *
 *      n        integer;  ( n > 0 )                                  *
 *               Dimension of  h and eivec,                           *
 *               length of the real parts vector  wr and of the       *
 *               imaginary parts vector  wi of the eigenvalues.       *
 *      low      integer;                                             *
 *      high     integer;  see balance                                *
 *      h        Square_Matrix;                                       *
 *               upper Hessenberg matrix as output of Elmhes          *
 *               (destroyed in the process).                          *
 *                                                                    *
 *   Output parameters:                                               *
 *   -----------------                                                *
 *      eivec    Square_Matrix;  (only if vec = Tue)                  *
 *               Matrix, which for vec = True contains the            *
 *               eigenvectors as follows:                             *
 *               For real eigebvalues the corresponding column        *
 *               contains the corresponding eigenvactor, while for    *
 *               complex eigenvalues the corresponding column contains*
 *               the real part of the eigenvactor with its imaginary  *
 *               part is stored in the subsequent column of eivec.    *
 *               The eigenvactor for the complex conjugate eigenvactor*
 *               is given by the complex conjugate eigenvactor.       *
 *      wr       Real_Vector;                                         *
 *               Real part of the n eigenvalues.                      *
 *      wi       Real_Vector;                                         *
 *               Imaginary parts of the eigenvalues                   *
 *      cnt      Integer_Vector;                                      *
 *               vector of iterations used for each eigenvalue.       *
 *               For a complex conjugate eigenvalue pair the second   *
 *               entry is negative.                                   *
 *                                                                    *
 *   Return code:                                                     *
 *   ------------                                                     *
 *      =   0    all ok                                               *
 *      = 4xx    Iteration maximum exceeded when computing evalue xx  *
 *      =  99    zero  matrix                                         *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   functions in use:                                                *
 *   ----------------                                                 *
 *                                                                    *
 *      hqrvec():  reverse transform for eigenvectors                 *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Constants used:     MACH_EPS, MAXIT                              *
 *   --------------                                                   *
 *                                                                    *
 *********************************************************************}
Label 10, 12, 15, 20, 30;
Var
  i, j, ll: integer;
  na, en, iter, k, l, m: integer;
  p, q, r, s, t, w, x, y, z: REAL_AR;
  debug:Boolean; fp:TEXT;
Begin
  debug:=False;
  if debug then
  begin
    Assign(fp,'hqr2.lst'); Rewrite(fp);
    WriteMat(fp,'Matrix h in begin hqr2:',h,n);
    Writeln(fp,' low=',low,' high=',high);
  end;

  p:=ZERO; q:=ZERO; r:=ZERO; 
  for i := 0 to n-1 do
    if (i < low) or (i > high) then
    begin
      wr[i] := h[i,i];
      wi[i] := ZERO;
      cnt[i] := 0
    end;

  en := high;
  t := ZERO;

  while en >= low do
  begin
    iter := 0;
    na := en - 1;

    while(True) do
    begin
      ll:=999;                          
      for l := en downto low+1 do               {search for small   }
      begin                                     {subdiagonal element}
        if (debug) then
          WriteLn(fp,'l=',l,' a=',ABS(h[l,l-1]),' b=',MACH_EPS*(ABS(h[l-1,l-1])+ABS(h[l,l])) );
        if ABS(h[l,l-1]) <= MACH_EPS * (ABS(h[l-1,l-1]) + ABS(h[l,l])) then
        begin
          ll:=l;     {save current index}
          goto 10    {exit l loop}
        end
      end;
10:   if ll<>999 then l:=ll else l:=0;   {restore l}
      if debug then
      begin
        Writeln(fp,' iter=',iter,' l=',l,' en=',en,' na=',na);
        WriteMat(fp,'Matrix h in hqr2 label 10:',h,n)
      end;

      x := h[en,en];
      if l = en then                            {found one evalue}
      begin
        wr[en] := x + t;
        h[en,en] := x + t;
        wi[en] := ZERO;
        cnt[en] := iter;
        Dec(en);
        goto 15;                                {exit from loop while(True) }
      end;

      y := h[na,na];
      w := h[en,na] * h[na,en];

      if l = na then                            {found two evalues}
      begin
        p := (y - x) * 0.5;
        q := p * p + w;
        z := SQRT(ABS(q));
        x := x + t;
        h[en,en] := x + t;
        h[na,na] := y + t;
        cnt[en] := -iter;
        cnt[na] := iter;
        if q >= ZERO then
        begin                                    {real eigenvalues}
          if p<ZERO then z:=p-z else z:=p+z;
          wr[na] := x + z;
          wr[en] := x - w / z;
          s := w - w / z;
          wi[na] := ZERO;
          wi[en] := ZERO;
          x := h[en,na];
          r := SQRT (x * x + z * z);

          if (vec) then
          begin
            p := x / r;
            q := z / r;
            for j := na to n-1 do
            begin
              z := h[na,j];
              h[na,j] := q * z + p * h[en,j];
              h[en,j] := q * h[en,j] - p * z
            end;

            for i := 0 to en do
            begin
              z := h[i,na];
              h[i,na] := q * z + p * h[i,en];
              h[i,en] := q * h[i,en] - p * z
            end;

            for i := low to high do
            begin
              z := eivec[i,na];
              eivec[i,na] := q * z + p * eivec[i,en];
              eivec[i,en] := q * eivec[i,en] - p * z
            end
          end {if (vec) }
        end {if q >= ZERO}
        else                                  {pair of complex}
        begin                                 {conjugate evalues}
          wr[na] := x + p;
          wr[en] := x + p;
          wi[na] :=   z;
          wi[en] := - z
        end;

        en := en - 2;
        goto 15;                              {exit while(True) }
      end; {if l = na}

      if iter >= MAXIT then
      begin
        cnt[en] := MAXIT + 1;
        rc:=en;
        writeln(' stop at iter >= MAXIT.');
        goto 20                               {MAXIT Iterations: close debugging file if open}
      end;                                    {and return }

      if (iter <> 0) and (iter MOD 10 = 0) then
      begin
        t := t + x;
        for i := low to en do h[i,i] := h[i,i] - x;
        s := ABS(h[en,na]) + ABS(h[na,en-2]);
        x := 0.75 * s; y:=x;
        w := -0.4375 * s * s
      end;

      Inc(iter);

      for m := en - 2 downto l do
      begin
        z := h[m,m];
        r := x - z;
        s := y - z;
        p := ( r * s - w ) / h[m+1][m] + h[m][m+1];
        q := h[m + 1][m + 1] - z - r - s;
        r := h[m + 2][m + 1];
        s := ABS (p) + ABS (q) + ABS (r);
        p := p / s;
        q := q / s;
        r := r / s;
        if m = l then goto 12;
        if ( ABS (h[m,m-1]) * (ABS (q) + ABS (r)) <= MACH_EPS * ABS (p)
                 * ( ABS (h[m-1][m-1]) + ABS (z) + ABS (h[m+1][m+1])) ) then
          goto 12                {exit m loop}
      end;

12:   for i := m + 2 to en do h[i,i-2] := ZERO;
      for i := m + 3 to en do h[i,i-3] := ZERO;

      for k := m to na do
      begin
        if k <> m then           {double  QR step, for rows l to en }
        begin                    {and columns m to en               }
          p := h[k,k-1];
          q := h[k+1,k-1];
          if k<>na then r := h[k+2,k-1] else r := ZERO;
          x := ABS (p) + ABS (q) + ABS (r);
          if x = ZERO then goto 30;                  {next k}
          p := p / x;
          q := q / x;
          r := r / x
        end;
        s := SQRT (p * p + q * q + r * r);
        if p < ZERO then s := -s;

        if k <> m then
          h[k,k-1] := -s * x
        else if l <> m then
          h[k,k-1] := -h[k,k-1];
        p := p + s;
        x := p / s;
        y := q / s;
        z := r / s;
        q := q / p;
        r := r / p;

        for j := k to n-1 do                   {modify rows}
        begin
          p := h[k,j] + q * h[k+1,j];
          if k <> na then
          begin
            p := p + r * h[k+2,j];
            h[k+2,j] := h[k+2,j] - p * z
          end;
          h[k+1,j] := h[k+1,j] - p * y;
          h[k,j]   := h[k,j] - p * x
        end;

        if k+3 < en then j:=k+3 else j:=en;
        for i := 0 to j do                     {modify columns}
        begin
          p := x * h[i,k] + y * h[i,k+1];
          if k <> na then
          begin
            p := p + z * h[i,k+2];
            h[i,k+2] := h[i,k+2] - p * r
          end;
          h[i,k+1] := h[i,k+1] - p * q;
          h[i,k]   := h[i,k] - p
        end;

        if (vec) then     {if eigenvectors are needed ..................}
        begin
          for i := low to high do
          begin
            p := x * eivec[i,k] + y * eivec[i,k+1];
            if k <> na then
            begin
              p := p + z * eivec[i,k+2];
              eivec[i,k+2] := eivec[i,k+2] - p * r
            end;
            eivec[i,k+1] := eivec[i,k+1] - p * q;
            eivec[i,k]   := eivec[i,k] - p
          end
        end;
30:   end; {of k loop}

    end;  {while(True)}

15:end; {while en >= low                          All evalues found    }

  if (vec) then                             {transform evectors back  }
    hqrvec (n, low, high, h, wr, wi, eivec,rc)
  else
    rc:=0;

20: if debug then Close(fp);
End; {hqr2}


PROCEDURE norm_1            {normalize eigenvectors to have one norm 1 .}
                  (n:    integer;       {Dimension of matrix ...........}
               VAR v:    Square_Matrix; {Matrix with eigenvektors ......}
                   wi:   Real_Vector;   {Imaginary parts of evalues ....}
               VAR rc:   integer        {return code ...................}
                  );
 {*********************************************************************
 *                                                                    *
 *  norm_1 normalizes the one norm of the column vectors in v.        *
 *  (special attention to complex vectors in v  is given)             *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Input parameters:                                                *
 *   ----------------                                                 *
 *      n        integer; ( n > 0 )                                   *
 *               Dimension of matrix v                                *
 *      v        Square_Matrix;                                       *
 *               Matrix of (not normalized) eigenvectors              *
 *      wi       Real_Vector;                                         *
 *               Imaginary parts of the eigenvalues                   *
 *                                                                    *
 *   Output parameter:                                                *
 *   ----------------                                                 *
 *      v        Square_Matrix;                                       *
 *               Matrix with normalized eigenvectors                  *
 *                                                                    *
 *   Return code:                                                     *
 *   -----------                                                      *
 *      = 0      all ok                                               *
 *      = 1      n < 1                                                *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   functions used:                                                  *
 *   --------------                                                   *
 *      REAL_AR  comabs():  complex absolute value                    *
 *               comdiv():  complex division                          *
 *                                                                    *
 *********************************************************************}
Label 10;
Var
  i,j: integer;
  maxi, tr, ti: REAL_AR;
Begin
  if n < 1 then
  begin
    rc := 1;
    goto 10
  end;

  for j := 0 to n-1 do
  begin
    if wi[j] = ZERO then
    begin
      maxi := v[0,j];
      for i := 1 to n-1 do
        if ABS(v[i,j]) > ABS(maxi) then maxi := v[i,j];

      if maxi <> ZERO then
      begin
        maxi := ONE / maxi;
        for i := 0 to n-1 do v[i,j] := v[i,j] * maxi
      end
    end
    else
    begin
      tr := v[0][j];
      ti := v[0][j+1];
      for i := 1 to n-1 do
        if Comabs(v[i,j], v[i,j+1]) > Comabs(tr, ti) then
        begin
          tr := v[i,j];
          ti := v[i,j+1]
        end;

      if (tr <> ZERO) or (ti <> ZERO) then
        for i := 0 to n-1 do
          Comdiv(v[i,j], v[i,j+1], tr, ti, v[i,j], v[i,j+1],rc);

      Inc(j)                                           {raise j by two}
    end
  end; {of j loop}
  rc := 0;
10: End;

{utility procedure used only for debugging file}
PROCEDURE WriteVec(VAR fp:TEXT;caption:string;V:Integer_Vector;n:Integer);
Var i: Integer;
Begin
  Writeln(fp,' ',caption);
  For i:=0 to n-1 do write(fp,V[i]:6);
  Writeln(fp)
End;

PROCEDURE eigen (         {Compute all evalues/evectors of a matrix ..}
           vec:    Boolean;         {switch for computing evectors ...}
           ev_norm:Boolean;         {normalize Eigenvectors? .........}
           n:      integer;         {size of matrix ..................}
           mat:    Square_Matrix;   {input matrix ....................}
       VAR eivec:  Square_Matrix;   {Eigenvectors ....................}
       VAR valre:  Real_Vector;     {real parts of eigenvalues .......}
       VAR valim:  Real_Vector;     {imaginary parts of eigenvalues ..}
       VAR cnt:    Integer_Vector;  {Iteration counter ...............}
       VAR rc:     Integer          {return code .....................}
          );
 {*********************************************************************
 *                                                                    *
 * The function  eigen  determines all eigenvalues and (if desired)   *
 * all eigenvectors of a real square  n * n  matrix via the QR method *
 * in the version of  Martin, Parlett, Peters, Reinsch and Wilkinson. *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Litterature:                                                     *
 *   -----------                                                      *
 *      1) Peters, Wilkinson: Eigenvectors of real and complex        *
 *         matrices by LR and QR triangularisations,                  *
 *         Num. Math. 16, p.184-204, (1970); [PETE70]; contribution   *
 *         II/15, p. 372 - 395 in [WILK71].                           *
 *      2) Martin, Wilkinson: Similarity reductions of a general      *
 *         matrix to Hessenberg form, Num. Math. 12, p. 349-368,(1968)*
 *         [MART 68]; contribution II,13, p. 339 - 358 in [WILK71].   *
 *      3) Parlett, Reinsch: Balancing a matrix for calculations of   *
 *         eigenvalues and eigenvectors, Num. Math. 13, p. 293-304,   *
 *         (1969); [PARL69]; contribution II/11, p.315 - 326 in       *
 *         [WILK71].                                                  *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Control parameters:                                              *
 *   ------------------                                               *
 *      vec      boolean;                                             *
 *               call for eigen :                                     *
 *      = False  compute eigenvalues only                             *
 *      = True   compute all eigenvalues and eigenvectors             *
 *      ortho    boolean flag that shows if transformation of mat to  *
 *               Hessenberg form shall be done orthogonally by        *
 *               `orthes' (flag True) or elementarily by `elmhes'     *
 *               (flag False). The Householder matrices used in       *
 *               orthogonal transformation have the advantage of      *
 *               preserving the symmetry of input matrices.           *
 *               (only elmhes is implemented here).                   *
 *      ev_norm  boolean flag that shows if Eigenvectors shall be     *
 *               normalized (flag True) or not (flag False).          *
 *                                                                    *
 *   Input parameters:                                                *
 *   ----------------                                                 *
 *      n        integer; ( n > 0 )                                   *
 *               size of matrix, number of eigenvalues                *
 *      mat      Square_Matrix;                                       *
 *               input matrix                                         *
 *                                                                    *
 *   Output parameters:                                               *
 *   -----------------                                                *
 *      eivec    Square_Matrix;     (only if vec = True)              *
 *               matrix, if  vec = True  that holds the eigenvectors  *
 *               thus :                                               *
 *               If the jth eigenvalue of the matrix is real then the *
 *               jth column is the corresponding real eigenvector;    *
 *               if the jth eigenvalue is complex then the jth column *
 *               of eivec contains the real part of the eigenvector   *
 *               while its imaginary part is in column j+1.           *
 *               (the j+1st eigenvector is the complex conjugate      *
 *               vector.)                                             *
 *      valre    Real_Vector;                                         *
 *               Real parts of the eigenvalues.                       *
 *      valim    Real_Vector;                                         *
 *               Imaginary parts of the eigenvalues                   *
 *      cnt      Integer_Vector;                                      *
 *               vector containing the number of iterations for each  *
 *               eigenvalue. (for a complex conjugate pair the second *
 *               entry is negative).                                  *
 *                                                                    *
 *   Return code:                                                     *
 *   ------------                                                     *
 *      =   0    all ok                                               *
 *      =   1    n < 1 or other invalid input parameter               *
 *      =   2    insufficient memory                                  *
 *      = 10x    error x from balance()                               *
 *      = 20x    error x from elmh()                                  *
 *      = 30x    error x from elmtrans()   (for vec = Trie only)      *
 *      = 4xx    error xx from hqr2()                                 *
 *      = 50x    error x from balback()    (for vec = True only)      *
 *      = 60x    error x from norm_1()     (for vec = True only)      *
 *                                                                    *
 * ------------------------------------------------------------------ *
 *                                                                    *
 *   Procedures used:                                                 *
 *   ---------------                                                  *
 *                                                                    *
 *   balance (): Balancing of an  n x n  matrix                       *
 *   elmhes ():  Transformation to upper Hessenberg form              *
 *   elmtrans(): intialize eigenvectors                               *
 *   hqr2 ():    compute eigenvalues/eigenvectors                     *
 *   balback (): Reverse balancing to obtain eigenvectors             *
 *   norm_1 ():  Normalize eigenvectors                               *
 *                                                                    *
 *********************************************************************}
Label 10;
Var
  i, low, high: Integer;
  d, scale: Real_Vector;
  perm: Integer_Vector;
  fp: TEXT;  {debugging file}
  debug: Boolean;
Begin
  debug:=False;
  if n < 1 then                              {case n < 1 .............}
  begin
    rc := 1;
    goto 10
  end;                                        

  for i := 0 to n-1 do
  begin
    cnt[i] := 0;
    d[i]:=0.0
  end;

  if n = 1 then                              {case n = 1 .............}
  begin
    eivec[0,0] := ONE;
    valre[0]    := mat[0,0];
    valim[0]    := ZERO;
    rc := 0;
    goto 10
  end;
                                              {balance mat for nearly}
  Balance(n, mat, scale, low, high);          {equal row and column  }
                                              {one norms             }

  Elmhes(n, low, high, mat, perm);            {reduce mat to upper   }
                                              {Hessenberg form       }

  if (vec) then                               {initialize eivec      }
    elmtrans (n, low, high, mat, perm, eivec);

  if debug then
  begin
    Assign(fp,'eigen.lst'); Rewrite(fp);
    WriteMat(fp,' Matrix A after Elmtrans:',mat,n);
    WriteVec(fp,' Vector perm after Elmhes:',perm,n);
    WriteMat(fp,' Matrix Eivec after Elmhes:',eivec,n)
  end;

  hqr2 (vec, n, low, high, mat,               {execute Francis QR    }
             valre, valim, eivec, cnt,rc);    {algorithm to obtain   }

  if rc<>0 then
  begin
    rc:=2;
    goto 10
  end;

  if (vec) then
  begin
    balback (n, low, high, scale, eivec);     {reverse balancing if  }
                                              {eigenvaectors are to  }
                                              {be determined         }
    if (ev_norm) then
      norm_1 (n, eivec, valim, rc);           {normalize eigenvectors}

    if rc<>0 then
    begin
      rc:=3;
      goto 10
    end
  end;

  rc := 0;
  if debug then Close(fp);
10: End;

END.
{ ------------------------ END feigen.pas --------------------------- }