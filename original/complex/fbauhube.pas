{ -------------------------- UNIT fbauhube.pas ----------------------------
* Ref.: "Numerical algorithms with C, By Gisela Engeln-Muellges and       *
*        Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].                  *
*                                                                         *
*                                       TPW Release By J-P Moreau, Paris. *
*                                               (www.jpmoreau.fr)         *
--------------------------------------------------------------------------}
UNIT FBauhube;

INTERFACE

Uses WinCrt;

Const
  ITERMAX = 1000;                    { Maximal number of function
                                       evaluations per root           }
  MACH_EPS = 1E-15;                  { Small real number              }
  EPSROOT = 31.622E6;                { SQRT(MACH_EPS)                 } 
  EPS = 64.0 * MACH_EPS;             { Accuracy for functional value  }
  BETA = 8.0 * EPS;
  QR = 0.1;                          { Real and imaginary parts of the}
  QI = 0.9;                          { spiralization constant         }

  NMAX = 100;

Type VEC = Array[0..NMAX] of Double;


  Function bauhub (         {Bauhuber's method for complex polynomials}
            real0,                 { Are the coefficients real ? .....}
            scale,                 { Scaling ? .......................}
            n: integer;            { degree of polynomial ............}
            ar,                    { Real parts of coefficients ......}
            ai: VEC;               { Imaginary parts, coefficients ...}
        Var rootr,                 { Real parts of roots .............}
            rooti,                 { Imaginary parts of roots ........}
            absf: VEC              { Absolute value of function values}
           ):Integer;

IMPLEMENTATION

  Procedure scpoly (n: Integer;      { length of vector .............}
                Var ar,              { Real part of the vector ......}
                    ai: VEC;         { Imaginary part of the vector .}
                Var scal: Double);   { Scaling factor ...............}
                    Forward;

  Function bauroot (n,                   { largest degree ...............}
                    iu: Integer;         { lowest degree ................}
                    ar,                  { Real parts of the coefficients}
                    ai: VEC;             { Imaginary parts, coefficients }
                Var x0r,                 { Real part of the root ........}
                    x0i:Double):Integer; { Imaginary part of the root ...}
                    Forward;

  Procedure chorner (n,              { highest degree in polynomial .}
                     iu: Integer;    { lowest degree in polynomial ..}
                     ar,             { Real parts of coefficients ...}
                     ai: VEC;        { Imaginary parts, coefficients }
                     xr,             { Real part of x ...............}
                     xi: Double;     { Imaginary part of x ..........}
                 Var pr,             { Real part of function value ..}
                     pi,             { Imaginary part of function v. }
                     p1r,            { Real part of first derivative }
                     p1i,            { Imaginary part, first deriv. .}
                     p2r,            { Real part, second derivative .}
                     p2i,            { Imaginary part, second deriv. }
                     rf1: Double);   { Error estimate for 1st deriv. }
                     Forward;

  Procedure polydiv (n,              { maximal degree ...............}
                     iu: Integer;    { minimal degree ...............}
                 Var ar,             { Real parts of coefficients ...}
                     ai: VEC;        { Imaginary parts, coefficients }
                     x0r,            { Real part of x ...............}
                     x0i: Double);   { Imaginary part of x ..........}
                     Forward;

 Function comabs        { Complex absolute value ....................}
              (
               ar,                { Real part .......................}
               ai: Double         { Imaginary part ..................}
              ): Double;
 {====================================================================*
 *                                                                    *
 *  Complex absolute value of   a                                     *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      ar,ai    REAL   ar, ai;                                       *
 *               Real, imaginary parts of  a                          *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      Absolute value of a                                           *
 *                                                                    *
 *   Macros used :    SQRT, ABS, SWAP                                 *
 *   =============                                                    *
 *                                                                    *
 *====================================================================}
 Var t: Double;
 Begin
 
  if (ar = 0.0) and (ai = 0.0) then
  begin
    comabs:=0.0;
    exit
  end;

  ar := ABS (ar);
  ai := ABS (ai);

  if ai > ar then                                  {Switch  ai and ar}
  begin
    t:=ai; ar:=ai; ai:=t
  end;

  if ai=0.0 then
    comabs:=ar
  else
    comabs:=ar * SQRT(1.0 + Sqr(ai/ar))
 End;

 Function comdiv         { Complex division ..........................}
           (
            ar,                    { Real part of numerator ..........}
            ai,                    { Imaginary part of numerator .....}
            br,                    { Real part of denominator ........}
            bi: Double;            { Imaginary part of denominator ...}
        Var cr,                    { Real part of quotient ...........}
            ci: Double             { Imaginary part of quotient ......}
           ): Integer;
 {====================================================================*
 *                                                                    *
 *  Complex division  c = a / b                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      ar,ai    REAL   ar, ai;                                       *
 *               Real, imaginary parts of numerator                   *
 *      br,bi    REAL   br, bi;                                       *
 *               Real, imaginary parts of denominator                 *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      cr,ci    REAL   *cr, *ci;                                     *
 *               Real , imaginary parts of the quotient               *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      ok                                                   *
 *      = 1      division by 0                                        *
 *                                                                    *
 *   Macro used :     ABS                                             *
 *   ============                                                     *
 *                                                                    *
 *====================================================================}
Label fin;
Var
  tmp: Double;
Begin
  if (br = 0.0) and (bi = 0.0) then
  begin
    comdiv:=1;
    goto fin
  end;

  if ABS (br) > ABS (bi) then
  begin
    tmp  := bi / br;
    br   := tmp * bi + br;
    cr   := (ar + tmp * ai) / br;
    ci   := (ai - tmp * ar) / br
  end
  else
  begin
    tmp  := br / bi;
    bi   := tmp * br + bi;
    cr   := (tmp * ar + ai) / bi;
    ci   := (tmp * ai - ar) / bi
  end;

  comdiv:=0;
fin: End;

 Function Max(x,y: Double): Double;
 Begin
   If x>=y then
     Max:=x
   Else
     Max:=y
 End;

 {calculate y^x}
 FUNCTION Power(y,x: Double): Double;
 BEGIN
   IF x<0 THEN EXIT;
   Power:=Exp(x*Ln(y))
 END;


  Function bauhub (         {Bauhuber's method for complex polynomials}
            real0,                 { Are the coefficients real ? .....}
            scale,                 { Scaling ? .......................}
            n: integer;            { degree of polynomial ............}
            ar,                    { Real parts of coefficients ......}
            ai: VEC;               { Imaginary parts, coefficients ...}
        Var rootr,                 { Real parts of roots .............}
            rooti,                 { Imaginary parts of roots ........}
            absf: VEC              { Absolute value of function values}
           ):Integer;
 {====================================================================*
 *                                                                    *
 *  bauhub uses Bauhuber's Method to find all real or complex roots   *
 *  of a polynomial of degree n :                                     *
 *                                               n-1           n      *
 *      P(x) = a[0] + a[1] * x + ... + a[n-1] * x    + a[n] * x ,     *
 *                                                                    *
 *  where a[i], i=0, ..., n, can be complex.                          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Application:                                                     *
 *   ===========                                                      *
 *      Find roots of arbitrary polynomials with complex coefficients.*
 *      If the polynomial roots are ill-condi=ditioned, i.e., if small*
 *      changes in the coefficients lead to large changes in the roots*
 *      the polynomial should not be scaled. Otherwise scaling helps  *
 *      with stability and performance.                               *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      real0    integer;                                             *
 *        = 0    Polynomial coefficients are complex                  *
 *       <> 0    Polynomial coefficients are real                     *
 *      scale    integer;                                             *
 *        = 0    no scaling                                           *
 *       <> 0    Scaling, see procedure scpoly()                      *
 *      n        integer;                                             *
 *               degree of the polynomial (must be >= 1)              *
 *      ar, ai   type VEC;                                            *
 *               Real and imaginary parts of the polynomial           *
 *               coefficients (ar[0], ..., ar[n])                     *
 *                                                                    *
 *   Output parameters:                                               *
 *   =================                                                *
 *      rootr    type VEC;   (Vector of length  n+1 )                 *
 *               rootr[0],..,rootr[n-1] are the real parts of the n   *
 *               roots                                                *
 *      rooti    type VEC;   (Vector of length n+1 )                  *
 *               rooti[0],..,rooti[n-1] are the imaginary parts       *
 *               of the roots                                         *
 *      absf     type VEC;                                            *
 *               absf[0],..,absf[n-1] are the magnitudes of the       *
 *               polynomial values at the computed roots              *
 *                                                                    *
 *   Return value:                                                    *
 *   ============                                                     *
 *      = 0      all is ok                                            *
 *      = 1      n < 1 or invalid input parameter                     *
 *      = 2      ar[n] = 0.0 and ai[n] = 0.0                          *
 *      = 3      Iteration maximum ITERMAX exceeded                   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions or Procedures used:                                    *
 *   ============================                                     *
 *     bauroot():   determines one root of the polynomial             *
 *     scpoly():    Scales the polynomial                             *
 *     chorner():   Evaluates the polynomial                          *
 *     polydiv():   factors off a root                                *
 *     comabs():    magnitude of a complex number                     *
 *                                                                    *
 *====================================================================}
Label fin; {return label}
Var
  i, res: Integer;
  x0r, x0i, tempr, tempi, t1, t2, t3, t4, t5: Double;
  scalefak: Double;
Begin

  scalefak := 1.0; 
  if n < 1 then     {n is to small!}
  begin
    bauhub:=1;
    goto fin
  end;

  if (ar[n] = 0.0) and (ai[n] = 0.0) then  { Leading coefficient must}
  begin                                    { differ from zero        }
    bauhub:=2;
    goto fin
  end;                              

  for i := 0 to n do              { store given coefficients in root }
  begin
    rootr[i] := ar[i];
    rooti[i] := ai[i];
    if i < n then absf[i] := 0.0
  end;

  scalefak := 1.0;
  if scale <> 0 then                  { Scale polynomial, if desired }
    scpoly (n, rootr, rooti, scalefak);

  x0r := 0.0;                                       { Starting value }
  x0i := 0.0;

  for i := 0 to n-1 do
  begin                                       { compute the ith root }
    res := bauroot (n, i, rootr, rooti, x0r, x0i);

    rootr[i] := scalefak * x0r;                         { store root }
    rooti[i] := scalefak * x0i;

    if res<>0 then
    begin
      bauhub := 3;                      { Iteration maximum reached? }
      goto fin
    end;

    { Polynomial value of input polynomial at (rootr[i], rooti[i])   }

    chorner (n, 0, ar, ai, rootr[i], rooti[i],tempr, tempi, t1, t2, t3, t4, t5);

    absf[i] := comabs(tempr, tempi);            { store error        }

    polydiv (n, i, rootr, rooti, x0r, x0i);     { reduce degree      }

    if real0<>0 then                    { New starting value         }
      x0i := -x0i                       { depending on real x..      }
    else
    begin
      x0r := 0.0;
      x0i := 0.0
    end
  end;

  bauhub := 0;                                         { normal exit }

fin: End; {bauhub}


  Procedure scpoly (n: Integer;      { length of vector .............}
                Var ar,              { Real part of the vector ......}
                    ai: VEC;         { Imaginary part of the vector .}
                Var scal: Double);   { Scaling factor ...............}
 {====================================================================*
 *                                                                    *
 *  scalpoly scales the polynomial P :                                *
 *                                               n-1           n      *
 *      P(x) = a[0] + a[1] * x + ... + a[n-1] * x    + a[n] * x ,     *
 *                                                                    *
 *  where all a[i], i=0, ..., n, can be complex.                      *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Eingabeparameter:                                                *
 *   ================                                                 *
 *      n        integer;                                             *
 *               degree of the polynomila (>=1)                       *
 *      ar, ai   type VEC;                                            *
 *               Real and imaginary parts of the coefficients         *
 *               a[0],..,a[n]                                         *
 *                                                                    *
 *   Output parameters:                                               *
 *   =================                                                *
 *      ar, ai   type VEC;                                            *
 *               Real and imaginary parts of the coefficients         *
 *               a[0],..,a[n] of the scaled polynomial.               *
 *      scal     REAL                                                 *
 *               Scaling factor                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *      comabs:  Modulus of a complex number                          *
 *      Power :  y^x                                                  *
 *      Max   :  Maximum of two real numbers                          *
 *                                                                    *
 *====================================================================}
Var
  p, pot: Double;
  i: Integer;
Begin

  scal := 0.0;
                              { scal =                               }
  p := comabs (ar[n], ai[n]); {              a[i]  1/(n-i)           }
  for i := 0 to n-1 do        {   max{ cabs( ---- )       , i=0..n-1 } 
                              {              a[n]                    }
    if (ar[i] <> 0.0) or (ai[i] <> 0.0) then
    begin
      ai[i] := ai[i] / p;
      ar[i] := ar[i] / p;

      pot := Power(comabs (ar[i],ai[i]), 1.0/(n-i));
      scal := Max(scal, pot)
    end;

  ar[n] := ar[n]/p;                    { Absolute value of a[n] = 1  }
  ai[n] := ai[n]/p;

  if scal = 0.0 then scal := 1.0;

  p:=1.0;
  for i := n-1 Downto 0 do
  begin
    p := p*scal;               {                    n-i              }
    ar[i]:=ar[i]/p;            { a[i] = a[i] / (scal    ), i=0..n-1  }
    ai[i]:=ai[i]/p             {                                     }
  end
End;


  Procedure chorner (n,              { highest degree in polynomial .}
                     iu: Integer;    { lowest degree in polynomial ..}
                     ar,             { Real parts of coefficients ...}
                     ai: VEC;        { Imaginary parts, coefficients }
                     xr,             { Real part of x ...............}
                     xi: Double;     { Imaginary part of x ..........}
                 Var pr,             { Real part of function value ..}
                     pi,             { Imaginary part of function v. }
                     p1r,            { Real part of first derivative }
                     p1i,            { Imaginary part, first deriv. .}
                     p2r,            { Real part, second derivative .}
                     p2i,            { Imaginary part, second deriv. }
                     rf1: Double);   { Error estimate for 1st deriv. }
 {====================================================================*
 *                                                                    *
 *  Horner scheme for polynomial with complex coefficients.           *
 *  We compute :                                                      *
 *    1. Polynomial value of the polynomial P (complex) of degree     *
 *       n - iu,                                                      *
 *    2. value of first derivative at x,                              *
 *    3. value of 2nd derivative at x,                                *
 *    4. an error estimate for the first derivative.                  *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;                                               *
 *               Maximal degree of the polynomial ( >= 1 )            *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficients of the  *
 *               polynomial with the coefficients a[iu], ..., a[n]    *
 *      x0r,x0i  REAL   x0r, x0i;                                     *
 *               Real and imaginary parts of the point of evaluation  *
 *                                                                    *
 *   Ausgabeparameter:                                                *
 *   ================                                                 *
 *      pr, pi   REAL   *pr, *pi;                                     *
 *               Real and imaginary part of the polynomial            *
 *      p1r, p1i REAL   *p1r, *p1i;                                   *
 *               Real and imaginary parts of the 1st derivative there *
 *      p2r, p2i REAL   *p2r, *p2i;                                   *
 *               Real and imaginary parts of the 2nd derivative       *
 *      rf1      REAL   *rf1;                                         *
 *               Error estimate for the first derivative              *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *     comabs():    modulus of a complex number                       *
 *     pow ():      Power function                                    *
 *====================================================================*
 *                                                                    *
 *   Constant used:  EPS                                              *
 *   =============                                                    *
 *                                                                    *
 *====================================================================}
Label continue, 10;
Var
    i,i1,j: Integer;
    temp, temp1: Double;
Begin

  p2r := ar[n];
  p2i := ai[n];

  pr := p2r;
  p1r := p2r;
  pi := p2i;
  p1i := p2i;

  rf1 := comabs(pr, pi);
  i1 := n - iu;

  j:=n-iu;
  for i := n - 1 Downto iu do
  begin
    temp := pr;                        { Polynomial value (pr,pi)    }
    pr := pr * xr - pi * xi + ar[i];
    pi := pi * xr + temp * xi + ai[i];
    if i = iu then goto continue;                      { exit i loop }

    temp := p1r;                       { 1st derivative (p1r,p1i)    }
    p1r := p1r * xr - p1i * xi;
    p1i := p1i * xr + temp * xi;

    temp := comabs(p1r, p1i);          { Error estimate for the 1st  }
    p1r := p1r + pr;                   { derivative of P             }
    p1i := p1i + pi;

    temp1 := comabs (pr, pi);
    temp := Max(temp, temp1);

    if temp > rf1 then
    begin
      rf1 := temp;
      i1 := j - 1
    end;

    if i - iu <= 1 then goto 10;                    {go to next i, j }

    temp := p2r;                         { 2nd derivative (p2r,p2i)  }
    p2r := p2r * xr - p2i * xi + p1r;
    p2i := p2i * xr + temp * xi + p1i;
10: Dec(j)
  end;

continue: temp := comabs(xr, xi);

  if temp <> 0.0 then
    rf1 := rf1 * Power(temp, 1.0*i1) * (i1 + 1)
  else
    rf1 := comabs(p1r, p1i);

  rf1 := rf1 * EPS;

  p2r := p2r + p2r;
  p2i := p2i + p2i

End; {chorner}

Procedure quadsolv       { Complex quadratic equation ................}
             (
               ar,                  {second degree coefficient .......}
               ai,
               br,                  {linear coefficient ..............}
               bi,   
               cr,                  {polynomial constant .............}
               ci: Double;
           Var tr,                  {solution ........................}
               ti: Double
             );
 {====================================================================*
 *                                                                    *
 *  Compute the least magnitude solution of the quadratic equation    *
 *  a * t**2 + b * t + c = 0. Here a, b, c and t are complex.         *
 *                                         2                          *
 *  Formeula used: t = 2c / (-b +/- sqrt (b  - 4ac)).                 *
 *  This formula is valid for a=0 .                                   *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *  Input parameters:                                                 *
 *  ================                                                  *
 *      ar, ai   coefficient of t**2             REAL   ar, ai;       *
 *      br, bi   coefficient of t                REAL   br, bi;       *
 *      cr, ci   constant term                   REAL   cr, ci;       *
 *                                                                    *
 *  Output parameter:                                                 *
 *  ================                                                  *
 *      tr, ti   complex solution of minimal magnitude                *
 *                                               REAL   *tr, *ti;     *
 *                                                                    *
 *  Macro used :     SQRT                                             *
 *  ============                                                      *
 *                                                                    *
 *====================================================================}
Var
   pr, pi, qr, qi, h: Double;
Begin
  pr := br * br - bi * bi;
  pi := 2.0 * br * bi;                       {  p = b * b             }

  qr := ar * cr - ai * ci;
  qi := ar * ci + ai * cr;                   {  q = a * c             }

  pr := pr - 4.0 * qr;       
  pi := pi - 4.0 * qi;                       { p = b * b - 4 * a * c  }

  h  := SQRT (pr * pr + pi * pi);            { q = sqrt (p)           }

  qr := h + pr;
  if qr > 0.0 then
    qr := SQRT (qr * 0.5)
  else
    qr := 0.0;

  qi := h - pr;
  if qi > 0.0 then
    qi := SQRT (qi * 0.5)
  else
    qi := 0.0;

  if pi < 0.0 then qi := -qi;

  h := qr * br + qi * bi;      { p = -b +/- q, choose sign for large }
                               { magnitude  p                        }
  if h > 0.0 then
  begin
    qr := -qr;
    qi := -qi
  end;

  pr := qr - br;
  pi := qi - bi;
  h := pr * pr + pi * pi;                      { t = (2 * c) / p     }

  if h = 0.0 then
  begin
    tr := 0.0;
    ti := 0.0
  end
  else
  begin
    tr := 2.0 * (cr * pr + ci * pi) / h;
    ti := 2.0 * (ci * pr - cr * pi) / h
  end
End; {quadsolv}

  Function bauroot (n,                   { largest degree ...............}
                    iu: Integer;         { lowest degree ................}
                    ar,                  { Real parts of the coefficients}
                    ai: VEC;             { Imaginary parts, coefficients }
                Var x0r,                 { Real part of the root ........}
                    x0i:Double):Integer; { Imaginary part of the root ...}
 {====================================================================*
 *                                                                    *
 *  bauroot computes one root of the polynomial P of degree n-iu:     *
 *                                                 n-iu               *
 *      P(x) = a[iu] + a[iu+1] * x + ... + a[n] * x                   *
 *                                                                    *
 *  with complex  coefficients a[i], i=iu, ..., n.                    *
 *  This program uses Newton's method on the function P(x) / P'(x).   *
 *  The iteration is stabilized using spiralization and extrapolation.*
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;                                               *
 *               Maximal degree of the polynomial ( >= 1 )            *
 *      iu       int iu;                                              *
 *               Index for the constant term of the polynomial,       *
 *               n-iu is the degree of the polynomial with            *
 *               coefficients a[iu], ..., a[n]                        *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficients         *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      x0r,x0i  REAL   *x0r, x0i;                                    *
 *               Real and imaginary part of the computed root         *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all is ok                                            *
 *      = 1      Division by zero                                     *
 *      = 2      Iteration number ITERMAX exceeeded                   *
 *      = 3      Improper input                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used:                                                  *
 *   ==============                                                   *
 *     chorner():   Complex Horner scheme                             *
 *     comabs():    magnitude of a complex number                     *
 *     comdiv():    Complex division                                  *
 *     quadsolv():  solves quadratic equations                        *
 *                                                                    *
 *   Constants used : TRUE, FALSE, ITERMAX,                           *
 *   ===============  QR, QI, MACH_EPS, EPS, EPSROOT, BETA            *
 *                                                                    *
 *====================================================================}
Label fin,10;
Var
    rc, result, endit, iter, i: Integer;

    xoldr, xoldi, xnewr, xnewi, h, h1, h2, h3, h4, dzmax, dzmin,
       dxr, dxi, tempr, tempi, abs_pold, abs_pnew, abs_p1new, temp,
       ss, u, v, bdze, pr, pi, p1r, p1i, p2r, p2i: Double;
Begin

  result:=2; endit:=0; iter:=0; i:=0;
  xoldr:=0.0; xoldi:=0.0; dxr:=0.0; dxi:=0.0; bdze:=0.0;

  if n < 1 then
  begin
    bauroot:=3;
    goto fin
  end;

  if n - iu = 1 then                        { Polynomial of degree 1 }
  begin
    quadsolv (0.0, 0.0, ar[n], ai[n], ar[n-1], ai[n-1], x0r, x0i);
    bauroot:=0;
    goto fin
  end;

  if n - iu = 2 then                        { Polynomial of degree 2 }
  begin
    quadsolv (ar[n],ai[n], ar[n-1],ai[n-1], ar[n-2],ai[n-2], x0r,x0i);
    bauroot:=0;
    goto fin
  end;

  xnewr := x0r;  xnewi := x0i;
  endit := 0;

  chorner (n, iu, ar, ai, xnewr, xnewi,     { Evaluate polynomial   }
           pr, pi, p1r, p1i, p2r, p2i, ss);

  Inc(iter);

  abs_pnew := comabs (pr, pi);
  if abs_pnew < EPS then
  begin
    bauroot:=0;                             { Starting value is a   }
    goto fin                                { good approximation    }
  end;

  abs_pold := abs_pnew;
  dzmin := BETA * (EPSROOT + comabs (xnewr, xnewi));

  while iter < ITERMAX do                       { Bauhuber-Iteration }
  begin
    abs_p1new := comabs (p1r, p1i);

    if abs_pnew > abs_pold then             { Spiralization step     }
    begin
      i := 0;                                { dx = dx * q            }
      Inc(iter);
      temp := dxr;

      dxr := QR * dxr - QI * dxi;
      dxi := QR * dxi + QI * temp
    end
    else
    begin
      dzmax := 1.0 + comabs(xnewr, xnewi);
      h1 := p1r * p1r - p1i * p1i - pr * p2r + pi * p2i;
      h2 := 2.0 * p1r * p1i - pr * p2i - pi * p2r;
      if (abs_p1new > 10.0 * ss) and (comabs(h1, h2) > 100.0 * ss * ss) then
      begin
        Inc(i);
        if i > 2 then i := 2;
        tempr := pr * p1r - pi * p1i;
        tempi := pr * p1i + pi * p1r;

        rc := comdiv(-tempr, -tempi, h1, h2, dxr, dxi);
        if rc <> 0 then
        begin
          bauroot := 1;
          goto fin
        end;

        if comabs (dxr, dxi) > dzmax then
        begin
          temp := dzmax / comabs (dxr, dxi);       { Newton step     }
          dxr := dxr * temp;
          dxi := dxi * temp;
          i := 0
        end;

        if (i = 2) and (comabs(dxr, dxi) < dzmin / EPSROOT) and (comabs(dxr, dxi) > 0.0) then
        begin
          i := 0;                            { Extrapolation step    }
          rc := comdiv(xnewr - xoldr, xnewi - xoldi, dxr, dxi, h3, h4);
          if rc <> 0 then
          begin
            bauroot := 1;
            goto fin
          end;

          h3 := h3 + 1.0;
          h1 := h3 * h3 - h4 * h4;
          h2 := 2.0 * h3 * h4;
          rc := comdiv (dxr, dxi, h1, h2, h3, h4);
          if rc <> 0 then
          begin
            bauroot := 1;
            goto fin
          end;

          if comabs(h3, h4) < 50.0 * dzmin then
          begin
            dxr := dxr + h3;
            dxi := dxi + h4
          end
        end;

        xoldr := xnewr;
        xoldi := xnewi;
        abs_pold := abs_pnew
      end
      else
      begin
        i := 0;                            { Close to a saddle point }
        h := dzmax / abs_pnew;
        dxr := h * pr;
        dxi := h * pi;

        xoldr := xnewr;
        xoldi := xnewi;
        abs_pold := abs_pnew;

        Repeat

          chorner (n, iu, ar, ai, xnewr+dxr, xnewi+dxi, u, v, h, h1, h2, h3, h4);
          Inc(iter);

          dxr := dxr + dxr;
          dxi := dxi + dxi;                   { dx = dx * 2.0         }

        Until ABS(comabs(u,v)/abs_pnew-1.0) >= EPSROOT
      end
    end;

    if endit<>0 then
    begin
      if comabs(dxr, dxi) < 0.1 * bdze then
      begin
        xnewr := xnewr + dxr;
        xnewi := xnewi + dxi
      end;

      result := 0;
      goto 10                                 { stop iteration   }
    end
    else
    begin
      xnewr := xoldr + dxr;
      xnewi := xoldi + dxi;
      dzmin := BETA * (EPSROOT + comabs (xnewr, xnewi));

      chorner (n, iu, ar, ai, xnewr, xnewi, pr, pi, p1r, p1i, p2r, p2i, ss);

      Inc(iter);
      abs_pnew := comabs( pr, pi);

      if abs_pnew = 0.0 then
      begin
        result := 0;
        goto 10
      end;

      if (comabs(dxr, dxi) < dzmin) or (abs_pnew < EPS) then
      begin
        endit := 1;
        bdze := comabs(dxr, dxi)
      end
    end
  end;  {End Bauhuber iteration}

10: x0r := xnewr;
    x0i := xnewi;

  bauroot:=result;
fin: End; {bauroot}


  Procedure polydiv (n,              { maximal degree ...............}
                     iu: Integer;    { minimal degree ...............}
                 Var ar,             { Real parts of coefficients ...}
                     ai: VEC;        { Imaginary parts, coefficients }
                     x0r,            { Real part of x ...............}
                     x0i: Double);   { Imaginary part of x ..........}
 {====================================================================*
 *                                                                    *
 *  polydiv computes the coefficients of the polynomial Q, with       *
 *  P(x) = Q(x) * (x - x0), where x0 is a computed root of P.         *
 *  Both P and Q and x0 may be complex.                               *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;                                               *
 *               Highest degree of the polynomial ( >= 1 )            *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficienten of P   *
 *               of degree n-iu, with a[iu], ..., a[n]                *
 *      x0r,x0i  REAL   x0r, x0i;                                     *
 *               Real and imaginary parts of the root x0              *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      ar, ai   REAL   ar[], ai[];                                   *
 *               Real and imaginary parts of the coefficients         *
 *               ar[iu+1],..,ar[n] of the remainder polynomial Q      *
 *                                                                    *
 *====================================================================}
 Var
    i: Integer;
    temp: Double;
 Begin
   for i := n - 1 Downto iu+1 do
   begin
     temp := ar[i+1];
     ar[i] := ar[i] + temp * x0r - ai[i+1] * x0i;
     ai[i] := ai[i] + ai[i+1] * x0r + temp * x0i
   end
 End;

END.

{ ------------------------- END fbauhube.pas  ----------