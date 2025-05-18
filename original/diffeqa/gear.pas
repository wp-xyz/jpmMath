{***********************************************************************
*                                                                      *
* Solve a first order system of DEs using the implicit Gear method     *
* of order 4.                                                          *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Klaus Niederdrenk, 1.22.1996 (FORTRAN 77)      *
* Adaptation:           Juergen Dietel, Computing Center, RWTH Aachen  *
* Source:               FORTRAN 77 source code                         *
* Date:                 2.26.1996                                      *
*                                                                      *
*                       TPW Release By J-P Moreau, Paris.              *
*                       (www.jpmoreau.fr)                              *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C,  By Gisela Engeln-Muellges       *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
***********************************************************************}
UNIT GEAR;

INTERFACE

Uses WinCrt1, FGAUSS, UAWP, T_DGLS;

  Procedure gear4       {Gear method of 4th order for DESs of 1st order}
    (
  var x:      double;           {starting or end point ................}
      xend:   double;           {desired end point (> x) ..............}
      n:      integer;          {number of DEs ........................}
      y:      pVEC;             {initial value or solution ............}
      epsabs,                   {absolute error bound .................}
      epsrel: double;           {relative error bound .................}
  var h:      double;           {starting or final step size ..........}
      fmax:   integer;          {maximal number of calls of dgl() .....}
  var aufrufe,                  {actual number of calls of dgl() ......}
      fehler: integer           {error code ...........................}
    );

IMPLEMENTATION

{print a REAL square matrix with name (debug only) }
Procedure PrintMat(Name:string; n:integer; mat:pMAT);
var
  i, j: integer;
begin
  writeln(Name);
  for i := 0 to n-1 do
  begin
    for j := 0 to n-1 do
      write(' ',mat^[i,j]);
    writeln
  end
end;

Procedure gear4        {Gear  method of 4th order for DESs of 1st order}
    (
  var x:      double;           {starting or end point ................}
      xend:   double;           {desired end point (> x) ..............}
      n:      integer;          {number of DEs ........................}
      y:      pVEC;             {initial value or solution ............}
      epsabs,                   {absolute error bound .................}
      epsrel: double;           {relative error bound .................}
  var h:      double;           {starting or final step size ..........}
      fmax:   integer;          {maximal number of calls of dgl() .....}
  var aufrufe,                  {actual number of calls of dgl() ......}
      fehler: integer           {error code ...........................}
    );
{***********************************************************************
* Compute the value of the solution at xend, starting with the IVP.    *
* We use the implicit method of Gear of 4th order with step size       *
* control which is especially suited for stiff DEs.                    *
* The local step size control insures that the two error bounds are met*
* The number of function calls of the right hand side is limited by    *
* fmax. This function can be used inside a loop to find solutions at   *
* a specified point to arbitrary accuracy.                             *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        initial value x0                                            *
* xend     final value for the integration (xend > x0)                 *
* n        number of DEs                                               *
* dgl      Function to compute the right hand side f(x0,y0) for the    *
*          system of DEs                                               *
* y        [0..n-1] solution vector y0 of the system of DEs at x0      *
* epsabs   absolute error bound (>= 0); if zero, we only check for the *
*          relative error.                                             *
* epsrel   relative error bound (>= 0); if zero, we only check for the *
*          absolute error.                                             *
* h        given step size for first step                              *
* fmax     maximal number of calls of the right hand side of the system*
*                                                                      *
* Output parameters:                                                   *
* =================                                                    *
* x        final x-value of the integration; normally equal to xend    *
* h        final step size; keep for further calls                     *
* y        [0..n-1]-vector, the solution of the system of DEs at x     *
* aufrufe  counter of calls of dgl()                                   *
*                                                                      *
* Return values (fehler):                                              *
* ======================                                               *
* Error code:                                                          *
*   =   0:  all ok                                                     *
*   =   1:  Both error bounds too small                                *
*   =   2:  xend <= x0                                                 *
*   =   3:  Step size h <= 0                                           *
*   =   4:  n <= 0                                                     *
*   =   5:  # calls > fmax; we have not reached the desired xend;      *
*           repeat function call with actual values of x, y, h.        *
*   =   6:  Jacobi matrix is singular; x, y, h are the last values     *
*           that could be computed with accuracy                       *
*   =   7:  ran out of memory (not used here)                          *
*   =   8:  error when calling gauss() for the second time;            *
*           should not occur.                                          *
*                                                                      *
***********************************************************************}
Label fin;   {to exit procedure}
Type
  pMAT4 = ^MAT4;
  MAT4 = Array[0..4,0..NMAX] of double;
Var
  eps1,                   {MACH_EPS^0.75; used instead of MACH_EPS
                           to check whether we have reached xend in
                           order to avoid minuscule step sizes}
  eps2,                   {100 * MACH_EPS; for comparison with zero}
  hs: double;             {optimal step size for Jacobi matrix
                           approximation      }

  hilf    : pVEC;              {pointer to a [0..NMAX]-vector       }
  zj, zjp1: pMAT4;
  f,                           {pointers to [0..NMAX]-vectors       }
  ykp1:  pVEC;
  fs,                        {pointers to [0..NMAX,0..NMAX]-matrices}
  fsg:   pMAT;
  con:   pVEC;                 { [0..NMAX]-vector                   }
  perm:  pIVEC;                { pointer to [0..NMAX] permutation
                                 vector for Gauss elimination       }

  sg,                          {sign of xend                        }
  xe:    double;               {|xend| - eps2, carrying the sign
                                of xend                             }
  amEnde:integer;              {Flag, indicating that we shall reach
                                xend with the current step          }
  ymax,                        {Maximum norm of the current
                                approximate value of y              }
  dummy,                       {aux storage for h                   }
  xka,
  xke, hka, hk1, diff, eps, q, halt, quot1, quot2, quot3, quot4: double;

  done, nochmal: Boolean;
  aufrufe_awp,
  signdet,                     {sign of determinant in Gauss        }
  iter, i, k: integer;

  dum, dum1: pVEC;             {dummy vectors for gauss             }

Begin

  { ------------------------- Initialize -------------------------- }
  dummy :=ZERO; done:=False;
  eps1  := POW(MACH_EPS, 0.75);
  eps2  := 100 * MACH_EPS;
  hs    := 10.0 * SQRT(MACH_EPS);

  { ----------------- Allocate matrices & vectors ---------------- }
  New(hilf); New(zj); New(zjp1); New(f); New(ykp1); New(fs);
  New(fsg); New(con); New(perm); New(dum); New(dum1);

  if xend >= ZERO then sg := ONE
                  else sg := -ONE;
  xe       := (ONE - sg * eps2) * xend;
  fehler   := 0;
  aufrufe  := 1;
  amEnde   := 0;

  { --------- check input parameters -------- }
  ymax := norm_max(y, n);
  if (epsabs <= eps2 * ymax) and (epsrel <= eps2) then
    fehler := 1
  else if xe < x then
    fehler := 2
  else if h < eps2 * ABS(x) then
    fehler := 3
  else if n <= 0 then
    fehler := 4;
  if fehler>0 then goto fin;

  { ------------ first integration step ----------- }
  if x + h > xe then
  begin
    h      := xend - x;
    dummy  := h;
    amEnde := 1
  end;
  for i:=0 to n-1 do hilf^[i]:=y^[i];
  xka := x;
  xke := xka;
  hka := 0.25 * h;
  hk1 := hka;
  for k := 1 to 4 do
  begin
    xke := xke + hka;

    awp(xka, xke, n, hilf, epsabs, epsrel, hk1,6, fmax - aufrufe, aufrufe_awp, fehler);

    aufrufe := aufrufe + aufrufe_awp;

    if fehler<>0 then goto fin;
    for i:=0 to n-1 do zjp1^[k,i] := hilf^[i];
  end;

  dgl(x,y,f);
  Inc(aufrufe);

  { -------- Compute first Gear-Nordsieck approximation ------ }
  for i := 0 to n-1 do
  begin
    zj^[0,i] := y^[i];
    zj^[1,i] := h * f^[i];
    zj^[2,i] := ONE / 24.0 * (35.0 * y^[i] - 104.0 * zjp1^[1,i] +
                114.0 * zjp1^[2,i] - 56.0 * zjp1^[3,i] + 11.0 * zjp1^[4,i]);
    zj^[3,i] := ONE / 12.0 * (-5.0 * y^[i] + 18.0  * zjp1^[1,i] -
                24.0 * zjp1^[2,i] + 14.0 * zjp1^[3,i] - 3.0 * zjp1^[4,i]);
    zj^[4,i] := ONE / 24.0 * (y^[i] - 4.0 * zjp1^[1,i] + 6.0 * zjp1^[2,i] -
                4.0 * zjp1^[3,i] + zjp1^[4,i])
  end;

  { -------------- adjust step size ------------- }
  Repeat
  
    { --- use Newton method for an implicit approximation --- }

    for i:=0 to n-1 do
      ykp1^[i] := zj^[0,i] + zj^[1,i] + zj^[2,i] + zj^[3,i] + zj^[4,i];

    dgl(x+h,ykp1,f);

    for k:=0 to n-1 do
    begin
      
      {copy_vector(hilf, ykp1, n); }
      for i:=0 to n-1 do hilf^[i]:=ykp1^[i];
      hilf^[k] := hilf^[k] - hs;
      dgl(x+h,hilf,dum);
      for i:=0 to n-1 do fs^[k,i]:=dum^[i];
      for i:=0 to n-1 do
        fs^[k,i] := -h * 0.48 * (f^[i] - fs^[k,i]) / hs;
      fs^[k,k] := fs^[k,k] + ONE
    end;

    aufrufe := aufrufe + n + 1;
    for i:=0 to n-1 do
    begin
      con^[i] := ykp1^[i] - 0.48 * (zj^[1,i] + 2.0 * zj^[2,i] +
                3.0 * zj^[3,i] + 4.0 * zj^[4,i]);
      for k:=0 to n-1 do
        fsg^[k,i] := fs^[i,k]
    end;

    gauss(1, n, fsg, fsg, perm, dum, dum1, signdet,fehler);

    if fehler>0 then  {error in gauss ?}
    begin
      fehler := 6;
      goto fin
    end;

    for iter:=1 to 3 do
    begin
      for i:=0 to n-1 do
      begin
        hilf^[i] := - ykp1^[i];
        for k:=0 to n-1 do
          hilf^[i] := hilf^[i] + fs^[k,i] * ykp1^[k];
        hilf^[i] := h * 0.48 * f^[i] + hilf^[i] + con^[i]
      end;

      gauss(2, n, fsg, fsg, perm, hilf, ykp1, signdet,fehler);

      if fehler > 0 then
      begin
        fehler:=8;
        goto fin
      end;

      dgl(x+h,ykp1,f);

    end;
    aufrufe := aufrufe + 3;

    { ---- Compute corresponding Gear-Nordsieck approximation ---- }

    for i:=0 to n-1 do
      hilf^[i] := h * f^[i] - zj^[1,i] - 2.0 * zj^[2,i] -
                 3.0 * zj^[3,i] - 4.0 * zj^[4,i];

    for i:=0 to n-1 do
    begin
      zjp1^[0,i] := ykp1^[i];
      zjp1^[1,i] := h * f^[i];
      zjp1^[2,i] := zj^[2,i] + 3.0 * zj^[3,i] + 6.0 * zj^[4,i] + 0.7 * hilf^[i];
      zjp1^[3,i] := zj^[3,i] + 4.0 * zj^[4,i] + 0.2 * hilf^[i];
      zjp1^[4,i] := zj^[4,i] + 0.02 * hilf^[i]
    end;

    { --- decide whether to accept last step ---

    copy_vector(hilf, zjp1^[4], n);
    copy_vector(con, zj^[4], n);           }
    for i:=0 to n-1 do
    begin
      hilf^[i]:=zjp1^[4,i];
      con^[i]:=zj^[4,i]
    end;

    diff := dist_max(hilf, con, n);
    ymax := norm_max(ykp1, n);
    eps  := (epsabs + epsrel * ymax) / 6.0;
    q    := SQRT(SQRT(eps / diff)) / 1.2;

    if diff < eps then
    begin

      { --- accept last step; prepare for next one --- }

      x := x + h;
      {copy_vector(y, ykp1, n);}
      for i:=0 to n-1 do y^[i]:=ykp1^[i];

      { stop integration, if interval end xend has been reached
        or if there has been too many function calls  }

      nochmal := FALSE;
      repeat
      
        if amEnde<>0 then
        begin
          h := dummy;
          goto fin
        end
        else if aufrufe > fmax then
        begin
          fehler := 5;
          goto fin
        end;

        { --- adjust step size for next step --- }
        halt := h;
        h    := Min(q, 2.0) * h;

        if x + h >= xe then
        begin
          dummy  := h;
          h      := xend - x;
          amEnde := 1;

          { --- close enough to xend  :=> stop integration --- }
          if h < eps1 * ABS(xend) then nochmal := TRUE
        end
      until not nochmal;

      { ------ compute Gear-Nordsieck approximation ---
        ------ for the next step                    --- }
      quot1 := h / halt;
      quot2 := sqr(quot1);
      quot3 := quot2 * quot1;
      quot4 := quot3 * quot1;
      for i:=0 to n-1 do
      begin
        zj^[0,i] := zjp1^[0,i];
        zj^[1,i] := quot1 * zjp1^[1,i];
        zj^[2,i] := quot2 * zjp1^[2,i];
        zj^[3,i] := quot3 * zjp1^[3,i];
        zj^[4,i] := quot4 * zjp1^[4,i]
      end
    end
    else
    begin

      { ------ repeat last step with smaller step size;  ---
        ------ adjust Gear-Nordsieck approximation --------- }
      halt  := h;
      h     := Max(0.5, q) * h;
      quot1 := h / halt;
      quot2 := sqr(quot1);
      quot3 := quot2 * quot1;
      quot4 := quot3 * quot1;
      for i:=0 to n-1 do
      begin
        zj^[1,i] := quot1 * zj^[1,i];
        zj^[2,i] := quot2 * zj^[2,i];
        zj^[3,i] := quot3 * zj^[3,i];
        zj^[4,i] := quot4 * zj^[4,i];
      end;
      amEnde := 0
    end
  Until done;

{free memory}
fin: Dispose(hilf); Dispose(zj); Dispose(zjp1); Dispose(f);
     Dispose(ykp1); Dispose(fs); Dispose(fsg); Dispose(con);
     Dispose(perm); Dispose(dum); Dispose(dum1)

End;

(*
char *gear_fehlertext          /* find error message and error class  */
/*.IX{gear\unt fehlertext}*/
    (
     int       fehlercode,     /* error number from gear()       .....*/
     fehler_t  *fehlerart      /* severity of error from gear()  .....*/
    )                          /* error text .........................*/

/***********************************************************************
* Use error codes from gear() to find the appropriate error message    *
* and classify the type of error as: warning, fatal, no or unknown     *
* error.                                                               *
*                                                                      *
* global variable:                                                     *
* :=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=                                                     *
* fehler_t                                                             *
***********************************************************************/

{
  static char *fehlertext[] :=
    {
      "gear(): no error",
      "gear(): at least one error bound too small",
      "gear(): xend <:= x0",
      "gear(): h <:= 0",
      "gear(): n <:= 0",
      "gear(): too many calls of dgl()",
      "gear(): Jacobi matrix singular",
      "gear(): out of memory",
      "gear(): impossible: 2nd call of gauss() with error"
    };

  switch (fehlercode)
  {
    case 0:  *fehlerart := KEIN_FEHLER;
             break;
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:  *fehlerart := FATAL;
             break;
    default: *fehlerart := UNBEKANNT;
             break;
  }

  return fehlertext[fehlercode];
}  *)

END.
{ -------------------------- END  gear.pas 