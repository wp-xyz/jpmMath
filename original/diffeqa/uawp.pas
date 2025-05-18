{***********************************************************************
*                                                                      *
* Solve an ordinary system of first order differential equations using *
* -------------------------------------------------------------------- *
* automatic step size control                                          *
* ----------------------------                                         *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Klaus Niederdrenk (FORTRAN)                    *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 6.2.1992, 10.2.1995                            *
*                                                                      *
*                       TPW Release By J-P Moreau, Paris.              *
*                       (www.jpmoreau.fr)                              *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C,  By Gisela Engeln-Muellges       *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
***********************************************************************}
UNIT UAWP;

INTERFACE

Uses Wincrt1, FGauss, T_Dgls;

     Procedure awp        {1st order DESs with autom. step size control}
        (
      var x:    double;           {initial/final x value ..............}
          xend: double;           {desired end point ..................}
          n:    integer;          {number of DEs ......................}
          y:    pVEC;             {initial/final y value ..............}
          epsabs,                 {absolute error bound ...............}
          epsrel: double;         {relative error bound ...............}
      var h:    double;           {initial/final step size ............}
          methode,                {desired method (3, 6, 7) ...........}
          fmax: integer;          {maximal # of calls of  dgl() .......}
      var aufrufe,                {actual # of calls of  dgl() ........}
          fehler: integer         {error code .........................}
        );

     Procedure PrintVec(name:String; n:integer; V:pVEC);

     Function dist_max    {Maximum norm of a difference vector ........}
        (
          vector1,
          vector2: pVEC;
          n: integer
        ): Double;

     Function POW(y,x:double): double;

     Function norm_max    { Find the maximum norm of a REAL vector .........}
        (
          vektor:pVEC;                 { pointer to vector ...... }
          n:INTEGER                    { length of vector ....... }
        ) : Double;                    { Maximum norm ........... }

     Function Max(a,b:double): double;
     Function Min(a,b:double): double;


IMPLEMENTATION

Var

{***********************************************************************
*          Global variables for implementation section.                *
*                                                                      *
***********************************************************************}
   yhilf,             {pointers to {[0..n-1] aux. vectors for the   }
   k1,                {embedding formulas in ruku23() and engl45(). }
   k2,                 
   k3,
   k4,
   k5,
   k6: pVEC;

   k7,                  {more pointers to [0..n-1] aux. vectors for }
   g6,                  {embedding formula  in prdo45               }
   g7: pVEC;            {these pointers are initialized by awp.     }      

   steif1,            {Flag, that is set in prdo45() if its         }
                      {stiffness test (dominant eigenvalue)         }
                      {indicates so. Otherwise no changes.          }
   steifanz,          {counter for number of successive successes   }
                      {of stiffness test of Shampine and Hiebert in }
                      {prdo45().                                    }
   steif2: integer;   {Flag, set in prdo45(), when the stiffness    }
                      {test of  Shampine and Hiebert wa successful  }
                      {three times in a row; otherwise no changes   }

Procedure PrintVec(name:String; n:integer; V:pVEC);
Var i:integer;
Begin
  writeln(name,':');
  for i:=0 to n-1 do write(V^[i]);
  writeln
End;


Function dist_max         {Maximum norm of a difference vector ........}
        (
          vector1,
          vector2: pVEC;
          n: integer
        ): Double;
{***********************************************************************
* Compute the maximum norm of the difference of two [0..n-1] vectors   *
*                                                                      *
* global name used:                                                    *
* ================                                                     *
*   ZERO                                                               *
***********************************************************************}
Var
     abstand,              {reference value for computation of distance}
     hilf: Double;         {distance of two vector elements            }
     i: integer;
Begin
  abstand:=ZERO;
  for i:=n-1 Downto 0 do
  begin
    hilf:=ABS(vector1^[i]-vector2^[i]);
    if hilf > abstand then abstand := hilf
  end;
  dist_max := abstand
End;


Procedure ruku23          {Runge-Kutta embedding f. of 2nd, 3rd degree}
        (
          x:  Double;
          y:  pVEC;
          n:  Integer;
          h:  Double;
      Var y2: pVEC;
      Var y3: pVEC
        );
{***********************************************************************
* Compute 2nd and 3rd order approximates y2, y3 at x + h starting with *
* a solution y at x by using Runge-Kutta embedding formulas on the     *
* first order system of n differential equations   y' = f(x,y) , as    *
* supplied by  dgl().                                                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x    x-value of left end point                                       *
* y    y-values at x                                                   *
* n    number of differential equations                                *
* h    step size                                                       *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* y2   2nd order approximation for y at x + h                          *
* y3   3rd order approximation for y at x + h                          *
*                                                                      *
* External function:                                                   *
* =================                                                    *
* dgl  function that evaluates the right hand side of the system       *
*      y' = f(x,y)                                                     *
***********************************************************************}
Var
  i: integer;     {loop variable}
Begin
  dgl(x,y,k1);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + h * k1^[i];
  dgl(x+h,yhilf,k2);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + 0.25 * h * (k1^[i] + k2^[i]);
  dgl(x+0.5*h,yhilf,k3);
  for i:=0 to n-1 do
  begin
    y2^[i] := y^[i] + 0.5 * h * (k1^[i] + k2^[i]);
    y3^[i] := y^[i] + h / 6.0 * (k1^[i] + k2^[i] + 4.0 * k3^[i])
  end
End;


Procedure engl45           {Einbettungsforml von England 4. und 5. Ord.}
    (
      x:    double;                {starting point of integration .....}
      y:    pVEC;                  {initial value at x ................}
      n:    Integer;               {number of differential equations ..}
      h:    double;                {step size .........................}
  Var y4,                       {4th order approximation for y at x + h}
      y5:   pVEC                {5th order approximation for y at x + h}
    );
{***********************************************************************
* Compute 4th and 5th order approximates y4, y5 at x + h starting with *
* a solution y at x by using the England embedding formulas on the     *
* first order system of n differential equations   y' = f(x,y) , as    *
* supplied by  dgl().                                                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x    initial x-value                                                 *
* y    y-values at x, type pVEC                                        *
* n    number of differential equations                                *
* h    step size                                                       *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* y4   4th order approximation for y at x + h (pVEC)                   *
* y5   5th order approximation for y at x + h (pVEC)                   *
*                                                                      *
* External function:                                                   *
* =================                                                    *
* dgl  function that evaluates the right hand side of the system       *
*      y' = f(x,y)                                                     *
***********************************************************************}
Var
  i: integer;      {loop variable}
Begin
  dgl(x,y,k1);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + 0.5 * h * k1^[i];
  dgl(x+0.5*h,yhilf,k2);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + (0.25 * h * (k1^[i] + k2^[i]));
  dgl(x+0.5*h,yhilf,k3);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + h * (-k2^[i] + 2.0 * k3^[i]);
  dgl(x+h,yhilf,k4);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + h / 27.0 * (7.0 * k1^[i] + 10.0 * k2^[i] + k4^[i]);
  dgl(x+2.0/3.0*h,yhilf,k5);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + h / 625.0 * (28.0 * k1^[i] - 125.0 * k2^[i] +
                546.0 * k3^[i] + 54.0 * k4^[i] - 378.0 * k5^[i]);
  dgl(x+h/5.0,yhilf,k6);
  for i:=0 to n-1 do
  begin
    y4^[i] := y^[i] + h / 6.0 * (k1^[i] + 4.0 * k3^[i] + k4^[i]);
    y5^[i] := y^[i] + h / 336.0 * (14.0 * k1^[i] + 35.0 * k4^[i] +
             162.0 * k5^[i] + 125.0 * k6^[i])
  end
End;


Procedure prdo45  {embedding formulas of Prince-Dormand of 4./5. order}
    (
      x:    double;                 {starting point of integration ....}
      y:    pVEC;                   {initial value at x ...............}
      n:    integer;                {number of DEs ....................}
      h:    double;                 {step size ........................}
  Var y4,                           {solution of 4th order at x+h .....}
      y5:   pVEC                    {solution of 5th order at x+h .....}
    );
{***********************************************************************
* Compute 4th and 5th order approximates y4, y5 at x + h starting with *
* a solution y at x by using the Prince-Dormand embedding formulas on  *
* the first order system of n differential equations y' = f(x,y) , as  *
* supplied by  dgl().                                                  *
* Simultaneously we perform two tests for stiffness whose results are  *
* stored in steif1 and steif2.                                         *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x    initial x-value                                                 *
* y    y-values at x (pVEC)                                            *
* n    number of differential equations                                *
* h    step size                                                       *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* y4   4th order approximation for y at x + h                          *
* y5   5th order approximation for y at x + h                          *
*                                                                      *
* External function:                                                   *
* =================                                                    *
* dgl  function that evaluates the right hand side of the system       *
*      y' = f(x,y)                                                     *
***********************************************************************}
Var
  i,               {loop variable}
  steifa: integer; {Flag which is set if the second test for stiffness
                    (Shampine und Hiebert) is positive; otherwise the
                    flag is erased.   }
Begin

  dgl(x,y,k1);                         {coefficients}
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + 0.2 * h * k1^[i];
  dgl(x+0.2*h,yhilf,k2);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + 0.075 * h * (k1^[i] + 3.0 * k2^[i]);
  dgl(x+0.3*h,yhilf,k3);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + h / 45.0 * (44.0 * k1^[i] - 168.0 * k2^[i] + 160.0 * k3^[i]);
  dgl(x+0.8*h,yhilf,k4);
  for i:=0 to n-1 do
    yhilf^[i] := y^[i] + h / 6561.0 * (19372.0 * k1^[i] - 76080.0 * k2^[i]
                + 64448.0 * k3^[i] - 1908.0 * k4^[i]);
  dgl(x+(8.0/9.0)*h,yhilf,k5);
  for i:=0 to n-1 do
    g6^[i] := y^[i] + h / 167904.0 * (477901.0 * k1^[i] - 1806240.0 * k2^[i]
              + 1495424.0 * k3^[i] + 46746.0 * k4^[i] - 45927.0 * k5^[i]);
  dgl(x+h,g6,k6);
  for i:=0 to n-1 do
    g7^[i] := y^[i] + h / 142464.0 * (12985.0 * k1^[i] + 64000.0 * k3^[i]
              + 92750.0 * k4^[i] - 45927.0 * k5^[i] + 18656.0 * k6^[i]);
  dgl(x+h,g7,k7);
  for i:=0 to n-1 do
  begin
    y5^[i] := g7^[i];
    y4^[i] := y^[i] + h / 21369600.0 * (1921409.0 * k1^[i] + 9690880.0 * k3^[i]
              + 13122270.0 * k4^[i]- 5802111.0 * k5^[i] + 1902912.0 * k6^[i]
              + 534240.0 * k7^[i])
  end;

  {Test for stiffness via dominant eigenvalue}

  if dist_max(k7, k6, n) > 3.3 * dist_max(g7, g6, n) then  steif1 := 1;

  {one step in steffness test of Shampine & Hiebert}

  for i:=0 to n-1 do
  begin
    g6^[i] := h * (2.2 * k2^[i] + 0.13 * k4^[i] + 0.144 * k5^[i]);
    g7^[i] := h * (2.134 * k1^[i] + 0.24 * k3^[i] + 0.1 * k6^[i])
  end;

  if dist_max(g6, g7, n) < dist_max(y4, y5, n) then steifa := 1
  else  steifa := 0;

  if steifa > 0 then
  begin
    Inc(steifanz);
    if steifanz >= 3 then steif2 := 1
  end
  else
    steifanz := 0
End;

  Function POW(y,x:double): double;
  Begin
    IF x<0 THEN EXIT;
    POW:=Exp(x*Ln(y));
  End;

  Function norm_max    { Find the maximum norm of a REAL vector .........}
             (
              vektor:pVEC;                 { pointer to vector ...... }
              n:INTEGER                    { length of vector ....... }
             ) : Double;                   { Maximum norm ........... }
  {***********************************************************************
  *          Return the maximum norm of a [0..n-1] vector v              *
  *                                                                      *
  ***********************************************************************}
  VAR
    norm : REAL;                                             { local max }
    betrag : REAL;                            { magnitude of a component }
    i : INTEGER;
  Begin
    norm:=0.0;
    for i:=0 to n-1 do
    begin
      betrag:=abs(vektor^[i]);
      if betrag > norm then norm := betrag
    end;   
    norm_max := norm
  End;

  Function Min(a,b:double): double;
  begin
    if a<b then Min:=a else Min:=b
  end;

  Function Max(a,b:double): double;
  begin
    if a>b then Max:=a else Max:=b
  end;

Procedure awp            {1st order DESs with autom. step size control}
        (
      var x:    double;           {initial/final x value ..............}
          xend: double;           {desired end point ..................}
          n:    integer;          {number of DEs ......................}
          y:    pVEC;             {initial/final y value ..............}
          epsabs,                 {absolute error bound ...............}
          epsrel: double;         {relative error bound ...............}
      var h:    double;           {initial/final step size ............}
          methode,                {desired method (3, 6, 7) ...........}
          fmax: integer;          {maximal # of calls of  dgl() .......}
      var aufrufe,                {actual # of calls of  dgl() ........}
          fehler: integer         {error code .........................}
        );
{***********************************************************************
* Compute the solution y of a system of first order ordinary           *
* differential equations       y' = f(x,y)   at xend from the given    *
* initial data (x0, y0).                                               *
* We use automatic step size control internally so that the error of   *
* y (absolutely or relatively) lies within the given error bounds      *
* epsabs and epsrel.                                                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        initial value for x                                         *
* y        initial values for y (type pVEC)                            *
* n        number of differential equations                            *
* xend     end of integration; xend > x0                               *
* h        initial step size                                           *
* epsabs   absolute error bound; >= 0; if = 0 we only check the        *
*          relative error.                                             *
* epsrel   relative error bound; >= 0; if = 0 we check only the        *
*          absolute eror.                                              *
* fmax     max number of evaluations of right hand side in dgl()       *
* methode  chooses the method                                          *
*          = 3: Runge-Kutta method of 2nd/3rd order                    *
*          = 6: England formula of 4th/5th order                       *
*          = 7: Formula of Prince-Dormand of 4th/5th order             *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x        final x-value for iteration. If fehler = 0  we usually have *
*            x = xend.                                                 *
* y        final y-values for the solution at x                        *
* h        final step size used; leave for subsequent calls            *
* aufrufe  actual number of calls of dgl()                             *
*                                                                      *
* Return value (fehler):                                               *
* =====================                                                *
* = 0: all ok                                                          *
* = 1: both error bounds chosen too small for the given mach. constant *
* = 2: xend <= x0                                                      *
* = 3: h <= 0                                                          *
* = 4: n <= 0                                                          *
* = 5: more right hand side calls than allowed: aufrufe > fmax,        *
*      x and h contain the current values when stop occured.           *
* = 6: improper input for embedding formula                            *
* = 7: lack of available memory (not used here)                        *
* = 8: Computations completed, but the Prince Dormand formula stiff-   *
*      ness test indicates possible stiffness.                         *
* = 9: Computations completed, but both Prince Dormand formula stiff-  *
*      ness tests indicate possible stiffness. Use method for stiff    *
*      systems instead !                                               *
* =10: aufrufe > fmax, see error code 5; AND the Prince Dormand formula*
*      indicates stiffness; retry using a stiff DE solver !            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* ruku23, engl45, dist_max, yhilf, k1, k2, k3, k4, k5, k6, MACH_EPS,   *
* Min, Max, norm_max, POW, ZERO, ONE, TWO, NULL, prdo45, steif1,       *
* steif2, steifanz.                                                    *
*                                                                      *
* External function:                                                   *
* =================                                                    *
* dgl  function that evaluates the right hand side of the system       *
*      y' = f(x,y)       (see t_dgls.pas).                             *
***********************************************************************}
Label  fin;    {to exit procedure}
Const
       MACH_2 = 100.0 * MACH_EPS;          {machine constant related
                                            value used for break-off
                                            criteria as zero }
Var
    xend_h,          {|xend| - MACH_2, carrying same sign as xend}
    ymax,            {Maximum norm of newest approximation of max
                      order  }
    hhilf,           {aux storage for the latest value of h
                      produced by step size control. It is saved
                      here in order to avoid to return a `h' that
                      resulted from an arbitrary reduction at the
                      end of the interval.  }
    diff,            {distance of the two approximations from the
                      embedding formula  }
    s: double;       {indicates acceptance level for results from
                      embeding formula   }
    y_bad,           {approximate solution of low order          }
    y_good: pVEC;    {ditto of high order                        }
    mach_1: double;  {machine constant dependent variable which
                      avoids using too little steps near xend    }
    i:      integer; {Loop variable}
    amEnde,          {flag that shows if the end of the interval
                     can be reached with the actual step size    }
    fertig: Boolean; {flag indicating end of iterations          }

Begin

  New(yhilf); New(k1); New(k2); New(k3); New(k4); New(k5); New(k6);
  New(k7); New(g6); New(g7); New(y_bad); New(y_good);

  fehler   := 0;                      {initialize some variables}
  mach_1   := POW(MACH_EPS, 0.75);
  amEnde   := FALSE;
  fertig   := FALSE;
  steif1   := 0;
  steif2   := 0;
  steifanz := 0;
  aufrufe  := 1;
  ymax     := norm_max(y, n);

  if xend >= ZERO then xend_h := xend * (1.0 - MACH_2)
                  else xend_h := xend * (1.0 + MACH_2);

{ ----------------------- check inputs ---------------------- }
  if (epsabs <= MACH_2 * ymax) and (epsrel <= MACH_2) then
  begin
    fehler:=1;
    goto fin
  end;
  if xend_h < x then
  begin
    fehler:= 2;
    goto fin
  end;
  if h < MACH_2 * ABS(x) then
  begin
    fehler:=3;
    goto fin
  end;
  if n <= 0 then
  begin
    fehler:=4;
    goto fin
  end;
  if (methode <> 3) and (methode <> 6) and (methode <> 7) then
  begin
    fehler:=6;
    goto fin
  end;

  {*********************************************************************
  *                                                                    *
  *                       I t e r a t i o n s                          *
  *                                                                    *
  *********************************************************************}
  if x + h > xend_h then            
  begin                             {almost at end point ?             }
    hhilf  := h;                    {A shortened step might be         }
    h      := xend - x;             {enough.                           }
    amEnde := TRUE
  end;

  Repeat                           {solve DE system by integrating from
                                    x0 to xend by suitable steps       }
    Case methode of
      3: ruku23(x, y, n, h, y_bad, y_good);
      6: engl45(x, y, n, h, y_bad, y_good);
      7: prdo45(x, y, n, h, y_bad, y_good)
    End;

    aufrufe := aufrufe + methode;

    diff := dist_max(y_bad, y_good, n);

    if diff < MACH_2 then         {compute s}
      s := 2.0
    else
    begin
      ymax := norm_max(y_good, n);
      s    := SQRT(h * (epsabs + epsrel * ymax) / diff);
      if methode <> 3 then  s := SQRT(s)
    end;

    if s > 1.0 then               {integration acceptable? }
    begin
      for i:=0 to n-1 do          {accept highest order solution}
        y^[i] := y_good^[i];      {move x                       }

      x := x + h;

      if (amEnde) then            {at end of interval? }
      begin
        fertig := TRUE;           {stop iteration }
        if methode = 7 then
        begin
          if (steif1>0) or (steif2>0)  then fehler := 8;
          if (steif1>0) and (steif2>0) then fehler := 9
        end
      end
      else if aufrufe > fmax then {too many calls of   dgl()? }
      begin
        hhilf  := h;                {save actual step size}
        fehler := 5;                {report error and stop}
        fertig := TRUE;
        if (methode = 7) and ((steif1>0) or (steif2>0)) then fehler := 10
      end
      else                          {Integration was successful    }
      begin                         {not at the interval end?      }
        h := h * min(2.0, 0.98*s);  {increase step size for next
                                    {step properly, at most by
                                    {factor two. Value `0.98*s' is
                                    {used in order to avoid that
                                    {the theoretical value s is
                                    {exceeded by accidental
                                    {rounding errors.              }
        if x + h > xend_h then      {nearly reached xend?          }
        begin
          hhilf  := h;              {=> One further step with      }
          h      := xend - x;       {reduced step size might be    }
          amEnde := TRUE;           {enough.                       }

          if h < mach_1 * ABS(xend) then    {very close to xend ?  }
            fertig := TRUE                  {finish iteration      }
        end
      end
    end
    else                            {step unsuccessful?            }
    begin                           {before repeating this step:   }
      h := h * Max(0.5, 0.98*s);    {reduce step size properly, at }
                                    {most by factor 1/2 (for factor}
      amEnde := FALSE               {0.98: see above).             }
    end

  Until fertig;
 
  h := hhilf;        {return the latest step size computed by step
                      size control and  error code to the caller }

     {free memory}
fin: Dispose(yhilf); Dispose(k1); Dispose(k2); Dispose(k3); Dispose(k4);
     Dispose(k5); Dispose(k6); Dispose(k7); Dispose(g6); Dispose(g7);
     Dispose(y_bad); Dispose(y_good)

End;

END. {of unit}

{ -------------------------- END  uawp.pas --------