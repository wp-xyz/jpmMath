{***********************************************************************
*                                                                      *
* Solve an ordinary system of differential equations of first order    *
* -------------------------------------------------------------------- *
* using the predictor-corrector method of Adams-Bashforth-Moulton      *
* ----------------------------------------------------------------     *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Turbo C 2.0                                    *
* Computer:             IBM PS/2 70 with 80387                         *
* Author:               Jobst Hoffmann (FORTRAN)                       *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 7.9.1992                                       *
*                                                                      *
*                                  TPW Release By J-P Moreau, Paris.   *
*                                         (www.jpmoreau.fr)            *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges        *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
***********************************************************************}
UNIT AB_MOU;

INTERFACE

Uses WinCrt1, FGauss, T_dgls;

Function prae_korr {Predictor-corrector meth. for 1st order DE systems}
             (
          Var x: Double;             { initial/final x-value ........}
              y: pVEC;               { initial value/ solution ......}
              n: Integer;            { number of DEs ................}
              xend: Double;          { desired final x-value ........}
          Var h: Double;             { starting/final step size .....}
              epsabs,                { absolute error bound .........}
              epsrel: Double;        { relative error bound .........}
              fmax: Integer;         { maximal # of calls for dgl() .}
          Var aufrufe:Integer;       { actual # of calls of dgl() ...}
              hmax: Double;          { maximal step size ............}
              neu: boolean          { delete old data ? ............}
             ): Integer;             { error code ...................}

IMPLEMENTATION

Const  HALF = 0.5;

Type               

Hilfstyp = Record  { Rcord of auxiliary vectors .....................}
       f,          { buffder for the starting values:  f[1],...,f[4] }
                   { contain the newest entries from rk_start() and  }
                   { abm_schritt(), f[0] is also needed              }
       tmp,        { aux buffer for the vectors k[i] from one Runge- }
                   { Kutta step; used in abm_schritt() for a linear  }
                   { combination of nodes in correktor step          }
       ki_sum,     { Linear combination A[0]*k[0]+..+A[3]*k[3] in the}
                   { Runge-Kutta method and end result of            }
                   { rk_schritt()                                    }
       y1,         { in rk_start() solution for step size  3*h, or   }
                   { in abm_schritt() solution of the predictor step }
                   { (also used for aux. purposes in rk_schritt() )  }
       y2,         { newest approximation                            }
       diff:VECTOR { error estimate for computed solution            }
End;

{For debug only}
Procedure Write_Hilf(hilf:Hilfstyp; n:integer);
  procedure printvec(V:VECTOR; n:integer);
  var i:integer;
  begin
    for i:=0 to n-1 do write(' ',V[i]);
    writeln
  end;
Begin
  writeln(' Record hilf');
  With hilf do
  begin
    writeln(' Vector f:');
    printvec(f,5);
    writeln(' Vector tmp:');
    printvec(tmp,n);
    writeln(' Vector ki_sum:');
    printvec(ki_sum,4);
    writeln(' Vector y1:');
    printvec(y1,n);
    writeln(' Vector y2:');
    printvec(y2,n);
    writeln(' Vector diff:');
    printvec(diff,n)
  end
End;

Function Max(a,b:double): double;
Begin
  if a>b then Max := a
         else Max := b
End;

Function Min(a,b:double): double;
Begin
  if a<b then Min := a
         else Min := b
End;

Function norm_max    { Find the maximum norm of a REAL vector .........}
           (
            vektor:pVEC;                 { vector ................. }
            n:INTEGER                    { length of vector ....... }
           ): Double;                    { Maximum norm ........... }
{***********************************************************************
* Return the maximum norm of a [0..n-1] vector  v.                     *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
*  None                                                                *
***********************************************************************}
VAR
  norm : Double;                                           { local max }
  betrag : Double;                          { magnitude of a component }
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

Procedure Copy_Vector(a,b:pVec;n:integer);
Var i:integer;
Begin
  For i:=0 to n-1 do a^[i]:=b^[i]
End;

Procedure Init0_vector(vektor:pVEC; n: integer);
{***********************************************************************
* initialize the [0..n-1] vector vektor to be the zero vector          *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
*  ZERO                                                                *
***********************************************************************}
Var i:integer;
Begin
  for i:=0 to n-1 do vektor^[i] := ZERO
End;

Procedure Inc_vector(ziel,quelle:pVEC; factor:double; n:integer);
{***********************************************************************
* add factor * quelle to the vector ziel                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
*  None                                                                *
***********************************************************************}
Var i:integer;
Begin
  for i:= 0 to n-1 do
    ziel^[i] := ziel^[i] + factor * quelle^[i]
End;

Procedure Add_vector(summe,summand1,summand2:pVEC; factor:Double; n:integer);
{***********************************************************************
* add factor * summand2 to the [0..n-1] vector summand1, store result  *
* summe.                                                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
*  None                                                                *
***********************************************************************}
Var i:integer;
Begin
  for i:=0 to n-1 do
    summe^[i] := summand1^[i] + factor * summand2^[i]
End;

Procedure rk_schritt(
                      x0: Double;
                      y0:pVEC;
                      n: integer;
                  Var hilf: Hilfstyp;
                      h:Double
                    );
{***********************************************************************
* perform one Runge-Kutta integration                                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0    x-value from wchich to integrate                               *
* y0    initial value of solution at  x0                               *
* n     number of DEs                                                  *
* h     step size for current Runge-Kutta integration                  *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* hilf  stores the result in ki_sum. The entries in y1 and tmp have no *
*       meaning.                                                       *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* hilfstyp, inc_vector, add_vector, ZERO, ONE.                         *
***********************************************************************}
Var A:  Array[0..3] of Double;
    a1: Array[0..3] of Double;
    i, j: integer;
Begin
  { ---- Coefficients of the classical Runge-Kutta method ----------- }
  A[0] := ONE / 6.0; A[1] := ONE / 3.0;
  A[2] := ONE / 3.0; A[3] := ONE / 6.0;
  a1[0] := ZERO; a1[1] := HALF;
  a1[2] := HALF; a1[3] := ONE;                {use also as the b[j,s] }

  dgl(x0, y0, @hilf.tmp);         {k[0] <- right hand side at (x0,y0) }

  for i := 0 to n-1 do                        {initialize ki_sum with }
    hilf.ki_sum[i] := A[0] * hilf.tmp[i];     {first term             }

  for j := 1 to 3 do             {compute k[1]..k[3], add to ki_sum   }
  begin                          {multiplied by factor                }
    add_vector(@hilf.y1, y0,                     {y1 <- y0 +          }
               @hilf.tmp, a1[j] * h, n);         {a1[j] * h * k[j - 1]}
    dgl(x0 + a1[j] * h, @hilf.y1, @hilf.tmp);         {compute k[j]   }

    inc_vector(@hilf.ki_sum,
               @hilf.tmp, A[j], n)             {ki_sum += A[j] * k[j] }
  end
End;


Function rk_start (
                    x: Double;
                Var x0: Double;
                    y: pVEC;
                    n: integer;
                    xend: Double;
                Var h: Double; 
                    hmax: Double;
                    new_step: Boolean;
                Var methode, aufrufe: integer;
                Var hilf: Hilfstyp
                  ): Integer;
{***********************************************************************
* Find the starting values using the Runge-Kutta method as needed in   *
* prae_korr() for the Adams-Bashforth-Moulton method                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x         x-value for start of integration                           *
* y         Initial value for the solution at x                        *
* n         number of DE equations                                     *
* xend      x-value where solution is wanted; xend may be less than x  *
* h         Step size                                                  *
* hmax      maximal step size, must be positive                        *
* new_step  <> 0: check new step size for properness                   *
*            = 0: do not check new step size                           *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x0       x-value until which we have integrated                      *
* h        step size for next integration                              *
* methode  value 0; hence the error estimate in prae_korr() is per-    *
*          formed for the factor for Runge-Kutta values                *
* aufrufe  actual number of calls of dgl                               *
* hilf     Record, which contains info of the solution:                *
*          y1:         approximate solution needed for error estimation*
*                      (Step size  3*h)                                *
*          y2:         actual approximate solution                     *
*          f[2]..f[4]: starting values                                 *
*          The entries in  ki_sum and tmp represent auxiliary values.  *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all ok                                                          *
* = 1: new step size too small relative to machine constant            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* Hilfstyp, rk_schritt, inc_vector, add_vector, MACH_EPS, Min,         *
* copy_vector, THREE, EIGHT.                                           *
***********************************************************************}
Label fin;
Var
    j: integer;                                   { loop counter     }

Begin
  if (new_step) then            { check new step size for properness }
  begin
    h := Min(h, hmax);
    h := Min(ABS(h), ABS(xend - x)/3.0);
    if xend <= x then h := -h;
    if ABS(h) <= 8.0 * MACH_EPS * ABS(x) then
    begin
      rk_start := 1;
      goto fin
    end
  end;

  x0 := x;                                      { save initial value }

  copy_vector(@hilf.y2, y, n);                             { y2 <- y }

  for j := 2 to 4 do                      { three steps with steps h }
  begin
    rk_schritt(x, @hilf.y2, n, hilf, h);         { accumulate ki_sum }
    x := x + h;
    inc_vector(@hilf.y2, @hilf.ki_sum, h, n);     { y2 += h*ki_sum   }
    dgl(x, @hilf.y2, @hilf.f[j])        { copy remainder of starting }
  end;                                  { values in  f[2]..f[4]      }

  { after three steps with size h, we now perform one step with 3*h  }
  { result put into y1.                                              }

  rk_schritt(x0, y, n, hilf, 3.0 * h);            { compute ki_sum   }
  add_vector(@hilf.y1, y, @hilf.ki_sum, 3.0*h, n);

  x0      := x0 + 3.0 * h;
  methode :=  0;
  aufrufe := aufrufe + 13;   {13 calls of dgl for the starting values}

  rk_start := 0;
fin: End;


{ ---------------------------------------------------------------------}
Procedure abm_schritt (
                    Var x0: Double;
                        n:  Integer;
                        h: Double;
                    Var methode, aufrufe: Integer;
                    Var hilf: Hilfstyp
                      );
{***********************************************************************
* Perform one step of the Adams-Bashforth-Moulton method               *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x0    initial x-value                                                *
* n     number of DEs in system                                        *
* h     step size                                                      *
* hilf  Structure with the vectors f[1]..f[4] for the starting values  *
*                                                                      *
* Output parameter:                                                    *
* =================                                                    *
* x0       final x-value of the integration                            *
* methode  equal to 1; hencethe error estimation in prae_korr() is     *
*          performed with the factor for Adams-Bashforth-Moulton.      *
* aufrufe  current number of calls of dgl()                            *
* hilf     Structure with the following output :                       *
*          y1:   approximate solution used for the error estimation    *
*                (from predictor step)                                 *
*          y2:   approximate solution from corrector step              *
*          f[4]: new node for starting values                          *
*          f[0],...,f[3] are alterd and tmp is used as aux storage     *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* Hilfstyp, Init0_vector, inc_vector, add_vector, ONE.                 *
***********************************************************************}
Var
  prae: Array[0..3] of Double;
  korr: Array[0..3] of Double;
  tmp: Double;            {aux vector for cyclic exchanges of adresses}
                          {of starting vectors                        }
  j: integer;             {loop counter                               }
Begin
  { ---- the coefficients for  Adams-Bashforth-Moulton -------------- }
  prae[0] :=  -9.0; prae[1] := 37.0;
  prae[2] := -59.0; prae[3] := 55.0;                      {Predictor  }
  korr[0] := ONE;  korr[1] := -5.0;
  korr[2] := 19.0; korr[3] :=  9.0;                       {corrector  }

  tmp := hilf.f[0];                 { the starting values are expected}
  for j := 0 to 3 do                { reside in f[0], ..., f[3] , but }
    hilf.f[j] := hilf.f[j+1];       { are stored in f[1]..f[4], we    }
  hilf.f[4] := tmp;                 { rotate their pointer            }

  init0_vector(@hilf.y1, n);                { one predictor step      }
  for j := 0 to 3 do
    inc_vector(@hilf.y1, @hilf.f[j], prae[j], n);
  add_vector(@hilf.y1, @hilf.y2, @hilf.y1, h/24.0, n);

  x0 := x0 + h;                          { move on by h; compute rhs  }
  dgl(x0, @hilf.y1, @hilf.f[4]);         { of DE at (x0,y1) : this    }
                                         { yields the start for the   }
                                         { corrector                  }

  init0_vector(@hilf.tmp, n);                { one corrector step     }
  for j := 0 to 3 do
    inc_vector(@hilf.tmp, @hilf.f[j+1], korr[j], n);
  inc_vector(@hilf.y2, @hilf.tmp, h/24.0, n);           { y2: new appr}
                                                        { solution    }
  dgl(x0, @hilf.y2, @hilf.f[4]);          { f[4] comtains the node for}
                                          { the starting values       }
                                         { derived from the corrector }
  methode := 1;
  aufrufe := aufrufe + 2           {2 calls of dgl for one  A-B-M step}
End;

{ --------------------------------------------------------------------}
Function prae_korr {Predictor-corrector meth. for 1st order DE systems}
             (
          Var x: Double;             { initial/final x-value ........}
              y: pVEC;               { initial value/ solution ......}
              n: Integer;            { number of DEs ................}
              xend: Double;          { desired final x-value ........}
          Var h: Double;             { starting/final step size .....}
              epsabs,                { absolute error bound .........}
              epsrel: Double;        { relative error bound .........}
              fmax: Integer;         { maximal # of calls for dgl ...}
          Var aufrufe:Integer;       { actual # of calls of dgl .....}
              hmax: Double;          { maximal step size ............}
              neu: boolean           { delete old data ? ............}
             ): Integer;             { error code ...................}
{***********************************************************************
* Solve a system of ordinary differential equations                    *
*                                                                      *
*      y' = f(x,y)        with the initial condition y(x0) = y0        *
*                                                                      *
* using the predictor-corrector method of Adams-Bashforth-Moulton.     *
* The needed starting values are obtained by an initial run of Runge-  *
* Kutta with the same order as the  A-B-M method.                      *
* We work with automatic step size control with doubling/halfing of    *
* the step size according to the subsequent error estimates.           *
*                                                                      *
* REMARK :                                                             *
*   This procedure frees certain buffers only at the beginning of a    *
*   new run or when memory is not sufficient, or for improper input.   *
*   Hence we advise to terminate a connected run of this function with *
*   one call with improper input such as n < 0 to clear all memory.    *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        x-value where a corresponding y-value is known              *
* y        [0..n-1] vector with the initial y-data                     *
* dgl      pointer to a function f which evaluates the right hand side *
*          of the DE system   y' = f(x,y)                              *
* n        number of DEs                                               *
* xend     x-value for which we want to find the solution; may be less *
*          than x                                                      *
* h        step size for next step; this is usually determined inside  *
*          this function                                               *
* epsabs\  error bounds for absolute and relative errors; each >= 0    *
* epsrel/  We apply a mixed test :                                     *
*              |local error|  <=  |y| * epsrel + epsabs.               *
*          For  epsrel = 0, we test for absolute error only.           *
*          For  epsabs = 0, we test for relative error.                *
*          Each error bound is internally set to 10 * machine constant *
*          if it should be prescribed to be too small.                 *
* fmax     upper bound for the allowed number of evaluations of the    *
*          right hand side  dgl() od the DE system                     *
* hmax     maximal  step size; must be positive                        *
* neu      If TRUE, we start off using Runge-Kutta, even if the earlier*
*          integration counld be continued with another A-B-M step     *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x        final x-value of integration (if run was successful: x=xend)*
* y        [0..n-1] approximate solution vector at x                   *
* h        final step size; should be used for start of next call.     *
*          If changed, please reset flag as well.                      *
* aufrufe  counter for calls of  dgl()                                 *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* error code :                                                         *
*   = 0: no error, xend was reached                                    *
*   = 1: after fmax calls of dgl() we have not reached xend:           *
*        a new call with unchanged parameters might be successful      *
*        (or try raising the error bounds)                             *
*   = 2: the step size is less than 8 * machine constant.              *
*        Before subsequent calls increase h and the error bounds.      *
*   = 3: epsabs or epsrel negative, or both equal to zero.             *
*   = 4: xend = x                                                      *
*   = 5: hmax negative                                                 *
*   = 6: n <= 0                                                        *
*   = 7: lack of available memory                                      *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* rk_start, abm_schritt, init_praeko, MACH_EPS, max, boolean, EIGHT,   *
* norm_max, copy_vector, REAL, ZERO, ONE, TEN, FALSE, TRUE.            *
***********************************************************************}
Label fin;
Const CHECK: Boolean = True;      { check new step size in rk_start  }

  { Vector with factors for error estimation :                       }
  { guess[0] for Runge-Kutta, guess[1] for Adams-Bashforth-Moulton   }
Var
  guess: Array[0..1] of Double;

  x0,                         { aux storage for       x              }
  h_save: Double;             { aux storage for       h              }
  
  methode: Integer;           { Number of most recently used method  }
                              { (Runge-Kutta or Adams-Bashforth-M.)  }

  hilf: Hilfstyp;             { local auxiliary variables            }

  ynorm,                      { Maximum norm from hilf.y2            }
  diffnorm: Double;           { Maximum norm from hilf.diff          }
  i,ii: integer;              { loop counters                        }
  rks : integer;              { Result of rk_start                   }
  folgeaufruf: Boolean;       { Flag indicating we can continue      }
                              { an earlier call with A-B-M method    }

Begin

  folgeaufruf := False;

  guess[0] := ONE/80.0;
  guess[1] := -19.0/270.0;

  { ------------------- Check input data --------------------------- }
  if (epsabs < ZERO) or (epsrel < ZERO) or (epsabs + epsrel <= ZERO) then
  begin
    folgeaufruf := FALSE;
    prae_korr := 3;
    goto fin
  end;
  if xend = x then
  begin
    folgeaufruf := FALSE;
    prae_korr := 4;
    goto fin
  end;
  if hmax <= ZERO then
  begin
    folgeaufruf := FALSE;
    prae_korr := 5;
    goto fin
  end;
  if n <= 0 then
  begin
    folgeaufruf := FALSE;
    prae_korr := 6;
    goto fin
  end;

  { -------------- Prepare integration loop ------------------------- }
  epsrel := Max(epsrel, 10.0 * MACH_EPS);
  epsabs := Max(epsabs, 10.0 * MACH_EPS);

  aufrufe := 1;

  if (neu) then                      {  Force start with R-K step ?  }
    folgeaufruf := FALSE;


  if (folgeaufruf) then              { Use values of previous call ? }
    abm_schritt(x0, n, h,            { Perform A-B-M step            }
                methode, aufrufe, hilf)

  else                               { very first call ?             }
  begin
    dgl(x, y, @hilf.f[1]);         { store beginning of starting     }
    Inc(aufrufe);                  { values in hilf.f[1]             }
    h_save := h;                   { save starting step size         }
    rks :=  rk_start(x, x0, y, n,  { new step size too small ?       }
                xend, h, hmax, CHECK,
                methode, aufrufe, hilf);
    if rks <> 0 then
    begin
      h           := h_save;
      folgeaufruf := FALSE;
      prae_korr := 0;
      goto fin
    end;
    { ---------- starting values available in                -------- }
    { ---------- hilf.f[1], hilf.f[2], hilf.f[3], hilf.f[4]. -------- }
  end;

  { ---------------------- Integration  loop ------------------------ }
  ii:=0;
  While True do
  begin
    Inc(ii);
    if aufrufe > fmax then               { excessive function calls ? }
    begin
      x := x0;
      copy_vector(y, @hilf.y2, n);    { newest approximation for  y   }
      folgeaufruf := TRUE;
      prae_korr := 1;                { stop and report excessive calls}
      goto fin
    end;

    for i := 0 to n-1 do                       { errro estimation     }
      hilf.diff[i] := guess[methode] * (hilf.y2[i] - hilf.y1[i]);

    diffnorm := norm_max(@hilf.diff, n);
    ynorm    := norm_max(@hilf.y2,   n);

    if diffnorm >= epsrel * ynorm + epsabs then   {error too large ?  }
    begin
      h := h * HALF;                             {halve the step size }
                                                 {and  repeat         }
      if ABS(h) <= 8.0*MACH_EPS*ABS(x0) then     {step size too small?}
      begin
        folgeaufruf := FALSE;
        prae_korr := 2;
        goto fin
      end;

      rk_start(x, x0, y, n,                  { compute new starting   }
               xend, h, hmax, Not CHECK,     { values with  Runge-    }
               methode, aufrufe, hilf);      { Kutta                  }
    end
    else                              { error not excessive ?         }
    begin                             { step was successful, continue }
      x := x0;                        { on with the previous step size}
      for i := 0 to n-1 do            { add estimated error onto new  }
        hilf.y2[i] := hilf.y2[i] + hilf.diff[i];      { approximation }
      copy_vector(y, @hilf.y2, n);    { newest approximation for y    }
      dgl(x0, y, @hilf.f[4]);         { correct last node for the next}
                                      { A-B-M step                    }
      Inc(aufrufe);

      {accuracy excessive ? }
      if diffnorm <= 0.02 * (epsrel * ynorm + epsabs) then
      begin
        h      := h + h;                          { double step size }
        h_save := Max(h, h_save);

        if ((h > ZERO) and (x0 >= xend)) or         { reached xend ? }
           ((h < ZERO) and (x0 <= xend)) then
        begin
          folgeaufruf := FALSE;
          prae_korr := 0;                             { all done !   }
          goto fin
        end;

        { Continue integration with doubled step size.               }
        { First find a new set of starting values via Runge-Kutta.   }

        copy_vector(@hilf.f[1],   {use final value of starting value }
                    @hilf.f[4],   {as the first entry for the new set}
                    n);

        if rk_start(x, x0, y, n,         { new step size too small ? }
                     xend, h, hmax, CHECK,
                     methode, aufrufe, hilf) <> 0 then 
        begin
          h           := h_save;
          folgeaufruf := FALSE;
          prae_korr := 0;
          goto fin
        end
      end
      else                                  { last step successful ? }
      begin
        if ((h > ZERO) and (x0 >= xend)) or { has xend been reached ?}
           ((h < ZERO) and (x0 <= xend)) then
        begin
          folgeaufruf := FALSE;
          prae_korr := 0;
          goto fin
        end;

        if ((h > ZERO) and (x0 + h >= xend)) or    {sufficintly close }
           ((h < ZERO) and (x0 + h <= xend)) then     {to  xend ?     }
        begin
          { The next step would lead beyond xend. Hence we noe use
            xend - x0 as the step size, and start anew using R-K.     }

          h := xend - x0;
          h_save := h;
          copy_vector(@hilf.f[1],     { assign the final entry of the }
                      @hilf.f[4],     { old starting values to be the }
                      n);             { first value of the new one    }
          if rk_start(x, x0, y, n,          {new step size too small? }
                       xend, h, hmax, CHECK,
                       methode, aufrufe, hilf) <> 0 then
          begin
            h           := h_save;
            folgeaufruf := FALSE;
            prae_korr := 0;
            goto fin
          end
        end
        else
          abm_schritt(x0, n, h,             { perform one A-B-M step }
                      methode,
                      aufrufe, hilf)
      end
    end
  end; {while True}
fin:End;

End.

{ -------------------------- END ab_mou.pas ---------------------------}