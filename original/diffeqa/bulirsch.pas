{***********************************************************************
*                                                                      *
* Solve a system of first degree ordinary differential equations using *
* -------------------------------------------------------------------- *
* the extrapolation method of Bulirsch-Stoer-Gragg                     *
* -------------------------------------------------                    *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Jobst Hoffmann (FORTRAN)                       *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 3.11.1992                                      *
*                                                                      *
*                                  TPW Release By J-P Moreau, Paris.   *
*                                         (www.jpmoreau.fr)            *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C,  By Gisela Engeln-Muellges       *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
***********************************************************************}
UNIT Bulirsch;

Interface

Uses WinCrt1, Fgauss, T_dgls;


  Function bul_stoe   { Extrapolation method for 1st order DE systems ..}
      (
       Var x: Double;            { initial x-value/ final x-value .....}
           xend: Double;         { desired end point ..................}
           n: integer;           { number of DEs ......................}
           y:pVEC;               { initial y-value/ final y-value .....}
           epsabs,               { absolute error bound ...............}
           epsrel:Double;        { relative error bound ...............}
       Var h: Double;            { initial/final step size ............}
           hmax:Double;          { maximal step size ..................}
           neu:Boolean;          { use an outside given x ? ...........}
           fmax:integer;         { maximal # of calls of  dgl() .......}
       Var aufrufe:integer       { actual # of calls of dgl() .........}
      ): Integer;                { error code .........................}


Implementation

{ ------------------------------------------------------------------- }
Const

  {the Bulirsch sequence}
  bufol:array[0..11] of integer = (2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128);

Type

  pHilf  = ^Thilf;
  Thilf  = Array[0..11] of VECTOR;

Var

  ext_max: integer;                  {maximal order of the extrapolation}

{ ------------------------------------------------------------------- }

Function Max(a,b:Double):Double;
Begin
  if a>b then Max:=a
         else Max:=b
End;

Function IMin(a,b:integer): integer;
Begin
  if a<b then IMin:=a
         else IMin:=b
End;

Function Min(a,b:Double):Double;
Begin
  if a<b then Min:=a
         else Min:=b
End;

Procedure extrapol (
 Var row: integer;
     fhilf:pHilf;
     y:pVEC;
     n:integer;
     epsabs,
     epsrel: Double;
 Var x, h: Double;
     h0: Double;
 Var ahead: Boolean;
 Var index: integer
    );
{***********************************************************************
* Perform one extrapolation step for function bul_stoe                 *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* row     pointer to the arrays  bufol and fhilf                       *
* fhilf   pointer to Matrix for extrapolation values                   *
* y       y-value of solution at x                                     *
* n       number of equations in the system                            *
* epsabs  absolute error bound                                         *
* epsrel  relative error bound                                         *
* x       x-value                                                      *
* h       local step sizes                                             *
* h0                                                                   *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* ahead  *ahead  = 1:  Step is acceptable                              *
*        *ahead != 1:  Step not acceptable                             *
* index  pointer to the arrays bufol and fhilf                         *
* x      final x-value                                                 *
* h      final step size                                               *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
*        ext_max, bufol, TRUE, Min, Max, Copy_vector, Power            *
***********************************************************************}
Var
    column,         { column of extrapolation scheme                 }
    i: integer;     { loop counter                                   }
    diff_max,       { maximal difference of two columns in scheme    }
    fhilf_i1,       { column in extrapolation scheme                 }
    fhilf_i2,       { column to the right of fhilf_i1                }
    bufol_r,        { element in Bulirsch sequence                   }
    bufol_i,        { ditto                                          }
    y_max,          { maximal value in a column                      }
    help: Double;   { aux variable                                   }

    maxcol: integer;

    Function Min(a,b:integer): integer;
    begin
      if a<b then Min:=a else Min:=b
    end;

    Function Max(a,b:Double): Double;
    begin
      if a>b then Max:=a else Max:=b
    end;

Begin

  y_max := 0.0;
  maxcol := Min(row+1,ext_max);
  ahead:=False;

  for column := 2 to maxcol do
  begin
    if ahead then exit;
    index := Min(11 - column, row - column + 3);
    diff_max := 0.0;
    for i := 0 to n-1 do
    begin
      fhilf_i1 := fhilf^[index - 1, i];
      fhilf_i2 := fhilf^[index - 2, i];
      bufol_r  := 1.0*bufol[row];
      bufol_i  := 1.0*bufol[index - 2];
      fhilf^[index - 2, i] := fhilf_i1 + (fhilf_i1 - fhilf_i2)/(Sqr(bufol_r / bufol_i) - 1.0);
      fhilf_i1 := fhilf^[index - 1, i];
      fhilf_i2 := fhilf^[index - 2, i];
      y_max    := Max(y_max, ABS(fhilf_i2));
      diff_max := Max(diff_max, ABS(fhilf_i2 - fhilf_i1))
    end;
   
    if diff_max < epsrel * y_max + epsabs then    {Step acceptable ?}
    begin
      x := x + h0;
      {Copy_vector(y, fhilf[index - 2], n); }
      for i:=0 to n-1 do y^[i]:=fhilf^[index-2, i];
      help  := 1.0*(column - ext_max);
      {Step size for next step}
      h     := 0.9 * h * Exp(help*Ln(0.6));  {0.6^help}
      row   := -1;
      ahead := True
    end
  end
End;



{ ------------------------------------------------------------------- }
Function bul_stoe   { Extrapolation method for 1st order DE systems ..}
    (
     Var x: Double;            { initial x-value/ final x-value .....}
         xend: Double;         { desired end point ..................}
         n: integer;           { number of DEs ......................}
         y:pVEC;               { initial y-value/ final y-value .....}
         epsabs,               { absolute error bound ...............}
         epsrel:Double;        { relative error bound ...............}
     Var h: Double;            { initial/final step size ............}
         hmax:Double;          { maximal step size ..................}
         neu:Boolean;          { use an outside given x ? ...........}
         fmax:integer;         { maximal # of calls of  dgl() .......}
     Var aufrufe:integer       { actual # of calls of dgl() .........}
    ): Integer;                { error code .........................}
{***********************************************************************
* Compute the solution of an intial value problem for a first order    *
* system of ordinary differential equations                            *
*                                                                      *
*      y' = f(x,y)        with initial value y(x0) = y0                *
*                                                                      *
* by using the extrapolation method of Bulirsch-Stoer-Gragg.           *
* The maximal extrapolation order is set internally depending on the   *
* machine constant; the step size control follows the ideas in :       *
*     HALL, G.; WATT, J.M: Modern Numerical Methods for Ordinary       *
*                          Differential Equations, Oxford 1976,        *
*                          Clarendon Press,  [HALL76]                  *
* REMARK :                                                             *
* The dynamic allocations from this function are only set free for a   *
* new run or with improper input (fehler = TRUE or Return value >= 2). *
* Hence we advise the user to end work with this function by one call  *
* with an improper value such as n = -2 in order to free up all storage*
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        x-value                                                     *
* y        [0..n-1] vector with initial y-value at x                   *
* dgl      pointer to a function which evaluates the right hand side   *
*          odf the DE system   y' = f(x,y)                             *
* n        number of equations in the system                           *
* xend     desired final value for x, xend < x is allowed              *
* h        step size for next step                                     *
* hmax     maximal step size; hmax > 0                                 *
* epsabs\  error bounds for absolute and relative errors, both non     *
* epsrel/  negative.                                                   *
*          We test for the mixed error:                                *
*              |lokaler Fehler|  <=  |y| * epsrel + epsabs.            *
*          If epsrel = 0, we test the absolute error; if epsabs = 0 we *
*          test the relative error. epsrel will be set to 10 * machine *
*          constant if it has been chosen smaller than this.           *
* neu      Flag, indicating whether this function is called to compute *
*          on after one partial run that was cut short due to too many *
*          evaluations:                                                *
*          neu = FALSE: continue with old values                       *
*          neu = TRUE:  start new                                      *
* fmax     upper bound for the number of allowed function evaluations  *
*          of the right hand side via dgl()                            *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x        final x-value of integration; usually x = xend              *
* y        [0..n-1] y-value vector at x                                *
* h        final step size; should be used unchanged in the next call  *
*          if neu = TRUE; should be redefined for neu = TRUE.          *
* aufrufe  number of calls of dgl()                                    *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* = 0: all ok; after a new setting for xend, the function can be used  *
*      again for the same problem, or all inputs may be changed        *
* = 1: method has not reached xend after  MAX_FCT right hand side      *
*      evaluations; try another call with same parameters or change    *
*      error bounds.                                                   *
* = 2: step size less than 4 * machine constant.                       *
*      For a future call, increase step size and error bounds.         *
* = 3: epsabs or epsrel negative, or both zero.                        *
* = 4: xend = x                                                        *
* = 5: hmax negative.                                                  *
* = 6: n <= 0.                                                         *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* extrapol, ext_max, bufol, Double, TRUE, FALSE, MACH_EPS, Min, IMin,  *
* Max, boolean, copy_vector, ABS, Exp, Ln.                             *
***********************************************************************}
Label fin;
Var

  fhilf: pHilf;               { Index 0..7: Matrix for the           }
                              {             extrapolation scheme     }
                              { Index 8..11: aux vectors             }
  x0:Double;                  { x-value on call of dgl()             }
  row:integer;                { row of extrapolation scheme          }
  absh,                       { magnitude of step size               }
  absh0,                      { magnitude of aux step size  h0       }
  h0,                         { aux step size                        }
  hilf: Double;               { aux variable                         }
  i, j,                       { loop counters                        }
  count,                      { Loop variable for mid-point rule     }
  index: integer;             { Index for the vectors fhilf and bufol}
  ahead: boolean;             { Flag indicating acceptance of step   }

  tmp,tmp1:pVec;              { Auxiliary vectors of size n          }
  icount: integer;
Begin

  New(tmp); New(tmp1); New(fhilf);

  For i:=0 to 11 do
    For j:=0 to NMAX do
      fhilf^[i,j] := ZERO;

  h0 := 0.0;

  {check input}
  if (epsabs < ZERO) or (epsrel < ZERO) or (epsabs + epsrel <= ZERO) then
  begin
    bul_stoe := 3;
    goto fin
  end
  else if xend = x then
  begin
    bul_stoe := 4;
    goto fin
  end
  else if hmax <= ZERO then
  begin
    bul_stoe := 5;
    goto fin
  end;
  if n <= 0 then
  begin
    bul_stoe := 6;
    goto fin
  end;

  if (Not neu) then  {new call or repeat call due to excessive evaluations? }
  begin
    x := fhilf^[9][0];
    {copy_vector(y, fhilf[8], n);}
    For i:=0 to n-1 do y^[i]:=fhilf^[8, i]
  end
  else
    row:=-1;
 
  aufrufe := 0;
  ahead   := TRUE;
  ext_max := Round(-Ln(MACH_EPS) / Ln(2.0) / 7.5);

  epsrel  := Max(epsrel, 10.0 * MACH_EPS);

  icount:=0;

  {main loop}
  While (True) do
  begin
    Inc(icount);
    if (neu) then
    begin
      if (ahead) then                                       { new step ? }
      begin
        absh := ABS(h);
        absh := Min(absh, hmax);
        h0 := Min(absh, ABS(xend - x));
        if xend < x then h0:=-h0;
        absh0 := ABS(h0);

        if absh0 <= 4.0 * MACH_EPS * ABS(x) then
        begin
          bul_stoe := 0;
          goto fin                                          {normal exit}
        end;
        ahead := FALSE
      end;

      Repeat
        Inc(row);                    { find step size for extrapolation }
        h := h0 / bufol[row];
        x0 := x;                     { Euler step; save initial values  }

        {copy_vector(fhilf[8], y, n);}
        For i:=0 to n-1 do fhilf^[8, i] := y^[i];

        For i:=0 to n-1 do
        begin
          tmp^[i]:=fhilf^[8, i];
          tmp1^[i]:=fhilf^[11, i]
        end;

        dgl(x0, tmp, tmp1);

        For i:=0 to n-1 do
        begin
          fhilf^[8, i]:=tmp^[i];
          fhilf^[11, i]:=tmp1^[i]
        end;

        for i := 0 to n-1 do
          fhilf^[9,i] := fhilf^[8,i] + h * fhilf^[11,i];
        x0 := x0 + h;

        For i:=0 to n-1 do
        begin
          tmp^[i]:=fhilf^[9, i];
          tmp1^[i]:=fhilf^[11, i]
        end;

        dgl(x0, tmp, tmp1);

        For i:=0 to n-1 do
        begin
          fhilf^[9, i]:=tmp^[i];
          fhilf^[11, i]:=tmp1^[i]
        end;

        { use mid-point rule  }
        for count := 1 to bufol[row] - 1 do
        begin
          for i := 0 to n-1 do
            fhilf^[10, i] := fhilf^[8, i] + 2.0 * h * fhilf^[11, i];

          x0 := x0 + h;

          For i:=0 to n-1 do
          begin
            tmp^[i]:=fhilf^[10, i];
            tmp1^[i]:=fhilf^[11, i]
          end;

          dgl(x0, tmp, tmp1);

          For i:=0 to n-1 do
          begin
            fhilf^[10, i]:=tmp^[i];
            fhilf^[11, i]:=tmp1^[i]
          end;

          for j := 8 to 10 do                    { store for next step }
          begin
            {copy_vector(fhilf[j], fhilf[j + 1], n); }
            for i:=0 to n-1 do fhilf^[j, i] := fhilf^[j+1, i]
          end
        end;

        For i:=0 to n-1 do
        begin
          tmp^[i]:=fhilf^[9, i];
          tmp1^[i]:=fhilf^[11, i]
        end;

        dgl(x0, tmp, tmp1);

        For i:=0 to n-1 do
        begin
          fhilf^[9, i]:=tmp^[i];
          fhilf^[11, i]:=tmp1^[i]
        end;

        for i := 0 to n-1 do            {stabilize with trapezoidal rule}
          fhilf^[row, i] := 0.5 * (fhilf^[9, i] + fhilf^[8, i] + h * fhilf^[11, i]);
        aufrufe := aufrufe + bufol[row] + 2
      Until row <> 0;

      {Extrapolation}
      extrapol(row, fhilf, y, n, epsabs, epsrel, x, h, h0, ahead, index);

      if aufrufe >= fmax then                        { too many calls ?}
      begin                                          { => stop         }
        fhilf^[9][0] := x;
        {copy_vector(fhilf[8], y, n); }
        for i:=0 to n-1 do fhilf^[8, i] := y^[i];
        x := x0;
        {copy_vector(y, fhilf[index - 2], n); }
        for i:=0 to n-1 do y^[i] := fhilf^[index-2, i];
        bul_stoe := 1;
        goto fin
      end
    end; {if neu}

    if (Not ahead) or (Not neu) then    { do we need to repeat step ? }
    begin
      neu := TRUE;                          { set flag for a new call }
      if row >= IMin(7, ext_max - 1) then

      { store differently since the extrapolation scheme has at most }
      { 11 rows, but we allow only 8 extrapolations                  }

        for j := 0 to 6 do
        begin
          {copy_vector(fhilf[j], fhilf[j + 1], n); }
          for i:=0 to n-1 do fhilf^[j, i]:=fhilf^[j+1, i]
        end;

      { accuracy could not be reached from the whole extrapolation   }
      { table; repeat step with smaller step size                    }

      if row >= ext_max + 2 then
      begin
        hilf := 1.0*(row - ext_max + 1);
        h0   := 0.9 * h * Exp(hilf*Ln(0.6));                {0.6^hilf}
        if ABS(h0) <= 4.0 * MACH_EPS * ABS(x0) then   {step too small}
        begin
          bul_stoe :=  2;
          goto fin
        end;
        row := -1
      end
    end
  end;  {while true}
fin:   Dispose(tmp); Dispose(tmp1); Dispose(fhilf)
End;

END.

{ -------------------------- END bulirsch.pas ----------------------