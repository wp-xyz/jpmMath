{ ----------------------- UNIT fpegasus.pas ------------------------- }
UNIT FPEGASUS;

{BIBLI 11}

INTERFACE

Uses Fonction;

Const
       ITERMAX = 300;                 { Maximal number of iterations }

       ABSERR = 0.0;                  { Admissable absolute error    }
       MACH_EPS = 1E-15;              { Small real number            }
       RELERR = 4.0 * MACH_EPS;       { Admissable relative error    }
       FCTERR = MACH_EPS*MACH_EPS;    { Maximal function errror      }

  Function pegasus        { Pegasus Method    .........................}
            (
             x1: Double;            { Starting value 1 ................}
         Var x2: Double;            { Starting value 2 / solution .....}
         Var f2: Double;            { Function value at x2 ............}
         Var iter: Integer          { Number of iterations ............}
            ): Integer;


IMPLEMENTATION


   Function pegasus       { Pegasus Method    .........................}
            (
             x1: Double;            { Starting value 1 ................}
         Var x2: Double;            { Starting value 2 / solution .....}
         Var f2: Double;            { Function value at x2 ............}
         Var iter: Integer          { Number of iterations ............}
            ): Integer;
 {====================================================================*
 *                                                                    *
 *  pegasus computes one zero of the continuous function fct,         *
 *  provided that the two starting values x1 and x2 satisfy:          *
 *                fct(x1) * fct(x2) <= 0  .                           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Applications:                                                    *
 *   ============                                                     *
 *      Determine one root of the continuous function fct, if an      *
 *      inclusion interval [x1, x2] is known for the root.            *
 *      fct is defined in unit Fonction.pas.                          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      x1,x2    Double;                                              *
 *               Starting values with fct(x1) * fct(x2) <= 0.         *
 *                                                                    *
 *   Output parameters:                                               *
 *   =================                                                *
 *      x2       Double;                                              *
 *               Computed approximation for the root of fct           *
 *                                                                    *
 *      f2       Double;                                              *
 *               Functional value at the approximate root, this       *
 *               must be nearly zero.                                 *
 *      iter     integer;                                             *
 *               Number of iterations that were performed.            *
 *                                                                    *
 *   Return values:                                                   *
 *   =============                                                    *
 *      = -1     No inclusion: fct(x2) * fct(x1) > 0                  *
 *      =  0     Root has been found with ABS(f2) < FCTERR            *
 *      =  1     Break-off with                                       *
 *                   ABS(xnew-xold) < ABSERR + xnew * RELERR,         *
 *               check functional value                               *
 *      =  2     Iteration limit reached                              *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants used:   ABSERR, RELERR, MACH_EPS, EPSROOT, ITERMAX     *
 *   ==============                                                   *
 *                                                                    *
 *====================================================================}
Label 10, fin;
Var
    f1, x3, f3, s12: Double;
    rc: integer;
Begin

  rc := 2;
  iter := 0;                           { Initialize iteration counter }

  f1 := fct(x1);                         { Function values at x1, x2 }
  f2 := fct(x2);

  if f1 * f2 > 0.0 then
  begin
    pegasus := -1;                            { No inclusion -> Error  }
    goto fin
  end;

  if f1 * f2 = 0.0 then                 { One starting value is a root }
  begin
    if f1 = 0.0  then
    begin
      x2 := x1;
      f2 := 0.0
    end;
    pegasus:=0;
    goto fin
  end;

  while iter <= ITERMAX do              { Pegasus iteration          }
  begin

    Inc(iter);

    s12 := (f2 - f1) / (x2 - x1);            { Secant slope          }

    x3  := x2 - f2 / s12;                    { new approximation     }
    f3  := fct(x3);

    if f2 * f3 <= 0.0 then                  { new inclusion interval }
    begin
      x1 := x2;
      f1 := f2
    end
    else
      f1 := f1 * f2 / (f2 + f3);

    x2 := x3;
    f2 := f3;

    if ABS(f2) < FCTERR then                     { Root found        }
    begin
      rc := 0;
      goto 10
    end;
                                    { Break-off with small step size }
    if ABS(x2 - x1) <= ABS(x2) * RELERR + ABSERR then
    begin
      rc := 1;
      goto 10
    end
  end;

10: if ABS(f1) < ABS(f2) then     { Choose approximate root with    }
  begin                           { least magnitude function value  }
    x2 := x1;
    f2 := f1
  end;

  pegasus := rc;

fin: End;

END.

{ ---------------------- END fpegasus.pas ------------------------ }