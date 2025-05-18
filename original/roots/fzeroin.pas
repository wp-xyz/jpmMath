{ ------------------------ UNIT fzeroin.pas -------------------------- }
UNIT FZEROIN;
{***********************************************************************
*                                                                      *
* Compute a zero of a continuous real valued function with the         *
* ------------------------------------------------------------         *
* Zeroin method.                                                       *
* --------------                                                       *
*                                                                      *
* exported funktion:                                                   *
*   - zeroin():  Zeroin- method for real functions                     *
*                                                                      *
* Programming language: ANSI C                                         *
* Compiler:             Borland C++ 2.0                                *
* Computer:             IBM PS/2 70 mit 80387                          *
* Author:               Siegmar Neuner                                 *
* Modifications:        Juergen Dietel, Rechenzentrum, RWTH Aachen     *
* Source:               FORTRAN source code                            *
* Date:                 11. 27. 1992                                   *
* -------------------------------------------------------------------- *
* Ref.: "Numerical Algorithms with C By G. Engeln-Mueller and F. Uhlig,*
*        Springer-Verlag, 1996" [BIBLI 11].                            *
*                                                                      *
*                                 TPW version by J-P Moreau, Paris.    *
*                                         (www.jpmoreau.fr)            *
***********************************************************************}
INTERFACE


Const
      FOUR = 4.0; HALF = 0.5; ONE = 1.0; THREE = 3.0; TWO = 2.0;
      ZERO = 0.0;

      MACH_EPS = 1e-12;          { machine smallest real number       }


Procedure zeroin                 { Find roots with the Zeroin method  }
          (
       VAR abserr: REAL;         { absolute error bound ..............}
       VAR relerr: REAL;         { relative error bound ..............}
           fmax:   integer;      { maximal number of calls for fkt() .}
           protnam:String;       { Name of the log file ..............}
           a:      REAL;         { [a,b] = inclusion interval ........}
       VAR b:      REAL;         { right endpoint or zero ............}
       VAR fb:     REAL;         { Function value at the root b ......}
       VAR fanz:   integer;      { number of actual function calls ...}
       VAR rc:     integer       { error code ........................}
          );


IMPLEMENTATION

{test function F(x) }
Function Fkt(x:REAL): REAL;
Begin
  Fkt := 5*x - exp(x)
End;


Function Sign(a,b:REAL): REAL;
Begin
  if b>=ZERO then Sign := ABS(a)
  else Sign:=-ABS(a)
End;


Procedure zeroin                 { Find roots with the Zeroin method  }
          (
       VAR abserr: REAL;         { absolute error bound ..............}
       VAR relerr: REAL;         { relative error bound ..............}
           fmax:   integer;      { maximal number of calls for fkt() .}
           protnam:String;       { Name of the log file ..............}
           a:      REAL;         { [a,b] = inclusion interval ........}
       VAR b:      REAL;         { right endpoint or zero ............}
       VAR fb:     REAL;         { Function value at the root b ......}
       VAR fanz:   integer;      { number of actual function calls ...}
       VAR rc:     integer       { error code ........................}
          );
{***********************************************************************
* Given a real valued function f on an interval [a,b] with             *
* f(a) * f(b) < 0, we compute a root of f in [a,b].                    *
* The Zeroin method combines bisection and secant methods with inverse *
* quadratic interpolation.                                             *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* abserr\  Error bounds with   abserr >= 0  and  relerr >= 0. Their    *
* relerr/  sum must exceed zero. For break-off we use the test         *
*              |xm|  <=  0.5 * (abserr + |b| * relerr),                *
*          where xm denotes half the length of the final inclusion     *
*          interval. For relerr = 0 we test only for absolute error,   *
*          while for  abserr = 0, we test the relative error only.     *
*          abserr and relerr are used as put in only when both exceed  *
*          four times the machine constant. Otherwise we set them to   *
*          this value.                                                 *
* fmax     upper limit of calls of function                            *
* protnam  Name of a file used for intermediate results. If the pointer*
*          is set to zero, we do not use this file.                    *
* a,b      end points of the interval, that includes a root            *
*                                                                      *
* Output parameters:                                                   *
* =================                                                    *
* abserr\  actual error bounds used                                    *
* relerr/                                                              *
* b        approximate root                                            *
* fb       value of f at the root b                                    *
* fanz     number of calls of  fkt()                                   *
*                                                                      *
* Error code rc:                                                       *
* =============                                                        *
* = -2: abserr or relerr is negative, or both are zero, or fmax < 1.   *
* = -1: The necessary assumption  f(a) * f(b) < 0  is violated.        *
* =  0: Desired accuracy has been reache :                             *
*           |xm|  <=  0.5 * (abserr + |b| * relerr).                   *
* =  1: b on output is a root with  fkt(b) = 0.                        *
* =  2: Either a or b is a root at input.                              *
* =  3: After fmax calls of Fkt the algorithm has not found a root     *
* =  4: Error when opening intermediate result file                    *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, REALFCT, FOUR, MACH_EPS, ZERO, FABS, HALF, NULL, bi, TWO, ONE, *
* THREE, FILE, fprintf                                                 *
***********************************************************************}
Label  50, fin;
Var
       fa,                 { Function value fkt(a)                    }
       fc,                 { Function value fkt(c)                    }
       eps,                { minimal value for error bounds abserr    }
                           { and relerr                               }
       tol1,               { auxiliary variable for mixed error bound }
       xm,                 { half of the current interval length      }
       c,                  { value inside [a,b]                       }
       d,                  { Distance to the nearest approximate root }
       e,                  { previous value of d                      }
       p,
       q,
       r,
       s,
       tmp: REAL;          { auxiliary variable to check inclusion    }
                           {   f(a) * f(b) < 0                        }
       fehler:integer;     { error code of this function              }
       prodat:TEXT;        { intermediate data file                   }

Begin
  { ------------------ initialize variables ------------------------- }
  c:=ZERO; d:=ZERO; e:=ZERO; rc:=0; fehler:=0;
  eps := FOUR * MACH_EPS;

  fa   := Fkt(a);                    { evaluate fkt() at the end      }
  fb   := Fkt(b);                    { points a and b                 }
  fanz := 2;

  { ----------- check   f(a) * f(b) < 0  ---------------------------- }
  tmp := fa * fb;
  if tmp > ZERO then
    rc := -1
  else if tmp = ZERO then
    rc := 2;
  if rc<>0 then goto fin;

  { ----------- check usability of given error bounds --------------- }
   if (abserr < ZERO) or (relerr < ZERO) or (abserr + relerr <= ZERO) or (fmax < 1) then
   begin
    rc := -2;
    goto fin
  end;

  if (relerr = ZERO) and (abserr < eps) then
    abserr := eps
  else if (abserr = ZERO) and (relerr < eps) then
    relerr := eps
  else
  begin
    if abserr < eps then abserr := eps;
    if relerr < eps then relerr := eps
  end;

  if length(protnam)<>0 then
  begin
    Assign(prodat,protnam);
    Rewrite(prodat)
  end
  else
  begin
    rc := 4;
    goto fin
  end;

  fc:=fb;
  while (TRUE) do                         { start iteration }
  begin

    if fb * (fc / ABS(fc)) > ZERO then   { no inclusion of a root   }
    begin
      c  := a;                            { between b and c ?        }
      fc := fa;                           { alter c so that b and c  }
      e  := b - a;                        { include the root of f    }
      d  := e
    end;


    if ABS(fc) < ABS(fb) then        { If fc has the smaller modulus }
    begin
      a   := b;                       { interchange interval end      }
      b   := c;                       { points.                       }
      c   := a;
      fa  := fb;
      fb  := fc;
      fc  := fa
    end;

    writeln(prodat, ' a  = ', a:20:14, ' b  = ', b:20:14, ' c  = ', c:20:14);
    writeln(prodat, ' fa = ', fa:20:14, ' fb = ', fb:20:14, ' fc = ', fc:20:14);

    tol1 := HALF * (abserr + relerr * ABS(b));
    xm   := HALF * (c - b);

    writeln(prodat, ' xm =', xm,' tol1 = ', tol1);

    if ABS(xm) <= tol1 then          { reached desired accuracy ? }
    begin
      fehler := 0;
      goto 50                        { normal exit }
    end;

    if fb = ZERO then                { Is the best approximate root b }
    begin                            { a precise root of f ?          }
      fehler := 1;
      goto 50                        { other normal exit, b is a root }
    end;

    r := ZERO;

    if (ABS(e) < tol1) or (ABS(fa) <= ABS(fb)) then
    begin
      e := xm;
      d := e;
      writeln(prodat, ' Bisection');
    end
    else
    begin
      if a <> c then         { if a is not equal to c :              }
      begin
        q := fa / fc;        { With a, b and c we have 3 points for  }
        r := fb / fc;        { an inverse quadratic interpolation    }
        s := fb / fa;
        p := s * (TWO * xm * q * (q - r) - (b - a) * (r - ONE));
        q := (q - ONE) * (r - ONE) * (s - ONE)
      end
      else                   { If a equals  c :                      }
      begin
        s := fb / fa;        { Use the secant method or linear       }
        p := TWO * xm * s;   { interpolation                         }
        q := ONE - s
      end;

      if p > ZERO then       { alter the sign of  p/q for the        }
        q := -q              { subsequent division                   }
      else
        p := ABS(p);

      if ((TWO * p  >=  THREE * xm * q - ABS(tol1 * q)) or (p >= ABS(HALF * e * q))) then
      begin
        e := xm; d := e;
        writeln(prodat, ' Bisection');
      end
      else
      begin
        e := d;       { Compute the quotient p/q for both iterations }
        d := p / q;   { which will be used to modify b               }
        if r = ZERO then
          writeln(prodat, ' Secant method')
        else
          writeln(prodat, ' Inverse quadratic interpolation')
      end
    end;

    a  := b;         { store the best approximation b and its        }
    fa := fb;        { function value fb in a and fa.                }

    if ABS(d) > tol1 then                       { d large enough?    }
    begin
      b := b + d;                               { add d to b         }
      writeln(prodat, ' Difference d from new b: ', d:20:14)
    end
    else                                        { d too small?       }
    begin                                       { add tol1 to b      }
      b := b + Sign(tol1, xm);
      if xm < ZERO then
        writeln(prodat, ' Subtract error bound: d = ', -tol1:20:14)
      else
        writeln(prodat, ' Add error bound: d = ', tol1:20:14)
    end;

    fb := fkt(b);                      { compute function value at b }
    Inc(fanz);                         {      up evaluation counter  }

    writeln(prodat, ' b = ',b:20:14,' fb= ', fb:20:14,' Number of functional evaluations = ', fanz);

    if fanz > fmax then                { too many function evaluations? }
    begin
      fehler := 3;
      goto 50;
    end
  end;                                                  { end iteration }

50: close(prodat);

  rc := fehler;
fin: End;

END.
{ ------------------------ END fzeroin.pas --------------------------- }