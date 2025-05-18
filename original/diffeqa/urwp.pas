{***********************************************************************
*                                                                      *
* Solve a two point boundary problem of first order with the shooting  *
* -------------------------------------------------------------------  *
* method                                                               *
* ------                                                               *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Klaus Niederdrenk (FORTRAN)                    *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 6.2.1992; 10.30.1995                           *
*                                                                      *
*                                                                      *
*                                  TPW Release By J-P Moreau, Paris.   *
*                                          (www.jpmoreau.fr)           *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C,  By Gisela Engeln-Muellges       *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
***********************************************************************}
UNIT Urwp;

INTERFACE

Uses WinCrt1, Fgauss, Uawp, Ab_mou, Bulirsch;

Type ar_real = DOUBLE;

Function rwp  {Shooting method for boundary value problem of 1st order}
       (
         a,                      {left end point .....................}
         b,                      {right end point ....................}
         h: ar_real;             {starting step size .................}
         y_start: pVec;          {initial approximation or solution ..}
                                 {initial value problem  y(a) ........}
         n: integer;             {number of differntial equations ....}
         awpnumm: integer;       {Number of desired IVP solver .......}
         epsawp,                 {error bound for initial value       }
                                 {problem ............................}
         epsrb: ar_real;         {error bound for boundary value      }
                                 {problem ............................}
         fmax,                   {maximal number of calls of dgl() ...}
         itmax: integer;         {maximal number of Newton iterations }
     var act_iter: integer       {actual number of Newton steps ......}
       ): integer;               {error code .........................}


IMPLEMENTATION

{print a vector of length n and its name}
Procedure zeig(s:String; v:pVec; n:integer);
Var i:integer;
Begin
  writeln(s);
  for i := 0 to n-1 do write(' ', v^[i]);
  writeln
End;

Procedure randbed(ya:pVec; yb:pVec; r:pVec);
Begin
  r^[0] := ya^[0];
  r^[1] := yb^[0]
End;

Function Power(y,x:ar_real):ar_real;
Begin
  Power:=0.0;
  if y<0.0 then exit
  else Power:=Exp(x*Ln(y))
End;

Procedure Copy_Vector(a,b:pVec;n:integer);
Var i:integer;
Begin
  For i:=0 to n-1 do a^[i]:=b^[i]
End;

{ ------------------------------------------------------------------- }
Function rwp  {Shooting method for boundary value problem of 1st order}
       (
         a,                      {left end point .....................}
         b,                      {right end point ....................}
         h: ar_real;             {starting step size .................}
         y_start: pVec;          {initial approximation or solution ..}
                                 {initial value problem  y(a) ........}
         n: integer;             {number of differntial equations ....}
         awpnumm: integer;       {Number of desired IVP solver .......}
         epsawp,                 {error bound for initial value       }
                                 {problem ............................}
         epsrb: ar_real;         {error bound for boundary value      }
                                 {problem ............................}
         fmax,                   {maximal number of calls of dgl() ...}
         itmax: integer;         {maximal number of Newton iterations }
     var act_iter: integer       {actual number of Newton steps ......}
       ): integer;               {error code .........................}

{***********************************************************************
* This function solves the first order boundary value problem          *
*                                                                      *
*     y' = F(x,y),    a <= x <= b,    R(y(a),y(b)) = 0.                *
*                                                                      *
* It uses the shooting method starting with an approximate y_start for *
* the associated initial value y(a), from which it constructs an       *
* approximate solution using as initial value solver the function awp()*
* The nonlinear system arising in the shooting method is solved using  *
* Newton's method.                                                     *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* a        left end point                                              *
* b        right end point  (b > a)                                    *
* h        suitable step size for solving the associated initial value *
*          problem for the shooting method approximately               *
* y_start  [0..n-1] vector, initial approximation for y(a)             *
* n        number of differential equations                            *
* dgl      pointer to the function which evaluates the right hand side *
*          of the differential equation  y' = f(x,y)                   *
* rand     pointer to the function that evaluates the boundary cond.   *
* awpnumm  label for the desired IVP solver in the shooting method:    *
*          = 1: Runge-Kutta embedding formula of 4/5th order; England  *
*               formula from awp().                                    *
*          = 2: Predictor-corrector method of order 4 by Adams-        *
*               Bashforth-Moulton (from prae_korr())                   *
*          = 3: Runge-Kutta embedding formula of 7/8th order (from     *
*               einb_rk())                                             *
*          = 4: Extrapolation method of  Bulirsch-Stoer (from          *
*               bul_stoe())                                            *
*          = 5: implicit Runge-Kutta-Gauss method (from implruku())    *
* epsawp   desired accuracy for the solution of the associated initial *
*          value problem                                               *
* epsrb    desired accuracy for which the approximation y_start for    *
*          the initial value  y(a) of a solution of the boundary value *
*          problem should satisfy the boundary condition R.            *
* fmax     upper bound for number of calls of dgl() in the system of   *
*          differential equations when solving the associated initial  *
*          value problem                                               *
* itmax    upper bound for number of Newton iterations when solving    *
*          the nonlinear systems in the shooting method                *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* y_start   [0..n-1] approximation for the initial value y(a) of a     *
*           solution  y of the boundary problem                        *
* act_iter  number of Newton iterations performed                      *
*                                                                      *
* Return value :                                                       *
* ==============                                                       *
* error code                                                           *
* = 0: all is ok                                                       *
* = 1: at least one of the accuracy bounds  epsawp, epsrb too small    *
* = 2: b <= a                                                          *
* = 3: h <= 0                                                          *
* = 4: n <= 0                                                          *
* = 5: improper input for awpnumm                                      *
* = 6: maximal number of allowed function evaluations exceeded         *
* = 7: act_iter > itmax. number of allowd Newton steps exceeded without*
*      finding a suitable value for  y_start.                          *
* = 8: The Jacobi matrix for Newton iterations is singular.            *
* = 9: lack of memory space                                            *
* >  9: error in one of the IVP solvers at first node                  *
* > 19: error in one of the IVP solvers at second node                 *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* ar_real, dglsysfnk, rndbedfnk, max, MACH_EPS, copy_vector, FABS,     *
* vminit, vmalloc, vmcomplete, vmfree, VEKTOR, VVEKTOR, MATRIX, awp,   *
* gauss, POW, ZERO, ONE, HALF, prae_korr, einb_rk, bul_stoe, implruku  *
***********************************************************************}
Label fin;
Const
      SINGU=1;             { Name for the return value of gauss(),    }
                           { which indicates that the matrix is       }
                           { singular                                 }
      ENGL45=6;            { awp() shall use the England embedding    }
                           { formula of 4th and 5th order             }

      FORMEL11=11;         { Macros for calling einb_rk(); to improve }
      MIT_HULL=1;          { readability                              }
      VON_VORNE=TRUE;
      NOT_SAVE=FALSE;

      MMAX=12;             { for calling implruku()                   }
      STUETZ='stuetz';
      AUS='';
      PROT='';

      HALF = 0.5;

Var
        Debug:Boolean;     { if True, print debugging infos            }
        yk,                { [0..n-1] vector with value of solution    }
                           { at right end point                        }
        yaj,               { [0..n-1] vector of modified y_start used  }
                           { for forming the Jacobi matrix             }
        r,                 { [0..n-1] vector with left boundary        }
                           { condition, which serves as the right hand }
                           { side for the linear system in the Newton  }
                           { method.                                   }
        rj,                { [0..n-1] vector with left end boundary    }
                           { condition of the modified initial value   }
                           { problem.                                  }
        d: pVec;           { [0..n-1] vector, the solution of the      }
                           { linear system.                            }
        amat: pMat;        { [0..n-1,0..n-1] array of the Jacobi       }
                           { matrix for the  Newton step.              }
        yk2,               { [0..n-1] vector with the values of the    }
                           { solution of the DE at the left endpoint   }
                           { for implruku().                           }
        g: pVec;           { [0..n-1] vector of weights for implruku() }
        xk,                { left end point for  awp()                 }
        hk,                { desired step size for awp()               }
        epsabs,            { absolute accuracy desired in  awp()       }
        epsrel,            { relative accuracy in awp()                }
        epsrel2,           { relative error bound for implruku()       }
        epsrel3,           { same as `epsrel2', but possibly altered   }
                           { during last call of  implruku()           }
        rmax,              { Maximum norm of the left boundary         }
                           { condition                                 }
        delta,             { Aux variable for  Jacobi matrix           }
        mach1,             { accuracy bounds depending on computer     }
        mach2: ar_real;    {                                           }
        aufrufe:integer;   { number of function evaluations in awp()   }
        pivot: pIVec;      { Permutation vector for Gaussian elimin.   }
        i,                 { loop counter                              }
        jacobi,            { counter for columns of Jacobi matrix      }
        mark,              { error code from gauss()                   }
        sign_det,          { sign of the permutation in gauss()        }
        fehler: integer;   { error code in  awp()                      }

Begin

  { --------- allocate storage for aux vectors ------------------------- }
  New(yk); New(yaj); New(r); New(rj); New(d); New(amat); New(yk2); New(g);
  New(pivot);

  mach1  := Power(MACH_EPS, 0.75);
  mach2  := 100.0 * MACH_EPS;

  Debug:=False;

  if epsrb < mach2 then rwp:=1;               {check input}
  if b <= a then rwp:=2;
  if h <= mach2 * ABS(a) then rwp:=3;
  if n <= 0 then rwp:=4;
  if (awpnumm < 1) or (awpnumm > 5) then rwp:=5;

  act_iter := 0;
  epsabs    := HALF * epsawp;
  epsrel    := epsabs;

  for i := 0 to n-1 do g^[i]:=ONE;   {weight vector for implruku(), set}
                                     {up as all ones.                  }

  if epsrel < 1E-10 then epsrel2:=1e-10 else epsrel2 := epsrel;

  {*********************************************************************
  * If  y_start proves to be a sufficiently good approximation for y(a)*
  * in the following loop, we report the success and check the         *
  * remaining input data in awp().                                     *
  *********************************************************************}

  While True do
  begin

    if Debug then zeig('Y_start:', y_start, n);

    copy_vector(yk, y_start, n);
    xk := a;
    hk := h;

    Case awpnumm of
      1: awp(xk, b, n, yk, epsabs, epsrel, hk, ENGL45, fmax, aufrufe, fehler);
      2: fehler := prae_korr(xk, yk, n, b, hk, epsabs,epsrel, fmax, aufrufe, ONE,VON_VORNE);
      3: fehler := bul_stoe(xk, b, n, yk, epsabs, epsrel,hk, ONE, VON_VORNE, fmax, aufrufe)
    end;

    if fehler <> 0 then
    begin
      rwp := 10+fehler;      {error in IVP solver ? }
      goto fin
    end;

    if Debug then zeig('Yk:', yk, n);

    randbed(y_start, yk, r);

    rmax:=0.0;
    for i := 0 to n-1 do rmax := Max(rmax, ABS(r^[i]));
    if rmax < epsrb then              {boundary condition satisfied ? }
    begin                                     {report success         }
      rwp := 0;
      goto fin
    end;

    Inc(act_iter);
    if act_iter > itmax then      {approximation not good enough      }
    begin                         {after  itmax Newton steps ?        }
      rwp := 6;                   {report failure                     }
      goto fin
    end;

    {find a better approximation  y_start for y(a) by using the
     Newton method with a Jacobi matrix that is approximated by one-
     sides difference quotients    }

    for jacobi := 0 to n-1 do
    begin
      copy_vector(yk,  y_start, n);
      copy_vector(yaj, y_start, n);

      if ABS(yk^[jacobi]) < mach2 then
      begin
        yk^[jacobi] := yk^[jacobi] + mach1;
        delta       := ONE / mach1
      end
      else
      begin
        yk^[jacobi] := yk^[jacobi] * (ONE + mach1);
        delta       := ONE / (mach1 * yk^[jacobi])
      end;

      yaj^[jacobi] := yk^[jacobi];
      xk := a;
      hk := h;

      if Debug then zeig('Yaj:', yaj, n);

      Case awpnumm of
        1: awp(xk, b, n, yk, epsabs, epsrel, hk, ENGL45, fmax, aufrufe, fehler);
        2: fehler := prae_korr(xk, yk, n, b, hk, epsabs, epsrel, fmax, aufrufe, ONE, VON_VORNE);
        3: fehler := bul_stoe(xk, b, n, yk, epsabs, epsrel, hk, ONE, VON_VORNE, fmax, aufrufe)
      End;

      if fehler <> 0 then                      {error in IVP solver ? }
      begin
        rwp := 20 + fehler;
        goto fin
      end;

      if Debug then zeig('Yk:', yk, n);

      randbed(yaj, yk, rj);

      for i := 0 to n-1 do
        amat^[i][jacobi] := (rj^[i] - r^[i]) * delta

    end;

    gauss(0, n, amat, amat, pivot, r, d, sign_det, mark);

    if mark = SINGU then                   {Jacobi matrix singular? }
    begin                                  {return error            }
      rwp := 8;
      goto fin
    end;                             

    if Debug then zeig('D:', d, n);

    for i := 0 to n-1 do                           {correct y_start }
      y_start^[i] := y_start^[i] - d^[i];

  end; {while True}
  {free memory}
fin:Dispose(yk); Dispose(yaj); Dispose(r); Dispose(rj); Dispose(d); Dispose(amat);
    Dispose(yk2); Dispose(g); Dispose(pivot);
End;

END.

{ ---------------------------- END rwp.pas -----------