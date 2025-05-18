{***********************************************************************
* Collection of Pascal subroutines to Integrate a System of Ordinary   * 
* Differential Equations By the Runge-Kutta-Fehlberg method (in double *
* precision).                                                          *
* -------------------------------------------------------------------- *
* REFERENCE:     H A Watts and L F Shampine,                           *
*                Sandia Laboratories,                                  *
*                Albuquerque, New Mexico.                              *
*                                                                      *
* See demonstration program trkf45.pas.                                *
*                                                                      *
*                                       TPW Release 1.1 By J-P Moreau. *
*                                            (www.jpmoreau.fr)         *
* -------------------------------------------------------------------- *
* Release 1.1 (03/15/05): added example #3 in procedure f().           *
***********************************************************************}
UNIT URKF45;

INTERFACE

Const NEQ = 5;                    {Maximum number of equations}
      EPSILON = 2.22e-16;

Type
     VEC   = Array[1..NEQ] of double;

Var
     num: Integer;                      {# of test}
     w1,w3,w4,w5,w6,w7:VEC;             {work space}
     w2,w8,w9: double;                  {work space}
     iw1,iw2,iw3,iw4,iw5: Integer;      {work space}

     {Visible from calling program}
     Procedure rkf45 (neqn:Integer; var y:VEC; t:double; tout:double;
                      relerr, abserr:double; var iflag:Integer);


IMPLEMENTATION

{User defined system of diff. equations}
Procedure f(t:double; y:VEC; var yp:VEC);
{--------------------------------------------------------------------
!
!  F evaluates the derivative for the ODE (TESTS #1 and #2).
!
!-------------------------------------------------------------------}
Begin
  Case num of
    1: yp[1] := 0.25 * y[1] * ( 1.0 - y[1] / 20.0);
    2: begin yp[1] := y[2];  yp[2] := - y[1] end;
    3: begin
         yp[1] := y[2];
         yp[2] := y[3];
         yp[3] := y[4];
         yp[4] := y[5];
         yp[5] := (45.0 * y[3] * y[4] * y[5] - 40.0 * y[4] * y[4] * y[4]) / (9.0 * y[3] * y[3])
       end
  end
end;


Procedure fehl(neqn:integer; Var y:VEC; t,h:double; yp:VEC;
                             Var f1, f2, f3, f4, f5, s:VEC);
{************************************************************************
!
!  FEHL takes one Fehlberg fourth-fifth order step (double precision).
!
!  Discussion:
!
!    FEHL integrates a system of NEQN first order ordinary differential
!    equations of the form
!      dY(i)/dT = F(T,Y(1),---,Y(NEQN))
!    where the initial values Y and the initial derivatives
!    YP are specified at the starting point T.
!
!    FEHL advances the solution over the fixed step H and returns
!    the fifth order (sixth order accurate locally) solution
!    approximation at T+H in array S.
!
!    The formulas have been grouped to control loss of significance.
!    FEHL should be called with an H not smaller than 13 units of
!    roundoff in T so that the various independent arguments can be
!    distinguished.
!
!  Author:
!
!    H A Watts and L F Shampine,
!    Sandia Laboratories,
!    Albuquerque, New Mexico.
!
!    RKF45 is primarily designed to solve non-stiff and mildly stiff
!    differential equations when derivative evaluations are inexpensive.
!
!  Parameters:
!
!    Input, external F, a subroutine of the form
!      Procedure f(t:double; y:VEC; var yp:VEC);
!    to evaluate the derivatives.
!      YP[I] := dY[I] / dT;
!
!    Input, integer NEQN, the number of equations to be integrated.
!
!    Input, double Y(NEQN), the current value of the dependent variable.
!
!    Input, double T, the current value of the independent variable.
!
!    Input, double H, the step size to take.
!
!    Input, double YP(NEQN), the current value of the derivative of the
!    dependent variable.
!
!    Output, double F1(NEQN), double F2(NEQN), double F3(NEQN), double F4(NEQN),
!    double F5(NEQN) are arrays of dimension NEQN which are needed for
!    internal storage.
!
!    Output, double S(NEQN), the computed estimate of the solution at T+H.
!******************************************************************************}
Var
    ch: double;
    i: integer;
Begin

  ch := h / 4.0;

  For i:=1 to neqn do
    f5[i] := y[i] + ch * yp[i];

  f(t + ch, f5, f1);

  ch := 3.0 * h / 32.0;

  For i:=1 to neqn do
    f5[i] := y[i] + ch * (yp[i] + 3.0 * f1[i]);
  f(t + 3.0 * h / 8.0, f5, f2);

  ch := h / 2197.0;

  For i:=1 to neqn do
    f5[i] := y[i] + ch * (1932.0 * yp[i] + (7296.0 * f2[i] - 7200.0 * f1[i]));

  f(t + 12.0 * h / 13.0, f5, f3);

  ch := h / 4104.0;

  For i:=1 to neqn do
    f5[i] := y[i] + ch * ((8341.0 * yp[i] - 845.0 * f3[i]) + (29440.0 * f2[i]
             - 32832.0 * f1[i]));

  f(t + h, f5, f4);

  ch := h / 20520.0;

  For i:=1 to neqn do
    f1[i] := y[i] + ch * ((-6080.0 * yp[i] + (9295.0 * f3[i] - 5643.0 * f4[i]))
             + (41040.0 * f1[i] - 28352.0 * f2[i]));

  f(t + h / 2.0, f1, f5);

{  Ready to compute the approximate solution at T+H.}

  ch := h / 7618050.0;

  For i:=1 to neqn do
    s[i] := y[i] + ch * ((902880.0 * yp[i] + (3855735.0 * f3[i] - 1371249.0 * f4[i]))
            + (3953664.0 * f2[i] + 277020.0 * f5[i]))

end; {fehl}


Procedure rkfs(neqn:integer; var y:VEC; var t:double; tout, relerr, abserr:double;
               var iflag:integer; var yp:VEC; var h:double; var f1, f2, f3, f4,f5:VEC;
               var savre, savae:Double; var nfe, kop, init, jflag, kflag:Integer);
               Forward;

Procedure rkf45 (neqn:Integer; var y:VEC; t:double; tout:double;
                 relerr, abserr:double; var iflag:Integer);
{*******************************************************************************
!
!  RKF45 carries out the Runge-Kutta-Fehlberg method (double precision).
!
!  Author:
!
!    H A Watts and L F Shampine,
!    Sandia Laboratories,
!    Albuquerque, New Mexico.
!
!    RKF45 is primarily designed to solve non-stiff and mildly stiff
!    differential equations when derivative evaluations are inexpensive.
!
!  Abstract:
!
!    RKF45 integrates a system of NEQN first order ordinary differential
!    equations of the form:
!
!      dY(i)/dT := F(T,Y(1),Y(2),...,Y(NEQN))
!
!    where the Y(1:NEQN) are given at T.
!
!    Typically the subroutine is used to integrate from T to TOUT but it
!    can be used as a one-step integrator to advance the solution a
!    single step in the direction of TOUT.  On return, the parameters in
!    the call list are set for continuing the integration.  The user has
!    only to call RKF45 again (and perhaps define a new value for TOUT).
!    Actually, RKF45 is an interfacing routine which calls subroutine
!    RKFS for the solution.  RKFS in turn calls subroutine FEHL which
!    computes an approximate solution over one step.
!
!  Reference:
!
!    E. Fehlberg,
!    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
!    NASA Technical Report R-315.
!
!    L F Shampine, H A Watts, S Davenport,
!    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
!    Sandia Laboratories Report SAND75-0182,
!    To appear in SIAM Review.
!
!  Parameters:
!
!    Input, external F, a subroutine of the form
!    Procedure f(t:double; y:VEC; var yp:VEC);
!    to evaluate the derivatives YP[I] := dY[I] / dT;
!
!    Input, integer NEQN, the number of equations to be integrated.
!
!    Input/output, double Y(NEQN), the solution vector at T.
!
!    Input/output, double T, the independent variable.
!
!    Input, double TOUT, the output point at which solution is desired.
!
!    Input, double RELERR, ABSERR, the relative and absolute error tolerances
!    for the local error test.  At each step the code requires:
!      abs(local error) <:= relerr*abs(y) + abserr
!    for each component of the local error and solution vectors
!
!    Output, integer IFLAG, indicator for status of integration.
!
!    Workspace, double W1(NEQN) to W9(NEQN), arrays to hold information internal
!    to RKF45 which is necessary for subsequent calls.
!
!    Workspace, integers IW1 to IW5, arrays used to hold information internal
!    to RKF45 which is necessary for subsequent calls.
!
!    Note: Workspaces are declared here statically in the INTERFACE section
!          of the present unit. 
!
!
!  first call
!
!    The user initializes the following parameters:
!
!      neqn -- number of equations to be integrated.  (neqn >:= 1)
!
!      y(neqn) -- vector of initial conditions
!
!      t -- starting point of integration , must be a variable
!
!      tout -- output point at which solution is desired.
!            t:=tout is allowed on the first call only, in which case
!            rkf45 returns with iflag:=2 if continuation is possible.
!
!      relerr,abserr -- relative and absolute local error tolerances
!            which must be non-negative. relerr must be a variable while
!            abserr may be a constant. the code should normally not be
!            used with relative error control smaller than about 1.e-8 .
!            to avoid limiting precision difficulties the code requires
!            relerr to be larger than an internally computed relative
!            error parameter which is machine dependent. in particular,
!            pure absolute error is not permitted. if a smaller than
!            allowable value of relerr is attempted, rkf45 increases
!            relerr appropriately and returns control to the user before
!            continuing the integration.
!
!      iflag -- +1,-1  indicator to initialize the code for each new
!            problem. normal input is +1. the user should set iflag:=-1
!            only when one-step integrator control is essential. in this
!            case, rkf45 attempts to advance the solution a single step
!            in the direction of tout each time it is called. since this
!            mode of operation results in extra computing overhead, it
!            should be avoided unless needed.
!
!
!  output from rkf45
!
!      y(neqn)   -- solution at t
!      t         -- last point reached in integration.
!      iflag = 2 -- integration reached tout. indicates successful retur
!                   and is the normal mode for continuing integration.
!            =-2 -- a single successful step in the direction of tout
!                   has been taken. normal mode for continuing
!                   integration one step at a time.
!            = 3 -- integration was not completed because relative error
!                   tolerance was too small. relerr has been increased
!                   appropriately for continuing.
!            = 4 -- integration was not completed because more than
!                   3000 derivative evaluations were needed. this
!                   is approximately 500 steps.
!            = 5 -- integration was not completed because solution
!                   vanished making a pure relative error test
!                   impossible. must use non-zero abserr to continue.
!                   using the one-step integration mode for one step
!                   is a good way to proceed.
!            = 6 -- integration was not completed because requested
!                   accuracy could not be achieved using smallest
!                   allowable stepsize. user must increase the error
!                   tolerance before continued integration can be
!                   attempted.
!            = 7 -- it is likely that rkf45 is inefficient for solving
!                   this problem. too much output is restricting the
!                   natural stepsize choice. use the one-step integrator
!                   mode.
!            = 8 -- invalid input parameters
!                   this indicator occurs if any of the following is
!                   satisfied -   neqn <= 0
!                                 t:=tout  and  iflag <> +1 or -1
!                                 relerr or abserr < 0.
!                                 iflag = 0  or < -2  or > 8
!      w1..w9, iw1..iw5 -- information which is usually of no interest
!                   to the user but necessary for subsequent calls.
!                   w1(neqn) contains the first derivatives of the
!                   solution vector y at t. w2 contains the stepsize
!                   h to be attempted on the next step. iw1 contains the
!                   derivative evaluation counter, etc.
!
!  subsequent calls:
!
!    RKF45 returns with all information needed to continue
!    the integration. if the integration reached tout, the user need onl
!    define a new tout and call RKF45 again.  In the one-step integrator
!    mode (iflag:-2) the user must keep in mind that each step taken is
!    in the direction of the current tout.  Upon reaching tout (indicated
!    by changing iflag to 2),the user must then define a new tout and
!    reset iflag to -2 to continue in the one-step integrator mode.
!
!    If the integration was not completed but the user still wants to
!    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3
!    the relerr parameter has been adjusted appropriately for continuing
!    the integration. in the case of iflag=4 the function counter will
!    be reset to 0 and another 3000 function evaluations are allowed.
!
!    However,in the case iflag=5, the user must first alter the error
!    criterion to use a positive value of abserr before integration can
!    proceed. if he does not, execution is terminated.
!
!    Also,in the case iflag=6, it is necessary for the user to reset
!    iflag to 2 (or -2 when the one-step integration mode is being used)
!    as well as increasing either abserr,relerr or both before the
!    integration can be continued. if this is not done, execution will
!    be terminated. the occurrence of iflag=6 indicates a trouble spot
!    (solution is changing rapidly,singularity may be present) and it
!    often is inadvisable to continue.
!
!    If iflag=7 is encountered, the user should use the one-step
!    integration mode with the stepsize determined by the code.
!    If the user insists upon continuing the integration with RKF45,
!    he must reset iflag to 2 before calling RKF45 again. otherwise,
!    execution will be terminated.
!
!    If iflag=8 is obtained, integration can not be continued unless
!    the invalid input parameters are corrected.
!
!    The arrays or scalars w1..w9 and integers iw1 to iw5 contain information
!    required for subsequent integration, and should not be altered.
!***************************************************************************}
Begin

  rkfs(neqn, y, t, tout, relerr, abserr, iflag, w1, w2, w3, w4, w5, w6,
       w7, w8, w9, iw1, iw2, iw3, iw4, iw5);

end;

{Emulation of Fortran Intrisic Functions}
Function ISign(a,b:Integer): Integer;
Begin
  if b<0 then ISign:=-ABS(a)
         else ISign:=ABS(a)
End;

Function Max(a,b:double): double;
Begin
  if a>=b then Max:=a
  else Max := b
End;

Function Min(a,b:double): double;
Begin
  if a<=b then Min:=a
  else Min := b
End;

Function Power(y,x: double): double;
Begin
  IF x<0 THEN EXIT;
  Power:=Exp(x*Ln(y))
End;

Function Sign(a,b:Double): Double;
Begin
  if b<0 then Sign:=-ABS(a)
         else Sign:=ABS(a)
End;

Procedure rkfs(neqn:integer; var y:VEC; var t:double; tout, relerr, abserr:double;
               var iflag:integer; var yp:VEC; var h:double; var f1, f2, f3, f4,f5:VEC;
               var savre, savae:Double; var nfe, kop, init, jflag, kflag:Integer);
{***************************************************************************************
!
!  RKFS implements the Runge-Kutta-Fehlberg method (double precision).
!
!  Discussion:
!
!    RKFS integrates a system of first order ordinary differential
!    equations as described in the comments for RKF45.
!
!    The arrays yp, f1, f2, f3, f4, and f5 (of dimension neqn) and
!    the variables h, savre, savae, nfe, kop, init, jflag and kflag are used
!    internally by the code and appear in the call list to eliminate
!    local retention of variables between calls.  Accordingly, they
!    should not be altered.  Items of possible interest are
!
!      YP  - the derivative of the solution vector at T;
!      H   - an appropriate stepsize to be used for the next step;
!      NFE - the number of derivative function evaluations.
!
!    The expense is controlled by restricting the number
!    of function evaluations to be approximately MAXNFE.
!    As set, this corresponds to about 500 steps.
!
!    REMIN is the minimum acceptable value of RELERR.  Attempts
!    to obtain higher accuracy with this subroutine are usually
!    very expensive and often unsuccessful.
!******************************************************************************}
Label  25,40,45,50,60,65,80,100,200,260, return;
Const  remin = 1e-12;
       maxnfe = 3000;
Var
    a,ae,dt,ee,eeoet,eps,esttol,et,hmin,rer,s,scale: double;
    tol,toln,ypk:double;
    hfaild,output:boolean;
    i,k,mflag: integer;
Begin

{  Check the input parameters } 

  eps := epsilon;

  if neqn < 1 then
  begin
    iflag := 8;
    goto return
  end;

  if relerr < 0.0 then
  begin
    iflag := 8;
    goto return
  end;

  if abserr < 0.0 then
  begin
    iflag := 8;
    goto return
  end;

  mflag := abs(iflag);

  if (abs(iflag) < 1) or (abs(iflag) > 8 ) then
  begin
    iflag := 8;
    goto return
  end;

{  Is this the first call? }

  if mflag = 1 then goto 50;

{  Check continuation possibilities }

  if (t=tout) and (kflag <> 3) then
  begin
    iflag := 8;
    goto return
  end;

  if mflag <> 2 then goto 25;

{  iflag = +2 or -2 }

  if kflag = 3 then goto 45;
  if init = 0 then goto 45;
  if kflag = 4 then goto 40;

  if (kflag = 5) and (abserr = 0.0) then Halt;

  if (kflag = 6) and (relerr <= savre) and (abserr <= savae) then Halt;

  goto 50;

{  iflag = 3,4,5,6,7 or 8 }

25:if iflag = 3 then goto 45;
   if iflag = 4 then goto 40;
   if (iflag = 5) and (abserr > 0.0) then goto 45;

{  Integration cannot be continued since user did not respond to
   the instructions pertaining to iflag:=5,6,7 or 8 }

  Halt;

{  Reset function evaluation counter }

40:nfe := 0;
  if mflag = 2 then goto 50;

{  Reset flag value from previous call }

45:iflag := jflag;

  if kflag = 3 then  mflag := abs(iflag);

{  Save input iflag and set continuation flag for subsequent input checking }

50:jflag := iflag;
   kflag := 0;

{  Save relerr and abserr for checking input on subsequent calls }

  savre := relerr;
  savae := abserr;

{  Restrict relative error tolerance to be at least as large as
   2*eps+remin to avoid limiting precision difficulties arising
   from impossible accuracy requests }

  rer := 2.0 * epsilon + remin;

{  The relative error tolerance is too small }

  if relerr < rer then
  begin
    relerr := rer;
    iflag := 3;
    kflag := 3;
    goto return
  end;

  dt := tout - t;

  if mflag = 1 then goto 60;
  if init = 0 then goto 65;
  goto 80;

{  Initialization:
     set initialization completion indicator, init
     set indicator for too many output points, kop
     evaluate initial derivatives
     set counter for function evaluations, nfe
     evaluate initial derivatives
     set counter for function evaluations, nfe
     estimate starting stepsize  }

60:init := 0;
  kop := 0;
  a := t;
  f(a, y, yp);
  nfe := 1;

  if t = tout then
  begin
    iflag := 2;
    goto return
  end;

65:init := 1;
  h := abs(dt);
  toln := 0.0;
  for k := 1 to neqn do
  begin
    tol := relerr * abs(y[k]) + abserr;
    if tol > 0.0 then
    begin
      toln := tol;
      ypk := abs(yp[k]);
      if ypk * Power(h,5.0) > tol then
        h := Power((tol/ypk),0.2)
    end
  end;

  if toln <= 0.0 then h := 0.0;

  h := Max(h, 26.0 * eps * Max(abs(t), abs(dt)));
  jflag :=  ISign(2, iflag);

{  Set stepsize for integration in the direction from T to TOUT }

80: h := Sign(h,dt);

{  Test to see if RKF45 is being severely impacted by too many output points }

  if abs(h) >= 2.0 * abs(dt) then kop := kop + 1;

{  Unnecessary frequency of output }

  if kop = 100 then
  begin
    kop := 0;
    iflag := 7;
    goto return
  end;

{  If too close to output point, extrapolate and return }

  if abs(dt) <= 26.0 * eps * abs(t) then
  begin
    for i:=1 to neqn do
      y[i] := y[i] + dt * yp[i];
    a := tout;
    f(a, y, yp);
    nfe := nfe + 1;
    t := tout;
    iflag := 2;
    goto return
  end;

{  Initialize output point indicator }

  output := false;

{  To avoid premature underflow in the error tolerance function,
   scale the error tolerances }

  scale := 2.0 / relerr;
  ae := scale * abserr;

{  Step by step integration }

100: hfaild := false;

{  Set smallest allowable stepsize }

  hmin := 26.0 * eps * abs(t);

{  Adjust stepsize if necessary to hit the output point.
   Look ahead two steps to avoid drastic changes in the stepsize and
   thus lessen the impact of output points on the code }

  dt := tout - t;
  if abs(dt) >= 2.0 * abs(h) then goto 200;

{  The next successful step will complete the integration to the output point }

  if abs(dt) <= abs(h) then
  begin
    output := true;
    h := dt;
    goto 200
  end;

  h := 0.5 * dt;    {reduce step

!     Core integrator for taking a single step
!
!     The tolerances have been scaled to avoid premature underflow in
!     computing the error tolerance function ET.
!     To avoid problems with zero crossings, relative error is measured
!     using the average of the magnitudes of the solution at the
!     beginning and end of a step.
!     The error estimate formula has been grouped to control loss of
!     significance.
!
!     To distinguish the various arguments, H is not permitted
!     to become smaller than 26 units of roundoff in T.
!     Practical limits on the change in the stepsize are enforced to
!     smooth the stepsize selection process and to avoid excessive
!     chattering on problems having discontinuities.
!     To prevent unnecessary failures, the code uses 9/10 the stepsize
!     it estimates will succeed.
!
!     After a step failure, the stepsize is not allowed to increase for
!     the next attempted step.  This makes the code more efficient on
!     problems having discontinuities and more effective in general
!     since local extrapolation is being used and extra caution seems
!     warranted.
!
!     Test number of derivative function evaluations.
!     If okay, try to advance the integration from T to T+H.

  Too many function calls }

200:if nfe > maxnfe then
  begin
    iflag := 4;
    kflag := 4;
    goto return
  end;

{  Advance an approximate solution over one step of length H }

  fehl(neqn, y, t, h, yp, f1, f2, f3, f4, f5, f1);
  nfe := nfe + 5;

{  Compute and test allowable tolerances versus local error estimates
   and remove scaling of tolerances.  Note that relative error is
   measured with respect to the average of the magnitudes of the
   solution at the beginning and end of the step }

  eeoet := 0.0;

  for k := 1 to neqn do
  begin

    et := abs(y[k]) + abs(f1[k]) + ae;

    if et <= 0.0 then
    begin
      iflag := 5;
      goto return
    end;

    ee := abs((-2090.0 * yp[k] + (21970.0E+00 * f3[k] - 15048.0 * f4[k]))
          + (22528.0 * f2[k] - 27360.0 * f5[k]));

    eeoet := Max(eeoet, ee/et)

  end;

  esttol := abs(h) * eeoet * scale / 752400.0;

  if esttol <= 1.0 then goto 260;

{  Unsuccessful step.  Reduce the stepsize, try again.
   The decrease is limited to a factor of 1/10 }

  hfaild := true;
  output := false;

  if esttol < 59049.0 then
    s := 0.9 / Power(esttol, 0.2)
  else
    s := 0.1;

  h := s * h;

  if abs(h) < hmin then
  begin
    iflag := 6;
    kflag := 6;
    goto return
  end
  else
    goto 200;

{  Successful step.  Store solution at T+H and evaluate derivatives there }

260: t := t + h;
  for i:=1 to neqn do  y[i] := f1[i];
  a := t;
  f(a, y, yp);
  nfe := nfe + 1;

{  Choose next stepsize.  The increase is limited to a factor of 5.
   If step failure has just occurred, next stepsize is not allowed to increase }

  if esttol > 0.0001889568 then
    s := 0.9 / Power(esttol, 0.2)
  else
    s := 5.0;

  if (hfaild) then  s := Min(s, 1.0);

  h := Sign(Max(s * abs(h), hmin), h);

{  End of core integrator }

{  Should we take another step? }

  if (output) then
  begin
    t := tout;
    iflag := 2
  end;

  if iflag > 0 then goto 100;

{ Integration successfully completed }

{  ne-step mode }
  iflag := - 2;

return: End;

END.

{end of file rkf45.pas}