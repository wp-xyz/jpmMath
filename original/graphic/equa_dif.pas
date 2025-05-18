{****************************************************
*             UNIT  Equa_dif.pas                    *
* ------------------------------------------------- *
*  Integration procedure odeint1 used by program    *
*  planets.pas                                      *
* ------------------------------------------------- *
*             Author: Jean-Pierre Dumont, France    *
*                      (www.jpmoreau.fr)            *
****************************************************}
Unit Equa_dif;

INTERFACE
Uses WinCrtMy, Type_def;

Const Nvar_max = 32;             {Cf. Planets.pas}
      Ncorps = 8;
      DtSave = 0.02;

TYPE Vecteur_etat = ARRAY [1..Nvar_max] OF REAL_AR;
     Pointers     = ARRAY [0..Nvar_max] OF INTEGER;
     Colonne      = ARRAY [1..Nvar_max] OF REAL_AR;

var ipt : Pointers;
    Masse,q1d,q2d,q1,q2 : Colonne;

PROCEDURE odeint1(VAR ystart: Vecteur_etat;
                      nvar  : INTEGER;
                      t1,t2,eps,h1,hmin: REAL_AR;
		  VAR nok,nbad: INTEGER);

IMPLEMENTATION


PROCEDURE Derivs (t: REAL_AR;
             VAR XT: vecteur_etat;
	     VAR XF: vecteur_etat);
{=======================================================================
            VERSION 1.0 - 12/27/1993
-----------------------------------------------------------------------}
VAR xi,xdi,xj,
    yi,ydi,yj,
    dxij,dyij,rij,
    xsm,ysm             : REAL_AR;
    i,ipti,iptj,j       : INTEGER;

BEGIN
FOR i := 1 TO Ncorps DO
  BEGIN
{Calculate intermediate scalar values}
    ipti  := ipt[i];
    xdi   := XT[ipti+1];
    ydi   := XT[ipti+2];
    xi    := XT[ipti+3];
    yi    := XT[ipti+4];
{Calculate second members}
    XF[ipti+3] := xdi;
    XF[ipti+4] := ydi;
    xsm := 0;
    ysm := 0;
    IF (i <> 1) THEN FOR j := 1 TO i-1 DO
      BEGIN
        iptj := ipt[j];
        xj   := XT[iptj+3];
        yj   := XT[iptj+4];
        dxij := xj-xi;
        dyij := yj-yi;
        rij  := sqr(dxij)+sqr(dyij);
        rij  := Masse[j]/(rij*sqrt(rij));
        xsm  := xsm + dxij*rij;
        ysm  := ysm + dyij*rij;
      END; {IF}
    IF (i <> Ncorps) THEN FOR j := i+1 TO Ncorps DO
      BEGIN
        iptj := ipt[j];
        xj   := XT[iptj+3];
        yj   := XT[iptj+4];
        dxij := xj-xi;
        dyij := yj-yi;
        rij  := sqr(dxij)+sqr(dyij);
        rij  := Masse[j]/(rij*sqrt(rij));
        xsm  := xsm + dxij*rij;
        ysm  := ysm + dyij*rij;
      END; {IF}
  XF[ipti+1] := xsm;
  XF[ipti+2] := ysm;
  END; {i}
END;   {Derivs}

{================================================================}
PROCEDURE RK4(VAR y,dydt : Vecteur_etat;
                        n: INTEGER;
                      t,h: REAL_AR;
                VAR yout : Vecteur_etat);
{=================================================================
       Version 1.1 - 12/31/1991
==================================================================
| This procedure integrates the 1st order differential system:
|      dy/dt = F(y,t)
| by a fourth order Runge-Kutta method to advance the solution
| on interval h of the independant variable t:
 ------------------------------------------------------------------
| INPUTS:
| y     = State vector at begin of integration step
| dydt  = its derivative at the same point.
| n     = number of equations of system
| t     = i.s. value at begin of time step
| h     = time step
 -------------------------------------------------------------------
| OUTPUT:
| yout  = State vector at end of integration step.
 -------------------------------------------------------------------
| Programs using procedure RK4 must provide access to:
|
| PROCEDURE Derivs(t:REAL_AR;y:Vecteur_etat;VAR dydt:Vecteur_etat);
| which returns the values of derivatives dydt at point t, knowing
| both t and the values of functions y.
|
/ Also define the type:
|
| TYPE Vecteur_etat = ARRAY [1..Nvar_max] OF REAL_AR;
| in calling program, where Nvar_max is the maximal number
| of equations of the system:  n <= Nvar_max.
 --------------------------------------------------------------------}
  VAR i          : INTEGER;
      th,hh,h6   : REAL_AR;
      dym,dyt,yt : Vecteur_etat;

  BEGIN
       hh := 0.5* h;
       h6 :=h/6;
       th := t+ hh;
       FOR i := 1 TO n DO yt[i] := y[i] + hh * dydt[i];
       derivs(th,yt,dyt);
       FOR i := 1 TO n DO yt[i] := y[i] + hh * dyt[i];
       derivs(th,yt,dym);
       FOR i := 1 TO n DO BEGIN
           yt[i] := y[i]+ h*dym[i];
           dym[i]:= dyt[i] + dym[i];
           END;
       derivs(t+h,yt,dyt);
       FOR i := 1 TO n DO yout[i] := y[i]+h6*(dydt[i]+dyt[i]+2.0*dym[i])
  END; { RK4 }

{//////////////////////////////////////////////////////////////////////}
PROCEDURE rkqc(VAR y,dydt: Vecteur_etat;
                        n: INTEGER;
                    VAR t: REAL_AR;
                 htry,eps: REAL_AR;
                VAR yscal: Vecteur_etat;
           VAR hdid,hnext: REAL_AR);
{=======================================================================
                     Version 1.2 - 03/26/1993
========================================================================
| Runge-Kutta integration step with control of truncation local error 
| to obtain a required precision and adjust time step consequently
----------------------------------------------------------------------
| INPUTS:
| y      = State vector of size n
| dydt   = its derivative at begin value of independant variable, t
| n      = number of equations of system
| t      = begin value of independant variable
| htry   = time step proposed as a try
| eps    = precision requirement:
|          Max (ycalc[i] - yvrai[i])/yscal[i] < eps
| yscal  = normalization vector of solution.
----------------------------------------------------------------------
| OUTPUTS:
|  y     = end state vector
|  t     = i.s. end value
|  hdid  = actual time step  
|  hnext = time step advised for the next integration step
----------------------------------------------------------------------
| Programs using procedure RKQC must provide access to:
|
| PROCEDURE Derivs(t:REAL_AR;y:Vecteur_etat;VAR dydt:Vecteur_etat);
| which returns the values of derivatives dydt at point t, knowing
| both t and the values of functions y.
|
/ Also define the type:
|
| TYPE Vecteur_etat = ARRAY [1..Nvar_max] OF REAL_AR;
| in calling program, where Nvar_max is the maximal number
| of equations of the system:  n <= Nvar_max.
=====================================================================}
LABEL 1;

CONST
    pgrow=-0.20;
    pshrnk=-0.25;
    fcor=0.06666666;   {1/15 }
    un = 1.0;
    safety=0.9;
    errcon=6E-4;
    tiny= 1E-20;    

VAR
    i,kount              : INTEGER;
    tsav,hh,h,temp,errmax: REAL_AR;
    dysav,ysav,ytemp     : Vecteur_etat;

BEGIN
         tsav:= t;          {Save begin time}
         FOR i:=1 TO n DO BEGIN
             ysav[i] := y[i];
             dysav[i]:= dydt[i];
             END;
         h:= htry;        {define increment for a try value}
    1:   hh := 0.5*h;     {take 2 half time steps}
         rk4(ysav,dysav,n,tsav,hh,ytemp);
         t:= tsav + hh;
         derivs(t,ytemp,dydt);
         rk4(ytemp,dydt,n,t,hh,y);
         t:= tsav + h;
         IF (t = tsav) THEN BEGIN
            writeln('Pause in RKQC procedure');
            writeln('Increment too small of independant variable');
            writeLn('Press any key to continue...');
            Readkey
            END;
         rk4(ysav,dysav,n,tsav,h,ytemp);
         errmax := 0;   {Evaluate error}
         temp :=0;
         FOR i := 1 TO n DO BEGIN
             ytemp[i] := y[i] - ytemp[i];   {ytemp = estimated error}
             IF (yscal[i]>tiny) THEN temp := abs(ytemp[i]/yscal[i]);
             IF ( errmax < temp) THEN errmax := temp;
         END; { i }
         errmax:= errmax/eps;    {real error / requirement}
         IF (errmax > un) THEN    {Error too big, reduce h}
            BEGIN
            h := safety*h*exp(pshrnk*ln(errmax));
            GOTO 1; {start again}
            END
            ELSE BEGIN       {the step has been a success!}
                 hdid := h;  {Calculate next time step}
                 IF (errmax > errcon) THEN
                    hnext:=safety*h*exp(pgrow*ln(errmax))
                    ELSE hnext:= 4.0*h;
                 END;
         FOR i := 1 TO n DO y[i]:=y[i]+ytemp[i]*fcor;
END;  {rkqc}

{//////////////////////////////////////////////////////////////////////}
PROCEDURE odeint1(VAR ystart: Vecteur_etat;
                      nvar  : INTEGER;
                      t1,t2,eps,h1,hmin: REAL_AR;
                  VAR nok,nbad: INTEGER);
{================================================================
                  Version 1.2 - 12/31/1991
=================================================================
|  This procedure integrates the 1st order differential system:
|       dy/dt = F(y,t)
|  where y and F are vectors of size nvar, between times
|  t1 and t2.
|
|  INPUTS:
|  ystart= begin coordinates vector and speeds
|  nvar  = number of equations
|  t1    = begin integration time
|  t2    = end integration time
|          t2 may be > or < t1
|  eps   = absolute required precision for solution
|  h1    = time increment proposed at beginning (try value)
|  hmin  = minimum time increment
|
|  OUTPUTS:
|  ystart= end coordinates vector and speeds
|  nok   = number of unchanged time steps
|  nbad  = number of modified time steps
====================================================================
| Programs using procedure RKQC must provide access to:
|
| PROCEDURE Derivs(t:REAL_AR;y:Vecteur_etat;VAR dydt:Vecteur_etat);
| which returns the values of derivatives dydt at point t, knowing
| both t and the values of functions y.
|
/ Also define the type:
|
| TYPE Vecteur_etat = ARRAY [1..Nvar_max] OF REAL_AR;
| in calling program, where Nvar_max is the maximal number
| of equations of the system:  n <= Nvar_max.
|
| Also define:        VAR  dtsave: REAL_AR;
|
| Initialize:         dtsave
/
====================================================================}

LABEL 99;

CONST maxstp= 10000;
      two   =    2 ;
      zero  =    0 ;
      tiny  = 1E-20;  

VAR   nstp,i : INTEGER;
      tsav,t,hnext,hdid,h  : REAL_AR;
      yscal,y,dydt         : vecteur_etat;

BEGIN
     t:= t1;
     IF ( t2 > t1) THEN h :=   abs(h1)
                   ELSE h := - abs(h1);
     nok:= 0;
     nbad :=0;
     FOR i := 1 TO nvar DO y[i] := ystart[i];
     tsav := t - dtsave*two;
     FOR nstp :=1 TO maxstp DO
         BEGIN
         derivs(t,y,dydt);
         FOR i:=1 TO nvar DO yscal[i]:=abs(y[i])+abs(dydt[i]*h); {+tiny;}
         IF (((t+h-t2)*(t+h-t1)) > zero) THEN h := t2 - t;
         rkqc(y,dydt,nvar,t,h,eps,yscal,hdid,hnext);
         IF (hdid = h ) THEN nok  := nok +1
                        ELSE nbad := nbad +1;
         IF (((t-t2)*(t2-t1)) >= zero) THEN
            BEGIN
            FOR i := 1 TO nvar DO ystart[i] := y[i];
            GOTO 99; { c'est fini }
            END;
         IF (abs(hnext) < hmin) THEN
            BEGIN
            writeLn(' Time step too small!');
            nok :=-1 ; {error flag}
            Readkey;
            GOTO 99;
            END;
         h := hnext;
         END; {nstp}
	 writeLn('Pause in ODEINT1 procedure- too many time steps!');
         Readkey;
99  : END; { odeint1 }
{//////////////////////////////////////////////////////////////////////}
END.

{end of file equa_dif.pas}