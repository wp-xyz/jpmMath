{************************************************************************
  *     PLANETS.PAS        Version 1.2 (TPW et VGA) 01/12/1994            *
  * --------------------------------------------------------------------- *
  * This program calculates the motions of N bodies submitted to their    *
  * gravity attractions. The integration of motions is made by a Runge-   *
  * Kutta method with controled time step to satisfy a required precision.*
  *                                                                       *
  *                                        TPW version By J-P Moreau      *
  *                                   (from a DOS program By J-P Dumont)  *
  *                                            (www.jpmoreau.fr)          *
  * --------------------------------------------------------------------- *
  *                                         Print option not implemented. *
  ************************************************************************}
  PROGRAM Planets;
  USES
          WinCrtMy,WinTypes,WObjects,WinProcs,Strings,         {   TPW    }
          Type_Def, CrtGr2D,                                   {J-P Moreau}
          Equa_dif;                                            {J-P Dumont}

  LABEL fin;

  CONST

  NcorpsMx   =       8;                          {Maximal number of bodies}
  Nddl       =       4;             {Number of degrees of freedom per body}
  NsmbMx     =      32;                  {Maximal number of second members}
  stab       =    0.05;

  VAR

  Cx,Cy,Echl,Xratio,Yratio,dt,dtmIN,eps,h1,Tmax,time,t2: REAL_AR;

  ystart                   : Vecteur_etat;              {see Equa_dif.pas}
  XRayon,YRayon            : REAL_AR;
  NstepMax                 : LONGINT;
  iStep,Nsmb,nstep,nok,nbad: INTEGER;
  Done                     : BOOLEAN;

  {//////////////////////////////////////////////////////////////////////}
  PROCEDURE initpointers ( VAR ipt : Pointers);
  {=======================================================================
  Calculate pointers to store state variables in Table
  -----------------------------------------------------------------------}
  VAR i  : INTEGER;

  BEGIN
    ipt[0]  := - Nddl;
    FOR i:= 1 TO Ncorps DO
    BEGIN
      ipt[i]  := ipt[i-1]  + Nddl;
    END;
  END;

  {----------------------------------------------------------------------}
  PROCEDURE ReadData ( VAR Masse,q1,q2,q1d,q2d : Colonne;
                       VAR Tmax,dtsave : REAL_AR;
                       VAR NCorps      : INTEGER);
  {=======================================================================
  VERSION 1.0 - 12/27/1993
  -----------------------------------------------------------------------}
  VAR i    : INTEGER;
  sign : REAL_AR;

  BEGIN

    Tmax     := 75;
    Dtsave   := 0.002;
    Ncorps   := 8;        {number of planets (maxi 8) }
    sign     := 1.0;

    Masse[1] := 10000;
    q1[1]    := 0;
    q2[1]    := 0;
    q1d[1]   := 0;
    q2d[1]   := 0;

    FOR i := 2 TO NCorps DO
    BEGIN
      sign     := - sign;
      Masse[i] := 1.0*(9.0 - i);
      q1[i]    := 1.0*(i-1.0);
      q2[i]    := 0;
      { Factor 1.0 circular trajectory}
      { Factor 1.414 parabolic trajectory}
      q2d[i]   := sign*1.38*sqrt(Masse[1]/q1[i]);
      q1d[i]   := 0;
      {centrifugal force = gravity force}
      q1[i]    := q1[1]+sign*q1[i];
    END;
  END; {ReadData}

  {----------------------------------------------------------------------}
  PROCEDURE InitialConditions
  (VAR XN : Vecteur_etat; q1,q2,q1d,q2d: Colonne);

  VAR i,ipti : INTEGER;

  BEGIN
    FOR i := 1 TO Ncorps DO
    BEGIN
      ipti := ipt[i];
      XN[ipti+1]:= q1d[i];
      XN[ipti+2]:= q2d[i];
      XN[ipti+3]:= q1[i];
      XN[ipti+4]:= q2[i];
    END;
  END;  {InitialConditions}

  {------------------------------------------------------------------------}
  PROCEDURE Display( XN : Vecteur_Etat);

  VAR
       xw,yw  : REAL_AR;
       i,ipti : INTEGER;
  BEGIN
    FOR i := 1 TO Ncorps DO
    BEGIN
      ipti := ipt[i];
      xw := XN[ipti+3];  {xi}
      yw := XN[ipti+4];  {yi}
      {Draw central star at first time}
      IF i=1 THEN
      BEGIN
        IF NOT Done THEN
        BEGIN
          Circle1(CrtDc,0,0,Xrayon,TRUE);
          Done:=TRUE
        END
      END
      {draw planets around central star}
      ELSE
      	Circle1(CrtDc,xw,yw,Yrayon,TRUE)
    END; {i loop}
  END;


  {----------------------------------------------------------------------}
  BEGIN   {main program}

    WinCrtInit('EIGHT PLANETS');
    Done:=FALSE;

    {define graph window in physical coordinates}
    Fenetre(-40,40,-40,40);
    {define graph window in pixels}
    Cloture(10,MaxX-25,80,MaxY-10);
    {draw a rectangle around graph zonz} 
    Bordure(CrtDC);

    ReadData(Masse,q1,q2,q1d,q2d,Tmax,dtsave,Ncorps);

    Nsmb := Nddl*Ncorps;

    NstepMax := Trunc(TMax/dtSave) + 1;

    initpointers (ipt);

    {Define useful constants}
    XRayon:= 1.0;                                {Radius of central star}
    YRayon:= 0.5;                                     {Radius of planets}

    eps   := 1.0E-12;                       {Required relative precision}
    h1    := 0.01*dtSave;
    dt    := h1;
    dtmIN := 0.001*h1;                      {Smallest time step allowed}
    time  := 0;
    iStep := 1;

    InitialConditions(ystart,q1,q2,q1d,q2d);

    {main integration loop}
    Repeat
      WHILE (time <= Tmax) AND (istep < nstepMax) AND (NOT KeyPressed) DO
      BEGIN
        INC(iStep);
        t2:=time+dtsave;
        odeint1(ystart,nsmb,time,t2,eps,dt,dtmIN,nok,nbad);
        IF (nok < 0 ) OR (KeyPressed) THEN
        BEGIN
          DEC(iStep);
       	  GOTO fin;          {no storage}
        END;
        time:=t2;
        {display planets once every four times}
        IF (Round(t2*100) MOD 4)=0 THEN Display (ystart)
      END;
      fin:  nstep:=iStep;        {actual number of time steps}
      SortieGraphique;
    Until rep='n';
    DoneWinCrt
  END.
  {---------------------------  E N D  ---------------------------------}

{end of file planets.pas}