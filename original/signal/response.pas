{**************************************************************
*   Program to demonstrate the Acceleration Shock Spectrum    *
*                                                             *
*                      Pascal version by J-P Moreau, Paris    *
* ----------------------------------------------------------- *
* SAMPLE RUN:                                                 *
*                                                             *
* Input data file (tfft.dat):                                 *
*                                                             *
* (The test signal contains 3 frequencies: 50, 250, 500 hz)   *
*                                                             *
* 1024                                                        *
* 0.00000000000000E+0000   0.00000000000000E+0000             *
* 1.95503421309917E-0004   3.53914399999776E+0001             *
* 3.91006842619834E-0004   5.95684899999760E+0001             *
* 5.86510263929974E-0004   6.54621699999552E+0001             *
* 7.82013685239669E-0004   5.24038399999845E+0001             *
* ...                      ...                                *
* 1.98826979472187E-0001   2.77372500000183E-0001             *
* 1.99022482893497E-0001  -2.43361500000174E+0000             *
* 1.99217986314807E-0001  -4.84236799999780E+0000             *
* 1.99413489736116E-0001  -6.02247899999929E+0000             *
* 1.99608993157426E-0001  -5.45615399999951E+0000             *
* 1.99804496578736E-0001  -3.22824200000105E+0000             *
* 2.00000000000045E-0001  -2.96010699999982E-0003             *
*                                                             *
* Output file (response.lst):                                 *
*                                                             *
*       Time       Mass Response                              *
*       (s)           (m/s2)                                  *
*  --------------------------------                           *
*      0.0000         0.000000                                *
*      0.0002         1.086857                                *
*      0.0004         6.069040                                *
*      0.0006        16.052865                                *
*      0.0008        29.932272                                *
*      0.0010        44.695332                                *
*      0.0012        56.432559                                *
*      0.0014        61.679698                                *
*      0.0016        58.610420                                *
*      0.0018        47.653752                                *
*      0.0020        31.323971                                *
*      0.0022        13.338342                                *
*      0.0023        -2.644980                                *
*      --/--          --/--                                   *
*      0.1990        -2.561815                                *
*      0.1992        -8.139107                                *
*      0.1994       -13.283557                                *
*      0.1996       -17.590142                                *
*      0.1998       -20.583757                                *
*      0.2000       -21.790780                                *
*                                                             *
**************************************************************}
PROGRAM Mass_Response;
Uses Wincrt,Type_def;

Var
        f_in,f_out: TEXT;
        i,j,ndata,nfreq: INTEGER;
        dt,dzeta,f,t,temp,tbegin,tend: REAL_AR;
        signal,response: RV;      {pointers to real vectors}


  {***********************************************************
  * This procedure calculates the acceleration of the sismic *
  * mass of an elementary oscillator (1 degree of freedom),  *
  * the basis of which is submitted to a given acceleration, *
  * x"(t).                                                   *
  *  The input signal x" is digitalized with a constant time *
  * step, Sampling_Incr. The table signal contains N accele- *
  * ration values.                                           *
  * -------------------------------------------------------- *
  * INPUTS:                                                  *
  *         Frequency....: eigen frequency of oscillator     *
  *         signal.......: input acceleration signal x"(t)   *
  *         Sampling_Incr: time step of signal (constant)    *
  *         Dzeta........: reduced damping factor            *
  *         N............: number of points of signal        * 
  * OUTPUT:                                                  *        
  *         x_y..........: table containing acceleration     *
  *                        response of oscillator (N points) *
  *                                                          *
  *                 Pascal version by J-P Moreau/J-P Dumont  *
  ***********************************************************}
  PROCEDURE _1dof_Oscillator_Response
  ( frequency:REAL_AR; signal: RV; VAR x_y: RV;
    Sampling_Incr,Dzeta: REAL_AR; N: INTEGER );

  VAR
  omega,Q,dQ2,delta,
  p0,p1,q1,q2,sq2,arg,cosA,
  yn,ynm1,ynm2,xn,xnm1          :REAL_AR;
  i                             :INTEGER;
  BEGIN                { Calculer la réponse }
    Omega:=2*Pi*frequency;
    IF dzeta<1E-6 THEN dzeta:=1E-6;
    Q:=1.0/(2.0*Dzeta);
    dQ2:=2*Q*Q;
    delta:=sqrt(2.0*dQ2-1.0);
    p0:=omega*Sampling_Incr/Q;
    q2:=Exp(-p0);
    sq2:=Sqrt(q2);
    arg:=0.5*p0*Delta;
    cosA:=cos(arg);
    q1:=-2.0*sq2*cosA;
    p1:=p0*sq2*((dQ2-1.0)*sin(arg)/delta-cosA);
    ynm1:=0.0;
    xn:=signal^[1];
    yn:=0.0;
    x_y^[1]:=yn;
    FOR i:=2 TO N DO
    BEGIN
      ynm2:=ynm1;ynm1:=yn;xnm1:=xn;
      xn:=signal^[i];
      yn:=p0*xn+p1*xnm1-q1*ynm1-q2*ynm2;
      x_y^[i]:=yn;
    END;
  END;  {_1dof_Oscillator_Response}


{main program}
BEGIN
  {open input and output file}
  Assign(f_in,'tfft.dat'); Reset(f_in);
  Assign(f_out,'response.lst'); Rewrite(f_out);
  {read number of input signal points in input file}
  read(f_in,ndata);
  {Dynamic allocation of vectors}
  New(signal); New(response);
  {read ndata couples T(i), Y(i) in input data file}
  for i:=1 to ndata do
  begin
    readln(f_in,temp,signal^[i]);
    if i=1 then tbegin:=temp;
    if i=ndata then tend:=temp
  end;
  close(f_in);
  {calculate sampling increment dt of input signal}
  dt := (tend-tbegin)/(ndata-1);
  {input damping factor dzeta}
  dzeta:=0.05;
  f:=250.0;

  {call elementary oscillator response procedure for frequency f}
  _1dof_Oscillator_Response(f,signal,response,dt,dzeta,ndata);
  
  {print response(i) to output file}
  t:=tbegin;
  writeln(f_out);
  writeln(f_out,'       Time       Mass Response   ');
  writeln(f_out,'       (s)           (m/s2)       ');
  writeln(f_out,'  --------------------------------');
  for i:=1 to ndata do
  begin
    writeln(f_out,'      ',t:6:4,'     ',response^[i]:12:6);
    t:=t+dt
  end;

  close(f_out);
  writeln;
  writeln(' Results in file response.lst.');
  Readkey;
  {free memory}
  Dispose(signal); Dispose(response);
  {close program main window and exit}
  DoneWinCrt

END.

{End of file response.pas}