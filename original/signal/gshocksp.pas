{**************************************************************
*   Program to demonstrate the Acceleration Shock Spectrum    *
*                                                             *
*                      Pascal version by J-P Moreau, Paris    *
*                             (with graphic options)          *
*                               (www.jpmoreau.fr)             *
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
* Output file (tshocksp.lst):                                 *
*                                                             *
*     Frequency     Shock Spectrum      Shock Spectrum        *
*       (Hz)           negative            positive           *
* -------------------------------------------------------     *
*        10.00        -3.037656           3.602233            *
*        20.00        -6.422721           9.166994            *
*        30.00       -16.325201          13.674868            *
*        40.00       -27.368246          27.100563            *
*        50.00       -41.964013          42.437765            *
*        60.00       -27.828617          26.667987            *
*        ...         ...                 ...                  *
*       740.00      -109.784192         112.831489            *
*       750.00      -109.242479         110.265918            *
*       760.00      -108.338595         108.819857            *
*       770.00      -107.086906         108.976530            *
*       780.00      -105.503821         110.103584            *
*       790.00      -103.607643         110.830962            *
*       800.00      -101.418411         109.903625            *
*                                                             *
* Note: the shock spectrum shows 3 peaks around 50, 250 and   *
*       500 hz.                                               *
**************************************************************}
Uses WinCrtMy,Type_def,Graph_2D;

Var
        f_in,f_out : TEXT;
        i,j,ndata,nfreq : INTEGER;
        df,dt,dzeta,f,fbegin,fend,temp,tbegin,tend : REAL_AR;
        mini,maxi : REAL_AR;
        signal,response,SP1,SP2 : RV;      {pointers to real vectors}


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

  {Graph of input signal F(t) on request}
  Procedure View_Input_Signal;
  Begin
    ClrScr;
    CourbeXY(CrtDC,ndata,10,Signal,tbegin,tend);
    Legendes(CrtDC,'INPUT SIGNAL','Time (s)','Acc. m/s2');
    SortieGraphique
  End;

  {Graph of calculated shock spectrum (negative and positive) on request}
  Procedure View_Shock_Spectrum;
  Var y1,y2,ymin,ymax:REAL_AR;
  Begin
    ClrScr;
    {Get minimum of negative spectrum in y1}
    MinMax(nfreq,SP1,ymin,ymax); y1:=ymin;
    {Get maximum of positive spectrum in y2}
    MinMax(nfreq,SP2,ymin,ymax); y2:=ymax;
    {Initialize virtual graphic window n° 10 with physical ranges in x and y}
    InitFenetre(CrtDC,10,0,fend,y1,y2);
    {Draw SP1 and SP2 in virtual graphic window}
    TracerXY(CrtDC,nfreq,SP1,fbegin,fend);
    TracerXY(CrtDC,nfreq,SP2,fbegin,fend);
    {Write graph title, names of axis x and y}
    Legendes(CrtDC,'  Shock Spectrum (Q=10)','Freq. (hz)','Acc. m/s2');
    {exit graph menu (print option not active here) }
    SortieGraphique;
  End;

{main program}
BEGIN
  {open main window application with title}
  WinCrtInit('Acceleration Shock Spectrum');
  {open input and output file}
  Assign(f_in,'tfft.dat'); Reset(f_in);
  Assign(f_out,'tshocksp.lst'); Rewrite(f_out);
  {read number of input signal points in input file}
  read(f_in,ndata);
  {Dynamic allocation of vectors}
  New(signal); New(response); New(SP1); New(SP2);
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

  writeln;
  write(' Data read over. Do you want to view input signal (y/n) ? ');
  if readkey='y' then View_Input_Signal;

  {input begin, end frequencies and frequency step df}
  fbegin:=10.0; fend:=800.0; df:=10.0;
  {calculate number of frequencies nfreq}
  nfreq:=1+Round((fend-fbegin)/df);
  {input damping factor dzeta}
  dzeta:=0.05;
  f:=fbegin;
  {frequency loop}
  for i:=1 to nfreq do
  begin
    {call elementary oscillator response procedure for frequency f}
    _1dof_Oscillator_Response(f,signal,response,dt,dzeta,ndata);
    {seek minimum and maximum values of response}
    mini:=response^[1]; maxi:=mini;
    for j:=2 to ndata do
    begin
      if response^[j]<mini then mini:=response^[j];
      if response^[j]>maxi then maxi:=response^[j]
    end;
    {store mini in SP1[i] and maxi in SP2[i] }
    SP1^[i]:=mini; SP2^[i]:=maxi;
    f:=f+df
  end;

  {print shock spectrum (negative and positive) to output file}
  f:=fbegin;
  writeln(f_out,'      Frequency     Shock Spectrum      Shock Spectrum   ');
  writeln(f_out,'        (Hz)           negative            positive      ');
  writeln(f_out,'  -------------------------------------------------------');
  for i:=1 to nfreq do
  begin
    writeln(f_out,'      ',f:8:2,'     ',SP1^[i]:12:6,'       ',SP2^[i]:12:6);
    f:=f+df;
  end;

  close(f_out); Clrscr;
  writeln;
  writeln(' Shock spectrum done. Numerical results in file tshocksp.lst.');
  write(' Do you want to view the shock spectrum (y/n) ? ');
  if readkey='y' then View_Shock_Spectrum;
  {free memory}
  Dispose(signal); Dispose(response); Dispose(SP1); Dispose(SP2);
  {close program main window and exit}
  DoneWinCrt

END.

{End of file gshocksp.pas}