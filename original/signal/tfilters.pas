{**************************************************************
*        Program to demonstrate the Butterworth filter        *
*           (removing frequencies greater then Fc)            *
*                                                             *
*                      Pascal Version By J-P Moreau, Paris    *
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
* Output file (tfilter.lst):                                  *
*                                                             *
*        Time        Input Signal      Filtered Signal        *
*  -------------------------------------------------------    *
*      0.0000         0.000000           0.000000             *
*      0.0002        35.391440           0.158005             *
*      0.0004        59.568490           1.278083             *
*      0.0006        65.462170           5.003994             *
*      0.0008        52.403840          12.809453             *
*      0.0010        26.229620          24.336649             *
*      ...           ...                ...                   *
*      0.1988         0.277373          -2.224546             *
*      0.1990        -2.433615          -0.802179             *
*      0.1992        -4.842368           0.588941             *
*      0.1994        -6.022479           1.320860             *
*      0.1996        -5.456154           1.036090             *
*      0.1998        -3.228242          -0.213007             *
*      0.2000        -0.002960          -1.978450             *
*                                                             *
* Note: The cut off frequency Fc is 500 hz, the filter order  *
*       is 4 (A Fourier analysis of the filtered signal shows *
*       that frequencies greater than Fc are satisfactorily   *
*       removed but with a slight time shift of the signal).  *
**************************************************************}
PROGRAM Butterworth_Filter;
Uses Wincrt,Type_def;

Const  MACH_EPS = 1e-12;

Type

    Filter_Coef = ARRAY[1..5,1..10] OF REAL_AR;
    Memory_Coef = ARRAY[1..2,1..10] OF REAL_AR;

Var
        f_in,f_out : TEXT;
        i,order,ndata : INTEGER;
        dt,Fc,t,temp,tbegin,tend : REAL_AR;
        signal,filtered : RV;      {pointers to real vectors}

        C        : Filter_Coef;
        D        : Memory_Coef;
        NSections: INTEGER;
        Tg,Xdc   : REAL_AR;

{Guard against divide by zero}
Function Tan(x: real_ar) : real_ar;
VAR Cx : real_ar;
begin
  Cx:=Cos(x);
  if (ABS(Cx)<MACH_EPS) then
    Tan := Sin(x)/MACH_EPS
  else
    Tan:=Sin(x)/Cx
End;


{**********************************************************************
*          Filtering a signal F(t) by Butterworth method              *
*             (removing frequencies greater then Fc)                  *
* ------------------------------------------------------------------- *  
* Calling mode:   Filter(Xs,Xd,NSections,C,D);                        *
* ------------------------------------------------------------------- *
* INPUTS:                                                             *
* -------                                                             *
*        Xd.......:  current value of input signal (double)           *
*        NSections: Number of required 2nd order sections (integer)   *
*                   = n/2     for n even                              *
*                   = (n+1)/2 for n odd                               *
*        n........: order of filter (1 to 4)                          *
*        C........: Table[1..5,1..NSections] of filter coefficients   *
*                   calculated previously by BUTTERWORTH procedure    *
*        D........: Table[1..2,1..NSections] of coefficients defining *
*                   the filter memory, initialized by INIT procedure. *
* ------------------------------------------------------------------- *
* OUTPUTS:                                                            *
* -------                                                             *
*        D........: Table updated after the call to Filter procedure  *
*        Xs.......: current value of filtered signal (double)         *
* ------------------------------------------------------------------- *
* Référence                                                           *
* ---------                                                           *
*  "Lawrence R.Rabiner et Bernard Gold                                *
*   Theory and application of digital processing.                     *
*   Prentice Hall Inc., EnglewoodclIFfs,NEW JERSEY,1975. [BIBLI 15]"  *
*                                                                     *
*                                Pascal Version By J-P Moreau, Paris  *
*               from Fortran Version By J-P Dumont / Tuan Dang Trong  *
**********************************************************************}
  PROCEDURE Filter ( VAR Xs    : real_ar;
                     Xd        : real_ar;
                     NSections : INTEGER;
                     C         : Filter_Coef;
                     VAR D     : Memory_Coef );
  VAR x,y,err : REAL_AR;
  i       : INTEGER;
  BEGIN
    x:=Xd;
    FOR i:=1 TO NSections DO
    BEGIN
      err:=x+C[1][i]*D[1][i]+C[2][i]*D[2][i];
      y:=C[5][i]*(Err +C[3][i]*D[1][i]+C[4][i]*D[2][i]);
      D[2][i]:=D[1][i];
      D[1][i]:=Err;
      x:=y;
    END;
    Xs:=x;
  END;

{*************************************************************************
*                       INIT FILTER PROCEDURE                            *
* ---------------------------------------------------------------------- *
* The filter response is initialized to stationnary value for a constant *
* input signal value.                                                    *
*                                                                        *
* Calling mode:   INIT(Xdc,C,NSections,D);                               *
* ---------------------------------------------------------------------- *
* INPUTS:                                                                *
* ------                                                                 *
*        Xdc......: constant input value (double)                        *
*        C........: Table[1..5,1..NSections] of filter coefficients      *
*                   calculated previously by BUTTERWORTH procedure       *
*        NSections: Number of required 2nd order sections (integer)      *
*                   = n/2     for n even                                 *
*                   = (n+1)/2 for n odd                                  *
*        n........: order of filter (1 to 4)                             *
* ---------------------------------------------------------------------- *
* OUTPUTS:                                                               *
* -------                                                                *
*        D........: Table[1..2,1..NSections] of coefficients defining    *
*                   the filter memory, initialized by INIT procedure.    *
*************************************************************************}
  PROCEDURE Init( Xdc: REAL_AR;
                  C: Filter_Coef;
                  NSections : INTEGER;
                  VAR D: Memory_Coef );
  VAR   dc,Csum : REAL_AR;
  i,j     : INTEGER;
  BEGIN
    dc:=Xdc;
    FOR j:=1 TO NSections DO
    BEGIN
      D[2][j]:=dc/(1-C[1][j]-C[2][j]);
      D[1][j]:=D[2][j];
      Csum:=0;
      FOR i:=1 TO 4 DO Csum:=Csum+ C[i][j];
      dc:=C[5][j]*(dc+D[2][j]*Csum);
    END; {j}
  END; {Init}
  

{*********************************************************************
*          Calculates the Butterworth filter coefficients            *
* ------------------------------------------------------------------ *  
*  Calling mode:   Butterworth(Fc,Ts,n,C,NSections,Tg);              *
* ------------------------------------------------------------------ *
*  INPUTS:                                                           *
*  ------                                                            *
*         Fc.......: Cut off frequency                               *
*         Ts.......: Sampling time of input signal                   *
*         n........: Order of filter                                 *
* ------------------------------------------------------------------ *
*  OUTPUTS:                                                          *
*  -------                                                           *
*         C........: Table[1..5,1..NSections] of filter coefficients *
*                    calculated previously by BUTTERWORTH procedure  *
*         NSections: Number of required 2nd order sections (integer) *
*                    = n/2     for n even                            *
*                    = (n+1)/2 for n odd                             *
*         Tg.......: Group delay in seconds                          *
*********************************************************************}
  PROCEDURE Butterworth( Fc,Ts : REAL_AR;
                         n     : INTEGER;
                         VAR C        : Filter_Coef;
                         VAR NSections: INTEGER;
                         VAR Tg       : REAL_AR );

  CONST Zero = 0;
  ONE  = 1;
  TWO  = 2;
  HALF = 0.5;
  
  VAR   Ns2,i,ModN                   : INTEGER;
  Arg,Rep,Omega,OmegaSq,Temp,W0,W1,m : REAL_AR;
  
  BEGIN
    Arg:=Pi*Ts*Fc;
    If abs(arg) > 2.0*Pi THEN
    BEGIN
      m:=Int(arg/2.0/Pi);
      arg:=arg-(m*2.0*pi)
    END;
    Omega:= Tan(Arg);   { Tan est défini dans Math.Pas }
    OmegaSq:=Omega*Omega;
    Modn:=(n MOD 2);
    If (Modn =0) THEN temp:=HALF
    ELSE temp:=Zero;
    Ns2:=(n DIV 2 );
    NSections:=Ns2+Modn;
    Tg:=zero;
    If (n>1) THEN FOR i:=1 TO Ns2 DO
    BEGIN
      Rep:=Omega*Cos(Pi*(i-temp)/n);
      Tg:=Tg+Ts*Rep/OmegaSq;
      W0:=TWO*Rep;
      W1:=ONE +W0+OmegaSq;
      C[1][i]:=-TWO*(OmegaSq-ONE)/W1;
      C[2][i]:=-(ONE-W0+OmegaSq)/W1;
      C[3][i]:=TWO;
      C[4][i]:=ONE;
      C[5][i]:=OmegaSq/W1;
    END;
    If (temp = Zero) THEN
    BEGIN
      C[1][Nsections]:=(ONE-Omega)/(ONE+Omega);
      C[2][NSections]:= zero;
      C[3][NSections]:= ONE;
      C[4][NSections]:= zero;
      C[5][NSections]:= Omega/(ONE+Omega);
      Tg:= Tg+Ts/(TWO*Omega);
    END
  END; {Butterworth}


{main program}
BEGIN
  {open input and output file}
  Assign(f_in,'tfft.dat'); Reset(f_in);
  Assign(f_out,'tfilter.lst'); Rewrite(f_out);
  {read number of input signal points in input file}
  read(f_in,ndata);
  {Dynamic allocation of vectors}
  New(signal); New(filtered);
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
  {input cut off frequencys and order of filter (1 to 4) }
  Fc:=500.0; order:=4; Xdc:=0; {no constant value}

  {call filtering procedures}
  {1. Calculate the filter coefficients}
  Butterworth(Fc,dt,order,C,NSections,Tg);

  {2. Initialize filter memory}
  Init(Xdc,C,NSections,D);

  {3. Recursively call Butterworth filter}
  FOR i:=1 TO ndata DO
  begin
    Filter(temp,signal^[i],NSections,C,D);
    filtered^[i]:=temp
  end;

  {print input and filtered signals to output file}
  t:=tbegin;
  writeln(f_out,'        Time        Input Signal      Filtered Signal   ');
  writeln(f_out,'  -------------------------------------------------------');
  for i:=1 to ndata do
  begin
    writeln(f_out,'      ',t:6:4,'     ',signal^[i]:12:6,'       ',filtered^[i]:12:6);
    t:=t+dt;
  end;

  close(f_out);
  writeln;
  writeln(' Results in file tfilter.lst.');
  Readkey;
  {free memory}
  Dispose(signal); Dispose(filtered);
  {close program main window and exit}
  DoneWinCrt

END.

{End of file tfilters.pas}