{****************************************************************
* NUMERICAL DECONVOLUTION OF A SUSPENDED CAPTOR'S MASS ACCELE-  *
* RATION RESPONSE SIGNAL  TO OBTAIN SPEED AT BASIS AND COMPUTE  *
* THE ACCELERATION SHOCK SPECTRUM OF THIS SPEED WITH OR WITHOUT *
* SIGNAL FILTERING OF HIGH FREQUENCIES                          *
* ------------------------------------------------------------- *
* The shock captor is considered as a damped elementary oscil-  *
* lator of one degree of freedom (M,K) with a resonance fre-    *
* quency F0 and a relative damping factor dzeta. Q is the       *
* quality factor used to compute the shock speCtrum.            *
*                                                               *
*                        Pascal version by J-P Moreau, Paris    *
*                      (in double precision with graph options) *
*                                 (www.jpmoreau.fr)             *
*                                                               *
* SAMPLE RUN:                                                   *
*                                                               *
* Read data from file signal.dat                                *
* Write results to file vsignal.lst.                            *
****************************************************************}  
PROGRAM Deconvolution;
Uses WinCrtMy,Type_def,Graph_2D;

Const
       ZERO = 0.0;
       ONE  = 1.0;
       TWO  = 2.0;
Var
       ACC,VIT,S : RV;
       i,ndata,nspec : INTEGER;
       fp_in,fp_out : TEXT;
       df,dzeta,f0,finf,fmax,fr,fsup,q,temp,ts:REAL_AR;
       tbegin,tend : REAL_AR;
       Title : STRING[80];
       input,output : STRING[12];
       nom : STRING[8];


   {------------------------------------------------
    Read    n : number of points of input signal
            t : sampling duration of signal in sec.
          A(i): Table with ndata acceleration values
    ------------------------------------------------}
    PROCEDURE READATA(VAR A:RV;VAR T:REAL_AR; VAR N:Integer);
    Var i:Integer;
    Begin                                                 
      READLN(fp_in,N);                                                             
      WRITELN(fp_out,' ',N,' POINTS READ.');
      READLN(fp_in,T);                                                              
      WRITELN(fp_out,' SAMPLING DURATION: ',T,' S');
      For I:=1 to N do
      begin
        READLN(fp_in,temp,A^[I]);
        if I=1 then tbegin:=temp;
        if I=N then tend:=temp;
      end
    End;

    {--------------------------------------------------
      DECONVOLUTION SUBROUTINE

      Inputs:
             frequency  resonance frequency of captor
             signal     input mass acceleration
         Sampling_Incr  sampling time of signal
             dzeta      relative damping factor
             Ndata      number of points of ACC(I)
      Output:
             x_y        speed of basis by deconvolution
    ---------------------------------------------------}
    PROCEDURE DECON(frequency: REAL_AR; signal: RV; VAR x_y  : RV;
                      Sampling_Incr,Dzeta: REAL_AR; Ndata: INTEGER );
    VAR
      Omega,Omega2,Q,OmegaQ,OmegaT,OmegaQT,
      q0,q1,q2,T,T2sur6,expA,yn,ypn,ynm1,
      xn,xpn,xnm1,xpnm1,xppn,xppnm1 : REAL_AR;
      i: INTEGER;
    BEGIN                     {Calculate speed of basis}
      Omega:=2*Pi*frequency;
      Omega2:=Omega*Omega;
      Q:=1/(2*Dzeta);
      T:=Sampling_Incr;
      T2sur6 :=T*T/6;
      OmegaT :=Omega*T;
      OmegaQ :=Omega*Q;
      OmegaQT:=Omega*Q*T;
      expA   :=exp(-omegaQT);
      q2:=OmegaQ*expA;
      q0:=(1-expA)/(Omega*OmegaT)+T/2;
      q1:=T+Q/Omega*expA-q0;
      yn:=0;xn:=0;
      xpn:=0;xppn:=Signal^[1];
      x_y^[1]:=xpn;
      FOR i:=2 TO Ndata DO
      BEGIN
        ynm1:=yn;xnm1:=xn;xpnm1:=xpn;xppnm1:=xppn;
        xppn:=signal^[i];
        ypn:=xpnm1+q0*xppn+q1*xppnm1+q2*(xnm1-ynm1);
        xpn:=xpnm1+0.5*T*(xppn+xppnm1);
        xn:=xnm1+T*xpnm1+T2sur6*(2*xppnm1+xppn);
        yn:=xn+xppn/Omega2+(xpn-ypn)/OmegaQ;
        x_y^[i]:=ypn
      END
    END;

    {************************************************************
    * This procedure calculates the acceleration of the sismic  *
    * mass of an elementary oscillator (1 degree of freedom),   *
    * the basis of which is submitted to a given speed, VIT(t). *
    * The input signal VIT is digitalized with a constant time  *
    * step, T. The table VIT(I) contains NDATA speed values.    *
    * --------------------------------------------------------- *
    * INPUTS:                                                   *
    *         FREQU........: eigen frequency of oscillator      *
    *         VIT..........: input speed signal of basis VIT(t) *
    *         T............: time step of signal (constant)     *
    *         DZETA........: reduced damping factor             *
    *         NDATA........: number of points of signal         *
    * OUTPUT:                                                   *
    *         ACC..........: table containing acceleration res- *
    *                        ponse of oscillator (NDATA points) *
    *                                                           *
    *                      Pascal version by J-P Moreau, Paris  *
    ************************************************************}
    PROCEDURE OSCIL(FREQU:REAL_AR;VIT:RV;VAR ACC:RV;T,DZETA:REAL_AR;NDATA:INTEGER);
    Var arg,cosarg,delta,expa,gomega,omega,omega2,sinarg : REAL_AR;
        q0,q1,q2,qsi,r0,r1,r2,r3,xpn,xpnp1,ypn,yppn,ypnp1,yppnp1 : REAL_AR;
        i:INTEGER;
    Begin                             
      OMEGA:=TWO*PI*FREQU;
      OMEGA2:=OMEGA*OMEGA;
      DELTA:=SQRT(ONE-DZETA*DZETA);
      EXPA:=EXP(-OMEGA*DZETA*T);
      GOMEGA:=OMEGA*DELTA;
      ARG:=GOMEGA*T;
      SINARG:=SIN(ARG);
      COSARG:=COS(ARG);
      QSI:=DZETA/DELTA;
      Q0:=(COSARG-QSI*SINARG)*EXPA;
      Q1:=OMEGA2*SINARG/GOMEGA*EXPA;
      Q2:=(ONE+(QSI*SINARG-COSARG)*EXPA)/T;
      R1:=SINARG/GOMEGA*EXPA;
      R0:=COSARG*EXPA+DZETA*OMEGA*R1;
      R2:=(ONE-DZETA*OMEGA*T)*R1/T-COSARG*EXPA;
      R3:=ONE-R1/T;
      XPNP1:=VIT^[1];
      YPPNP1:=ZERO;
      YPNP1:=ZERO;
      ACC^[1]:=Q2*XPNP1;
                                                                                
      For I:=2 to NDATA do
      begin
        YPN:=YPNP1;
        YPPN:=YPPNP1;
        XPN:=XPNP1;
        XPNP1:=VIT^[I];
        YPPNP1:=Q0*YPPN+Q1*(XPN-YPN)+Q2*(XPNP1-XPN);
        ACC^[I]:=YPPNP1;
        YPNP1:=R0*YPN+R1*YPPN+R2*XPN+R3*XPNP1
      end                                                                          
    End;                                                                       
                                                                                
   {-----------------------------------------------------------
     SEEK MINIMUM AND MAXIMUM OF A TABLE C(I)

     Inputs:
                C(I)    Table with NDATA values of acceleration
                NDATA   Number of points of ACC(I)
     Outputs:
                CMIN    minimum value of table C
                CMAX    maximum value of table C
    -----------------------------------------------------------}
    PROCEDURE EXTREM(C:RV; VAR CMIN:REAL_AR; VAR CMAX:REAL_AR;NDATA:INTEGER);
    Var I:INTEGER;
    Begin                                      
      CMIN:=ZERO;
      CMAX:=ZERO;
      For I:=1 to NDATA do
      begin
        IF C^[I] > CMAX then CMAX:=C^[I];
        IF C^[I] < CMIN then CMIN:=C^[I]
      end                                                                    
    End;                                                                       


   {-------------------------------------------------------------
     ACCELERATION SHOCK SPECTRUM OF A SPEED

     Inputs:
                VIT(I)  speed of basis from deconvolution
                TS      sampling time of speed
                FMIN    begin frequency of spectrum (<>0)
                FMAX    end frequency of spectrum (<=Nyquist)
                DZETA   relative damping factor
                N       number of points of VIT(I)
                NFREQ   number of frequencies of spectrum (<=400)
     Output:
                S(I)    shock spectrum (negative & positive)
    -------------------------------------------------------------}
    PROCEDURE SPEC(VIT:RV;DZETA,FMIN,FMAX:REAL_AR;N:INTEGER;T:REAL_AR;NFREQ:INTEGER;VAR S:RV);
    Var deltaf,f,ymin,ymax:REAL_AR; I:Integer;
    Begin                      
      DELTAF:=(FMAX-FMIN)/(NFREQ-1);
      {the algorithm fails if FMIN=0}
      IF FMIN < DELTAF then FMIN:=DELTAF;
      For I:=1 to NFREQ do
      begin
        F:=FMIN+(I-1)*DELTAF;
        OSCIL(F,VIT,ACC,T,DZETA,N);                                           
        EXTREM(ACC,YMIN,YMAX,N);                                              
        S^[2*I-1]:=YMAX;
        S^[2*I]:=YMIN
      end
    End;

  {Graph of input acceleration ACC(t) on request}
  Procedure View_Input_Signal;
  Begin
    ClrScr;
    CourbeXY(CrtDC,ndata,10,ACC,tbegin,tend);
    Legendes(CrtDC,'MASS INPUT ACCELERATION','Time (s)','Acc. m/s2');
    SortieGraphique
  End;


  {Graph of speed VIT(t) on request}
  Procedure View_Speed;
  Begin
    ClrScr;
    CourbeXY(CrtDC,ndata,10,VIT,tbegin,tend);
    Legendes(CrtDC,'BASIS SPEED','Time (s)','Speed m/s');
    SortieGraphique
  End;

  {Graph of shock spectrum S(freq) on request}
  Procedure View_Shock_Spectrum;
  Begin
    ClrScr;
    CourbeXY(CrtDC,nspec,10,S,ZERO,fsup);
    Legendes(CrtDC,'SHOCK SPECTRUM','Freq. (hz)','Acc. m/s2');
    SortieGraphique
  End;

{main program}
BEGIN
    {open main window application with title}
    WinCrtInit('Deconvolution of a signal F(t)');
    New(ACC); New(VIT); New(S);                                                
    writeln;
    writeln(' DECONVOLUTION OF A SIGNAL F(T)');
    writeln;
    write(' Input data file name (without .dat): '); readln(nom);
    writeln;
      
    OUTPUT:='v'+nom+'.lst';
    INPUT:=NOM+'.dat';

    Assign(fp_out,output); Rewrite(fp_out);
    Assign(fp_in,input); Reset(fp_in);

{ read title (1 line maxi 80 CharaCters) }
    READLN(fp_in,Title);
    WRITELN(fp_out,' DECONVOLUTION OF SIGNAL:');
    WRITELN(fp_out,' ',Title);
{ print resume of program's functions}

    WRITELN(fp_out);
    WRITELN(fp_out,' THIS PROGRAM ALLOWS TO:');
    WRITELN(fp_out);
    WRITELN(fp_out,'  - READ MEASUREMENTS IN A DATA FILE,');
    WRITELN(fp_out,'  - FILTER THESE MEASUREMENTS (OPTIONAL),');
    WRITELN(fp_out,'  - RECONSTRUCT THE SPEED AT BASIS BY DECONVOLUTION,');
    WRITELN(fp_out,'  - CALCULATE THE ACCELERATION SHOCK SPECTRUM');
    WRITELN(fp_out);
                                                                                
{ Frequency of captor }                                                         
    READLN(fp_in, F0);
    WRITELN(fp_out,' FREQUENCY OF CAPTOR: ',F0:6:1,' HERTZ');
                                                                                
{ rod damping }                                                 
    READLN(fp_in, DZETA);
    WRITELN(fp_out,' DAMPING OF CAPTOR: ',DZETA:5:2);
                                                                                
{ Value of Q for the shock spectrum }                                   
    READLN(fp_in, Q);                                                                    
    WRITELN(fp_out,' Q VALUE FOR THE SHOCK SPECTRUM: ',Q:4:0);
                                                                                
{ Number of frequencies for the shock spectrum }                       
    READLN(fp_in, NSPEC);
    WRITELN(fp_out,' NUMBER OF FREQUENCIES FOR SPECTRUM: ',NSPEC);
                                                                                
{ FREQUENCIES MINI & MAXI OF SPECTRUM }                                   
    READLN(fp_in, FINF);                                                              
    READLN(fp_in, FSUP);
    WRITELN(fp_out,' BEGIN FREQUENCY OF SPECTRUM: ',FINF:6:1);
    WRITELN(fp_out,' END FREQUENCY OF SPECTRUM  : ',FSUP:6:1);
    WRITELN(fp_out);
                                                                                
{ READ MEASUREMENTS }                                                          
    READATA(ACC,TS,NDATA);
    CLOSE(fp_in);

    Write(' Do you want to view the input signal (y/n) ? ');
    if readkey='y' then View_Input_Signal;

    Clrscr;
    Writeln;
    Writeln(' Computing basis speed...');

    FMAX:=0.5/TS;
    WRITELN(fp_out,' NYQUIST FREQUENCY: ',FMAX:6:1);
                                                                               
{ CALL DECONVOLUTION }                                                    
    DECON(f0,ACC,VIT,Ts,Dzeta,Ndata);
    IF FSUP > FMAX then FSUP:=FMAX;

{ print speed results to output file }
    WRITELN(fp_out);
    WRITELN(fp_out);
    WRITELN(fp_out,'   3RD COLUMN: SPEED OBTAINED BY DECONVOLUTION FROM');
    WRITELN(fp_out,'   ONEFILTERED MASS ACCELERATION.');
    WRITELN(fp_out,'   4TH COLUMN: SPEED OBTAINED BY DECONVOLUTION FROM');
    WRITELN(fp_out,'   FILTERED MASS ACCELERATION IF IT EXISTS.');
    WRITELN(fp_out);

    WRITELN(fp_out,'    TIME (S)     ACC M/S2      BASIS SPEED M/S  ');
    WRITELN(fp_out,'    ---------   -----------   ----------------- ');
    For I:=1 to NDATA do
      WRITELN(fp_out, (I-1)*TS:10:4,'  ',ACC^[I]:15:6,VIT^[I]:15:6);

{ Compute shock spectrum }
    SPEC(VIT,DZETA,FINF,FSUP,NDATA,TS,NSPEC,S);

{ print spectrum results to output file }
    DF:=(FSUP-FINF)/(NSPEC-1);
    WRITELN(fp_out);
    WRITELN(fp_out);
    WRITELN(fp_out,' THE ACCELERATION SHOCK SPECTRUM  IS COMPUTED FROM');
    WRITELN(fp_out,' THE FILTERED ACCELERATION IF IT EXISTS, ELSE FROM');
    WRITELN(fp_out,' THE ONEFILTERED ACCELERATION.');
    WRITELN(fp_out);
    WRITELN(fp_out,' FREQUENCY STEP OF SHOCK SPECTRUM:');
    WRITELN(fp_out,'   DF = ',DF:6:2);
    WRITELN(fp_out);
    WRITELN(fp_out,'   FREQUENCY HZ            SPECTRUM +             SPECTRUM -  IN M/S2');
    WRITELN(fp_out,'   ------------           -------------          -------------       ');
    For I:=1 to NSPEC do
    begin
      FR:=FINF+(I-1)*DF;                                                          
      WRITELN(fp_out, FR,S^[2*I-1],S^[2*I])
    end;
    WRITELN(fp_out);
    WRITELN(fp_out,' End of file ',output);
    CLOSE(fp_out);
    writeln(' Results in file ',output);
    Write(' Do you want to view the speed (y/n) ? ');
    if readkey='y' then View_Speed;
    Clrscr; writeln;
    Write(' Do you want to view the shock spectrum (y/n) ? ');
    if readkey='y' then View_Shock_Spectrum;
    Dispose(ACC); Dispose(VIT); Dispose(S);
    DoneWinCrt

END.
                                                                                      
{End of file gdeconv.pas}