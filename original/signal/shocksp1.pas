    {************************************************************
    * Calculate the acceleration shock spectrum of a speed V(t) *
    * --------------------------------------------------------- *
    * SAMPLE RUN:                                               *
    * Input speed: in file speed.dat                            *
    * Output shock spectrum in file shocksp1.txt                *
    *                                                           *
    *                      Pascal Version By J-P Moreau, Paris. *
	*                               (www.jpmoreau.fr)           *
    ************************************************************}
    Program Shocksp1;

    Uses Wincrt, Type_def;

    Var
        VIT, S: RV;
        I, NDATA, NFREQ: Integer;
        DF, dzeta, F, FMAX, FMIN, TS: REAL;
        fp_in, fp_out: TEXT;


    Procedure READATA(VIT:RV; Var TS:REAL; Var NDATA:Integer);
    {------------------------------------------------
    'Read ndata: number of points of input signal
    '        ts: sampling duration of signal in sec.
    '    VIT(i): Table with ndata speed values
    '------------------------------------------------}
    Var TEMP: REAL; I: Integer;
    Begin
        Readln(fp_in, NDATA);
        Writeln(fp_out, ' ', NDATA, ' POINTS READ.');
        Readln(fp_in, TS);
        Writeln(fp_out, ' SAMPLING DURATION: ', TS);
        For I := 1 To NDATA do
            Readln(fp_in, TEMP, VIT^[I]);
    End;

    Procedure OSCIL(FREQU: REAL; VIT: RV; ACC: RV; T, DZETA: REAL; NDATA: Integer);
        {************************************************************
        * This subroutine calculates the acceleration of the sismic *
        * mass of an elementary oscillator (1 degree of freedom),   *
        * the basis of which is submitted to a given speed, VIT(t). *
        * The input signal VIT is digitalized with a constant time  *
        * step, TS. The table VIT(I) contains NDATA speed values.   *
        * --------------------------------------------------------- *
        * INPUTS:                                                   *
        *         FREQU........: eigen frequency of oscillator      *
        *         VIT..........: input speed signal of basis VIT(t) *
        *         TS...........: time step of signal (constant)     *
        *         DZETA........: reduced damping factor             *
        *         NDATA........: number of points of signal         *
        * OUTPUT:                                                   *
        *         ACC..........: table containing acceleration res- *
        *                        ponse of oscillator (NDATA points) *
        *                                                           *
        *             VB 2008 Express Version By J-P Moreau, Paris. *
        ************************************************************}
    Var
        ARG, COSARG, DELTA, EXPA, GOMEGA, OMEGA, OMEGA2, SINARG: Double;
        Q0, Q1, Q2, QSI, R0, R1, R2, R3, XPNP1, YPNP1, YPPNP1: Double;
        XPN, YPN, YPPN: Double; I: Integer;
    Begin
        OMEGA := 2 * PI * FREQU;
        OMEGA2 := OMEGA * OMEGA;
        DELTA := Sqrt(1.0 - DZETA * DZETA);
        EXPA := Exp(-OMEGA * DZETA * T);
        GOMEGA := OMEGA * DELTA;
        ARG := GOMEGA * T;
        SINARG := Sin(ARG);
        COSARG := Cos(ARG);
        QSI := DZETA / DELTA;
        Q0 := (COSARG - QSI * SINARG) * EXPA;
        Q1 := OMEGA2 * SINARG / GOMEGA * EXPA;
        Q2 := (1.0 + (QSI * SINARG - COSARG) * EXPA) / T;
        R1 := SINARG / GOMEGA * EXPA;
        R0 := COSARG * EXPA + DZETA * OMEGA * R1;
        R2 := (1 - DZETA * OMEGA * T) * R1 / T - COSARG * EXPA;
        R3 := 1 - R1 / T;
        XPNP1 := VIT^[1];
        YPPNP1 := 0.0;
        YPNP1 := 0.0;
        ACC^[1] := Q2 * XPNP1;
        For I := 2 To NDATA do
        begin
            YPN := YPNP1;
            YPPN := YPPNP1;
            XPN := XPNP1;
            XPNP1 := VIT^[I];
            YPPNP1 := Q0 * YPPN + Q1 * (XPN - YPN) + Q2 * (XPNP1 - XPN);
            ACC^[I] := YPPNP1;
            YPNP1 := R0 * YPN + R1 * YPPN + R2 * XPN + R3 * XPNP1
        end

    End; {OSCIL}


    Procedure EXTREM(ACC: RV; Var YMIN, YMAX: REAL; NDATA: Integer);
    {-----------------------------------------------------------
    ' SEEK MINIMUM AND MAXIMUM OF A TABLE ACC(I)
    '
    ' Inputs:
    '            ACC(I)  Table with NDATA values of acceleration
    '            NDATA   Number of points of ACC(I)
    ' Outputs:
    '            YMIN    minimum value of table ACC
    '            YMAX    maximum value of table ACC
    '----------------------------------------------------------}
    Var I: Integer;
    Begin
        YMIN := 0.0; YMAX := 0.0;
        For I := 1 To NDATA do
        begin
            If ACC^[I] > YMAX Then YMAX := ACC^[I];
            If ACC^[I] < YMIN Then YMIN := ACC^[I]
        end

    End; {EXTREM}


    Procedure SPEC(VIT: RV; TS, FMIN, FMAX, DZETA: REAL; NDATA, NFREQ: Integer; S: RV);
    {-------------------------------------------------------------
    ' ACCELERATION SHOCK SPECTRUM OF A SPEED
    '
    ' Inputs:
    '            VIT(I)  speed of basis (from deconvolution or else)
    '            TS      sampling time of speed
    '            FMIN    begin frequency of spectrum (<>0)
    '            FMAX    end frequency of spectrum (<=Nyquist)
    '            DZETA   relative damping factor
    '            NDATA   number of points of VIT(I)
    '            NFREQ   number of frequencies of spectrum (<=1024)
    ' Output:
    '            S(I)    shock spectrum (negative & positive)
    '-------------------------------------------------------------}
    Var ACC: RV; DELTAF, FR, YMIN, YMAX: REAL; II: Integer;
    Begin
        New(ACC);
        DELTAF := (FMAX - FMIN) / (NFREQ - 1);
        { The algorithm fails if FMIN=0 }
        If FMIN < DELTAF Then FMIN := DELTAF;
        For II := 1 To NFREQ do
        begin
            FR := FMIN + (II - 1) * DELTAF;
            OSCIL(FR, VIT, ACC, TS, DZETA, NDATA);
            EXTREM(ACC, YMIN, YMAX, NDATA);
            S^[2 * II - 1] := YMAX;
            S^[2 * II] := YMIN
        end;
        Dispose(ACC)
    End;


    {main program}
    BEGIN

        New(VIT); New(S);

        Writeln;
        Writeln(' Computing...');

        Assign(fp_in, 'speed.dat'); Reset(fp_in);
        Assign(fp_out, 'shocksp1.txt'); Rewrite(fp_out);

        Writeln(fp_out);
        Writeln(fp_out, ' Acceleration Shock Spectrum of a speed(T) at Basis');
        Writeln(fp_out);

        dzeta := 0.05;    {damping coefficient

        READ INPUT SPEED }
        READATA(VIT, TS, NDATA);

        Close(fp_in);

        FMAX := 0.5 / TS;
        Writeln(fp_out, ' NYQUIST FREQUENCY:', FMAX);

        FMIN := 10.0; FMAX := 1000.0; NFREQ := 991;

        {Calculate shock spectrum S }
        SPEC(VIT, TS, FMIN, FMAX, dzeta, NDATA, NFREQ, S);

        {Write spectrum to output file }
        DF := (FMAX - FMIN) / (1.0*(NFREQ - 1));
        Writeln(fp_out);
        Writeln(fp_out);
        Writeln(fp_out, ' THE ACCELERATION SHOCK SPECTRUM  IS COMPUTED FROM');
        Writeln(fp_out, ' THE FILTERED BASIS SPEED IF IT EXISTS, ELSE FROM');
        Writeln(fp_out, ' THE UNFILTERED BASIS SPEED.');
        Writeln(fp_out);
        Writeln(fp_out, ' FREQUENCY STEP OF SHOCK SPECTRUM: ', DF);
        Writeln(fp_out);
        Writeln(fp_out);
        Writeln(fp_out, '         FREQUENCY HZ              SPECTRUM +          SPECTRUM -  IN M/S2');
        Writeln(fp_out, '         ------------             -----------          -------------------');

        For I := 1 To NFREQ do
        begin
            F := FMIN + (I - 1) * DF;
            Write(fp_out, '   ', F);
            Write(fp_out, '  ', S^[2 * I - 1]);
            Writeln(fp_out, '  ', S^[2 * I])
        end;

        {Closing section}
        Writeln(fp_out);
        Writeln(fp_out, ' End of file shockspl.txt.');
        Writeln(fp_out);
        Close(fp_out);

        Dispose(VIT); Dispose(S);

        ClrScr;
        Writeln;
        Writeln(' Results in file shockspl.txt...');
        Writeln;

        ReadKey;
        DoneWinCrt

    END.

{End of file shocksp1.pas}