{*****************************************************************
* Calculate the acceleration of the sismic mass of an elementary *
* oscillator (1 degree of freedom), the basis of which is submit-*
* ted to a given speed, VIT(t). The input signal VIT is digitali-*
* zed with a constant time step, TS. The table VIT(I) contains   *
* NDATA speed values.                                            *
*                                                                *
*                           Pascal Version By J-P Moreau, Paris. *
* -------------------------------------------------------------- *
* SAMPLE RUN:                                                    *
*                                                                *
* Read data from file speed.dat                                  *
* Write results to file response.txt.                            *
*****************************************************************}
Program Response1;
Uses WinCrt,Type_def;

Const
       ZERO = 0.0;
       ONE  = 1.0;
       TWO  = 2.0;
Var        
       ACC, VIT: RV;
       I, NDATA: Integer;
       dzeta, F0, FMAX, T, TBEGIN, TEND, TS: double;

       fp_in,fp_out: TEXT;


    Procedure READATA(VIT:RV; Var TS:Double; Var NDATA:Integer);
        {------------------------------------------------
        'Read ndata: number of points of input signal
        '        ts: sampling duration of signal in sec.
        '    ACC(i): Table with ndata acceleration values
        '-----------------------------------------------} 
    var TEMP: double; I:integer;
    Begin
        Readln(fp_in, NDATA);
        Writeln(fp_out, ' ' ,NDATA, ' POINTS READ.');
        Readln(fp_in, TS);
        Writeln(fp_out, ' SAMPLING DURATION: ', TS);
        For I := 1 To NDATA do
            Readln(fp_in, TEMP, VIT^[I]);
    End;


    Procedure OSCIL(FREQU:double; VIT, ACC:RV; TS, DZETA: double; NDATA: Integer);
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
        *                     Pascal Version By J-P Moreau, Paris.  *
        ************************************************************}
    Var
        ARG, DELTA, EXPA, GOMEGA, OMEGA, OMEGA2, T: double;
        COSARG, Q0, Q1, Q2, QSI, R0, R1, R2, R3, SINARG: double;
        XPN, XPNP1, YPN, YPNP1, YPPN, YPPNP1: double;
        I:integer;
    Begin
        T := TS;
        OMEGA := 2 * PI * FREQU;
        OMEGA2 := OMEGA * OMEGA;
        DELTA := Sqrt(1 - DZETA * DZETA);
        EXPA := Exp(-OMEGA * DZETA * T);
        GOMEGA := OMEGA * DELTA;
        ARG := GOMEGA * T;
        SINARG := Sin(ARG);
        COSARG := Cos(ARG);
        QSI := DZETA / DELTA;
        Q0 := (COSARG - QSI * SINARG) * EXPA;
        Q1 := OMEGA2 * SINARG / GOMEGA * EXPA;
        Q2 := (ONE + (QSI * SINARG - COSARG) * EXPA) / T;
        R1 := SINARG / GOMEGA * EXPA;
        R0 := COSARG * EXPA + DZETA * OMEGA * R1;
        R2 := (1 - DZETA * OMEGA * T) * R1 / T - COSARG * EXPA;
        R3 := ONE - R1 / T;
        XPNP1 := VIT^[1];
        YPPNP1 := ZERO;
        YPNP1 := ZERO;
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

    End;

{main program}
BEGIN

        New(ACC); New(VIT);

        Assign(fp_out,'response.txt'); Rewrite(fp_out);
        Assign(fp_in,'speed.dat'); Reset(fp_in);

        Writeln(fp_out);
        Writeln(fp_out, ' Mass Response to speed(T) at Basis');
        Writeln(fp_out);

        F0 := 250.0; dzeta := 0.05;

        {Frequency of oscillator}
        Writeln(fp_out, ' FREQUENCY OF OSCILLATOR: ', F0, ' HERTZ');
        Writeln(fp_out);

        {oscillator damping}
        Writeln(fp_out, ' DAMPING OF OSCILLATOR: ', DZETA);
        Writeln(fp_out);

        {READ INPUT SPEED}
        READATA(VIT, TS, NDATA);

        Close(fp_in);

        TBEGIN := 0.0; TEND := TBEGIN + (TS * (NDATA - 1));

        FMAX := 0.5 / TS;
        Writeln(fp_out, ' NYQUIST FREQUENCY: ', FMAX);

        {Compute mass response}
        OSCIL(F0, VIT, ACC, TS, dzeta, NDATA);

        Writeln(fp_out);
        Writeln(fp_out, '    TIME (S)    MASS ACC. (M/S2) ');
        Writeln(fp_out, '  ------------  -----------------');

        T := TBEGIN;

        For I := 1 To NDATA do
        begin
            Writeln(fp_out, ' ', T, '  ', ACC^[I]);
            T := T + TS
        end;

        Close(fp_out);


        Dispose(ACC); Dispose(VIT);

        writeln;
        writeln(' Results in file response.txt...');
        readkey;

        DoneWinCrt

END.

{end of file repons1.pas}