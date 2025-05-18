{**************************************************************
* Program to demonstrate the Fast Fourier Transform Procedure *
*      (frequency analysis of a discrete signal F(t))         *
*                                                             *
*                      Pascal version by J-P Moreau, Paris    *
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
* Output file (tfft.lst):                                     *
*                                                             *
*        Frequency                Value                       *
* --------------------------------------------------          *
*  0.00000000000000E+0000  2.71689845900710E-0001             *
*  4.99511718749888E+0000  5.46335766118318E-0001             *
*  9.99023437499775E+0000  5.55559620102031E-0001             *
*  1.49853515624966E+0001  5.72241896169544E-0001             *
*  1.99804687499955E+0001  5.98828307688564E-0001             *
*  2.49755859374944E+0001  6.40060507907037E-0001             *
*  ...                     ...                                *
*  2.52253417968695E+0003  2.05196435017854E-0002             *
*  2.52752929687445E+0003  2.05175809733191E-0002             *
*  2.53252441406195E+0003  2.05156130592457E-0002             *
*  2.53751953124945E+0003  2.05133281918108E-0002             *
*  2.54251464843694E+0003  2.05124075002914E-0002             *
*  2.54750976562444E+0003  2.05113721630994E-0002             *
*  2.55250488281194E+0003  2.05108424250966E-0002             *
*                                                             *
**************************************************************}
PROGRAM Test_FFT;
Uses Wincrt,Type_def;

Var
        f_in,f_out : TEXT;
        i,n1,ndata,p : INTEGER;
        df,dt,f,temp,tbegin,tend : REAL_AR;
        signal : RV;        {pointer to real vector}
        csignal: CV;        {pointer to complex vector}


  PROCEDURE FFT(VAR A:CV; M:INTEGER);
  {**************************************************************
  *             FAST FOURIER TRANSFORM PROCEDURE                *
  * ----------------------------------------------------------- *
  * This procedure calculates the fast Fourier transform of a   *
  * real function sampled in N points 0,1,....,N-1. N must be a *
  * power of two (2 power M).                                   *
  *  T being the sampling duration, the maximum frequency in    *
  * the signal can't be greater than fc = 1/(2T). The resulting *
  * specter H(k) is discrete and contains the frequencies:      *
  *         fn = k/(NT) with k = -N/2,.. -1,0,1,2,..,N/2.       *
  *         H(0) corresponds to null fréquency.                 *
  *         H(N/2) corresponds to fc frequency.                 * 
  * ----------------------------------------------------------- *
  * INPUTS:                                                     *
  *                                                             *
  *        A(i)  complex vector of size N, the real part of     *
  *              which contains the N sampled points of real    *
  *              signal to analyse (time spacing is constant).  *
  *                                                             *
  *          M   integer such as N=2^M                          *
  *                                                             *
  * OUTPUTS:                                                    *
  *                                                             *
  *        A(i)  complex vector of size N, the vector modulus   *
  *              contain the frequencies of input signal, the   *
  *              vector angles contain the corresponding phases *
  *              (not used here).                               *
  *                                                             *
  *                                     J-P Moreau/J-P Dumont   *
  **************************************************************}
  VAR U,W,T                   :Complex;
  N,NV2,NM1,J,I,IP,K,L,LE,LE1 :INTEGER;
  Phi,temp                    :REAL_AR;
  
  BEGIN
    N:=(2 SHL (M-1));
    NV2:=(N SHR 1);
    NM1:=N-1;
    J:=1;
    FOR I:=1 TO NM1 DO
    BEGIN
      if (I<J) then
      BEGIN
        T:=A^[J];
        A^[J]:=A^[I];
        A^[I]:=T;
      END; (* if *)
      K:=NV2;
      if K<J then
      REPEAT
        J:=J-K;
        K:=(K SHR 1);
      UNTIL K>=J;
      J:=J+K;
    END;  (* Do *)
    LE:=1;
    FOR L:=1 TO M DO
    BEGIN
      LE1:=LE;
      LE:=(LE SHL 1);
      U[1]:=1.0;
      U[2]:=0.0;
      Phi:= Pi/LE1;
      W[1]:=Cos(Phi);
      W[2]:=Sin(Phi);
      FOR J:=1 TO LE1 DO
      BEGIN
        I:=J-LE;
        WHILE I< N-LE DO
        BEGIN
          I:=I+LE;
          IP:=I+LE1;
          T[1]:=A^[ip][1]*U[1]-A^[ip][2]*U[2];
          T[2]:=A^[ip][1]*U[2]+A^[ip][2]*U[1];
          A^[ip][1]:=A^[i][1]-T[1];
          A^[ip][2]:=A^[i][2]-T[2];
          A^[i][1]:=A^[i][1]+T[1];
          A^[i][2]:=A^[i][2]+T[2];
        END; (* WHILE *)
        temp:=U[1];
        U[1]:=W[1]*U[1]-W[2]*U[2];
        U[2]:=W[1]*U[2]+W[2]*temp;
      END;
    END;
    FOR I:=1 TO N DO
    BEGIN
      A^[I][1]:=A^[I][1]/N;
      A^[I][2]:=A^[I][2]/N;
    END
  END;

{main program}
BEGIN
  {open input and output file}
  Assign(f_in,'tfft.dat'); Reset(f_in);
  Assign(f_out,'tfft.lst'); Rewrite(f_out);
  {read number of input signal points in input file}
  read(f_in,ndata);
  {take nearest power of two}
  if  ndata > 2048 then  begin ndata:=2048; p:=11 end;
  if (ndata < 2048) and (ndata>1023) then begin ndata:=1024; p:=10 end;
  if (ndata < 1024) and (ndata>511) then  begin ndata:=512; p:=9 end;
  if (ndata <  512) and (ndata>255) then  begin ndata:=256; p:=8 end;
  if (ndata <  256) and (ndata>127) then  begin ndata:=128; p:=7 end;
  if (ndata <  128) and (ndata> 63) then  begin ndata:=64; p:=6 end;
  if (ndata <   64) and (ndata> 31) then  begin ndata:=32; p:=5 end;
  if  ndata <   32 then
  begin
    writeln(f_out,' Error: number of points too small (<32) !');
    close(f_out);
    writeln(' Results in file tfft.lst (error).');
    Readkey;
    DoneWinCrt
  end;
  {Dynamic allocation of vectors}
  New(signal); New(csignal);
  {read ndata couples T(i), Y(i) in input data file}
  for i:=1 to ndata do
  begin
    readln(f_in,temp,signal^[i]);
    if i=1 then tbegin:=temp;
    if i=ndata then tend:=temp
  end;
  {Put input signal in real part of complex vector csignal}
  for i:=1 to ndata do
  begin
    csignal^[i][1]:=signal^[i];
    csignal^[i][2]:=0.0
  end;
  {call FFT procedure}
  FFT(csignal,p);
  {get frequencies in signal real vector} 
  n1:=ndata DIV 2;
  for i:=1 to n1 do
  begin
    temp:=sqrt(csignal^[i][1]*csignal^[i][1]+csignal^[i][2]*csignal^[i][2]);
    if (i>1) then temp:=2*temp;
    signal^[i]:=temp
  end;
  {calculat sampling time range dt}
  dt:=(tend-tbegin)/(ndata-1);
  {calculate frequency step df}
  df:=1.0/(ndata*dt);
  {print frequency specter to output file}
  f:=0.0;
  writeln(f_out,'        Frequency                Value            ');
  writeln(f_out,'--------------------------------------------------');
  for i:=1 to n1 do
  begin
    writeln(f_out,' ',f,' ',signal^[i]);
    f:=f+df;
  end;

  close(f_out);
  writeln;
  writeln(' Results in file tfft.lst.');
  Readkey;
  DoneWinCrt

END.

{End of file tfft.pas}