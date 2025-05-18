{*******************************************************************
* SMOOTHING AN ARRAY OF N ORDINATES Y's (ASCENDING ORDER ABCISSAS) *
* ---------------------------------------------------------------- *
* Description:                                                     *
* This program uses the procedure SMOOFT for smoothing an array of *
* given ordinates (y's) that are in order of increasing abscissas  *
* (x's), but without using the abscissas themselves supposed to be *
* equally spaced. It first removes any linear trend, then uses a   *
* Fast Fourier Transform procedure (REALFT) to low-pass filter the *
* data. The linear trend is reinserted at the end. One user-speci- *
* fied parameter, EPS, enters "the amount of smoothing", given as  *
* the number of points over which the data should be smoothed.     *
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
* Input data file contains:                                        *
* 1024                                                             *
* 0.00000000000000E+0000 7.50000000000000E-0001                    *
* 9.21288168207468E-0003 7.77368637910513E-0001                    *
* 1.84257633641494E-0002 8.34466556277221E-0001                    *
* 2.76386450462383E-0002 9.03071871110114E-0001                    *
* 3.68515267282987E-0002 9.92958153417021E-0001                    *
* 4.60644084103592E-0002 1.09195646826811E+0000                    *
* 5.52772900924197E-0002 1.15230452277865E+0000                    *
* 6.44901717745370E-0002 1.06763022290215E+0000                    *
* 7.37030534565974E-0002 1.34541171127239E+0000                    *
* 8.29159351386579E-0002 1.48611048393104E+0000                    *
* 9.21288168207184E-0002 1.09349703210864E+0000                    *
* 1.01341698502779E-0001 1.72386602840743E+0000                    *
* 1.10554580184839E-0001 1.14317464708984E+0000                    *
* ---------------------- ----------------------                    *
* 9.37871355209791E+0000 2.43969819122867E+0001                    *
* 9.38792643377383E+0000 2.42468007203424E+0001                    *
* 9.39713931544975E+0000 2.42436619192304E+0001                    *
* 9.40635219712567E+0000 2.42829449073179E+0001                    *
* 9.41556507880159E+0000 2.42980085689633E+0001                    *
* 9.42477796047751E+0000 2.43119449022633E+0001                    *
*                                                                  *
* Output file contains (here EPS=30):                              *
*                                                                  *
*       Time          Y        Smoothed Y                          *
* ----------------------------------------                         *
*     0.000000     0.750000     0.788489                           *
*     0.009213     0.777369     0.816559                           *
*     0.018426     0.834467     0.846303                           *
*     0.027639     0.903072     0.877999                           *
*     0.036852     0.992958     0.911780                           *
*     0.046064     1.091956     0.947587                           *
*     0.055277     1.152305     0.985167                           *
*     0.064490     1.067630     1.024100                           *
*     0.073703     1.345412     1.063860                           *
*     0.082916     1.486110     1.103897                           *
*     0.092129     1.093497     1.143733                           *
*     0.101342     1.723866     1.183046                           *
*     0.110555     1.143175     1.221739                           *
*     --------     --------     --------                           *
*     9.378714    24.396982    24.223682                           *
*     9.387926    24.246801    24.246966                           *
*     9.397139    24.243662    24.271064                           *
*     9.406352    24.282945    24.295721                           *
*     9.415565    24.298009    24.320887                           *
*     9.424778    24.311945    24.346705                           *
*                                                                  *
* On request, a graph showing both curves is drawn.                *
*                                                                  *
* ---------------------------------------------------------------- *
* Reference:  "Numerical Recipes By W.H. Press, B. P. Flannery,    *
*              S.A. Teukolsky and W.T. Vetterling, Cambridge       *
*              University Press, 1986" [BIBLI 08].                 *
*                                                                  *
*                              TPW Release By J-P Moreau, Paris.   *
*                                      (www.jpmoreau.fr)           *
*******************************************************************}
PROGRAM TEST_SMOOFT;

Uses WinCrtMy, Type_def, Graph_2d;

Var
      i, ndata: Integer;
      signal, ysave: RV;            {see type_def.pas}
      dt,t,tbegin,temp,tend: REAL;

      fp_in,fp_out: TEXT;
      rep:char;

      Procedure realft(var data:RV; n,isign:Integer); Forward;


Function MAX(a,b:REAL):REAL;
Begin
  if A>=b then MAX:=a
          else MAX:=b
End;


Procedure SMOOFT(Var Y:RV; N,PTS: Integer);
{---------------------------------------------------------------
  Smooths an array Y of length N, with a window whose full width
  is of order PTS neighboring points, a user supplied value. 
  Array Y is modified.
!--------------------------------------------------------------}
Label 1, return;
Var
    J,K,M,MO2,NMIN: Integer;
    CNST,FAC,RN1,Y1,YN: REAL;
Begin
  
  M:=2;
  NMIN:=N+Round(2.0*PTS);
1:IF M < NMIN THEN
  begin
    M:=2*M;
    GOTO 1
  end;
  IF M > Size then
  begin
    writeln(' Sample too big.');
    goto return
  end;
  CNST:=(PTS/M)*(PTS/M);
  Y1:=Y^[1];
  YN:=Y^[N];
  RN1:=1.0/(N-1);
  For J:=1 to N do      {Remove linear trend}
    Y^[J]:=Y^[J]-RN1*(Y1*(N-J)+YN*(J-1));
  IF N+1 <= M THEN
    For J:=N+1 to M do Y^[J]:=0.0;
  MO2:=M div 2;
  REALFT(Y,MO2,1);      {Fourier transform}
  Y^[1]:=Y^[1]/MO2;
  FAC:=1.0;
  For J:=1 to MO2-1 do
  begin
    K:=2*J+1;
    IF FAC <> 0.0 THEN
    begin
      FAC:=MAX(00.,(1.0-CNST*J*J)/MO2);
      Y^[K]:=FAC*Y^[K];
      Y^[K+1]:=FAC*Y^[K+1]
    end
    else
    begin
      Y^[K]:=0.0;
      Y^[K+1]:=0.0
    end
  end;
  FAC:=MAX(0.0,(1.0-0.25*PTS*PTS)/MO2);  {Last point}
  Y^[2]:=FAC*Y^[2];
  REALFT(Y,MO2,-1);       {Inverse Fourier transform}
  For J:=1 to N do        {Restore linear trend}
    Y^[J]:=RN1*(Y1*(N-J)+YN*(J-1))+Y^[J];
return: End;


Procedure four1(Var data:RV; nn,isign:Integer); 
{-------------------------------------------------------------------------------------------- 
 Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or replaces 
 data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as -1. 
 data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn 
 MUST be an integer power of 2 (this is not checked for here). 
!-------------------------------------------------------------------------------------------}
Label 1,2;
Var
    i,istep,j,m,mmax,n: Integer; 
    tempi,tempr: REAL; 
    theta,wi,wpi,wpr,wr,wtemp: Double; {Double precision for the trigonometric recurrences}
Begin     
  n:=2*nn;
  j:=1;
  i:=1;
  while i<=n do   {This is the bit-reversal section of the routine}
  begin
    if j > i then
    begin 
      tempr:=data^[j];      {Exchange the two complex numbers}
      tempi:=data^[j+1];
      data^[j]:=data^[i];
      data^[j+1]:=data^[i+1];
      data^[i]:=tempr;
      data^[i+1]:=tempi
    end; 
    m:=nn;
1:  if (m >=2) and (j > m) then
    begin 
      j:=j-m;
      m:=m div 2;
      goto 1 
    end; 
    j:=j+m;
    Inc(i,2)
  end; 
  mmax:=2;                   {Here begins the Danielson-Lanczos section of the routine}
2:if n > mmax then                                  {Outer loop executed log2 nn times}
  begin
    istep:=2*mmax;
    theta:=6.28318530717959/(isign*mmax); {Initialize for the trigonometric recurrence}
    wpr:=-2.0*sin(0.5*theta)*sin(0.5*theta);
    wpi:=sin(theta);
    wr:=1.0;
    wi:=0.0;
    m:=1;
    While m<=mmax do                              {Here are the two nested inner loops}
    begin
      i:=m;
      While i<=n do
      begin
        j:=i+mmax;                            {This is the Danielson-Lanczos formula: }
        tempr:=wr*data^[j]-wi*data^[j+1];
        tempi:=wr*data^[j+1]+wi*data^[j];
        data^[j]:=data^[i]-tempr;
        data^[j+1]:=data^[i+1]-tempi;
        data^[i]:=data^[i]+tempr;
        data^[i+1]:=data^[i+1]+tempi;
        Inc(i,istep)
      end;
      wtemp:=wr;  {Trigonometric recurrence}
      wr:=wr*wpr-wi*wpi+wr;
      wi:=wi*wpr+wtemp*wpi+wi;
      Inc(m,2)
    end; 
    mmax:=istep;
    goto 2   {Not yet done}
  end        {All done} 
End; 

Procedure realft(Var data:RV; n,isign:Integer);
{--------------------------------------------------------------------------------------------
 USES four1
 Calculates the Fourier transform of a set of n real-valued data points. Replaces this data 
 (which is stored in array data(1:n)) by the positive frequency half of its complex Fourier 
 transform. The real-valued first and last components of the complex transform are returned 
 as elements data(1) and data(2), respectively. n must be a power of 2. This routine 
 also calculates the inverse transform of a complex data array if it is the transform of real 
 data. (Result in this case must be multiplied by 2/n.) 
--------------------------------------------------------------------------------------------}
Var
    i,i1,i2,i3,i4,n2p3: Integer; 
    c1,c2,h1i,h1r,h2i,h2r,wis,wrs: REAL;
    theta,wi,wpi,wpr,wr,wtemp: Double;    {Double precision for the trigonometric recurrences}
Begin     
  theta:=PI/(n div 2);                    {Initialize the recurrence}
  c1:=0.5;
  if isign =  1 then
  begin 
    c2:=-0.5;
    four1(data,n div 2,1);                {The forward transform is here}
  end   
  else
  begin
    c2:=0.5;                              {Otherwise set up for an inverse transform}
    theta:=-theta
  end;
  wpr:=-2.0*sin(0.5*theta)*sin(0.5*theta);
  wpi:=sin(theta);
  wr:=1.0+wpr;
  wi:=wpi;
  n2p3:=n+3;
  For i:=2 to n div 4 do                  {Case i=1 done separately below}
  begin
    i1:=2*i-1;
    i2:=i1+1;
    i3:=n2p3-i2;
    i4:=i3+1;
    wrs:=wr;
    wis:=wi;
    h1r:=c1*(data^[i1]+data^[i3]);    {The two separate transforms are separated out of data}
    h1i:=c1*(data^[i2]-data^[i4]);
    h2r:=-c2*(data^[i2]+data^[i4]);
    h2i:=c2*(data^[i1]-data^[i3]);
    data^[i1]:=h1r+wrs*h2r-wis*h2i;  {Here they are recombined to form the true transform
                                     of the original real data  } 
    data^[i2]:=h1i+wrs*h2i+wis*h2r;
    data^[i3]:=h1r-wrs*h2r+wis*h2i;
    data^[i4]:=-h1i+wrs*h2i+wis*h2r;
    wtemp:=wr;                      {The recurrence}
    wr:=wr*wpr-wi*wpi+wr;
    wi:=wi*wpr+wtemp*wpi+wi
  end;

  if isign = 1 then
  begin 
    h1r:=data^[1];
    data^[1]:=h1r+data^[2];
    data^[2]:=h1r-data^[2]     {Squeeze the first and last data together to get}
  end                          {them all within the original array}
  else
  begin
    h1r:=data^[1];
    data^[1]:=c1*(h1r+data^[2]);
    data^[2]:=c1*(h1r-data^[2]);
    four1(data,n div 2,-1)  {This is the inverse transform for the case isign:-1}
  end 
End; 


{main program}
BEGIN

  WinCrtInit(' SMOOTH');
  New(signal); New(ysave);

  {open input and output files}
  Assign(fp_in,'smooth.dat'); Reset(fp_in);
  Assign(fp_out,'smooth.lst'); Rewrite(fp_out);

  {read number of input signal points in input file}
  readln(fp_in, ndata);

  {take nearest power of two}
  if (ndata > 2048) then ndata:=2048;
  if (ndata<2048) and (ndata>1023) then ndata:=1024;
  if (ndata<1024) and (ndata>511) then ndata:=512;
  if (ndata<512) and (ndata>255) then ndata:=256;
  if (ndata<256) and (ndata>127) then ndata:=128;
  if (ndata<128) and (ndata>63) then ndata:=64;
  if (ndata<64) and (ndata>31) then ndata:=32;
  if  (ndata<32) then
  begin
    writeln(fp_out,' Error: number of points too small (<32) !');
    close(fp_out);
    writeln(' Results in file smooth.lst (error).');
    Halt(0)
  end;

  {read ndata couples T(i), Y(i) in input data file}
  For i:=1 to ndata do
  begin
    readln(fp_in, temp, signal^[i]);
    if i = 1 then  tbegin:=temp;
    if i = ndata then tend:=temp
  end;
  close(fp_in);

  For i:=1 to ndata do ysave^[i]:=signal^[i];

  SMOOFT(signal,ndata,30);

  dt:=(tend-tbegin)/(ndata-1);
  t:=tbegin-dt;
  writeln(fp_out,'      Time          Y        Smoothed Y ');
  writeln(fp_out,'----------------------------------------');
  For i:=1 to ndata do
  begin
    t:=t+dt;
    writeln(fp_out, '   ', t:10:6, '   ', ysave^[i]:10:6, '   ', signal^[i]:10:6)
  end;

  close(fp_out);
  writeln;
  writeln(' Results in file smooth.lst.');
  writeln;
  write(' Do you want a graph of signals (y/n) ? '); readln(rep);

  if rep='y' then
  begin
    clrscr;
    {manual scaling}
    Ech_auto:=False;              {see Graph_2d.pas}
    X_mini:=0.0; X_maxi:=tend;
    Y_mini:=0.0; Y_maxi:=20.0;
    CourbeXY(CrtDc,ndata,5,ysave,tbegin,tend);
    Legendes(CrtDC,'Input Data','X','Y');
    CourbeXY(CrtDc,ndata,3,signal,tbegin,tend);
    Legendes(CrtDC,'Smoothed Data','X','Y');
    SortieGraphique
  end;

  Dispose(signal); Dispose(ysave);
  DoneWinCrt

END.

{end of file smooth.pas}