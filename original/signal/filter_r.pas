{****************************************************************
*   Numerical Bandpass Filter (see demo program test_fil.pas)   *
* ------------------------------------------------------------- *
* References:                                                   * 
*                                                               *
* http://en.wikipedia.org/wiki/Digital_biquad_filter            *  
* http://www.musicdsp.org/archive.php?classid=3#225             *  
* http://www.musicdsp.org/showone.php?id=197                    *
* http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt           *
* http://www.musicdsp.org/archive.php?classid=3#225             *
* http://www.musicdsp.org/showArchiveComment.php?ArchiveID=225  *
*                                                               *
*                    Turbo-pascal Release By J-P Moreau, Paris. *
*                              (www.jpmoreau.fr)                *
****************************************************************} 
UNIT Filter_r;

interface

Uses WinCrt;

const
  kLowPass=0;       {-LowPass}
  kHighPass=1;      {-HiPass}
  kBandPassCSG=2;   {-BandPass CSG}
  kBandPassCZPG=3;  {-BandPass CZPG}
  kNotch=4;         {-Notch}
  kAll=5;           {-AllPass}
  kPeaking=6;       {-Peaking}
  kLowShelf=7;      {-LowShelf}
  kHighShelf=8;     {-HiShelf}

  NMAX = 2047;

type

  psingle = ^Tab;
  Tab = array[0..NMAX] of single;


  {Butterworth all-purpose filter}
  TRbjEqFilter = Object
    b0a0,b1a0,b2a0,a1a0,a2a0:single;
    in1,in2,ou1,ou2:single;
    fSampleRate:single;
    fMaxBlockSize:integer;
    fFilterType:integer;
    fFreq,fQ,fDBGain:single;
    fQIsBandWidth:boolean;

    out1: psingle;
    
    constructor create(SampleRate:single;MaxBlockSize:integer);
    procedure SetQ(NewQ:single);
    procedure CalcFilterCoeff;
    procedure CalcFilterCoeffs(pFilterType:integer;pFreq,pQ,pDBGain:single;pQIsBandWidth:boolean);
    procedure Process(Input:psingle; sampleframes:integer);
    function  Process1(input:single):single;
  end;

implementation

FUNCTION sinh(x:single):single;
VAR expx: single;
Begin
  expx := exp(x);
  sinh := 0.5*(expx-1.0/expx)
End;

FUNCTION Power(y,x:single):single;
begin
  IF x<0 THEN EXIT;
  Power := Exp(x*Ln(y))
end;

constructor TRbjEqFilter.create(SampleRate:single;MaxBlockSize:integer);
begin
  fMaxBlockSize:=MaxBlockSize;
  fSampleRate:=SampleRate;

  New(out1);

  fFilterType:=0;
  fFreq:=500;
  fQ:=0.3;
  fDBGain:=0;
  fQIsBandWidth:=true;

  in1:=0;
  in2:=0;
  ou1:=0;
  ou2:=0;
end;

procedure TRbjEqFilter.SetQ(NewQ:single);
begin
  fQ:=(1-NewQ)*0.98;
end;

procedure TRbjEqFilter.CalcFilterCoeffs(pFilterType:integer;pFreq,pQ,pDBGain:single;pQIsBandWidth:boolean);
begin
  fFilterType:=pFilterType;
  fFreq:=pFreq;
  fQ:=pQ;
  fDBGain:=pDBGain;
  fQIsBandWidth:=pQIsBandWidth;

  CalcFilterCoeff;
end;

procedure TRbjEqFilter.CalcFilterCoeff;
var
  alpha,a0,a1,a2,b0,b1,b2:single;
  A,beta,omega,tsin,tcos:single;
begin
  {peaking, LowShelf or HiShelf}
  if fFilterType>=6 then
  begin
    A:=power(10.0,(fDBGain/40.0));
    omega:=2*pi*fFreq/fSampleRate;
    tsin:=sin(omega);
    tcos:=cos(omega);

    if fQIsBandWidth then
      alpha:=tsin*sinh(Ln(2.0)/2.0*fQ*omega/tsin)
    else
      alpha:=tsin/(2.0*fQ);

    beta:=sqrt(A)/fQ;

    {peaking}
    if fFilterType=6 then
    begin
      b0:=1.0+alpha*A;
      b1:=-2.0*tcos;
      b2:=1.0-alpha*A;
      a0:=1.0+alpha/A;
      a1:=-2.0*tcos;
      a2:=1.0-alpha/A;
    end else
    {lowshelf}
    if fFilterType=7 then
    begin
      b0:=(A*((A+1.0)-(A-1.0)*tcos+beta*tsin));
      b1:=(2.0*A*((A-1.0)-(A+1.0)*tcos));
      b2:=(A*((A+1.0)-(A-1.0)*tcos-beta*tsin));
      a0:=((A+1.0)+(A-1.0)*tcos+beta*tsin);
      a1:=(-2.0*((A-1.0)+(A+1.0)*tcos));
      a2:=((A+1.0)+(A-1.0)*tcos-beta*tsin);
    end;
    {hishelf}
    if fFilterType=8 then
    begin
      b0:=(A*((A+1.0)+(A-1.0)*tcos+beta*tsin));
      b1:=(-2.0*A*((A-1.0)+(A+1.0)*tcos));
      b2:=(A*((A+1.0)+(A-1.0)*tcos-beta*tsin));
      a0:=((A+1.0)-(A-1.0)*tcos+beta*tsin);
      a1:=(2.0*((A-1.0)-(A+1.0)*tcos));
      a2:=((A+1.0)-(A-1.0)*tcos-beta*tsin);
    end;
  end else  {other filter types}
  begin
    omega:=2*pi*fFreq/fSampleRate;
    tsin:=sin(omega);
    tcos:=cos(omega);

    if fQIsBandWidth then
      alpha:=tsin*sinh(Ln(2)/2*fQ*omega/tsin)
    else
      alpha:=tsin/(2*fQ);
    {lowpass}
    if fFilterType=0 then
    begin
      b0:=(1-tcos)/2;
      b1:=1-tcos;
      b2:=(1-tcos)/2;
      a0:=1+alpha;
      a1:=-2*tcos;
      a2:=1-alpha;
    end else {hipass}
    if fFilterType=1 then
    begin
      b0:=(1+tcos)/2;
      b1:=-(1+tcos);
      b2:=(1+tcos)/2;
      a0:=1+alpha;
      a1:=-2*tcos;
      a2:=1-alpha;
    end else {bandpass CSG}
    if fFilterType=2 then
    begin
      b0:=tsin/2;
      b1:=0;
      b2:=-tsin/2;
      a0:=1+alpha;
      a1:=-1*tcos;
      a2:=1-alpha;
    end else {bandpass CZPG}
    if fFilterType=3 then
    begin
      b0:=alpha;
      b1:=0.0;
      b2:=-alpha;
      a0:=1.0+alpha;
      a1:=-2.0*tcos;
      a2:=1.0-alpha;
    end else  {notch}
    if fFilterType=4 then
    begin
      b0:=1.0;
      b1:=-2.0*tcos;
      b2:=1.0;
      a0:=1.0+alpha;
      a1:=-2.0*tcos;
      a2:=1.0-alpha;
    end else   {allpass}
    if fFilterType=5 then
    begin
      b0:=1.0-alpha;
      b1:=-2.0*tcos;
      b2:=1.0+alpha;
      a0:=1.0+alpha;
      a1:=-2.0*tcos;
      a2:=1.0-alpha;
    end;
  end;

  b0a0:=b0/a0;
  b1a0:=b1/a0;
  b2a0:=b2/a0;
  a1a0:=a1/a0;
  a2a0:=a2/a0;

end; {CalcFilterCoeff}


function TRbjEqFilter.Process1(input:single):single;
var
  LastOut:single;
begin
  {filter}
  LastOut:= b0a0*input + b1a0*in1 + b2a0*in2 - a1a0*ou1 - a2a0*ou2;

  {push in/out buffers}
  in2:=in1;
  in1:=input;
  ou2:=ou1;
  ou1:=LastOut;

  {return output}
  Process1:=LastOut;
end;

{
Note:
use Process1(input:single):single;
for per sample processing
use Process(Input:psingle;sampleframes:integer);
for block processing. The input is a pointer to
the start of an array of single which contains
the audio data.
i.e.
RBJFilter.Process(@WaveData[0],256);
}

procedure TRbjEqFilter.Process(Input:psingle; sampleframes:integer);
var
  i:integer;
  LastOut:single;
begin
  New(Out1);
  for i:=0 to SampleFrames-1 do
  begin
    {filter}
    LastOut:= b0a0*(input^[i])+ b1a0*in1 + b2a0*in2 - a1a0*ou1 - a2a0*ou2;
    { LastOut:=input^;
      push in/out buffers }
    in2:=in1;
    in1:=input^[i];
    ou2:=ou1;
    ou1:=LastOut;

    Out1^[i]:=LastOut

  end;
  Dispose(Out1)
end;

end.

{end of file filter_r.pas}
 