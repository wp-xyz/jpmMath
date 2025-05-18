{**********************************************************************
*   Demonstration of a bandpass numerical filter using Unit filter_r  *
* ------------------------------------------------------------------- *
* SAMPLE RUN:                                                         *
* Input file 'input.txt' contains:                                    *
*                                                                     *
* Input Signal                                                        *
* 200                                                                 *
* 0.00000000000000E+0000   0.00000000000000E+0000                     *
* 2.51256700736135E-0003   1.40504705905914E+0000                     *
* 5.02513401472271E-0003   1.39680230617523E+0000                     *
* 7.53770102208762E-0003   3.99903088808060E-0001                     *
* 1.00502680294454E-0002   0.00000000000000E+0000                     *
* 1.25628350368032E-0002   7.07106769084930E-0001                     *
* 1.50754020441752E-0002   1.26007354259491E+0000                     *
* 1.75879690515330E-0002   4.31350797414780E-0001                     *
* 2.01005360588908E-0002  -1.17557048797607E+0000                     *
* --/--                                                               *
* 4.77387731398721E-0001   1.00000000000000E+0000                     *
* 4.79900298405937E-0001   1.84206306934357E+0000                     *
* 4.82412865413608E-0001   1.17557048797607E+0000                     *
* 4.84925432420823E-0001  -4.31350797414780E-0001                     *
* 4.87437999428039E-0001  -1.26007354259491E+0000                     *
* 4.89950566435709E-0001  -7.07106769084930E-0001                     *
* 4.92463133442925E-0001   1.78893358460108E-0018                     *
* 4.94975700450141E-0001  -3.99903088808060E-0001                     *
* 4.97488267457811E-0001  -1.39680230617523E+0000                     *
* 5.00000834465027E-0001  -1.40504705905914E+0000                     *
*                                                                     *
* Output file 'output.txt' contains (lowpass - Fc=60):                *
*                                                                     *
* Filtered Signal (lowpass)                                           *
* 200                                                                 *
* 0.00000000000000E+0000   0.00000000000000E+0000                     * 
* 2.51256697811186E-0003   1.23316287994385E-0001                     *
* 5.02513395622373E-0003   4.30956214666367E-0001                     *
* 7.53770116716623E-0003   6.37625455856323E-0001                     *
* 1.00502679124475E-0002   5.75907468795776E-0001                     *
* 1.25628346577287E-0002   4.80039566755295E-0001                     *
* --/--                                                               *
* 4.87439036369324E-0001   3.87427806854248E-0001                     *
* 4.89951610565186E-0001  -1.84587817639112E-0002                     *
* 4.92464184761047E-0001  -1.86481356620789E-0001                     *
* 4.94976758956909E-0001  -1.93247571587563E-0001                     *
* 4.97489333152771E-0001  -3.17190080881119E-0001                     *
* 5.00001907348633E-0001  -5.91049015522003E-0001                     *
*                                                                     *
* Output file 'output1.txt' contains (highpass - Fc=60):              *
*                                                                     *
* Filtered Signal (highpass)                                          *
* 200                                                                 *
* 0.00000000000000E+0000   1.34853705763817E-0001                     * 
* 2.51256697811186E-0003  -2.01714970171452E-0002                     *
* 5.02513395622373E-0003  -4.67873126268387E-0001                     *
* 7.53770116716623E-0003  -5.71431457996368E-0001                     *
* 1.00502679124475E-0002  -1.53637155890465E-0001                     *
* 1.25628346577287E-0002   2.12560698390007E-0001                     *
* --/--                                                               *
* 4.87439036369324E-0001  -4.87992614507675E-0002                     *
* 4.89951610565186E-0001   3.67969155311584E-0001                     *
* 4.92464184761047E-0001   2.29071646928787E-0001                     *
* 4.94976758956909E-0001  -2.04981118440628E-0001                     *
* 4.97489333152771E-0001  -2.70452022552490E-0001                     *
* 5.00001907348633E-0001   1.68433189392090E-0001                     *
* ------------------------------------------------------------------- *
* References:                                                         *
*                                                                     *
*  http://en.wikipedia.org/wiki/Digital_biquad_filter                 *  
*  http://www.musicdsp.org/archive.php?classid=3#225                  *  
*  http://www.musicdsp.org/showone.php?id=197                         *
*  http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt                *
*  http://www.musicdsp.org/archive.php?classid=3#225                  *
*  http://www.musicdsp.org/showArchiveComment.php?ArchiveID=225       *
*                                                                     *
*                          Turbo-pascal Release By J-P Moreau, Paris. *
*                                     (www.jpmoreau.fr)               *
***********************************************************************
Note: Here, only the lowpass and highpass options are demonstrated.
------------------------------------------------------------------    }
                
Program Test_filter;

Uses WinCrt {Borland unit},  filter_r;


Var

  f1, f2: Text;
  lRGJFilter: TRbjEqFilter;  {see unit filter_r}
  i: integer;
  lIn: psingle;              {input table of single}

  title: String;

  kSamples: integer;         {number of points}
  kSamplesPerSec: single;    {sample rate}

  t,ti,tf,dt: single;        {t = current time}
                             {ti=starting time}
                             {tf=ending time}
  Fc: single;                {cut-off frequency}

BEGIN

  New(lIn);                  {allocate memory}

  ClrScr;
  Writeln;

  Fc := 60.0; 

  {read data from input text file}
  Assign(f1,'input.txt'); Reset(f1);
  Readln(f1,title);
  Readln(f1,kSamples);
  ReadLn(f1, ti, lIn^[0]);
  for i := 1 to kSamples-2 do
    ReadLn(f1, t, lIn^[i]);
  ReadLn(f1, tf, lIn^[kSamples-1]);
  Close(f1);

  dt := (tf-ti)/(kSamples-1);              {time increment}
  kSamplesPerSec := 1.0*kSamples/(tf-ti);  {points/sec}

  {apply filter and save results}
  lRGJFilter.Create(kSamplesPerSec, 0);

  lRGJFilter.CalcFilterCoeffs(kLowPass,Fc,0.3,0, False);

  Assign(f2,'output.txt'); ReWrite(f2);
  writeln(f2,'Filtered Signal (lowpass)');
  writeln(f2,kSamples);

  t:=ti;
  for i := 0 to kSamples-1 do
  begin
    WriteLn(f2, t, '  ', lRGJFilter.Process1(lIn^[i]));
    t := t + dt
  end;

  Close(f2);

  writeln(' Lowpass results in file output.txt...');

  lRGJFilter.CalcFilterCoeffs(kHighPass,Fc,0.3,0, False);

  Assign(f2,'output1.txt'); ReWrite(f2);
  writeln(f2,'Filtered Signal (highpass)');
  writeln(f2,kSamples);

  t:=ti;
  for i := 0 to kSamples-1 do
  begin
    WriteLn(f2, t, '  ', lRGJFilter.Process1(lIn^[i]));
    t := t + dt
  end;

  Close(f2);

  writeln;
  writeln(' Highpass results in file output1.txt...');

  Readkey;

  Dispose(lIn);   {free memory}

  DoneWinCrt

END.

{end of file test_fil.pas}