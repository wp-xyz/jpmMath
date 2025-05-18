{*********************************************************
*    This program calculates the Fourier coefficients    *
*   of a periodic function F(x) by using the procedure   *
*   AnalyticFourierHn of unit Fourier.                   *
* ------------------------------------------------------ *
* SAMPLE RUN:                                            *
*                                                        *
* Example 1: F(x) has the following shape:               *
*                                                        *
*  1 ---- 2  pi/4                                        *
*    !  !                                                *
*    !  !                                                *
* ------!-------------------------------------------> x  *
*  -pi  !    0         ! pi                              *
*       !              !                                 *
*     3 ---------------- 4 -pi/4                         *
*                                                        *
* Calculate the Fourier coefficients of a periodic       *
* function F(x):                                         *
*                                                        *
* Begin x of period: -3.1415927                          *
* End of x period: 3.1415927                             *
*                                                        *
* Lowest  harmonic: 0                                    *
* Highest harmonic: 10                                   *
*                                                        *
* a0 = -0.392684481                                      *
* b0 =  0.000000000                                      *
*                                                        *
* a1 = -0.500000000                                      *
* b1 = -0.500029201                                      *
*                                                        *
* a2 = -0.000029201                                      *
* b2 = 0.5000000000                                      *
* -----------------                                      *
* a9 = -0.055555556                                      *
* b9 = -0.055584756                                      *
*                                                        *
* a10 = -0.000029201                                     *
* b10 =  0.100000000                                     *
*                                                        *
* Note that for this kind of sharp angles function, the  *
*      method for discreet functions is far better.      *
*                                                        *
* Example 2: F1(x) equals 0 from -pi to 0 and sin(x)     *
*            from 0 to pi.                               *
*                                                        *
* Calculate the Fourier coefficients of a periodic       *
* function F(x):                                         *
*                                                        *
* Begin x of period: -3.1415927                          *
* End of x period: 3.1415927                             *
*                                                        *
* Lowest  harmonic: 0                                    *
* Highest harmonic: 10                                   *
*                                                        *
* a0 =  0.318309881                                      *
* b0 =  0.000000000                                      *
*                                                        *
* a1 =  0.000000012                                      *
* b1 =  0.499999996                                      *
*                                                        *
* a2 = -0.212206596                                      *
* b2 =  0.000000010                                      *
* -----------------                                      *
* a9 = -0.000000000                                      *
* b9 = -0.000000002                                      *
*                                                        *
* a10 = -0.006430503                                     *
* b10 =  0.000000001                                     *
*                                                        *
*                     TPW version By J-P Moreau, Paris.  *
*                            (www.jpmoreau.fr)           *
*********************************************************}
Program Test_AnalyticFourier;
Uses WinCrt, Fourier;

Var
    count,h1,h2,i: Integer;
    x1,x2: Double;
    ai,bi: Double;

BEGIN
  writeln;
  writeln(' Calculate the Fourier coefficients of a periodic function F(x):');
  writeln('    (Function F(x) must be defined in file fourier.pas).');
  writeln;
  write(' Begin x of period: '); readln(x1);
  write(' End x of period: '); readln(x2);
  writeln;
  write(' Lowest  harmonic: '); readln(h1);
  write(' Highest harmonic: '); readln(h2);
  count:=1;
  for i:=h1 to h2 do
  begin
    AnalyticFourierHn(x1,x2,i,ai,bi);
    writeln;
    writeln('  a',i,' = ',ai:12:9);
    writeln('  b',i,' = ',bi:12:9);
    if count MOD 5 =0  then
    begin
      count:=0;
      readkey
    end;
    Inc(count)
  end;
  writeln;
  ReadKey;
  DoneWinCrt
END.

{end of file analfour.pas}