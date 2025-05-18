{*********************************************************
* This program calculates the Fourier coefficients of a  *
* periodic discrete function F(x) by using the procedure *
* DiscreetFourierHn of unit Fourier.                     *
* ------------------------------------------------------ *
* SAMPLE RUN:                                            *
*                                                        *
* F(x) has the following shape:                          *
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
* discreet function:                                     *
*                                                        *
* Number of points:  4                                   *
*                                                        *
*  x1 = -3.1415927                                       *
*  y1 = 0.7853982                                        *
*                                                        *
*  x2 = -1.5707963                                       *
*  y2 = 0.7853982                                        *
*                                                        *
*  x3 = -1.5707963                                       *
*  y3 = -0.7853982                                       *
*                                                        *
*  x4 = 3.1415927                                        *
*  y4 = -0.7853982                                       *
*                                                        *
* Lowest  harmonic: 0                                    *
* Highest harmonic: 10                                   *
*                                                        *
* a0 = -0.392699087                                      *
* b0 =  0.000000000                                      *
*                                                        *
* a1 = -0.500000023                                      *
* b1 = -0.500000048                                      *
*                                                        *
* a2 = 0.0000000025                                      *
* b2 = 0.5000000023                                      *
* -----------------                                      *
* a9 = -0.055555558                                      *
* b9 = -0.055555583                                      *
*                                                        *
* a10 = -0.000000025                                     *
* b10 =  0.100000005                                     *
*                                                        *
*                     TPW version by J-P Moreau, Paris.  *
*                            (www.jpmoreau.fr)           *
*********************************************************}
Program Test_DiscreetFourier;
Uses WinCrt, Fourier;

Var
    count,h1,h2,i,npoints: Integer;
    Xi,Yi: pTab;
    ai,bi: REAL;

BEGIN
  New(Xi); New(Yi);
  writeln;
  writeln(' Calculate the Fourier coefficients of a periodic discrete function:');
  writeln;
  write(' Number of points: '); readln(npoints);
  writeln;
  for i:=1 to npoints do
  begin
    write('  x',i,'= '); readln(Xi^[i]);
    write('  y',i,'= '); readln(Yi^[i]);
    writeln
  end;
  write(' Lowest  harmonic: '); readln(h1);
  write(' Highest harmonic: '); readln(h2);
  count:=1;
  for i:=h1 to h2 do
  begin
    DiscreetFourierHn(npoints,Xi,Yi,i,ai,bi);
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
  Dispose(Xi); Dispose(Yi);
  DoneWinCrt
END.

{end of file discfour.pas}