{***********************************************************
*  This program calculates the statistical moments of a    *
*  distribution: Mean, Variance, Skewness, etc.            *
* -------------------------------------------------------- *
* Ref.: "Numerical Recipes, by W.H. Press, B.P. Flannery,  *
*        S.A. Teukolsky and T. Vetterling, Cambridge       *
*        University Press, 1986".                          *
*                                                          *
*                    Pascal Release By J-P Moreau, Paris.  *
* -------------------------------------------------------- *
* SAMPLE RUN:                                              *
*                                                          *
* Number of data: 5                                        *
*   1: 12                                                  *
*   2: 9                                                   *
*   3: 7                                                   *
*   4: 15                                                  *
*   5: 6                                                   *
*                                                          *
* Average ..........:   9.80000000000000E+0000             *
* Average  Deviation:   2.96000000000000E+0000             *
* Standard Deviation:   3.70135110466435E+0000             *
* Variance .........:   13.7000000000000E+0000             *
* Skewness .........:   0.29154869588874E+0000             *
* Kurtosis .........:  -1.90779903031595E+0000             *
*                                                          *  
***********************************************************}
Program TMOMENT;
Uses WinCrt;

Const NMAX = 1024;
Type  pVEC = ^VEC;
       VEC = Array[1..NMAX] of double;
Var
    i,ndata: integer;
    Y: pVEC;
    average,avdev,stdev,variance,skewness,kurtosis: Double;


{***********************************************************
* Given an array of Data of length n, this routine returns * 
* its mean ave, average deviation adev, standard deviation *
* sdev, variance var, skewness skew, and kurtosis curt.    *
***********************************************************}
Procedure Moment(data:pVEC; n:Integer; Var ave,adev,sdev,var0,skew,curt:Double);
Label return;
Var j: integer; p,s: Double;
Begin
  if n <= 1 then
  begin
    writeln(' N must be at least 2!');
    goto return
  end;
  s:=0.0;
  For j:=1 to n do  s:=s+data^[j];
{ calculate mean }
  ave:=s/n;
  adev:=0.0;
  var0:=0.0;
  skew:=0.0;
  curt:=0.0;
  For j:=1 to n do
  begin
    s:=data^[j]-ave;
    adev:=adev+abs(s);
    p:=s*s;
    var0:=var0+p;
    p:=p*s;
    skew:=skew+p;
    p:=p*s;
    curt:=curt+p
  end;
  adev:=adev/n;
  var0:=var0/(n-1);
  sdev:=sqrt(var0);
  if var0 <> 0.0 then
  begin
    skew:=skew/(n*sdev*sdev*sdev);
    curt:=curt/(n*var0*var0)-3.0
  end
  else
    writeln(' No skew or kurtosis when zero variance.');
return: End;


{main program}
BEGIN

  New(Y);
  writeln;
  write(' Number of data: '); readln(ndata);

  For i:=1 to ndata do
  begin
    write('  ',i,': '); readln(Y^[i])
  end;

  Moment(Y,ndata,average,avdev,stdev,variance,skewness,kurtosis);

{ print results }
  writeln;
  writeln(' Average ..........: ', average);
  writeln(' Average  Deviation: ', avdev);
  writeln(' Standard Deviation: ', stdev);
  writeln(' Variance .........: ', variance);
  writeln(' Skewness .........: ', skewness);
  writeln(' Kurtosis .........: ', kurtosis);
  writeln;
  ReadKey;
  DoneWinCrt

END.

{end of file tmoment.pas}