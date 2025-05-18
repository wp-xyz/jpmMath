{*******************************************************************
* This Program tests Module STATS for basic statistical Functions  *
* ---------------------------------------------------------------- *
* SAMPLE RUN:                                                      *
*                                                                  *
*  Population mean and standard deviation:   0.075133   0.879506   *
*                                                                  *
*  Sample mean and standard deviation:   0.075133   0.883904       *
*                                                                  *
*  For full regression:                                            *
*               regression coefficients:   1.006200   1.986324     *
*  standard error of estimate of y on x:   9.918043                *
*               correlation coefficient:   0.985044                *
*                                                                  *
*  For regression forced through (0,0):                            *
*               regression coefficients:   0.000000   2.006249     *
*  standard error of estimate of y on x:   9.935045                *
*               correlation coefficient:   0.984993                *
*                                                                  *
* ---------------------------------------------------------------- *
* Ref.: "Problem Solving with Fortran 90 By David R.Brooks,        *
*        Springer-Verlag New York, 1997".                          *
*                                                                  *
*                           Pascal Release By J-P Moreau, Paris.   *
*                                   (www.jpmoreau.fr)              *
*******************************************************************}
Program Tstats;

Uses WinCrt1, Stats;

Var
    x, y: pVec;
    avg, std_dev,
    a, b,          {for linear regression y=ax+b}
    s_yx,          {standard error of estimate of y on x}
    corr: Real;    {correlation coefficient}
    n: Integer;    {# of points in array}
    i: Integer;

Begin

  Randomize;

  n := 100;

  New(x); New(y);

{ Test basic statistics }
  NormalArray(x,n);

  NormalStats(x,n,'p',avg,std_dev);
  
  writeln;
  writeln(' Population mean and standard deviation: ',avg:10:6,' ',std_dev:10:6);
  
  NormalStats(x,n,'s',avg,std_dev);       

  writeln;
  writeln(' Sample mean and standard deviation: ',avg:10:6,' ',std_dev:10:6);
  writeln;

{ Test linear regression }
  NormalArray(y,n);

{ Create a linear relationship with 'noise' }
  For i:=1 to n do
  begin
    x^[i] := 1.0*i;
    y^[i] := 2.0*i + 10.0*y^[i]
  end;
  
{ Set a <> zero for full regression analysis }
  a:=1.0;
  LinearReg(x,y,n,'s',a,b,s_yx,corr);
  
  writeln(' For full regression:');
  writeln('              regression coefficients: ',a:10:6,' ',b:10:6);
  writeln(' standard error of estimate of y on x: ',s_yx:10:6);
  writeln('              correlation coefficient: ',corr:10:6);  	    
  writeln;

{ Set a = zero for full regression forced through (0,0) }
  a:=0.0;
  LinearReg(x,y,n,'s',a,b,s_yx,corr);
  
  writeln(' For regression forced through (0,0):');
  writeln('              regression coefficients: ',a:10:6,' ',b:10:6);
  writeln(' standard error of estimate of y on x: ',s_yx:10:6);
  writeln('              correlation coefficient: ',corr:10:6);
  writeln;

  Readkey;
  Dispose(x); Dispose(y);
  DoneWinCrt

End.

{end of file tsats.pas}