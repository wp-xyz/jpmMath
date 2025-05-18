{***********************************************************************
*     Test program: Numerical differentation according to  Romberg     *
* -------------------------------------------------------------------- *
* SAMPLE RUN:                                                          *
* (Find the first derivative of function f(x)=1/x, for x=0.12).        *
*                                                                      *
* Numerical differentation according to Romberg                        *
*                                                                      *
* Test function   f(x) = 1/x                                           *
*                                                                      *
* Put in x-value at which you want to evaluate derivative: 0.12        *
* Put in desired accuracy: 1e-10                                       *
* Maximal number of columns in Romberg scheme: 4                       *
* Starting step size: 0.005                                            *
*                                                                      *
* Results:                                                             *
*                                                                      *
* er_app=  8.89315288077341E-0011                                      *
* res   = -6.94444444444433E+0001                                      *
* nend  =  3                                                           *
* hend  =  6.25000000000000E-0004                                      *
*                                                                      *
* -------------------------------------------------------------------- *
* "Numerical Algorithms with C,  By Gisela Engeln-Muellges             *
*  and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].                 *
*                                                                      *
*                                TPW Release 1.0 By J-P Moreau, Paris. *
*                                         (www.jpmoreau.fr)            *
***********************************************************************}
Program Test_Difrom;

Uses WinCrt, Difrom;

Var
    error, n, nend: Integer;
    x,prec,h: Double;
    res,schaetz, hend: Double;


{print results}
Procedure aus(error:Integer; schaetz:Double; res: Double; nend:Integer; hend:Double);
Begin
  if error <> 1 then
  begin
    writeln;
    writeln(' Results:');
    writeln;
    writeln(' er_app= ', schaetz);
    writeln(' res   = ', res);
    writeln(' nend  = ', nend:2);
    writeln(' hend  = ', hend)
  end
End;

{main program}
BEGIN

  writeln;
  writeln(' Numerical differentation according to Romberg');
  writeln;
  writeln(' Test function   f(x) = 1/x');
  writeln;

  write(' Put in x-value at which you want to evaluate derivative: ');
  readln(x);

  write(' Put in desired accuracy: '); readln(prec);

  write(' Maximal number of columns in Romberg scheme: '); readln(n);

  write(' Starting step size: '); readln(h);

  difrom1(x,prec,n,h,res,schaetz,nend,hend,error);

  aus(error,schaetz,res,nend,hend);

  ReadKey;
  DoneWinCrt

END.

{end of file tdifrom.pas}