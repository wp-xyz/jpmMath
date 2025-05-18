{******************************************************************************
*       Test program for pegasus (6 functions defined in fonction.pas)        *
* --------------------------------------------------------------------------- *
* SAMPLE RUN:                                                                 *
*                                                                             *
* File tpegasus.lst contains:                                                 *
*                                                                             *
* --------------------------------------------------------------------------- *
* Pegasus method for searching real roots of real valued, nonlinear functions *
* --------------------------------------------------------------------------- *
*                                                                             *
* f(x) = 5 * x - exp (x)                                                      *
* Starting value x0 =  0.00000000000000E+0000                                 *
* Return code       = 1                                                       *
* Root              =  2.59171101819074E-0001                                 *
* Function value    =  8.53267109746092E-0017                                 *
* Iterations        = 7                                                       *
*                                                                             *
* f(x) = (((((x-6)*x+15)*x-20)*x+15)*x-6)*x+1                                 *
* Starting value x0 =  5.00000000000000E-0001                                 *
* Return code       = -1                                                      *
* Root              =  1.50000000000000E+0000                                 *
* Function value    =  1.56250000000000E-0002                                 *
* Iterations        = 0                                                       *
*                                                                             *
* f(x) = sin (x)                                                              *
* Starting value x0 =  3.00000000000000E+0000                                 *
* Return code       = 1                                                       *
* Root              =  3.14159265358979E+0000                                 *
* Function value    =  1.22514845490862E-0016                                 *
* Iterations        = 6                                                       *
*                                                                             *
* f(x) = 1 + sin (x)                                                          *
* Starting value x0 = -2.00000000000000E+0000                                 *
* Return code       = -1                                                      *
* Root              = -1.50000000000000E+0000                                 *
* Function value    =  2.50501339594557E-0003                                 *
* Iterations        = 0                                                       *
*                                                                             *
* f(x) = exp(x) -(1.0 + x + x*x*0.5)                                          *
* Starting value x0 =  2.00000000000000E+0000                                 *
* Return code       = -1                                                      *
* Root              =  3.00000000000000E+0000                                 *
* Function value    =  1.15855369231877E+0001                                 *
* Iterations        = 0                                                       *
*                                                                             *
* f(x) = (x-1.0)*(x-1.0)*( sin(PI*x) - log(abs(2.0*x/(x+1.0)))                *
* Starting value x0 =  2.00000000000000E+0000                                 *
* Return code       = -1                                                      *
* Root              =  3.00000000000000E+0000                                 *
* Function value    = -1.62186043243266E+0000                                 *
* Iterations        = 0                                                       *
* --------------------------------------------------------------------------- *
*                                                                             *
*---------------------------------------------------------------------------- *
* Ref.: "Numerical algorithms with C, By Gisela Engeln-Muellges and Frank     *
*        Uhlig, Springer-Verlag, 1996" [BIBLI 11].                            *
*                                                                             *
*                                           TPW Release By J-P Moreau, Paris. *
*                                                  (www.jpmoreau.fr)          *
******************************************************************************}
Program Test_Pegasus;

Uses WinCrt, Fonction, FPegasus;   

Var

    i, rc, iter: Integer;
    x1,x2, f: Double;
    texte: String;
    fp: TEXT;

BEGIN

  Assign(fp, 'tpegasus.lst');  Rewrite(fp);    {output file}

  Writeln(fp,' ---------------------------------------------------------------------------');
  Writeln(fp,' Pegasus method for searching real roots of real valued, nonlinear functions');
  Writeln(fp,' ---------------------------------------------------------------------------');

  for i := 1 to 6 do
  begin

    Case i of             {select suitable example}
      1: begin
           x1 := 0.0;
           x2 := 1.0;
           NumFunc := 1;
           texte := 'f(x) = 5 * x - exp (x)'
         end;
      2: begin
           x1 := 0.5;
	   x2 := 1.5;
           NumFunc := 2;
           texte := 'f(x) = (((((x-6)*x+15)*x-20)*x+15)*x-6)*x+1'
         end;
      3: begin
           x1 := 3.0;
           x2 := 3.5;
           NumFunc := 3;
           texte := 'f(x) = sin (x)'
         end;
      4: begin
           x1 := -2.0;
           x2 := -1.5;
           NumFunc := 4;
           texte := 'f(x) = 1 + sin (x)'
         end;
      5: begin
           x1 := 2.0;
           x2 := 3.0; 
           NumFunc := 5;
           texte := 'f(x) = exp(x) -(1.0 + x + x*x*0.5)'
         end;
      6: begin
           x1 := 2.0;
           x2 := 3.0;
           NumFunc := 6;
           texte := 'f(x) = (x-1.0)*(x-1.0)*( sin(PI*x) - log(abs(2.0*x/(x+1.0)))'
         end
    end;

    writeln(fp);
    writeln(fp,' ', texte);
    writeln(fp,' Starting value x0 = ',x1);

    rc := pegasus (x1, x2, f, iter);

    writeln(fp,' Return code       = ', rc);
    writeln(fp,' Root              = ', x2);
    writeln(fp,' Function value    = ', f);
    writeln(fp,' Iterations        = ', iter);
  end; 
  Writeln(fp,' ---------------------------------------------------------------------------');
  close(fp);

  writeln;
  writeln(' Results in tpegasus.lst.');
  ReadKey;
  DoneWinCrt

END.

{end of file tpegasus.pas}