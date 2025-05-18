{********************************************************************
*   Find a real root of a real function F(x) by the Zeroin method   *
* ----------------------------------------------------------------- *
* SAMPLE RUN:                                                       *
* (Find a root of function 5x - exp(x) between 0 and 1)             *
*                                                                   *
* The output file tzeroin.lst contains:                             *
*                                                                   *
* ----------------------------------------------------------------- *
* Zeroin method for real valued, nonlinear functions                *
*------------------------------------------------------------------ *
* Find a real root of function:                                     *
* F(x) = 5x - Exp(x)                                                *
* Starting value x0 =  0.00000000000000E+0000                       * 
* Return code       = 0                                             *
* Root              =     0.25917110181899                          *
* Function value    =    -0.00000000000030                          *
* Function calls    = 8                                             *
* Absolute error    =     0.00000000010000                          *
* Relative error    =     0.00000000100000                          *
*------------------------------------------------------------------ *
* Ref.: "Numerical Algorithms with C  By G. Engeln-Mueller and      *
*        F. Uhlig, Springer-Verlag, 1996" [BIBLI 11].               *
*                                                                   *
*                                TPW Version By J-P Moreau, Paris.  *
*                                        (www.jpmoreau.fr)          *
********************************************************************}
Uses WinCrt, Basis, Fzeroin;
          
VAR
    i,rc,nmax,niter: Integer;
    abserr, relerr, x1, x2, f: REAL;
    fp: TEXT;

BEGIN

  {the function Fkt(x) is defined in unit Fzeroin.pas}

  abserr:=100.0*MACH_EPS;
  relerr:=1000.0*MACH_EPS;

  Assign(fp, 'tzeroin.lst'); Rewrite(fp);    {open output file}

  WriteHead(fp,' Zeroin method for real valued, nonlinear functions'); 

  nmax := 100; x1:=ZERO; x2:=ONE;

  writeln(fp,' Find a real root of function:');
  writeln(fp,' F(x) = 5x - Exp(x)');
  writeln(fp,' Starting value x0 = ', x1);

  zeroin(abserr, relerr, nmax, 'tzeroin.log', x1, x2, f, niter, rc);

  writeln(fp,' Return code       = ', rc);
  writeln(fp,' Root              = ', x2:20:14);
  writeln(fp,' Function value    = ', f:20:14);
  writeln(fp,' Function calls    = ', niter);
  writeln(fp,' Absolute error    = ', abserr:20:14); 
  writeln(fp,' Relative error    = ', relerr:20:14); 

  WriteEnd(fp);
    
  close(fp);

  writeln;
  writeln(' Results in tzeroin.lst.');
  Readkey;
  DoneWinCrt

END.

{end of file tzeroin.pas}