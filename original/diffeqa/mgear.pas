{***********************************************************************
*    Solve a first order Stiff System of Differential Equations using  *
*    the implicit Gear method of order 4.                              *
* -------------------------------------------------------------------  *
* Mode of operation:                                                   *
* ==================                                                   *
* This program solves can solve two of the 11 examples of file t_dgls  *
* using the implicit Gear method of order 4 (see file gear.pas).       *
* To test other systems of DEs, please proceed as explained in file    *
* t_dgls.pas.                                                          *
*                                                                      *
*   Inputs:                                                            *
*   =======                                                            *
*   bspnummer  Number of DE system from t_dgls.pas                     *
*   epsabs     desired absolute error bound                            *
*   epsrel     desired relative error bound                            *
*   x0         left edge of integration                                *
*   y0[0]   \  known approximation for the solution at x0              *
* ..  .      >                                                         *
*   y0[n-1] /                                                          *
*   h          initial step size                                       *
*   xend       right endpoint of integration                           *
*   fmax       maximal number of calls of the right hand side          *
*                                                                      *
*   The size n of the DE system is passed on from t_dgls.pas.          *
* -------------------------------------------------------------------- *
* SAMPLE RUN                                                           *
*                                                                      *
* Example #1:                                                          *
* (Solve set of differential equations (n=2):                          *
*     f[0] = y[0] * y[1] + COS(x) - HALF * SIN(TWO * x);               *
*     f[1] = y[0] * y[0] + y[1] * y[1] - (ONE + SIN(x));               *
*  Find values of f(0), f(1) at x=1.5).                                *
*                                                                      *
* Input example number (0 to 11): 0                                    *
* abs. epsilon: 1e-6                                                   *
* rel. epsilon: 1e-8                                                   *
* x0: 0                                                                *
* y0[0]: 0.5                                                           *
* y0[1]: 0.5                                                           *
* initial step size h: 0.0001                                          *
* right edge xend: 1.5                                                 *
* maximal number of calls of right hand side: 6000                     *
*                                                                      *
* Input data:                                                          *
* -----------                                                          *
* Example  =  0                                                        *
* n        =  2                                                        *
* x0       =  0.000000000000000e+0000                                  *
* xend     =  1.500000000000000e+0000                                  *
* epsabs   =  1.000000000000000e-0006                                  *
* epsrel   =  1.000000000000000e-0008                                  *
* fmax     = 6000                                                      *
* h        =  1.000000000000000e-0004                                  *
* y0[0]    =  5.000000000000000e-0001                                  *
* y0[1]    =  5.000000000000000e-0001                                  *
*                                                                      *
* Output data:                                                         *
* ------------                                                         *
* error code from gear4: 0                                             *
* final local step size:  6.06783655110639E-0002                       *
* number of calls of right hand size: 360                              *
* Integration stopped at x =  1.50000000000000E+0000                   *
*                                                                      *
* approximate solution y1(x) =  1.23598612893280E+0000                 *
* approximate solution y2(x) = -1.04949617987249E-0001                 *
*                                                                      *
* Example #2:                                                          *
* (Solve set of differential equations (n=5):                          *
*   f[0] = y[1];                                                       *
*   f[1] = y[2];                                                       *
*   f[2] = y[3];                                                       *
*   f[3] = y[4];                                                       *
*   f[4] = ((REAL)45.0 * y[2] * y[3] * y[4] -                          *
          (REAL)40.0 * y[3] * y[3] * y[3]) / (NINE * y[2] * y[2]);     *
*  Find values of f(0), ..., f(4) at x=1.5).                           *
*                                                                      *
* Input example number (0 to 11): 3                                    *
* abs. epsilon: 1e-10                                                  *
* rel. epsilon: 1e-10                                                  *
* x0: 0                                                                *
* y0[0]: 1                                                             *
* y0[1]: 1                                                             *
* y0[2]: 1                                                             *
* y0[3]: 1                                                             *
* y0[4]: 1                                                             *
* initial step size h: 0.001                                           *
* right edge xend: 1.5                                                 *
* maximal number of calls of right hand side: 6000                     *
*                                                                      *
* Input data:                                                          *
* -----------                                                          *
* Example  = 3                                                         *
* n        = 5                                                         *
* x0       =  0.00000000000000E+0000                                   *
* xend     =  1.50000000000000E+0000                                   *
* epsabs   =  1.00000000000000E-0010                                   *
* epsrel   =  1.00000000000000E-0010                                   *
* fmax     = 6000                                                      *
* h        =  1.00000000000000E-0003                                   *
* y0[0]    =  1.00000000000000E+0000                                   *
* y0[1]    =  1.00000000000000E+0000                                   *
* y0[2]    =  1.00000000000000E+0000                                   *
* y0[3]    =  1.00000000000000E+0000                                   *
* y0[4]    =  1.00000000000000E+0000                                   *
*                                                                      *
* Output data:                                                         *
* ------------                                                         *
* error code from gear4: 0                                             *
* final local step size:  4.86347662773039E-0003                       *
* number of calls of right hand size: 3423                             *
* Integration stopped at x =  1.50000000000000E+0000                   *
*                                                                      *
* approximate solution y1(x) =  4.36396102990278E+0000                 *
* approximate solution y2(x) =  4.00000000763432E+0000                 *
* approximate solution y3(x) =  2.82842715661993E+0000                 *
* approximate solution y4(x) =  4.86163052645233E-0008                 *
* approximate solution y5(x) = -3.77123622295570E+0000                 *
*                                                                      *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C, by Gisela Engeln-Muellges        *
*        and Frank Uhlig, Springer-Verlag, 1996".                      *
*                                                                      *
*                                 TPW Release By J-P Moreau, Paris.    *
*                                         (www.jpmoreau.fr)            *
***********************************************************************}
Program MGEAR;
Uses WinCrt1, FGauss, Gear, T_Dgls, Uawp;

Label fin;   {to exit program}
Var 
  epsabs,                { absolute error bound                      }
  epsrel,                { relative error bound                      }
  x0:     double;        { left edge of integration interval         }
  y0,                    { pointer to [0..n-1]-vector: initial values}
  yex:    pVEC;          { pointer to [0..n-1]-vector: exact solution}
  h,                     { initial, final step size                  }
  xend:   double;        { right edge of integration interval        }
  fmax,                  { maximal number of calls of right side     }
                         { in gear4()                                }
  aufrufe,               { actual number of function calls           }
 {bspnummer                Number of the system of DEs from t_dgls,
  n                        number of DEs in system (see t_dlgs.pas)  }
  fehler,                { error code from umleiten, gear4           }
  i:      integer;       { loop counter                              }

Begin

  New(y0); New(yex);

{ -------------------- read input  -------------------- }
  writeln;
  write(' Input example number (0 to 11): ');
  readln(bspnummer);
  if Not(bspnummer in [0,1,3]) then
  begin
    writeln;
    writeln(' Example not registered !');
    goto fin
  end;

  write(' abs. epsilon: '); readln(epsabs);
  write(' rel. epsilon: '); readln(epsrel);
  write(' x0: '); readln(x0);

  dgl(x0,y0,yex);   {dummy call to get example's n value}

  for i := 0 to n-1 do
  begin
    write(' y0[',i,']: '); readln(y0^[i])
  end;

  write(' initial step size h: '); readln(h);
  write(' right edge xend: '); readln(xend);
  write(' maximal number of calls of right hand side: ');
  readln(fmax);

  { ----------------- print input data -------------------- }
  writeln;
  writeln(' ==============================================');
  writeln('   Solve a first order ordinary system of DEs  ');
  writeln(' using the implicit method of Gear of 4th order');
  writeln(' ==============================================');
  writeln;
  writeln(' System of DEs:');
  writeln(' ------------- ');
  for i:=0 to n-1 do writeln(dgltxt[i]);
  writeln(' Example n° ', bspnummer);
  writeln(' Number of DEs = ', n);
  writeln(' x0     = ', x0);
  writeln(' xend   = ', xend);
  writeln(' epsabs = ', epsabs);
  writeln(' epsrel = ', epsrel);
  writeln(' fmax   = ', fmax);
  writeln(' h      = ', h);
  for i:=0 to n-1 do writeln(' y0[',i,']  = ', y0^[i]);

  { ----------------- Solve system of DEs ------------------ }
    gear4(x0, xend, n, y0, epsabs, epsrel, h, fmax, aufrufe, fehler);

  if fehler <> 0 then  {something went wrong}
  begin 
    writeln(' Gear4: error n° ', 10+fehler);
    goto fin
  end;

  { ---------------------- print results ------------------- }
  writeln;
  writeln(' Output data:');
  writeln(' -----------');
  writeln(' error code from gear4: ', fehler); 
  writeln(' final local step size: ', h);
  writeln(' number of calls of right hand side: ', aufrufe);
  writeln(' Integration stopped at x = ', x0);

  for i:=0 to n-1 do
    writeln(' approximate solution y',i+1,'(x) = ', y0^[i]);

fin: Dispose(y0); Dispose(yex);

  readkey;
  DoneWinCrt

End.

{ -------------------------- END mgear.pas -