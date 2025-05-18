{***********************************************************************
* This program tests the procedure brown1 which is designed to solve   *
* a nonlinear system of equations                                      *
*             f0  (x0,x1,...,xn-1)  =  0                               *
*             f1  (x0,x1,...,xn-1)  =  0                               *
*             ...                                                      *
*             fn-1(x0,x1,...,xn-1)  =  0                               *
* using Brown's method (see unit Brown.pas).                           *
*                                                                      *
* Scope of program:                                                    *
* ================                                                     *
* The program reads the input from stdin and writes output onto file   *
* brownts1.lst.                                                        *
*                                                                      *
* After reading, the input is printed out for control purposes. Then   *
* Brown's method is performed and, barring premature errors, the end   *
* results are printed.                                                 *
*                                                                      *
* The system to solve is defined in unit Brown.pas.                    *
* -------------------------------------------------------------------- *
* SAMPLE RUN:                                                          *
*                                                                      *
* Desired precision : 1e-8                                             *
* Protocole (0 or 1): 0                                                *
* Maximum number of iterations : 25                                    *
* Components of starting vector: 1 1 1 1                               *
*                                                                      *
* The output file brownts1.lst contains:                               *
*                                                                      *
* -------------------------------------------------------------------- *
* Brown's method for nonlinear systems of equations                    *
* -------------------------------------------------------------------- *
* Example n° 3                                                         *
* System to be solved:                                                 *
* 10*x0 + x1 + x2 + x3 - 20 + sin(x0)^2 + cos(x1)^2 = 0                *
* x0 + 20*x1 + x2 + x3 - 48 + 1/(x[0]^6)            = 0                *
* (x0+x1)^2 + 30*x2 + x3 - 97 + Ln(x0) + Ln(x1+x2)  = 0                *
*  x0 + x1 + x2 + 40*x3 -166 + x0*x0                = 0                *
* Starting vector:                                                     *
*    1.000    1.000    1.000    1.000                                  *
* Error bound =  1.00000000000000E-0008                                *
* Maximal number of iterations = 25                                    *
*                                                                      *
* Intermediate results are not kept.                                   *
* Solution vector:                                                     *
*  1.040648  1.972398  2.745049  3.978974                              *
* Number of iterations: 4                                              *
* -------------------------------------------------------------------- *
* Ref.: "Numerical Algorithms with C By G. Engeln-Mueller and F. Uhlig,*
*        Springer-Verlag, 1996" [BIBLI 11].                            *
*                                                                      *
*                                  TPW Version By J-P Moreau, Paris.   *
*                                         (www.jpmoreau.fr)            *
***********************************************************************}
Program Test_Brown;
Uses WinCrt, Basis, Brown;

Label fin;
Var
     n,             { size of system                             }
     maxit,         { maximal number of iterations               }
     error,         { error code of procedure brown1             }
     itanz,         { number of iterations performed             }
     i,             { Loop variable                              }
     prot:Integer;  { protocole flag                             }
     eps: REAL;     { desired accuracy                           }
     x0,            { [0..n-1] starting vector                   }
     x1: VEC;       { [0..n-1] approximate solution              }
     fp: TEXT;      { Output text file                           }    
BEGIN  
  { -------------------- read input -------------------------------- }
  writeln;
  n := 4;           {size of system}
  write(' Desired precision : '); readln(eps);
  writeln;
  write(' Protocole (0 or 1): '); readln(prot);
  writeln;
  write(' Maximal number of iterations : '); readln(maxit);
  writeln;    
  write(' Componants of starting vector: ');
  for i := 0 to n-1 do read(x0[i]);
  writeln;
  { ------------ print input for checking purposes ----------------- }
  Assign(fp,'brownts1.lst'); Rewrite(fp);
  WriteHead(fp,' Brown''s method for nonlinear systems of equations');
  Writeln(fp,' Example n° 3'); 
  Writeln(fp,' System to be solved:');
  Writeln(fp,' 10*x0 + x1 + x2 + x3 - 20 + sin(x0)^2 + cos(x1)^2 = 0');
  Writeln(fp,' x0 + 20*x1 + x2 + x3 - 48 + 1/(x[0]^6)            = 0');
  Writeln(fp,' (x0+x1)^2 + 30*x2 + x3 - 97 + Ln(x0) + Ln(x1+x2)  = 0');
  Writeln(fp,'  x0 + x1 + x2 + 40*x3 -166 + x0*x0                = 0');
  Writeln(fp,' Starting vector:');
  for i := 0 to n-1 do write(fp,x0[i]:9:3);
  Writeln(fp);
  Writeln(fp,' Error bound = ', eps);
  Writeln(fp,' Maximal number of iterations = ', maxit);
  Writeln(fp);
  if prot<>0 then
    Writeln(fp,' Intermediate results are saved.')
  else
    Writeln(fp,' Intermediate results are not kept.');

  { ------------ solve nonlinear system ---------------------------- }
  Brown1(n, x0, eps, prot, maxit, x1, itanz, error);

  Case error of
    0: goto fin;
    1: writeln(' brown: too many steps');
    2: writeln(' brown: linearized system singular');
    3: writeln(' brown: lack of memory');
    4: writeln(' brown: wrong input parameter: n < 1 or maxit < 1');
    5: writeln(' brown: error calling function')
  End;

  { ---------------------- print solution --------------------------- }
fin:Writeln(fp,' Solution vector:');
  for i := 0 to n-1 do f_aff_reel(fp,x1[i]);
  Writeln(fp);
  Writeln(fp,' Number of iterations: ', itanz);   
  WriteEnd(fp);
  Close(fp);
  
  writeln(' Results in file brownts1.lst.');
  writeln;
  ReadKey; DoneWinCrt
END. 

{ TPW version by J-P Moreau, Paris }