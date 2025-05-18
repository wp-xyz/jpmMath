{***********************************************************************
* This program solves a nonlinear system of two variables:             *
*                                                                      *
*             f(x,y)  =  0                                             *
*             g(x,y)  =  0                                             *
*                                                                      *
* using the Non_linear_system procedure.                               *
* -------------------------------------------------------------------- *
* REFERENCE: "Mathematiques en Turbo-Pascal part 1, by M. Ducamp and   *
*             A. Reverchon, Editions EYROLLES, Paris, 1991" [BIBLI 03] *
* -------------------------------------------------------------------- *
* SAMPLE RUN:                                                          *
*                  (Solve the non linear system                        *
*                  x^2 + x + y^2 - 2         = 0                       *
*                  x^2 + y - y^2 - 1 + ln(x) = 0 )                     *
*                                                                      *
* Desired precision........... : 1e-6                                  *
* Maximal number of iterations : 100                                   *
* Componants of starting vector:  1  1                                 *
*                                                                      *
* Results in file nlinsyst.lst:                                        *
* (Error status: 0)                                                    *
*                                                                      *
* -------------------------------------------------------------------- *
*  Method for nonlinear systems of two equations                       *
* -------------------------------------------------------------------- *
*                                                                      *
*  System to be solved:                                                *
*  x^2 + x + y^2 - 2          =  0                                     *
*  x^2 + y - y^2 - 1 + ln(x)  =  0                                     *
*                                                                      *
*  Starting vector:                                                    *
*  1.00000000000000E+0000   1.00000000000000E+0000                     *
*                                                                      *
*  Error bound =  1.00000000000000E-0006                               *
*  Maximal number of iterations = 100                                  *
*                                                                      *
*                                                                      *
*  Solution vector:                                                    *
*  9.15554449245427E-0001  4.96191093452388E-0001                      *
*                                                                      *
*  Number of iterations: 6                                             *
* -------------------------------------------------------------------- *
*                                                                      *
***********************************************************************}
Program Test_NLSYST;
Uses WinCrt,Basis_R,Type_def;

Var

          n,             { size of system                            }
          maxit,         { maximal number of iterations              }
          it:INTEGER;            { number of iterations performed            }

          eps:REAL_AR;   { desired accuracy                          }

          x0,y0,         { starting point                            }
          x,y:REAL_AR;   { approximate solution                      }

          fp:TEXT;       { Output file                               }

          error:INTEGER; { error code: 0 = OK, 1 = error in evluating}
                         { f(x,y) or g(x,y),   2 = singular system   }

          Function f(x,y:REAL_AR): REAL_AR;
          Begin
            f:=x*x+x+y*y-2.0
          End;

          Function g(x,y:REAL_AR;VAR rc:INTEGER): REAL_AR;
          Begin
            rc:=0;
            if x>0.0 then
              g:=x*x+y-y*y-1.0+LN(x)
            else
              rc:=1
          End;


PROCEDURE Non_linear_system( x0,y0:REAL_AR;prec:REAL_AR;maxiter:INTEGER;
                             VAR x,y:REAL_AR; VAR iter:INTEGER);
Const h = 0.01;
Var   rc : INTEGER;
      a,b,c,d,t,m,n,p,q : REAL_AR;
Begin
  error:=1; iter:=0;
  x:=x0; y:=y0;
  Repeat
    Inc(iter);
    if iter > maxiter then exit;
    a:=f(x+h,y);
    b:=g(x+h,y,rc); if rc<>0 then exit;
    a:=(a-f(x-h,y))/2.0/h; 
    b:=(b-g(x-h,y,rc))/2.0/h; if rc<>0 then exit;
    c:=f(x,y+h);
    d:=g(x,y+h,rc); if rc<>0 then exit;
    c:=(c-f(x,y-h))/2.0/h; 
    d:=(d-g(x,y-h,rc))/2.0/h; if rc<>0 then exit;
    t:=a*d-b*c;
    if ABS(t)<1e-12 then
    begin
      error:=2;
      exit
    end;
    m:=f(x,y); 
    n:=g(x,y,rc); if rc<>0 then exit;
    p:=(m*d-n*c)/t;
    q:=(n*a-m*b)/t;
    x:=x-p; y:=y-q
  Until (ABS(p)+ABS(q)) < prec;
  error:=0
End;


{main program}
BEGIN            
  { -------------------- read input -------------------------------- }
  n:=2;   {size of system}
  writeln;
  write(' Desired precision..: '); readln(eps);
  write(' Maximal number of iterations : '); readln(maxit);
  write(' Componants of starting vector:');
  write('  '); read(x0); write('  '); readln(y0);

  { ------------ print input for checking purposes ----------------- }
  Assign(fp,'nlinsyst.lst'); Rewrite(fp);
  WriteHead(fp,' Method for nonlinear systems of two equations');
  writeln(fp);
  writeln(fp,' System to be solved:');
  writeln(fp,' x^2 + x + y^2 - 2          =  0');
  writeln(fp,' x^2 + y - y^2 - 1 + ln(x)  =  0');
  writeln(fp);
  writeln(fp,' Starting vector:');
  write(fp,' ',x0,'  ',y0);
  writeln(fp);
  writeln(fp);
  writeln(fp,' Error bound = ',eps);
  writeln(fp,' Maximal number of iterations = ',maxit);
  writeln(fp);

  { ------------ solve nonlinear system ---------------------------- }
  Non_linear_system(x0,y0,eps,maxit,x,y,it);

  { --------------------- print solution --------------------------- }
  writeln(fp);
  writeln(fp,' Solution vector:');
  writeln(fp,' ',x,' ',y);
  writeln(fp);
  writeln(fp,' Number of iterations: ', it);
  WriteEnd(fp);
  close(fp);

  writeln;
  writeln(' Results in file nlinsyst.lst.');
  writeln(' Error status: ',error);
  Readkey; DoneWinCrt;

END.

{ End of file nlinsyst.pas}