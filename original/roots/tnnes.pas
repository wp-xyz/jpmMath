{*******************************************************************************
*               RESOLUTION OF A SET OF NON-LINEAR EQUATIONS                    *
* ---------------------------------------------------------------------------- *
* COPYRIGHT 1991, BY R.S. BAIN                                                 *
* Originally wrapped by bain@luther.che.wisc.edu on Tue Mar 30 16:27:48 1993   *
*                                                                              *
* As of this writing, we do not know how to contact Rod Bain, the              * 
* author of NNES.                                                              *
*                                                                              *
*    From F90 Release By David M. Gay (dmg@bell-labs.com)                      *
*    Bell Labs, Murray Hill                                                    *
*    8 February 1999                                                           *
* ---------------------------------------------------------------------------- *
* SAMPLE RUNS:                                                                 *
*    Example #1                                                                *
*    Solve non-linear system of size n=2:                                      *
*    x1^2 + x1 + x2^2 - 2           = 0                                        *
*    x1^2 + x2 - x2^2 - 1 + log(x1) = 0                                        *
*                                                                              *
*    with initial guess; x1=1, x2=1                                            *
*    (The system is described in unit fcn.pas).                                *
*                                                                              *
* With default options (output=2), the output file EX1.OUT contains:           *
*                                                                              *
*                                                                              *
*   *-----------------------------------------------------------------------*  *
*   *-----------------------------------------------------------------------*  *
*   *                                                                       *  *
*   *                           N N E S                                     *  *
*   *                                                                       *  *
*   *       NONMONOTONIC NONLINEAR EQUATION SOLVER VERSION 1.05             *  *
*   *                                                                       *  *
*   *                 COPYRIGHT 1991, BY R.S. BAIN                          *  *
*   *                                                                       *  *
*   *-----------------------------------------------------------------------*  *
*   *-----------------------------------------------------------------------*  *
*                                                                              *
*                                                                              *
*   *-----------------------------------------------------------------------*  *
*   *-----------------------------------------------------------------------*  *
*   *                                                                       *  *
*   *    INITIAL ESTIMATES               INITIAL FUNCTION VALUES            *  *
*   *                                                                       *  *
*   *  X(  1) =        1.000             F(  1) =        1.000              *  *
*   *  X(  2) =        1.000             F(  2) =        0.000              *  *
*   *                                                                       *  *
*   *  INITIAL OBJECTIVE FUNCTION VALUE =      0.500                        *  *
*   *                                                                       *  *
*   *-----------------------------------------------------------------------*  *
*   *-----------------------------------------------------------------------*  *
*   *                                                                       *  *
*   *    CONVERGENCE REACHED; TERMINATION CODE: ...............     2       *  *
*   *                                                                       *  *
*   *                                                                       *  *
*   *          FINAL ESTIMATES                 FINAL FUNCTION VALUES        *  *
*   *                                                                       *  *
*   *    X(  1) =        0.915553            F(  1) =        0.000005       *  *
*   *    X(  2) =        0.496200            F(  2) =       -0.000004       *  *
*   *                                                                       *  *
*   *    FINAL OBJECTIVE FUNCTION VALUE:        0.000                       *  *
*   *                                                                       *  *
*   *                                                                       *  *
*   *    TOTAL NUMBER OF ITERATIONS: ..........................     4       *  *
*   *    TOTAL NUMBER OF LINE SEARCH FUNCTION EVALUATIONS: ....     5       *  *
*   *    TOTAL NUMBER OF EXPLICIT JACOBIAN EVALUATIONS: .......     4       *  *
*   *    TOTAL NUMBER OF FUNCTION EVALUATIONS: ................    13       *  *
*   *                                                                       *  *
*   *-----------------------------------------------------------------------*  *
*                                                                              *
*    Example #2                                                                *
*    Solve non-linear system of size n=4:                                      *
*    10*x1^2 +       + x3    + x4 - 20 + sin(x1)^2 + cos(x1)^2 = 0             *
*    x1      + 20*x2 + x3    + x4 - 48 + 1/x1^6                = 0             *
*    (x1+x2]^2       + 30*x3 + x4 - 97 + Ln(x1) + Ln(x2+x3)    = 0             *
*    x1      + x2    + x3 + 40*x4 -166 + x1*x1                 = 0             *
*                                                                              *
*    with initial guess; x1=1, x2=1, x3=1, x4=1                                *
*    (The system is described in unit fcn.pas).                                *
*                                                                              *
* With default options (output=4), the output file EX2.OUT contains at the end:*
*                                                                              *
*  *    -- / --                                                            *   *  
*  *    SUMMARY OF ITERATION RESULTS                                       *   *
*  *                                                                       *   *
*  *        UPDATED ESTIMATES               UPDATED FUNCTION VALUES        *   *
*  *                                                                       *   *
*  *      X(  1) =        1.041               F(  1) =       -0.000        *   *
*  *      X(  2) =        1.972               F(  2) =       -0.000        *   *
*  *      X(  3) =        2.745               F(  3) =       -0.000        *   *
*  *      X(  4) =        3.979               F(  4) =       -0.000        *   *
*  *                                                                       *   *
*  *      OBJECTIVE FUNCTION VALUE:        0.000                           *   *
*  *                                                                       *   *
*  *      STEP ACCEPTANCE CODE, ACPCOD:        12                          *   *
*  *                                                                       *   *
*  *    CONVERGENCE TESTING                                                *   *
*  *                                                                       *   *
*  *    MAXIMUM STEP SIZE:        0.000   STPTOL:       0.000              *   *
*  *                                                                       *   *
*  *    MAXIMUM ABSOLUTE FUNCTION:        0.000     FTOL:       0.000      *   *
*  *                                                                       *   *
*  *                                                                       *   *
*  *    CONVERGENCE REACHED; TERMINATION CODE: ...............     2       *   *
*  *                                                                       *   *
*  *                                                                       *   *
*  *          FINAL ESTIMATES                 FINAL FUNCTION VALUES        *   *
*  *                                                                       *   *
*  *    X(  1) =        1.040648            F(  1) =       -0.000000       *   *
*  *    X(  2) =        1.972398            F(  2) =       -0.000003       *   *
*  *    X(  3) =        2.745049            F(  3) =       -0.000001       *   *
*  *    X(  4) =        3.978974            F(  4) =       -0.000000       *   *
*  *                                                                       *   *
*  *    FINAL OBJECTIVE FUNCTION VALUE:        0.000                       *   *
*  *                                                                       *   *
*  *                                                                       *   *
*  *    TOTAL NUMBER OF ITERATIONS: ..........................     7       *   *
*  *    TOTAL NUMBER OF LINE SEARCH FUNCTION EVALUATIONS: ....     8       *   *
*  *    TOTAL NUMBER OF EXPLICIT JACOBIAN EVALUATIONS: .......     7       *   *
*  *    TOTAL NUMBER OF FUNCTION EVALUATIONS: ................    36       *   *
*  *                                                                       *   *
*  *    NUMBER OF STEPS ACCEPTED BY FUNCTION VALUE ONLY: .....     0       *   *
*  *    NUMBER OF STEPS ACCEPTED BY STEP SIZE VALUE ONLY: ....     0       *   *
*  *    NUMBER OF STEPS ACCEPTED BY EITHER CRITERION: ........     7       *   *
*  *                                                                       *   *
*  *-----------------------------------------------------------------------*   *
*                                     Pascal Release 1.1 By J-P Moreau, Paris  *
*                                        (See help file nneshelp.txt).         *
*                                              (www.jpmoreau.fr)               *
*                                                                              *
*  Release 1.1:  added example of size 4 (03/29/2007).                         * 
********************************************************************************
This program uses units utils.pas, fcn.pas, unnes.pas, unnes1.pas, unnes2.pas}
PROGRAM Test_Unnes;
Uses WinCrt{Borland unit}, Utils, Unnes, Unnes2{Moreau units}; 

Const
      mgll = 10;
      nunit=10; {not used here}

Var
      alpha, confac, delfac, delta, epsmch, etafac, fcn1new, fdtolj,
      ftol, lam0, mstpf, nsttol, omega, ratiof, sigma, stptol: double;

      acptcr, example, i, itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns, maxqns,
      minqns, n, narmij, niejev, njacch, njetot, output, qnupdm, stopcr, supprs,
      trmcod, trupdm: integer;

      a, h, jac, plee: pMat;
      boundl, boundu, delf, fsave, ftrack, fvec, fvecc, hhpi, rdiag, s, sbar,
      scalef, scalex, sn, ssdhat, strack, vhat, xc, xplus, xsave: pVec;

      absnew, cauchy, deuflh, geoms, linesr, newton, overch: boolean;

      help: string;

BEGIN

  example:=2;   {# example in FCN.PAS}

  {open output text file}
  if example=1 then
    Assign(fp_out,'EX1.OUT')
  else
    Assign(fp_out,'EX2.OUT');
  Rewrite(fp_out);

  {allocate memory}
  New(a); New(h); New(jac); New(plee);
  New(boundl); New(boundu); New(delf); New(fsave); New(ftrack); New(fvec);
  New(fvecc); New(hhpi); New(rdiag); New(s); New(sbar); New(scalef);
  New(scalex); New(sn); New(ssdhat); New(strack); New(vhat); New(xc);
  New(xplus); New(xsave);

  if example=1 then n:=2
  else if example=2 then n:=4
  else n:=6;

  {initial solution in vector xc}
  For i:=1 to n do xc^[i]:=one;

  setup(absnew, cauchy, deuflh, geoms, linesr, newton, overch, acptcr,
        itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns, maxqns,
        minqns, n, narmij, niejev, njacch, output, qnupdm, stopcr, supprs,
        trupdm, alpha, confac, delta, delfac, epsmch, etafac, fdtolj,
        ftol, lam0, mstpf, nsttol, omega, ratiof, sigma, stptol, boundl,
        boundu, scalef, scalex, help);

  {more conditions for example #3}
  if example = 3 then
  begin
    boundl^[3]:=zero;
    boundl^[4]:=zero;
    boundl^[6]:=zero;
    nsttol:=1E-16;
    stptol:=1E-16;
    bypass:=True
  end;

  nnes(absnew, cauchy, deuflh, geoms, linesr, newton, overch, acptcr,
       itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns, maxqns, mgll,
       minqns, n, narmij, niejev, njacch, njetot, nunit, output, qnupdm,
       stopcr, supprs, trmcod, trupdm, alpha, confac, delta, delfac,
       epsmch, etafac, fcn1new, fdtolj, ftol, lam0, mstpf, nsttol, omega,
       ratiof, sigma, stptol, a, boundl, boundu, delf, fsave, ftrack,
       fvec, fvecc, h, hhpi, jac, plee, rdiag, s, sbar, scalef, scalex,
       sn, ssdhat, strack, vhat, xc, xplus, xsave, help);

  writeln;
  if example=1 then
    writeln(' Results in file ex1.out.')
  else
    writeln(' Results in file ex2.out.');
  writeln(' Program terminated...');

  ReadKey;

  Close(fp_out);

  {free memory}
  Dispose(a); Dispose(h); Dispose(jac); Dispose(plee);
  Dispose(boundl); Dispose(boundu); Dispose(delf); Dispose(fsave); Dispose(ftrack);
  Dispose(fvec); Dispose(fvecc); Dispose(hhpi); Dispose(rdiag); Dispose(s); Dispose(sbar);
  Dispose(scalef); Dispose(scalex); Dispose(sn); Dispose(ssdhat); Dispose(strack);
  Dispose(vhat); Dispose(xc); Dispose(xplus); Dispose(xsave);

  {exit program}
  DoneWinCrt

END.

{end of file tnnes.pas}