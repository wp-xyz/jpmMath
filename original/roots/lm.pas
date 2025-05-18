{*********************************************************************
* PROCEDURES FOR THE LEAST-SQUARES SOLUTION OF M NONLINEAR EQUATIONS *
* IN N VARIABLES USING THE Levenberg-Marquardt ALGORITHM.            *
* ------------------------------------------------------------------ *
*  REFERENCE                                                         *
*  From F77 program By:                                              *
*  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.         *
*  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             *
*                                                                    *
*                            Pascal Release By J-P Moreau, Paris.    *
*                                     (www.jpmoreau.fr)              *
*********************************************************************}
UNIT LM;

INTERFACE


Const  SIZE = 25;

       zero = 0.0; zp25 = 0.25; zp5 = 0.5; one = 1.0; two = 2.0; five = 5.0; eight = 8.0;
       ten = 10.0; c13 = 13.0; c14 = 14.0; c29 = 29.0; c45 = 45.0;
       p5 = 0.5; p25 = 0.25; p05 = 0.05;

       epsmch = 2.25E-16;

Type
     pMat = ^MAT;
     MAT =  Array[1..SIZE,1..SIZE] of double;
     pVec = ^VEC;
     VEC = Array[1..SIZE] of double;
     pIVec = ^IVEC;
     IVEC = Array[1..SIZE] of integer;
Var
     nprob, nfev, njev: Integer;

     Procedure ssqfcn (m, n:integer; x, fvec:pVec; nprob:integer);
     Function  enorm(n:integer; x:pVec): double;
     Procedure lmdif1(m, n: integer; x, fvec: pVec; tol: double; Var info:integer; iwa: pIVec);


IMPLEMENTATION

     Procedure fdjac2(m, n:integer; x, fvec:pVec; fjac:pMat; var iflag:integer;
                      epsfcn: double); Forward;
     Procedure qrfac(m, n:integer; a:pMat; pivot:Boolean; ipvt:pIVec; rdiag, acnorm:pVec);
     Forward;
     Procedure lmpar(n:integer; r:pMat; ipvt:pIVec; diag, qtb:pVec; delta:double;
                     var par:double; x, sdiag: pVec); Forward;
     Procedure lmder(m, n:integer; x, fvec:pVec; fjac:pMat; ftol, xtol, gtol:double;
                     var maxfev:integer; mode, factor, nprint:integer; var info, nfev,
                     njev:integer; ipvt:pIVec); Forward;
     Procedure qrsolv(n:integer; r:pMat; ipvt:pIVec; diag, qtb, x, sdiag:pVec); Forward;


     {Utility functions}
     Function Dot_Product(n:integer; a,b:pVec):double;
     Var i:integer; sum:double;
     begin
       sum:=zero;
       For i:=1 to n do sum:=sum+a^[i]*b^[i];
       Dot_Product:=sum
     end;

     Function MAX(a,b:double): double;
     begin
       if a>=b then MAX:=a else MAX:=b
     end;

     Function MIN(a,b:double): double;
     begin
       if a<=b then MIN:=a else MIN:=b
     end;


Procedure ssqfcn (m, n:integer; x, fvec:pVec; nprob:integer);
{     ****************************************************************

!     Procedure SSQFCN

!     THIS PROCEDURE DEFINES THE FUNCTIONS OF THREE NONLINEAR
!     LEAST SQUARES PROBLEMS. 

!     THE PROCEDURE STATEMENT IS

!       PROCEDURE SSQFCN(M, N, X, FVEC, NPROB)

!     WHERE

!       M AND N ARE POSITIVE INTEGER INPUT VARIABLES. N MUST NOT
!         EXCEED M.

!       X IS AN INPUT ARRAY OF LENGTH N (OF TYPE pVec).

!       FVEC IS AN OUTPUT ARRAY OF LENGTH M WHICH CONTAINS THE NPROB
!         FUNCTION EVALUATED AT X (OF TYPE pVec).

!       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
!         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.

!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!     *************************************************************** }
Var
    i, iev, j, nm1: integer;

Begin
{ FUNCTION ROUTINE SELECTOR }
  Case nprob of

    1:begin     { Example #1 (n=2) }
{     X1^2 + X1 + X2^2 - 2           = 0
      X1^2 + X2 - X2^2 - 1 + log(X1) = 0 }
      fvec^[1] := x^[1]*x^[1] + x^[1] + x^[2]*x^[2] - two;
      fvec^[2] := x^[1]*x^[1] + x^[2] - x^[2]*x^[2] - one + ln(x^[1])
      end;
    2:begin     { Example #2 (n=4) }
{     10.0*x + x2 + x3  + x4 - 20.0 + Sqr(sin(x1)) + Sqr(cos(x2)) = 0
      x1 + 20.0*x2 + x3 + x4 - 48.0 + one/pow^6                   = 0
      Sqr(x1 + x2) + 30.0*x3 + x4 - 97.0 + log(x1) + log(x2+x3)   = 0
      x1     + x2  + x3 + 40.0*x4 - 166.0 + Sqr(x1)               = 0  }
      fvec^[1] := 10.0*x^[1] + x^[2] + x^[3] + x^[4] - 20.0 + Sqr(sin(x^[1])) + Sqr(cos(x^[2]));
      fvec^[2] := x^[1] + 20.0*x^[2] + x^[3] + x^[4] - 48.0 + one/(x^[1]*x^[1]*x^[1]*x^[1]*x^[1]*x^[1]);
      fvec^[3] := Sqr(x^[1] + x^[2]) + 30.0*x^[3] + x^[4] - 97.0 + ln(x^[1]) + ln(x^[2]+x^[3]);
      fvec^[4] := x^[1] + x^[2] + x^[3] + 40.0*x^[4] - 166.0 + Sqr(x^[1])
      end;
    3:begin     { Example #3 (n=6)  Stiff system }
{    Hiebert's 2nd Chemical Engineering Problem

     source: Hiebert; Sandia Technical Report #SAND80-0181
             Sandia National Laboratories, Albuquerque, NM (1980)

      X1 + X2 + X4 - .001 = 0
      X5 + X6 -55         = 0
      X1 + X2 + X3 + 2X5 + X6 - 110.001 = 0
      X1 - 0.1X2        = 0
      X1 - 10000 X3 X4  = 0
      X5 - 5.5e15 X3 X6 = 0
      solution: (8.264e-5, 8.264e-4, 9.091e-5, 9.091e-5, 55, 1.1e-10)  }
      fvec^[1] := x^[1] + x^[2] + x^[4] - 0.001;
      fvec^[2] := x^[5] + x^[6] - 55.0;
      fvec^[3] := x^[1] + x^[2] + x^[3] + two * x^[5] + x^[6] - 110.001;
      fvec^[4] := x^[1] - 0.1 * x^[2];
      fvec^[5] := x^[1] - 1.0E04 * x^[3] * x^[4];
      fvec^[6] := x^[5] - 5.5E15 * x^[3] * x^[6]
      end
  end

End; {ssqfcn}

Procedure fcn (m, n:integer; x, fvec:pVec; var iflag:integer);
{  ************************************************************

!  THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
!  CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE NONLINEAR
!  LEAST-SQUARES SOLVER. FCN SHOULD ONLY CALL THE TESTING
!  FUNCTION SUBROUTINE SSQFCN WITH THE APPROPRIATE VALUE OF
!  PROBLEM NUMBER (NPROB).

!  SUBPROGRAMS CALLED

!    USER SUPPLIED: SSQFCN

!  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!  ************************************************************ }
Begin
  ssqfcn (m, n, x, fvec, nprob);
  IF iflag = 1 then nfev := nfev + 1;
  IF iflag = 2 then njev := njev + 1
End; {fcn}

Procedure lmdif(m, n:integer; x, fvec:pVec; ftol, xtol:double; gtol:double;
                var maxfev:integer; epsfcn:double; mode:integer; factor:double;
                nprint:integer; var info, nfev:integer; fjac:pMat; ipvt:pIVec);
Forward;

Procedure lmdif1(m, n: integer; x, fvec: pVec; tol: double; Var info:integer; iwa: pIVec);
 
{ From F90 Code converted using TO_F90 by Alan Miller
!  **************************************************************************

!  procedure lmdif1

!  The purpose of lmdif1 is to minimize the sum of the squares of m nonlinear
!  functions in n variables by a modification of the Levenberg-Marquardt
!  algorithm.  This is done by using the more general least-squares
!  solver lmdif.  The user must provide a subroutine which calculates the
!  functions.  The jacobian is then calculated by a forward-difference
!  approximation.

!  the subroutine statement is

!    procedure lmdif1(m, n, x, fvec, tol, info, iwa)

!  where

!    fcn is the name of the user-supplied procedure which calculates
!      the functions (removed from list of arguments).
!      fcn should be written as follows:

!      Procedure fcn (m, n:integer; x, fvec:pVec; var iflag:integer);
!      begin
!      ----------
!      calculate the functions at x and return this vector in fvec.
!      ----------
!      end;

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of lmdif1.
!      In this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number of functions.

!    n is a positive integer input variable set to the number of variables.
!      n must not exceed m.

!    x is an array of length n.  On input x must contain an initial estimate
!      of the solution vector.  On output x contains the final estimate of
!      the solution vector.

!    fvec is an output array of length m which contains
!      the functions evaluated at the output x.

!    tol is a nonnegative input variable.  Termination occurs when the
!      algorithm estimates either that the relative error in the sum of
!      squares is at most tol or that the relative error between x and the
!      solution is at most tol.

!    info is an integer output variable.  If the user has terminated execution,
!      info is set to the (negative) value of iflag.  See description of fcn.
!      Otherwise, info is set as follows.

!      info = 0  improper input parameters.

!      info = 1  algorithm estimates that the relative error
!                in the sum of squares is at most tol.

!      info = 2  algorithm estimates that the relative error
!                between x and the solution is at most tol.

!      info = 3  conditions for info = 1 and info = 2 both hold.

!      info = 4  fvec is orthogonal to the columns of the
!                jacobian to machine precision.

!      info = 5  number of calls to fcn has reached or exceeded 200*(n+1).

!      info = 6  tol is too small. no further reduction in
!                the sum of squares is possible.

!      info = 7  tol is too small.  No further improvement in
!                the approximate solution x is possible.

!    iwa is an integer work array of length n.

!    wa is a work array of length lwa.

!    lwa is a positive integer input variable not less than m*n+5*n+m.

!  subprograms called

!    user-supplied ...... fcn

!    minpack-supplied ... lmdif

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  ************************************************************************* }
Label Return;
Var
  maxfev, mode, nprint: Integer;
  epsfcn, ftol, gtol, xtol: double;
  fjac: pMat;

Const
  factor = 100.0;
Begin
  info := 0;
  New(fjac);

{ check the input parameters for errors }
  IF (n <= 0) OR (m < n) OR (tol < zero) Then goto Return;

  maxfev := 200*(n + 1);
  ftol := tol;
  xtol := tol;
  gtol := zero;
  epsfcn := zero;
  mode := 1;
  nprint := 0;

  lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,
           mode, factor, nprint, info, nfev, fjac, iwa);

  IF info = 8 then info := 4;

Return: Dispose(fjac);

End; {lmdif1}



Procedure lmdif(m, n:integer; x, fvec:pVec; ftol, xtol:double; gtol:double;
                var maxfev:integer; epsfcn:double; mode:integer; factor:double;
                nprint:integer; var info, nfev:integer; fjac:pMat; ipvt:pIVec);
{  ***************************************************************************

!  procedure lmdif

!  The purpose of lmdif is to minimize the sum of the squares of m nonlinear
!  functions in n variables by a modification of the Levenberg-Marquardt
!  algorithm.  The user must provide a subroutine which calculates the
!  functions.  The jacobian is then calculated by a forward-difference
!  approximation.

!  the procedure statement is

!    procedure lmdif(m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,
!                    diag, mode, factor, nprint, info, nfev, fjac,
!                    ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)

!  where:

!    fcn is the name of the user-supplied procedure which calculates
!      the functions (removed from list of arguments).
!      fcn should be written as follows:

!      Procedure fcn (m, n:integer; x, fvec:pVec; var iflag:integer);
!      begin
!      ----------
!      calculate the functions at x and return this vector in fvec.
!      ----------
!      end;

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of lmdif1.
!      In this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number of functions.

!    n is a positive integer input variable set to the number of variables.
!      n must not exceed m.

!    x is an array of length n.  On input x must contain an initial estimate
!      of the solution vector.  On output x contains the final estimate of the
!      solution vector.

!    fvec is an output array of length m which contains
!      the functions evaluated at the output x.

!    ftol is a nonnegative input variable.  Termination occurs when both the
!      actual and predicted relative reductions in the sum of squares are at
!      most ftol.  Therefore, ftol measures the relative error desired
!      in the sum of squares.

!    xtol is a nonnegative input variable.  Termination occurs when the
!      relative error between two consecutive iterates is at most xtol.
!      Therefore, xtol measures the relative error desired in the approximate
!      solution.

!    gtol is a nonnegative input variable.  Termination occurs when the cosine
!      of the angle between fvec and any column of the jacobian is at most
!      gtol in absolute value.  Therefore, gtol measures the orthogonality
!      desired between the function vector and the columns of the jacobian.

!    maxfev is a positive integer input variable.  Termination occurs when the
!      number of calls to fcn is at least maxfev by the end of an iteration.

!    epsfcn is an input variable used in determining a suitable step length
!      for the forward-difference approximation.  This approximation assumes
!      that the relative errors in the functions are of the order of epsfcn.
!      If epsfcn is less than the machine precision, it is assumed that the
!      relative errors in the functions are of the order of the machine
!      precision.

!    diag is an array of length n.  If mode = 1 (see below), diag is
!      internally set.  If mode = 2, diag must contain positive entries that
!      serve as multiplicative scale factors for the variables.

!    mode is an integer input variable.  If mode = 1, the variables will be
!      scaled internally.  If mode = 2, the scaling is specified by the input
!      diag. other values of mode are equivalent to mode = 1.

!    factor is a positive input variable used in determining the initial step
!      bound.  This bound is set to the product of factor and the euclidean
!      norm of diag*x if nonzero, or else to factor itself.  In most cases
!      factor should lie in the interval (.1,100.). 100. is a generally
!      recommended value.

!    nprint is an integer input variable that enables controlled printing of
!      iterates if it is positive.  In this case, fcn is called with iflag := 0
!      at the beginning of the first iteration and every nprint iterations
!      thereafter and immediately prior to return, with x and fvec available
!      for printing.  If nprint is not positive, no special calls
!      of fcn with iflag := 0 are made.

!    info is an integer output variable.  If the user has terminated
!      execution, info is set to the (negative) value of iflag.
!      See description of fcn.  Otherwise, info is set as follows.

!      info = 0  improper input parameters.

!      info = 1  both actual and predicted relative reductions
!                in the sum of squares are at most ftol.

!      info = 2  relative error between two consecutive iterates <:= xtol.

!      info = 3  conditions for info := 1 and info := 2 both hold.

!      info = 4  the cosine of the angle between fvec and any column of
!                the Jacobian is at most gtol in absolute value.

!      info = 5  number of calls to fcn has reached or exceeded maxfev.

!      info = 6  ftol is too small. no further reduction in
!                the sum of squares is possible.

!      info = 7  xtol is too small. no further improvement in
!                the approximate solution x is possible.

!      info = 8  gtol is too small. fvec is orthogonal to the
!                columns of the jacobian to machine precision.

!    nfev is an integer output variable set to the number of calls to fcn.

!    fjac is an output m by n array. the upper n by n submatrix
!      of fjac contains an upper triangular matrix r with
!      diagonal elements of nonincreasing magnitude such that

!             t     t           t
!            p *(jac *jac)*p = r *r,

!      where p is a permutation matrix and jac is the final calculated
!      Jacobian.  Column j of p is column ipvt(j) (see below) of the
!      identity matrix. the lower trapezoidal part of fjac contains
!      information generated during the computation of r.

!    ldfjac is a positive integer input variable not less than m
!      which specifies the leading dimension of the array fjac.

!    ipvt is an integer output array of length n.  ipvt defines a permutation
!      matrix p such that jac*p = q*r, where jac is the final calculated
!      jacobian, q is orthogonal (not stored), and r is upper triangular
!      with diagonal elements of nonincreasing magnitude.
!      Column j of p is column ipvt(j) of the identity matrix.

!    qtf is an output array of length n which contains
!      the first n elements of the vector (q transpose)*fvec.

!    wa1, wa2, and wa3 are work arrays of length n.

!    wa4 is a work array of length m.

!  subprograms called

!    user-supplied ...... fcn

!    dpmpar, enorm, fdjac2, lmpar, max, min, qrfac

!    pascal-supplied ... abs, sqrt, mod

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  ************************************************************************* }
Label 20,30,40,60,80,120,170,200,260,290,300;
Var
  i, iflag, iter, j, l: Integer;
  actred, delta, dirder, fnorm, fnorm1, gnorm,
  par, pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm: double;
  diag, qtf, tmp, tmp1, wa1, wa2, wa3, wa4: pVec;

Const
  p1 = 0.1; p75 = 0.75; p0001 = 0.0001;

Begin

  New(diag); New(qtf); New(tmp); New(tmp1); New(wa1); New(wa2);
  New(wa3);  New(wa4);

  info := 0;
  iflag := 0;
  nfev := 0;

{ check the input parameters for errors }

  IF (n <= 0) OR (m < n) OR (ftol < zero) OR (xtol < zero) OR (gtol < zero)
        OR (maxfev <= 0) OR (factor <= zero) Then goto 300;

  IF mode <> 2 then GOTO 20;

  For j := 1 to n do
    IF diag^[j] <= zero then GOTO 300;


{ evaluate the function at the starting point and calculate its norm }

20: iflag := 1;
  fcn(m, n, x, fvec, iflag);
  nfev := 1;
  IF iflag < 0 then GOTO 300;

  fnorm := enorm(m, fvec);

{ initialize levenberg-marquardt parameter and iteration counter }

  par := zero;
  iter := 1;

{ beginning of the outer loop.

  calculate the jacobian matrix }

30:iflag := 2;
  fdjac2(m, n, x, fvec, fjac, iflag, epsfcn);
  nfev := nfev + n;
  IF iflag < 0 then GOTO 300;

{ If requested, call fcn to enable printing of iterates }

  IF nprint <= 0 then GOTO 40;
  iflag := 0;
  IF (iter-1 mod nprint) = 0 then fcn(m, n, x, fvec, iflag);
  IF iflag < 0 then GOTO 300;

{ Compute the qr factorization of the jacobian }

40: qrfac(m, n, fjac, true, ipvt, wa1, wa2);

{ On the first iteration and if mode is 1, scale according
  to the norms of the columns of the initial jacobian }

  IF iter <> 1 then GOTO 80;
  IF mode = 2 then GOTO 60;
  For j := 1 to n do
  begin
    diag^[j] := wa2^[j];
    IF wa2^[j] = zero then diag^[j] := one
  end;

{ On the first iteration, calculate the norm of the scaled x
  and initialize the step bound delta }

60: For j := 1 to n do wa3^[j] := diag^[j]*x^[j];
  xnorm := enorm(n, wa3);
  delta := factor*xnorm;
  IF delta = zero then delta := factor;

{ Form (q transpose)*fvec and store the first n components in qtf }

80: For j := 1 to m do wa4^[j] := fvec^[j];
  For  j := 1 to n do
  begin
    IF fjac^[j,j] = zero then GOTO 120;
    For i:=j to m do
    begin
      tmp^[i-j+1]:=fjac^[i,j];
      tmp1^[i-j+1]:=wa4^[i]
    end;   
    sum := DOT_PRODUCT(m-j+1, tmp, tmp1);
    temp := -sum/fjac^[j,j];
    For  i := j to m do
      wa4^[i] := wa4^[i] + fjac^[i,j]*temp;
120:fjac^[j,j] := wa1^[j];
    qtf^[j] := wa4^[j]
  end;

{ compute the norm of the scaled gradient }

  gnorm := zero;
  IF fnorm = zero then GOTO 170;
  For j := 1 to n do
  begin
    l := ipvt^[j];
    IF wa2^[l] = zero then exit;
    sum := zero;
    For  i := 1 to j do
      sum := sum + fjac^[i,j]*(qtf^[i]/fnorm);
    gnorm := MAX(gnorm, ABS(sum/wa2^[l]))
  end;

{ test for convergence of the gradient norm }

170:IF gnorm <= gtol then info := 4;
  IF info <> 0 then GOTO 300;

{ rescale if necessary }

  IF mode = 2 then GOTO 200;
  For  j := 1 to  n do
    diag^[j] := MAX(diag^[j], wa2^[j]);

{ beginning of the inner loop.

 determine the Levenberg-Marquardt parameter }

200: lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2);

{ store the direction p and x + p. calculate the norm of p }

  For j := 1 to n do
  begin
    wa1^[j] := -wa1^[j];
    wa2^[j] := x^[j] + wa1^[j];
    wa3^[j] := diag^[j]*wa1^[j]
  end;
  pnorm := enorm(n, wa3);

{ on the first iteration, adjust the initial step bound }

  IF iter = 1 then delta := MIN(delta, pnorm);

{ evaluate the function at x + p and calculate its norm }

  iflag := 1;
  fcn(m, n, wa2, wa4, iflag);
  nfev := nfev + 1;
  IF iflag < 0 then GOTO 300;
  fnorm1 := enorm(m, wa4);

{ compute the scaled actual reduction }

  actred := -one;
  IF p1*fnorm1 < fnorm then actred := one - Sqr(fnorm1/fnorm);

{ Compute the scaled predicted reduction and
  the scaled directional derivative }

  For j := 1 to n do
  begin
    wa3^[j] := zero;
    l := ipvt^[j];
    temp := wa1^[l];
    For i := 1 to j do
      wa3^[i] := wa3^[i] + fjac^[i,j]*temp
  end;
  temp1 := enorm(n,wa3)/fnorm;
  temp2 := (SQRT(par)*pnorm)/fnorm;
  prered := temp1*temp1 + temp2*temp2/p5;
  dirder := -(temp1*temp1 + temp2*temp2);

{ compute the ratio of the actual to the predicted reduction }

  ratio := zero;
  IF prered <> zero then  ratio := actred/prered;

{ update the step bound }

  IF ratio <= p25 THEN
  begin
    IF actred >= zero then temp := p5;
    IF actred < zero  then temp := p5*dirder/(dirder + p5*actred);
    IF (p1*fnorm1 >= fnorm) OR (temp < p1) then temp := p1;
    delta := temp*MIN(delta,pnorm/p1);
    par := par/temp
  end
  ELSE
  begin
    IF (par <> zero) AND (ratio < p75) then GOTO 260;
    delta := pnorm/p5;
    par := p5*par
  end;

{ test for successful iteration }

260: IF ratio < p0001 then GOTO 290;

{ successful iteration. update x, fvec, and their norms }

  For j := 1 to n do
  begin
    x^[j] := wa2^[j];
    wa2^[j] := diag^[j]*x^[j]
  end;
  For j:=1 to m do fvec^[j] := wa4^[j];
  xnorm := enorm(n, wa2);
  fnorm := fnorm1;
  iter := iter + 1;

{ tests for convergence }

290: IF (ABS(actred) <= ftol) AND (prered <= ftol) AND (p5*ratio <= one) then info := 1;
  IF delta <= xtol*xnorm then info := 2;
  IF (ABS(actred) <= ftol) AND (prered <= ftol) AND (p5*ratio <= one)
    AND (info = 2) then info := 3;
  IF info <> 0 then GOTO 300;

{ tests for termination and stringent tolerances }

  IF nfev >= maxfev then info := 5;
  IF (ABS(actred) <= epsmch) AND (prered <= epsmch) AND (p5*ratio <= one) then info := 6;
  IF delta <= epsmch*xnorm then info := 7;
  IF gnorm <= epsmch then info := 8;
  IF info <> 0 then GOTO 300;

{ end of the inner loop. repeat if iteration unsuccessful }

  IF ratio < p0001 then GOTO 200;

{ end of the outer loop }

  GOTO 30;

{ termination, either normal or user imposed }

300:IF iflag < 0 then info := iflag;
  iflag := 0;
  IF nprint > 0 then fcn(m, n, x, fvec, iflag);

  Dispose(diag); Dispose(qtf); Dispose(tmp); Dispose(tmp1); Dispose(wa1);
  Dispose(wa2);  Dispose(wa3); Dispose(wa4)

End; {lmdif}


Procedure lmder1(m, n:integer; x, fvec:pVec; fjac:pMat; tol:double;
                 var info:integer; ipvt:pIVec);
{  ***************************************************************

!  procedure lmder1

!  The purpose of lmder1 is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of the
!  levenberg-marquardt algorithm.  This is done by using the more
!  general least-squares solver lmder.  The user must provide a
!  subroutine which calculates the functions and the jacobian.

!  the procedure statement is

!    procedure lmder1(m, n, x, fvec, fjac, tol, info, ipvt)

!  where

!    fcn is the name of the user-supplied procedure which calculates
!      the functions (removed from list of arguments).
!      fcn should be written as follows:

!      Procedure fcn (m, n:integer; x, fvec:pVec; var iflag:integer);
!      begin
!      ----------
!      calculate the functions at x and return this vector in fvec.
!      ----------
!      end;

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of lmdif1.
!      In this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number of functions.

!    n is a positive integer input variable set to the number
!      of variables.  n must not exceed m.

!    x is an array of length n. on input x must contain
!      an initial estimate of the solution vector. on output x
!      contains the final estimate of the solution vector.

!    fvec is an output array of length m which contains
!      the functions evaluated at the output x.

!    fjac is an output m by n array. the upper n by n submatrix
!      of fjac contains an upper triangular matrix r with
!      diagonal elements of nonincreasing magnitude such that

!             t     t           t
!            p *(jac *jac)*p := r *r,

!      where p is a permutation matrix and jac is the final calculated
!      Jacobian.  Column j of p is column ipvt(j) (see below) of the
!      identity matrix.  The lower trapezoidal part of fjac contains
!      information generated during the computation of r.

!    ldfjac is a positive integer input variable not less than m
!      which specifies the leading dimension of the array fjac.

!    tol is a nonnegative input variable. termination occurs
!      when the algorithm estimates either that the relative
!      error in the sum of squares is at most tol or that
!      the relative error between x and the solution is at most tol.

!    info is an integer output variable.  If the user has terminated
!      execution, info is set to the (negative) value of iflag.
!      See description of fcn.  Otherwise, info is set as follows.

!      info := 0  improper input parameters.

!      info := 1  algorithm estimates that the relative error
!                in the sum of squares is at most tol.

!      info := 2  algorithm estimates that the relative error
!                between x and the solution is at most tol.

!      info := 3  conditions for info := 1 and info := 2 both hold.

!      info := 4  fvec is orthogonal to the columns of the
!                jacobian to machine precision.

!      info := 5  number of calls to fcn with iflag := 1 has reached 100*(n+1).

!      info := 6  tol is too small.  No further reduction in
!                the sum of squares is possible.

!      info := 7  tol is too small.  No further improvement in
!                the approximate solution x is possible.

!    ipvt is an integer output array of length n. ipvt
!      defines a permutation matrix p such that jac*p := q*r,
!      where jac is the final calculated jacobian, q is
!      orthogonal (not stored), and r is upper triangular
!      with diagonal elements of nonincreasing magnitude.
!      column j of p is column ipvt(j) of the identity matrix.

!    wa is a work array of length lwa.

!    lwa is a positive integer input variable not less than 5*n+m.

!  subprograms called

!    user-supplied ...... fcn

!    minpack-supplied ... lmder

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  ************************************************************************** }
Label 10;
Var
  maxfev, mode, njev, nprint: Integer;
  ftol, gtol, xtol: double;
Const
  factor = 100;
Begin
  info := 0;

{ check the input parameters for errors }

  IF (n <= 0) OR (m < n) OR (tol < zero) then GOTO 10;

{ call lmder }

  maxfev := 100*(n + 1);
  ftol := tol;
  xtol := tol;
  gtol := zero;
  mode := 1;
  nprint := 0;
  lmder(m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev,
           mode, factor, nprint, info, nfev, njev, ipvt);
  IF info = 8 then info := 4;

10: End; {lmder1}



Procedure lmder(m, n:integer; x, fvec:pVec; fjac:pMat; ftol, xtol, gtol:double;
                var maxfev:integer; mode, factor, nprint:integer; var info, nfev,
                njev:integer; ipvt:pIVec);
{  ****************************************************************************

!  procedure lmder

!  the purpose of lmder is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of
!  the levenberg-marquardt algorithm. the user must provide a
!  procedure which calculates the functions and the jacobian.

!  the procedure statement is

!    procedure lmder(m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!                    maxfev,diag,mode,factor,nprint,info,nfev,
!                    njev,ipvt,qtf,wa1,wa2,wa3,wa4)

!  where:

!    fcn is the name of the user-supplied procedure which calculates
!      the functions (removed from list of arguments).
!      fcn should be written as follows:

!      Procedure fcn (m, n:integer; x, fvec:pVec; var iflag:integer);
!      begin
!      ----------
!      calculate the functions at x and return this vector in fvec.
!      ----------
!      end;

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of lmdif1.
!      In this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number
!      of functions.

!    n is a positive integer input variable set to the number
!      of variables. n must not exceed m.

!    x is an array of length n. on input x must contain
!      an initial estimate of the solution vector. on output x
!      contains the final estimate of the solution vector.

!    fvec is an output array of length m which contains
!      the functions evaluated at the output x.

!    fjac is an output m by n array. the upper n by n submatrix
!      of fjac contains an upper triangular matrix r with
!      diagonal elements of nonincreasing magnitude such that

!             t     t           t
!            p *(jac *jac)*p := r *r

!      where p is a permutation matrix and jac is the final calculated
!      jacobian.  Column j of p is column ipvt(j) (see below) of the
!      identity matrix.  The lower trapezoidal part of fjac contains
!      information generated during the computation of r.

!    ldfjac is a positive integer input variable not less than m
!      which specifies the leading dimension of the array fjac.

!    ftol is a nonnegative input variable.  Termination occurs when both
!      the actual and predicted relative reductions in the sum of squares
!      are at most ftol.   Therefore, ftol measures the relative error
!      desired in the sum of squares.

!    xtol is a nonnegative input variable. termination
!      occurs when the relative error between two consecutive
!      iterates is at most xtol. therefore, xtol measures the
!      relative error desired in the approximate solution.

!    gtol is a nonnegative input variable.  Termination occurs when the
!      cosine of the angle between fvec and any column of the jacobian is
!      at most gtol in absolute value.  Therefore, gtol measures the
!      orthogonality desired between the function vector and the columns
!      of the jacobian.

!    maxfev is a positive integer input variable.  Termination occurs when
!      the number of calls to fcn with iflag := 1 has reached maxfev.

!    diag is an array of length n.  If mode := 1 (see below), diag is
!      internally set.  If mode := 2, diag must contain positive entries
!      that serve as multiplicative scale factors for the variables.

!    mode is an integer input variable.  if mode := 1, the
!      variables will be scaled internally.  if mode := 2,
!      the scaling is specified by the input diag.  other
!      values of mode are equivalent to mode := 1.

!    factor is a positive input variable used in determining the
!      initial step bound. this bound is set to the product of
!      factor and the euclidean norm of diag*x if nonzero, or else
!      to factor itself. in most cases factor should lie in the
!      interval (.1,100.).100. is a generally recommended value.

!    nprint is an integer input variable that enables controlled printing
!      of iterates if it is positive.  In this case, fcn is called with
!      iflag := 0 at the beginning of the first iteration and every nprint
!      iterations thereafter and immediately prior to return, with x, fvec,
!      and fjac available for printing.  fvec and fjac should not be
!      altered.  If nprint is not positive, no special calls of fcn with
!      iflag := 0 are made.

!    info is an integer output variable.  If the user has terminated
!      execution, info is set to the (negative) value of iflag.
!      See description of fcn.  Otherwise, info is set as follows.

!      info := 0  improper input parameters.

!      info := 1  both actual and predicted relative reductions
!                in the sum of squares are at most ftol.

!      info := 2  relative error between two consecutive iterates
!                is at most xtol.

!      info := 3  conditions for info := 1 and info := 2 both hold.

!      info := 4  the cosine of the angle between fvec and any column of
!                the jacobian is at most gtol in absolute value.

!      info := 5  number of calls to fcn with iflag := 1 has reached maxfev.

!      info := 6  ftol is too small.  No further reduction in
!                the sum of squares is possible.

!      info := 7  xtol is too small.  No further improvement in
!                the approximate solution x is possible.

!      info := 8  gtol is too small.  fvec is orthogonal to the
!                columns of the jacobian to machine precision.

!    nfev is an integer output variable set to the number of
!      calls to fcn with iflag = 1.

!    njev is an integer output variable set to the number of
!      calls to fcn with iflag = 2.

!    ipvt is an integer output array of length n.  ipvt
!      defines a permutation matrix p such that jac*p := q*r,
!      where jac is the final calculated jacobian, q is
!      orthogonal (not stored), and r is upper triangular
!      with diagonal elements of nonincreasing magnitude.
!      column j of p is column ipvt(j) of the identity matrix.

!    qtf is an output array of length n which contains
!      the first n elements of the vector (q transpose)*fvec.

!    wa1, wa2, and wa3 are work arrays of length n.

!    wa4 is a work array of length m.

!  subprograms called

!    user-supplied ...... fcn, ssqfcn

!    dpmpar,enorm,lmpar,qrfac,max,min

!    pascal-supplied ... ABS, SQRT, mod

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  ***************************************************************************** }
Label 20,30,40,60,80,120,170,200,260,290,300;
Var
  i, iflag, iter, j, l: Integer;
  actred, delta, dirder, fnorm, fnorm1, gnorm, par, pnorm, prered,
  ratio, sum, temp, temp1, temp2, xnorm: double;
  diag, qtf, tmp, tmp1, wa1, wa2, wa3, wa4: pVec;
Const
  p1 = 0.1; p75 = 0.75; p0001 = 0.0001;
Begin

  New(diag); New(qtf); New(tmp); New(tmp1); New(wa1); New(wa2); New(wa3); New(wa4);

  info := 0;
  iflag := 0;
  nfev := 0;
  njev := 0;

{ check the input parameters for errors }

  IF (n <= 0) OR (m < n) OR (ftol < zero) OR (xtol < zero) OR (gtol < zero)
    OR (maxfev <= 0) OR (factor <= zero) then GOTO 300;
  IF mode <> 2 then GOTO 20;
  For j := 1 to n do
    IF diag^[j] <= zero then GOTO 300;

{ evaluate the function at the starting point and calculate its norm }

20: iflag := 1;
  fcn(m, n, x, fvec, iflag);
  nfev := 1;
  IF iflag < 0 then GOTO 300;
  fnorm := enorm(m, fvec);

{ initialize levenberg-marquardt parameter and iteration counter }

  par := zero;
  iter := 1;

{ beginning of the outer loop.

  calculate the jacobian matrix }

30: iflag := 2;
  fcn(m, n, x, fvec, iflag);
  njev := njev + 1;
  IF iflag < 0 then GOTO 300;

{ if requested, call fcn to enable printing of iterates }

  IF nprint <= 0 then GOTO 40;
  iflag := 0;
  IF (iter-1) mod nprint = 0 then fcn(m, n, x, fvec, iflag);
  IF iflag < 0 then GOTO 300;

{ compute the qr factorization of the jacobian }

40: qrfac(m, n, fjac, true, ipvt, wa1, wa2);

{ on the first iteration and if mode is 1, scale according
  to the norms of the columns of the initial jacobian }

  IF iter <> 1 then GOTO 80;
  IF mode = 2 then GOTO 60;
  For j := 1 to n do
  begin
    diag^[j] := wa2^[j];
    IF wa2^[j] = zero then diag^[j] := one
  end;

{ on the first iteration, calculate the norm of the scaled x
  and initialize the step bound delta }

60:For j:=1 to n do wa3^[j] := diag^[j]*x^[j];
  xnorm := enorm(n,wa3);
  delta := factor*xnorm;
  IF delta = zero then delta := factor;

{ form (q transpose)*fvec and store the first n components in qtf }

80:For j:=1 to m do wa4^[j] := fvec^[j];
  For j := 1 to n do
  begin
    IF fjac^[j,j] = zero then goto 120;
    For i:=j to m do
    begin
      tmp^[i-j+1]:=fjac^[i,j];
      tmp1^[i-j+1]:=wa4^[i]
    end;
    sum := DOT_PRODUCT(m-j+1,tmp,tmp1);
    temp := -sum/fjac^[j,j];
    For i := j to m do wa4^[i] := wa4^[i] + fjac^[i,j]*temp;
120: fjac^[j,j] := wa1^[j];
    qtf^[j] := wa4^[j]
  end;

{ compute the norm of the scaled gradient }

  gnorm := zero;
  IF fnorm = zero then GOTO 170;
  For j := 1 to n do
  begin
    l := ipvt^[j];
    IF wa2^[l] = zero then exit;
    sum := zero;
    For i := 1 to j do sum := sum + fjac^[i,j]*(qtf^[i]/fnorm);
    gnorm := MAX(gnorm,ABS(sum/wa2^[l]))
  end;

{ test for convergence of the gradient norm }

170:IF gnorm <= gtol then info := 4;
  IF info <> 0 then GOTO 300;

{ rescale if necessary }

  IF mode = 2 then GOTO 200;
  For j := 1 to n do diag^[j] := MAX(diag^[j], wa2^[j]);

{ beginning of the inner loop.

  determine the levenberg-marquardt parameter }

200: lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2);

{ store the direction p and x + p. calculate the norm of p }

  For j := 1 to n do
  begin
    wa1^[j] := -wa1^[j];
    wa2^[j] := x^[j] + wa1^[j];
    wa3^[j] := diag^[j]*wa1^[j]
  end;
  pnorm := enorm(n, wa3);

{ on the first iteration, adjust the initial step bound }

  IF iter = 1 then  delta := MIN(delta,pnorm);

{ evaluate the function at x + p and calculate its norm }

  iflag := 1;
  fcn(m, n, wa2, wa4, iflag);
  nfev := nfev + 1;
  IF iflag < 0 then GOTO 300;
  fnorm1 := enorm(m, wa4);

{ compute the scaled actual reduction }

  actred := -one;
  IF p1*fnorm1 < fnorm then actred := one - Sqr(fnorm1/fnorm);

{ compute the scaled predicted reduction and
  the scaled directional derivative }

  For j := 1 to n do
  begin
    wa3^[j] := zero;
    l := ipvt^[j];
    temp := wa1^[l];
    For i:=1 to j do wa3^[i] := wa3^[i] + fjac^[i,j]*temp
  end;
  temp1 := enorm(n,wa3)/fnorm;
  temp2 := (SQRT(par)*pnorm)/fnorm;
  prered := temp1*temp1 + temp2*temp2/p5;
  dirder := -(temp1*temp1 + temp2*temp2);

{ compute the ratio of the actual to the predicted reduction }

  ratio := zero;
  IF prered <> zero then ratio := actred/prered;

{ update the step bound }

  IF ratio <= p25 THEN
  begin
    IF actred >= zero then temp := p5;
    IF actred < zero then temp := p5*dirder/(dirder + p5*actred);
    IF (p1*fnorm1 >= fnorm) OR (temp < p1) then temp := p1;
    delta := temp*MIN(delta, pnorm/p1);
    par := par/temp
  end
  ELSE
  begin
    IF (par <> zero) AND (ratio < p75) then GOTO 260;
    delta := pnorm/p5;
    par := p5*par
  end;

{ test for successful iteration }

260:IF ratio < p0001 then GOTO 290;

{ successful iteration. update x, fvec, and their norms }

  For j := 1 to n do
  begin
    x^[j] := wa2^[j];
    wa2^[j] := diag^[j]*x^[j]
  end;
  For i:=1 to m do fvec^[i] := wa4^[i];
  xnorm := enorm(n,wa2);
  fnorm := fnorm1;
  iter := iter + 1;

{ tests for convergence }

290:IF (ABS(actred) <= ftol) AND (prered <= ftol) AND (p5*ratio <= one) then info := 1;
  IF delta <= xtol*xnorm then info := 2;
  IF (ABS(actred) <= ftol) AND (prered <= ftol) AND (p5*ratio <= one) AND
     (info = 2) then info := 3;
  IF info <> 0 then GOTO 300;

{ tests for termination and stringent tolerances }

  IF nfev >= maxfev then info := 5;
  IF (ABS(actred) <= epsmch) AND (prered <= epsmch) AND (p5*ratio <= one) then info := 6;
  IF delta <= epsmch*xnorm then info := 7;
  IF gnorm <= epsmch then info := 8;
  IF info <> 0 then GOTO 300;

{ end of the inner loop. repeat if iteration unsuccessful }

  IF ratio < p0001 then GOTO 200;

{ end of the outer loop }

  GOTO 30;

{ termination, either normal or user imposed }

300: IF iflag < 0 then info := iflag;
  iflag := 0;
  IF nprint > 0 then fcn(m, n, x, fvec, iflag);

  Dispose(diag); Dispose(qtf); Dispose(tmp); Dispose(tmp1); Dispose(wa1);
  Dispose(wa2);  Dispose(wa3); Dispose(wa4)

End; {lmder}


Procedure lmpar(n:integer; r:pMat; ipvt:pIVec; diag, qtb:pVec; delta:double;
                var par:double; x, sdiag: pVec);
{  *************************************************************************

!  procedure lmpar

!  Given an m by n matrix a, an n by n nonsingular diagonal matrix d,
!  an m-vector b, and a positive number delta, the problem is to determine a
!  value for the parameter par such that if x solves the system

!        a*x = b ,     sqrt(par)*d*x = 0 ,

!  in the least squares sense, and dxnorm is the euclidean
!  norm of d*x, then either par is zero and

!        (dxnorm-delta) <= 0.1*delta ,

!  or par is positive and

!        abs(dxnorm-delta) <= 0.1*delta .

!  This procedure completes the solution of the problem if it is provided
!  with the necessary information from the r factorization, with column
!  qpivoting, of a.  That is, if a*p := q*r, where p is a permutation matrix,
!  q has orthogonal columns, and r is an upper triangular matrix with diagonal
!  elements of nonincreasing magnitude, then lmpar expects the full upper
!  triangle of r, the permutation matrix p, and the first n components of
!  (q transpose)*b.
!  On output lmpar also provides an upper triangular matrix s such that

!         t   t                   t
!        p *(a *a + par*d*d)*p = s *s .

!  s is employed within lmpar and may be of separate interest.

!  Only a few iterations are generally needed for convergence of the algorithm.
!  If, however, the limit of 10 iterations is reached, then the output par
!  will contain the best value obtained so far.

!  the procedure statement is

!    procedure lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag, wa1,wa2)

!  where

!    n is a positive integer input variable set to the order of r.

!    r is an n by n array. on input the full upper triangle
!      must contain the full upper triangle of the matrix r.
!      On output the full upper triangle is unaltered, and the
!      strict lower triangle contains the strict upper triangle
!      (transposed) of the upper triangular matrix s.

!    ldr is a positive integer input variable not less than n
!      which specifies the leading dimension of the array r.

!    ipvt is an integer input array of length n which defines the
!      permutation matrix p such that a*p := q*r. column j of p
!      is column ipvt(j) of the identity matrix.

!    diag is an input array of length n which must contain the
!      diagonal elements of the matrix d.

!    qtb is an input array of length n which must contain the first
!      n elements of the vector (q transpose)*b.

!    delta is a positive input variable which specifies an upper
!      bound on the euclidean norm of d*x.

!    par is a nonnegative variable. on input par contains an
!      initial estimate of the levenberg-marquardt parameter.
!      on output par contains the final estimate.

!    x is an output array of length n which contains the least
!      squares solution of the system a*x := b, sqrt(par)*d*x := 0,
!      for the output par.

!    sdiag is an output array of length n which contains the
!      diagonal elements of the upper triangular matrix s.

!    wa1 and wa2 are work arrays of length n.

!  subprograms called

!    dpmpar, enorm, max, min, qrsolv

!    pascal-supplied ... ABS, SQRT

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  ************************************************************* }
Label 120,150,220;
Var
  i, iter, j, jm1, jp1, k, l, nsing: Integer;
  dxnorm, dwarf, fp, gnorm, parc, parl, paru, sum, temp: double;
  tmp, wa1, wa2: pVec;
Const
  p1 = 0.1; p001 = 0.001;
Begin

  New(tmp); New(wa1); New(wa2);

{ dwarf is the smallest positive magnitude }
  dwarf := 2.5E-16;

{ compute and store in x the gauss-newton direction. if the
  jacobian is rank-deficient, obtain a least squares solution }

  nsing := n;
  For j := 1 to n do
  begin
    wa1^[j] := qtb^[j];
    IF (r^[j,j] = zero) AND (nsing = n) then nsing := j - 1;
    IF nsing < n then wa1^[j] := zero
  end;

  For k := 1 to nsing do
  begin
    j := nsing - k + 1;
    wa1^[j] := wa1^[j]/r^[j,j];
    temp := wa1^[j];
    jm1 := j - 1;
    For i:=1 to jm1 do
      wa1^[i] := wa1^[i] - r^[i,j]*temp
  end;
  For j := 1 to n do
  begin
    l := ipvt^[j];
    x^[l] := wa1^[j]
  end;

{ initialize the iteration counter.
  evaluate the function at the origin, and test
  for acceptance of the gauss-newton direction }

  iter := 0;
  For i:=1 to n do wa2^[i] := diag^[i]*x^[i];
  dxnorm := enorm(n, wa2);
  fp := dxnorm - delta;
  IF fp <= p1*delta then GOTO 220;

{ if the jacobian is not rank deficient, the newton
  step provides a lower bound, parl, for the zero of
  the function.  Otherwise set this bound to zero }

  parl := zero;
  IF nsing < n then GOTO 120;
  For j := 1 to n do
  begin
    l := ipvt^[j];
    wa1^[j] := diag^[l]*(wa2^[l]/dxnorm)
  end;
  For j := 1 to  n do
  begin
    For i:=1 to j-1 do tmp^[i]:=r^[i,j];
    sum := DOT_PRODUCT(j-1, tmp, wa1);
    wa1^[j] := (wa1^[j] - sum)/r^[j,j]
  end;
  temp := enorm(n,wa1);
  parl := ((fp/delta)/temp)/temp;

{ calculate an upper bound, paru, for the zero of the function }

120:For j := 1 to n do
  begin
    For i:=1 to j do tmp^[i]:=r^[i,j];
    sum := DOT_PRODUCT(j, tmp, qtb);
    l := ipvt^[j];
    wa1^[j] := sum/diag^[l]
  end;
  gnorm := enorm(n,wa1);
  paru := gnorm/delta;
  IF paru = zero then paru := dwarf/MIN(delta,p1);

{ if the input par lies outside of the interval (parl,paru),
  set par to the closer endpoint }

  par := MAX(par,parl);
  par := MIN(par,paru);
  IF par = zero then par := gnorm/dxnorm;

{ beginning of an iteration }

150: iter := iter + 1;

{ evaluate the function at the current value of par }

  IF par = zero then par := MAX(dwarf, p001*paru);
  temp := SQRT(par);
  For i:=1 to n do wa1^[i] := temp*diag^[i];

  qrsolv(n, r, ipvt, wa1, qtb, x, sdiag);

  For i:=1 to n do wa2^[i] := diag^[i]*x^[i];
  dxnorm := enorm(n, wa2);
  temp := fp;
  fp := dxnorm - delta;

{ if the function is small enough, accept the current value
  of par. also test for the exceptional cases where parl
  is zero or the number of iterations has reached 10 }

  IF ((ABS(fp) <= p1*delta) OR (parl = zero)) AND ((fp <= temp)
     AND (temp < zero) OR (iter = 10)) then GOTO 220;

{ compute the newton correction }

  For j := 1 to n do
  begin
    l := ipvt^[j];
    wa1^[j] := diag^[l]*(wa2^[l]/dxnorm)
  end;
  For j := 1 to n do
  begin
    wa1^[j] := wa1^[j]/sdiag^[j];
    temp := wa1^[j];
    jp1 := j + 1;
    For i:=jp1 to n do wa1^[i] := wa1^[i] - r^[i,j]*temp
  end;
  temp := enorm(n,wa1);
  parc := ((fp/delta)/temp)/temp;

{ depending on the sign of the function, update parl or paru }

  IF fp > zero then parl := MAX(parl,par);
  IF fp < zero then paru := MIN(paru,par);

{ compute an improved estimate for par }

  par := MAX(parl, par+parc);

{ end of an iteration }

  GOTO 150;

220:IF iter = 0 then par := zero;

  Dispose(tmp); Dispose(wa1); Dispose(wa2)

End; {lmpar}



Procedure qrfac(m, n:integer; a:pMat; pivot:Boolean; ipvt:pIVec; rdiag, acnorm:pVec);
{  *************************************************************************

!  procedure qrfac

!  This procedure uses Householder transformations with column pivoting
!  (optional) to compute a qr factorization of the m by n matrix a.
!  That is, qrfac determines an orthogonal matrix q, a permutation matrix p,
!  and an upper trapezoidal matrix r with diagonal elements of nonincreasing
!  magnitude, such that a*p = q*r.  The householder transformation for
!  column k, k = 1,2,...,min(m,n), is of the form

!                        t
!        i - (1/u(k))*u*u

!  where u has zeros in the first k-1 positions.  The form of this
!  transformation and the method of pivoting first appeared in the
!  corresponding linpack procedure.

!  the procedure statement is

!    procedure qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm, wa)

! N.B. 3 of these arguments have been omitted in this version.

!  where

!    m is a positive integer input variable set to the number of rows of a.

!    n is a positive integer input variable set to the number of columns of a.

!    a is an m by n array.  On input a contains the matrix for
!      which the qr factorization is to be computed.  On output
!      the strict upper trapezoidal part of a contains the strict
!      upper trapezoidal part of r, and the lower trapezoidal
!      part of a contains a factored form of q (the non-trivial
!      elements of the u vectors described above).

!    lda is a positive integer input variable not less than m
!      which specifies the leading dimension of the array a.

!    pivot is a logical input variable.  If pivot is set true,
!      then column pivoting is enforced.  If pivot is set false,
!      then no column pivoting is done.

!    ipvt is an integer output array of length lipvt.  ipvt
!      defines the permutation matrix p such that a*p = q*r.
!      Column j of p is column ipvt(j) of the identity matrix.
!      If pivot is false, ipvt is not referenced.

!    lipvt is a positive integer input variable.  If pivot is false,
!      then lipvt may be as small as 1.  If pivot is true, then
!      lipvt must be at least n.

!    rdiag is an output array of length n which contains the
!      diagonal elements of r.

!    acnorm is an output array of length n which contains the norms of the
!      corresponding columns of the input matrix a.
!      If this information is not needed, then acnorm can coincide with rdiag.

!    wa is a work array of length n.  If pivot is false, then wa
!      can coincide with rdiag.

!  procedures or functions called:

!    dpmpar, enorm, MAX, MIN

!    pascal-supplied ... SQRT

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  *************************************************************************** }
Label 40,45,50;
Var
  i, j, jp1, k, kmax, minmn: Integer;
  ajnorm, sum, temp: double;
  tmp, tmp1, wa: pVec;
Begin

  New(tmp); New(tmp1); New(wa);

{ compute the initial column norms and initialize several arrays }

  For j := 1 to n do
  begin
    For i:=1 to m do tmp^[i]:=a^[i,j];
    acnorm^[j] := enorm(m,tmp);
    rdiag^[j] := acnorm^[j];
    wa^[j] := rdiag^[j];
    IF (pivot) then ipvt^[j] := j
  end;

{ Reduce a to r with Householder transformations }

  if m<=n then minmn := m else minmn := n;
  For j := 1 to minmn do
  begin
    IF NOT pivot then GOTO 40;
  
{ Bring the column of largest norm into the pivot position }
  
    kmax := j;
    For k := j to n do
      IF rdiag^[k] > rdiag^[kmax] then kmax := k;
    IF kmax = j then GOTO 40;
    For i := 1 to m do
    begin
      temp := a^[i,j];
      a^[i,j] := a^[i,kmax];
      a^[i,kmax] := temp
    end;
    rdiag^[kmax] := rdiag^[j];
    wa^[kmax] := wa^[j];
    k := ipvt^[j];
    ipvt^[j] := ipvt^[kmax];
    ipvt^[kmax] := k;
  
{ Compute the Householder transformation to reduce the
  j-th column of a to a multiple of the j-th unit vector }
  
40: For i:=j to m do tmp^[i-j+1]:=a^[i,j];
    ajnorm := enorm(m-j+1, tmp);
    IF ajnorm = zero then goto 50;
    IF a^[j,j] < zero then ajnorm := -ajnorm;
    For i:=j to m do a^[i,j] := a^[i,j]/ajnorm;
    a^[j,j] := a^[j,j] + one;

{ Apply the transformation to the remaining columns and update the norms }
  
    jp1 := j + 1;
    For k := jp1 to n do
    begin
      For i:=j to m do
      begin
        tmp^[i-j+1]:=a^[i,j];
        tmp1^[i-j+1]:=a^[i,k]
      end;
      sum := DOT_PRODUCT(m-j+1,tmp,tmp1);
      temp := sum/a^[j,j];
      For i:=j to m do a^[i,k] := a^[i,k] - temp*a^[i,j];
      IF (NOT pivot) OR (rdiag^[k] = zero) then goto 50;
      temp := a^[j,k]/rdiag^[k];
      rdiag^[k] := rdiag^[k]*SQRT(MAX(zero, one-temp*temp));
      IF p05*Sqr(rdiag^[k]/wa^[k]) > epsmch then goto 45;
      For i:=jp1 to m do tmp^[i-j]:=a^[i,k];
      rdiag^[k] := enorm(m-j, tmp);
      wa^[k] := rdiag^[k];
45: end;
    rdiag^[j] := -ajnorm
  end; {j loop}

50: Dispose(tmp); Dispose(tmp1); Dispose(wa)

End; {qrfac}


Procedure qrsolv(n:integer; r:pMat; ipvt:pIVec; diag, qtb, x, sdiag:pVec);
{  *************************************************************************

!  procedure qrsolv

!  Given an m by n matrix a, an n by n diagonal matrix d, and an m-vector b,
!  the problem is to determine an x which solves the system

!        a*x = b ,     d*x = 0 ,

!  in the least squares sense.

!  This procedure completes the solution of the problem if it is provided
!  with the necessary information from the qr factorization, with column
!  pivoting, of a.  That is, if a*p = q*r, where p is a permutation matrix,
!  q has orthogonal columns, and r is an upper triangular matrix with diagonal
!  elements of nonincreasing magnitude, then qrsolv expects the full upper
!  triangle of r, the permutation matrix p, and the first n components of
!  (q transpose)*b.  The system a*x = b, d*x = 0, is then equivalent to

!               t       t
!        r*z = q *b ,  p *d*p*z = 0 ,

!  where x := p*z. if this system does not have full rank,
!  then a least squares solution is obtained.  On output qrsolv
!  also provides an upper triangular matrix s such that

!         t   t               t
!        p *(a *a + d*d)*p = s *s .

!  s is computed within qrsolv and may be of separate interest.

!  the procedure statement is

!    procedure qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)

! N.B. Arguments LDR and WA have been removed in this version.

!  where

!    n is a positive integer input variable set to the order of r.

!    r is an n by n array.  On input the full upper triangle must contain
!      the full upper triangle of the matrix r.
!      On output the full upper triangle is unaltered, and the strict lower
!      triangle contains the strict upper triangle (transposed) of the
!      upper triangular matrix s.

!    ldr is a positive integer input variable not less than n
!      which specifies the leading dimension of the array r.

!    ipvt is an integer input array of length n which defines the
!      permutation matrix p such that a*p = q*r.  Column j of p
!      is column ipvt(j) of the identity matrix.

!    diag is an input array of length n which must contain the
!      diagonal elements of the matrix d.

!    qtb is an input array of length n which must contain the first
!      n elements of the vector (q transpose)*b.

!    x is an output array of length n which contains the least
!      squares solution of the system a*x = b, d*x = 0.

!    sdiag is an output array of length n which contains the
!      diagonal elements of the upper triangular matrix s.

!    wa is a work array of length n.

!  subprograms called

!    pascal-supplied ... ABS, SQRT

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  ************************************************************************* }
Var
    i, j, k, kp1, l, nsing: Integer;
    COS, cotan, qtbpj, SIN, sum, TAN, temp: double;
    tmp, tmp1, wa: pVec;
Begin

  New(tmp); New(tmp1); New(wa);

{ Copy r and (q transpose)*b to preserve input and initialize s.
  In particular, save the diagonal elements of r in x }

  For j := 1 to n do
  begin
    For i:=j to n do r^[i,j] := r^[j,i];
    x^[j] := r^[j,j];
    wa^[j] := qtb^[j]
  end;

{ Eliminate the diagonal matrix d using a givens rotation }

  For j := 1 to n do
  begin

{ Prepare the row of d to be eliminated, locating the
  diagonal element using p from the qr factorization }
  
    l := ipvt^[j];
    IF diag^[l] = zero then exit;
    For i:=j to n do sdiag^[i] := zero;
    sdiag^[j] := diag^[l];
  
{ The transformations to eliminate the row of d modify only a single
  element of (q transpose)*b beyond the first n, which is initially zero }
  
    qtbpj := zero;
    For k := j to n do
    begin
    
{ Determine a givens rotation which eliminates the
  appropriate element in the current row of d }
    
      IF sdiag^[k] = zero then exit;
      IF ABS(r^[k,k]) < ABS(sdiag^[k]) THEN
      begin
        cotan := r^[k,k]/sdiag^[k];
        SIN := p5/SQRT(p25 + p25*cotan*cotan);
        COS := SIN*cotan
      end
      ELSE
      begin
        TAN := sdiag^[k]/r^[k,k];
        COS := p5/SQRT(p25 + p25*TAN*TAN);
        SIN := COS*TAN
      end;
    
{ Compute the modified diagonal element of r and
  the modified element of ((q transpose)*b,0) }
    
      r^[k,k] := COS*r^[k,k] + SIN*sdiag^[k];
      temp := COS*wa^[k] + SIN*qtbpj;
      qtbpj := -SIN*wa^[k] + COS*qtbpj;
      wa^[k] := temp;
    
{ Accumulate the tranformation in the row of s }
    
      kp1 := k + 1;
      For i := kp1 to n do
      begin
        temp := COS*r^[i,k] + SIN*sdiag^[i];
        sdiag^[i] := -SIN*r^[i,k] + COS*sdiag^[i];
        r^[i,k] := temp
      end
    end; {k loop}
  
{ Store the diagonal element of s and restore
  the corresponding diagonal element of r }
  
    sdiag^[j] := r^[j,j];
    r^[j,j] := x^[j]
  end; {j loop}

{ Solve the triangular system for z.  If the system is singular,
  then obtain a least squares solution }

  nsing := n;
  For j := 1 to n do
  begin
    IF (sdiag^[j] = zero) AND (nsing = n) then nsing := j - 1;
    IF nsing < n then wa^[j] := zero
  end;

  For  k := 1 to nsing do
  begin
    j := nsing - k + 1;
    For i:=j+1 to nsing do
    begin
      tmp^[i-j]:=r^[i,j];
      tmp1^[i-j]:=wa^[i]
    end;
    sum := DOT_PRODUCT(nsing-j, tmp, tmp1);
    wa^[j] := (wa^[j] - sum)/sdiag^[j]
  end;

{ Permute the components of z back to components of x }

  For j := 1 to n do
  begin
    l := ipvt^[j];
    x^[l] := wa^[j]
  end;

  Dispose(tmp); Dispose(tmp1); Dispose(wa)

End; {qrsolv}


FUNCTION enorm(n:integer; x:pVec): double;
Var i:word; temp:double;
Begin
  temp:=zero;
  For i:=1 to n do temp:=temp+Sqr(x^[i]);
  enorm:=SQRT(temp)
End;


Procedure fdjac2(m, n:integer; x, fvec:pVec; fjac:pMat; var iflag:integer;
                 epsfcn: double);
{  ***********************************************************************

!  procedure fdjac2

!  this procedure computes a forward-difference approximation
!  to the m by n jacobian matrix associated with a specified
!  problem of m functions in n variables.

!  the procedure statement is

!    procedure fdjac2(m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)

!  where:

!    fcn is the name of the user-supplied procedure which calculates
!      the functions (removed from list of arguments).
!      fcn should be written as follows:

!      Procedure fcn (m, n:integer; x, fvec:pVec; var iflag:integer);
!      begin
!      ----------
!      calculate the functions at x and return this vector in fvec.
!      ----------
!      end;

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of lmdif1.
!      In this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number of functions.

!    n is a positive integer input variable set to the number of variables.
!      n must not exceed m.

!    x is an input array of length n.

!    fvec is an input array of length m which must contain the
!      functions evaluated at x.

!    fjac is an output m by n array which contains the
!      approximation to the jacobian matrix evaluated at x.

!    ldfjac is a positive integer input variable not less than m
!      which specifies the leading dimension of the array fjac.

!    iflag is an integer variable which can be used to terminate
!      the execution of fdjac2.  see description of fcn.

!    epsfcn is an input variable used in determining a suitable step length
!      for the forward-difference approximation.  This approximation assumes
!      that the relative errors in the functions are of the order of epsfcn.
!      If epsfcn is less than the machine precision, it is assumed that the
!      relative errors in the functions are of the order of the machine
!      precision.

!    wa is a work array of length m.

!  subprograms called

!    user-supplied ...... fcn, ssqfcn

!    minpack-supplied ... dpmpa, max

!    pascal-supplied ... ABS, SQRT

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

! ********************************************************** }
Var
    i,j: Integer;
    eps, h, temp: double;
    wa: pVec;
Begin

  New(wa);

  eps := SQRT(MAX(epsfcn, epsmch));
  For j := 1 to n do
  begin
    temp := x^[j];
    h := eps*ABS(temp);
    IF h = zero then h := eps;
    x^[j] := temp + h;
    fcn(m, n, x, wa, iflag);
    IF iflag < 0 then EXIT;
    x^[j] := temp;
    For i:=1 to m do fjac^[i,j] := (wa^[i] - fvec^[i])/h
  end;

  Dispose(wa)

End; {fdjac2}


END. {UNIT Levenberg_Marquardt}