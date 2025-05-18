{****************************************************
*      Program to demonstrate the real domain       *
*               Zbrent subroutine                   *
* ------------------------------------------------- *
* Reference:  BORLAND MATHEMATICAL LIBRARY          *
*                                                   *
*                TPW version by J-P Moreau, Paris.  *
*                       (www.jpmoreau.fr)           *
* ------------------------------------------------- *
* Example:    Find a real root of f(x)=(x+1)^5      *
*                                                   *
* SAMPLE RUN:                                       *
*                                                   *
*  Input interval (X1,X2):                          *
*                                                   *
*        X1 = -2                                    *
*        X2 =  0                                    *
*                                                   *
*  Convergence criterion: 1e-10                     *
*  Maximum number of iterations: 10                 *
*                                                   *
*  The estimated root is:                           *
*                                                   *
*        X = -1.00000000000000E+000                 *
*                                                   *
*  The associated Y value is Y = 0.00000000000E+000 *
*                                                   *
*  The number of iterations was: 2                  *
*  The error code is: 0                             *
*                                                   *
****************************************************}
PROGRAM DEMO_ZBRENT;
Uses WinCrt;

Var
        e,r,x1,x2,yr  : DOUBLE;
        err,maxiter,k : INTEGER;


{Test function for Brent method}
FUNCTION AFunction(x:DOUBLE): DOUBLE;
Begin
  AFunction := (x+1)*(x+1)*(x+1)*(x+1)*(x+1)
End;

{TRUE if x1*x2 negative}
FUNCTION RootBracketed(x1,x2: DOUBLE): BOOLEAN;
VAR result: BOOLEAN;
BEGIN
  IF ((x1 > 0.0) AND (x2 > 0.0)) OR
    ((x1 < 0.0) AND (x2 < 0.0)) THEN
    result := FALSE
  ELSE
    result := TRUE;
  RootBracketed := result;
END;

{returns the minimum of two real numbers}
FUNCTION Minimum(x1,x2: DOUBLE):DOUBLE;
VAR result: DOUBLE;
BEGIN
  IF x1 < x2 THEN result := x1 ELSE result := x2;
  Minimum := result;
END;

{****************************************************
*              Brent Method Function                *
* ------------------------------------------------- *
* The purpose is to find a real root of a real      *
* function f(x) using Brent method.                 *
*                                                   *
* INPUTS:  x1,x2     : interval of root             *
*          Tolerance : desired accuracy for root    *
*          maxIter   : maximum number of iterations *
*                                                   *
* OUTPUTS: The function returns the root value      *
*          ValueAtRoot : value of f(root)           *
*          niter    : number of done iterations     *
*          error    : =0, all OK                    *
*                   : =1, no root found in interval *
*                   : =2, no more iterations !      *
****************************************************}  
FUNCTION BrentRoots(x1, x2, Tolerance: DOUBLE;
                    maxIterations: INTEGER;
		    VAR  valueAtRoot: DOUBLE;
                    VAR  niter : INTEGER;
		    VAR  error: INTEGER ): DOUBLE;
CONST FPP = 1.0E-11;
      nearzero = 1.0E-20;

VAR result, AA, BB, CC, DD, EE, FA, FB, FC, Tol1, PP, QQ, RR, SS, XM: DOUBLE;
    i,j,k: INTEGER;
    done: BOOLEAN;

BEGIN
  i := 0; done := FALSE;   error := 0;
  AA := x1;  BB := x2;  FA := AFunction(AA); FB := AFunction(BB);
  IF NOT(rootBracketed(FA,FB)) THEN
    error := 1
  ELSE
  BEGIN
    FC := FB;
    REPEAT
      IF NOT(RootBracketed(FC,FB)) THEN
      BEGIN
        CC := AA; FC := FA; DD := BB - AA; EE := DD;
      END;
      IF ABS(FC) < ABS(FB) THEN
      BEGIN
        AA := BB; BB := CC; CC := AA;
        FA := FB; FB := FC; FC := FA;
      END;
      Tol1 := 2.0 * FPP * ABS(BB) + 0.5 * Tolerance;
      XM := 0.5 * (CC-BB);
      IF (ABS(XM) <= Tol1) OR (ABS(FA) < nearzero) THEN
      BEGIN
        result := BB;
        done := TRUE;
        valueAtRoot := AFunction(result);
      END
      ELSE
      BEGIN
        IF (ABS(EE) >= Tol1) AND (ABS(FA) > ABS(FB)) THEN
        BEGIN
          SS := FB/ FA;
          IF ABS(AA - CC) < nearzero THEN
          BEGIN
            PP := 2.0 * XM * SS;
            QQ := 1.0 - SS;
          END
          ELSE
          BEGIN
            QQ := FA/FC;
            RR := FB /FC;
            PP := SS * (2.0 * XM * QQ * (QQ - RR) - (BB-AA) * (RR - 1.0));
            QQ := (QQ - 1.0) * (RR - 1.0) * (SS - 1.0);
          END;
          IF PP > nearzero THEN QQ := -QQ;
          PP := ABS(PP);
          IF (2.0 * PP) < Minimum(3.0*XM *QQ-ABS(Tol1 * QQ), ABS(EE * QQ)) THEN
          BEGIN
            EE := DD;  DD := PP/QQ;
          END
          ELSE
          BEGIN
            DD := XM;   EE := DD;
          END;
        END
        ELSE
        BEGIN
          DD := XM;
          EE := DD;
        END;
        AA := BB;
        FA := FB;
        IF ABS(DD) > Tol1 THEN
          BB := BB + DD
        ELSE
          IF xm > 0.0 THEN BB := BB + ABS(Tol1)
          ELSE BB := BB - ABS(Tol1);
        FB := AFunction(BB);
        Inc(i);
      END;
    UNTIL done OR (i = maxIterations);
    IF i = maxIterations THEN error := 2;
  END;
  BrentRoots := result;
  niter := i
END;


{main program}
BEGIN
  Clrscr;
  Writeln;
  Writeln(' Input interval (X1,X2):');
  Writeln;
  Write('       X1 = '); read(x1);
  Write('       X2 = '); read(x2);
  Writeln;
  Write(' Convergence criterion: '); read(e);
  Write(' Maximum number of iterations: '); read(maxiter);

  r := BrentRoots(x1,x2,e,maxiter,yr,k,err);

  Writeln;
  Writeln(' The estimated root is:');
  Writeln;
  Writeln('       X = ', r);
  Writeln;
  Writeln(' The associated Y value is Y = ', yr);
  Writeln;
  Writeln(' The number of iterations was: ',k);  
  Writeln(' The error code is: ',err);
  Writeln;
  Readkey; DoneWinCrt
END.

{End of file zbrent.pas}