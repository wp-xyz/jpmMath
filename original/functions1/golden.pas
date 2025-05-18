{*******************************************************
*        Find minimum of a real function Y=F(X)        *
*        using routine for Golden Section Search       *
* ---------------------------------------------------- *
* REFERENCE: "Numerical Recipes, The Art of Scientific *
*             Computing By W.H. Press, B.P. Flannery,  *
*             S.A. Teukolsky and W.T. Vetterling,      *
*             Cambridge University Press, 1986"        *
*             [BIBLI 08].                              *
* ---------------------------------------------------- *
* Sample run:                                          *
*                                                      *
* Find a minimum of function X*SIN(X)-2*COS(X)         *
* between X=4 and X=6.                                 *
*                                                      *
* X1=4  X2=5  X3=6  TOL=1e-6                           *
*                                                      *
* Function minimum is: -5.53452877400217E+0000         *
* for X=5.23293710661618E+0000                         *
*                                                      *
* ---------------------------------------------------- *
*                        Pascal version by J-P Moreau. *
*                             (www.jpmoreau.fr)        *
*******************************************************}
PROGRAM SEEK_MINIMUM;
Uses WinCrt;

Var
        X1,X2,X3,XMINI,YMINI,TOL : DOUBLE;

{Function to be analyzed}
FUNCTION F(X:DOUBLE):DOUBLE;
Begin
  F:=X*SIN(X)-2.0*COS(X)
End;

FUNCTION GOLDEN(AX,BX,CX,TOL:DOUBLE; VAR XMIN:DOUBLE):DOUBLE;
{Given a function F, and given a bracketing triplet of abscissas 
 AX, BX, CX (such that BX is between AX and CX, and F(BX) is less 
 than both F(AX) and F(CX)), this routine performs a golden section 
 search for the minimum, isolating it to a fractional precision of 
 about TOL. The abscissa of the minimum is returned as XMIN, and the minimum
 function value is returned as GOLDEN, the returned function value.}
VAR   f0,f1,f2,f3,x0,x1,x2,x3:DOUBLE;
      R,C:DOUBLE;
Begin
  R:=0.61803399; C:=1.0-R;    {golden ratios}
  X0:=AX;  {At any given time we will keep trace of 4 points: X0,X1,X2,X3}
  X3:=CX;
  IF ABS(CX-BX) > ABS(BX-AX) THEN
  begin
    X1:=BX; X2:=BX+C*(CX-BX)
  end
  ELSE
  begin
    X2:=BX; X1:=BX-C*(BX-AX)
  end;
  F1:=F(X1); F2:=F(X2);    {Initial function evaluations}
  {main loop}
  While ABS(X3-X0) > TOL*(ABS(X1)+ABS(X2)) Do
  begin
    IF F2 < F1 THEN
    begin
      X0:=X1; X1:=X2;
      X2:=R*X1+C*X3;
      F0:=F1; F1:=F2;
      F2:=F(X2)
    end
    ELSE
    begin
      X3:=X2; X2:=X1;
      X1:=R*X2+C*X0;
      F3:=F2; F2:=F1;
      F1:=F(X1)
    end;
  end;

  IF F1 < F2 THEN
  begin
    GOLDEN:=F1;
    XMIN:=X1
  end
  ELSE
  begin
    GOLDEN:=F2;
    XMIN:=X2
  end
END; {Function Golden}


{main}
BEGIN

  X1:=4.0; X2:=5.0; X3:=6.0;
  TOL:=1e-6;

  YMINI:=GOLDEN(X1,X2,X3,TOL,XMINI);

  writeln;
  writeln(' Function minimum is ',YMINI);
  writeln;
  writeln(' for X=',XMINI);
  writeln;

  Readkey; DoneWinCrt

END.

  
{end of file Golden.pas}