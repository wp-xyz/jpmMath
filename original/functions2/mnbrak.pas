{*******************************************************
*    Bracketing a minimum of a real function Y=F(X)    *
*             using MNBRAK subroutine                  *
* ---------------------------------------------------- *
* REFERENCE: "Numerical Recipes, The Art of Scientific *
*             Computing By W.H. Press, B.P. Flannery,  *
*             S.A. Teukolsky and W.T. Vetterling,      *
*             Cambridge University Press, 1986"        *
*             [BIBLI 08].                              *
* ---------------------------------------------------- *
* Sample run:                                          *
*                                                      *
* Find 3 points of function X*SIN(X)-2*COS(X)          *
* bracketing a minimum, given two initial points.      *
*                                                      *
* X1=4.0  X2=5.0                                       *
*                                                      *
* The three points are:                                *
*                                                      *
* X1=  4.000000  X2=  5.000000  X3=  6.618034          *
*                                                      *
* Corresponding function values:                       *
*                                                      *
* F1= -1.719923  F2= -5.361946  F3=  0.285940          *
*                                                      *
*                       Pascal version by J-P Moreau.  *
*                             (www.jpmoreau.fr)        *
*******************************************************}
PROGRAM TEST_MNBRAK;
Uses WinCrt;

VAR  X1, X2, X3, F1, F2, F3 : DOUBLE;


{Function to be analyzed}
FUNCTION FUNC(X:DOUBLE):DOUBLE;
Begin
  FUNC:=X*SIN(X)-2.0*COS(X)
End;  

FUNCTION MAX(a,b:DOUBLE):DOUBLE;
Begin
  if a>=b then MAX:=a else MAX:=b
End;

FUNCTION Sign(a,b:DOUBLE):DOUBLE;
Begin
  if b>=0 then Sign:=ABS(a)
          else Sign:=-ABS(a)
End;


PROCEDURE MNBRAK(AX,BX:DOUBLE; VAR CX,FA,FB,FC:DOUBLE);
{Given a function FUNC(X), and given distinct initial points AX and
 BX, this routine searches in the downhill direction (defined by the
 function as evaluated at the initial points) and returns new points
 AX, BX, CX which bracket a minimum of the function. Also returned
 are the function values at the three points, FA, FB and FC.}
CONST GOLD=1.618034; GLIMIT=100.0; TINY=1E-20;
{The first parameter is the default ratio by which successive intervals
 are magnified; the second is the maximum magnification allowed for
 a parabolic-fit step.}
LABEL 1;
VAR  dum,fu,q,r,u,ulim : DOUBLE;
Begin
FA:=FUNC(AX);
FB:=FUNC(BX);
IF FB > FA THEN
begin
  DUM:=AX;
  AX:=BX;
  BX:=DUM;
  DUM:=FB;
  FB:=FA;
  FA:=DUM
end;
CX:=BX+GOLD*(BX-AX);
FC:=FUNC(CX);
1:IF FB>=FC THEN
begin
  R:=(BX-AX)*(FB-FC);
  Q:=(BX-CX)*(FB-FA);
  U:=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R));
  ULIM:=BX+GLIMIT*(CX-BX);
  IF (BX-U)*(U-CX)>0 THEN
  begin
    FU:=FUNC(U);
    IF FU < FC THEN
    begin
      AX:=BX;
      FA:=FB;
      BX:=U;
      FB:=FU;
      GOTO 1
    end
    ELSE IF FU > FB THEN
    begin
      CX:=U;
      FC:=FU;
      GOTO 1
    end;
    U:=CX+GOLD*(CX-BX);
    FU:=FUNC(U)
  end
  ELSE IF (CX-U)*(U-ULIM) > 0 THEN
  begin
    FU:=FUNC(U);
    IF FU < FC THEN
    begin
      BX:=CX;
      CX:=U;
      U:=CX+GOLD*(CX-BX);
      FB:=FC;
      FC:=FU;
      FU:=FUNC(U)
    end
  end
  ELSE IF (U-ULIM)*(ULIM-CX)>=0 THEN
  begin
    U:=ULIM;
    FU:=FUNC(U)
  end
  ELSE
  begin
    U:=CX+GOLD*(CX-BX);
    FU:=FUNC(U)
  end;
  AX:=BX;
  BX:=CX;
  CX:=U;
  FA:=FB;
  FB:=FC;
  FC:=FU;
  GOTO 1
end;
End;

{main program}
BEGIN
  X1:=4.0; X2:=5.0;

  MNBRAK(X1,X2,X3,F1,F2,F3);

  writeln;
  writeln(' The three points are:');
  writeln;
  writeln(' X1:=',X1:10:6,'  X2:=',X2:10:6,'  X3:=',X3:10:6);
  writeln;
  writeln(' Corresponding function values:');
  writeln;
  writeln(' F1:=',F1:10:6,'  F2:=',F2:10:6,'  F3:=',F3:10:6);
  writeln;
  ReadKey; DoneWinCrt
END.

{End of file mnbrak.pas}