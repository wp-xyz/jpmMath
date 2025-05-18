{*******************************************************
*        Find minimum of a real function Y=F(X)        *
* ---------------------------------------------------- *
* REFERENCE: "Numerical Recipes, The Art of Scientific *
*             Computing by W.H. Press, B.P. Flannery,  *
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
* Function minimum is: -5.53452877399827E+0000         *
* for X=5.23293671903524E+0000                         *
*                                                      *
* ---------------------------------------------------- *
*                        Pascal version by J-P Moreau. *
*                             (www.jpmoreau.fr)        *
*******************************************************}
PROGRAM SEEK_MINIMUM;
Uses WinCrt;

CONST ITMAX=100; CGOLD=0.3819660; ZEPS=1.0E-10;
{Maximum allowed number of iterations; golden ration; and a small
 number which protects against trying to achieve fractional accuracy
 for a minimum that happens to be exactly zero.}

Var
        X1,X2,X3,XMINI,YMINI,TOL : DOUBLE;

{Function to be analyzed}
FUNCTION F(X:DOUBLE):DOUBLE;
Begin
  F:=X*SIN(X)-2.0*COS(X)
End;

FUNCTION MIN(a,b:DOUBLE):DOUBLE;
Begin
  if a<=b then MIN:=a else MIN:=b
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

FUNCTION BRENT(AX,BX,CX,TOL:DOUBLE;Var XMIN:DOUBLE):DOUBLE;
Label    1,2,3;
Var      a,b,d,e,fu,fv,fx,fw,u,v,w,x : DOUBLE;
         etemp,p,q,r,tol1,tol2,xm : DOUBLE;
         iter : INTEGER;
Begin
{Given a function F, and a bracketing triplet of abscissas
 AX,BX,CX (such that BX is between AX and CX, and F(BX) is less 
 than both F(AX) and F(CX)), this routine isolates the minimum 
 to a fractional precision of about TOL using Brent's method.
 The abscissa of the minimum is returned in XMIN, and the minimum
 function value is returned as BRENT, the returned function value.}
A:=MIN(AX,CX);
B:=MAX(AX,CX);
V:=BX;
W:=V;
X:=V;
D:=0.0;
E:=0.0;
FX:=F(X);
FV:=FX;
FW:=FX;
For ITER:=1 to ITMAX do                                  {main loop}
begin
  XM:=0.5*(A+B);
  TOL1:=TOL*ABS(X)+ZEPS;
  TOL2:=2.*TOL1;
  IF ABS(X-XM)<=(TOL2-0.5*(B-A)) THEN GOTO 3;   {Test for done here}
  IF ABS(E)>TOL1 THEN              {Construct a trial parabolic fit}
  begin
    R:=(X-W)*(FX-FV);
    Q:=(X-V)*(FX-FW);
    P:=(X-V)*Q-(X-W)*R;
    Q:=2.0*(Q-R);            {bug corrected 07/24/2006 (0.2 instead of 2.0) }
    IF Q>0 then  P:=-P;
    Q:=ABS(Q);
    ETEMP:=E;
    E:=D;
    IF (ABS(P)>=ABS(0.5*Q*ETEMP)) OR (P<=Q*(A-X)) OR (P>=Q*(B-X)) THEN GOTO 1;
{   The above conditions determine the acceptability of the 
    parabolic fit. Here it is o.k.:                         }
    D:=P/Q;
    U:=X+D;
    IF (U-A<TOL2) OR (B-U<TOL2) THEN D:=SIGN(TOL1,XM-X);
    GOTO 2;
  end;
1:IF X>=XM THEN
    E:=A-X
  ELSE
    E:=B-X;
  D:=CGOLD*E;
2:IF ABS(D)>=TOL1 THEN
    U:=X+D
  ELSE
    U:=X+SIGN(TOL1,D);
  FU:=F(U);    {This the one function evaluation per iteration}
  IF FU<=FX THEN
  begin
    IF U>=X THEN
      A:=X
    ELSE
      B:=X;
    V:=W;
    FV:=FW;
    W:=X;
    FW:=FX;
    X:=U;
    FX:=FU
  end
  ELSE
  begin
    IF U<X THEN 
      A:=U
    ELSE
      B:=U;
    IF (FU<=FW) OR (W=X) THEN
    begin
      V:=W;
      FV:=FW;
      W:=U;
      FW:=FU
    end
    ELSE IF (FU<=FV) OR (V=X) OR (V=W) THEN
    begin
      V:=U;
      FV:=FU
    end
  end
end;
writeln(' Brent exceed maximum iterations.');
3:XMIN:=X;   {exit section}
  BRENT:=FX
End;

{main}
BEGIN

  X1:=4.0; X2:=5.0; X3:=6.0;
  TOL:=1e-6;

  YMINI:=BRENT(X1,X2,X3,TOL,XMINI);

  writeln;
  writeln(' Function minimum is ',YMINI);
  writeln;
  writeln(' for X=',XMINI);
  writeln;

  Readkey; DoneWinCrt

END.

  
!end of file Brent.pas