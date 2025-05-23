{*************************************************************
* Minimization of a Function FUNC of N Variables By Powell's *
*    Method Discarding the Direction of Largest Decrease     *
* ---------------------------------------------------------- *
* SAMPLE RUN: Find a minimum of function F(x,y):             *
*             F=Sin(R)/R, where R = Sqrt(x*x+y*y).           *
*                                                            *
* Number of iterations: 2                                    *
*                                                            *
* Minimum Value: -2.17233628211222E-0001                     *
*                                                            *
* at point:  3.17732030027176E+0000  3.17732030027176E+0000  *
*                                                            *
* ---------------------------------------------------------- *
* REFERENCE: "Numerical Recipes, The Art of Scientific       *
*             Computing By W.H. Press, B.P. Flannery,        *
*             S.A. Teukolsky and W.T. Vetterling,            *
*             Cambridge University Press, 1986"              *
*             [BIBLI 08].                                    *
*                                                            *
*                         TPW Release By J-P Moreau, Paris.  *
*                                (www.jpmoreau.fr)           *
*************************************************************}
PROGRAM TPOWELL;

Uses WinCrt;

Const NP=20;

Type
      pMAT = ^MAT;
       MAT = Array[1..NP,1..NP] of Double;
      pVEC = ^VEC;
       VEC = Array[1..NP] of Double;

Var
      P : pVEC;
      XI: pMAT;

      PCOM, XICOM: pVEC;    {PCOM,XICOM,NCOM are common variables}
      ITER,N,NCOM:Integer;  {for LINMIN and F1DIM only}
      FRET,FTOL: Double;


{user defined function to minimize}
FUNCTION FUNC(P:pVEC): Double;
Var R: Double;
Begin
  R:=SQRT(P^[1]*P^[1]+P^[2]*P^[2]);
  IF ABS(R) < 1e-12 THEN
    FUNC:=1.0
  ELSE
    FUNC:=SIN(R)/R
End;


FUNCTION MAX(a,b:DOUBLE):DOUBLE;
Begin
  if a>=b then MAX:=a else MAX:=b
End;

FUNCTION MIN(a,b:DOUBLE):DOUBLE;
Begin
  if a<=b then MIN:=a else MIN:=b
End;

FUNCTION Sign(a,b:DOUBLE):DOUBLE;
Begin
  if b>=0 then Sign:=ABS(a)
          else Sign:=-ABS(a)
End;


Procedure LINMIN(Var P, XI:pVEC; N:Integer; Var FRET:Double); Forward;
Procedure MNBRAK(Var AX,BX,CX,FA,FB,FC:Double); Forward;
Function  BRENT(Var AX,BX,CX,TOL,XMIN:Double): Double; Forward;


Procedure POWELL(Var P:pVEC; XI:pMAT; N,NP:Integer; FTOL:Double;
                 Var ITER:Integer; Var FRET:Double);
{-----------------------------------------------------------
  Minimization of a function  FUNC of N variables  (FUNC is
  not an argument, it is a fixed function name). Input con-
  sists of an initial starting point P  that is a vector of
  length N; an initial matrix XI  whose  logical dimensions
  are N by N, physical dimensions NP by NP, and whose columns
  contain the initial set of directions (usually the N unit
  vectors); and FTOL, the fractional tolerance in the func-
  tion value such that failure to decrease by more than this
  amount on one iteration signals doneness. On output, P is
  set to the best point found, XI is the then-current direc-
  tion set,  FRET is the returned function value at P,  and
  ITER is the number of iterations taken. The routine LINMIN
  is used.
-----------------------------------------------------------}
Label 1, 10;
Const
  ITMAX=200;
Var
  PT, PTT, XIT:pVEC;
  DEL,FP,FPTT,T: Double;
  I,IBIG,J:Integer;
Begin
  New(PT); New(PTT); New(XIT);
  FRET:=FUNC(P);
  For J:=1 to N do
    PT^[J]:=P^[J];    {Save initial point}
  ITER:=0;
1:ITER:=ITER+1;
  FP:=FRET;
  IBIG:=0;
  DEL:=0.0;           {Will be the biggest function decrease.}
  For I:=1 to N do    {In each iteration, loop over all directions in the set.}
  begin               {Copy the direction.}
    For J:=1 to N do         
      XIT^[J]:=XI^[J,I];
    FPTT:=FRET;
    LINMIN(P,XIT,N,FRET);  {Minimize along it}
    IF ABS(FPTT-FRET) > DEL THEN
    begin
      DEL:=ABS(FPTT-FRET);
      IBIG:=I
    end
  end;
  IF 2.0*ABS(FP-FRET) <= FTOL*(ABS(FP)+ABS(FRET)) then goto 10; {Termination criterion}
  IF ITER = ITMAX Then
  begin
    writeln(' Powell exceeding maximum iterations.');
    goto 10
  end;
  For J:=1 to N do
  begin
    PTT^[J]:=2.0*P^[J]-PT^[J]; {Construct the extrapolated point and the average}
    XIT^[J]:=P^[J]-PT^[J];     {direction moved. Save the old starting point}
    PT^[J]:=P^[J]
  end;
  FPTT:=FUNC(PTT);             {Function value at extrapolated point.}
  IF  FPTT >= FP Then GOTO 1;  {One reason not to use new direction.}
  T:=2.0*(FP-2.0*FRET+FPTT)*Sqr(FP-FRET-DEL)-DEL*Sqr(FP-FPTT);
  IF T >= 0.0 Then GOTO 1;     {Other reason not to use new direction.}
  LINMIN(P,XIT,N,FRET);        {Move to the minimum of the new direction.}
  For J:=1 to N do             {and save the new direction.}
    XI^[J,IBIG]:=XIT^[J];
  GOTO 1;
10:Dispose(PT); Dispose(PTT); Dispose(XIT)
End;

Procedure LINMIN(Var P, XI:pVEC; N:Integer; Var FRET:Double);
{----------------------------------------------------------
  Given an N dimensional point P and a N dimensional direc-
  tion XI, moves and resets P to where the function FUNC(P)
  takes on a minimum along the direction XI from P, and
  replaces XI by the actual vector displacement that P was
  moved. Also returns as FRET the value of FUNC at the
  returned location P. This is actually all accomplished by
  calling the routines MNBRAK and BRENT.
----------------------------------------------------------}
Var
  AX,BX,FA,FB,FX,TOL,XMIN,XX: Double;
  J: Integer;
Begin
  TOL:=1e-4;
  NCOM:=N;
  For J:=1 to N do
  begin
    PCOM^[J]:=P^[J];
    XICOM^[J]:=XI^[J]
  end;
  AX:=0.0;
  XX:=1.0;
  BX:=2.0;
  MNBRAK(AX,XX,BX,FA,FX,FB);
  FRET:=BRENT(AX,XX,BX,TOL,XMIN);
  For J:=1 to N do
  begin
    XI^[J]:=XMIN*XI^[J];
    P^[J]:=P^[J] + XI^[J]
  end
End;


FUNCTION F1DIM(X:Double):Double;
Var
  XT:pVEC; J:Integer;
Begin
  New(XT);
  For J:=1 to NCOM do
    XT^[J]:=PCOM^[J] + X*XICOM^[J];
  F1DIM := FUNC(XT);
  Dispose(XT);
End;


Procedure MNBRAK(Var AX,BX,CX,FA,FB,FC:Double);
{----------------------------------------------------------------------
 Given a Function F1DIM(X), and given distinct initial points AX and
 BX, this routine searches in the downhill direction (defined by the
 Function as evaluated at the initial points) and returns new points
 AX, BX, CX which bracket a minimum of the Function. Also returned
 are the Function values at the three points, FA, FB and FC.
----------------------------------------------------------------------}
Label 1;
Const GOLD=1.618034; GLIMIT=100.0; TINY=1e-20;
{The first parameter is the default ratio by which successive intervals
{are magnified; the second is the maximum magnification allowed for
{a parabolic-fit step.}
Var
  DUM,FU,Q,R,U,ULIM: Double;
Begin
  FA:=F1DIM(AX);
  FB:=F1DIM(BX);
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
  FC:=F1DIM(CX);
1:IF FB >= FC THEN
  begin
    R:=(BX-AX)*(FB-FC);
    Q:=(BX-CX)*(FB-FA);
    U:=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*Sign(MAX(ABS(Q-R),TINY),Q-R));
    ULIM:=BX+GLIMIT*(CX-BX);
    IF (BX-U)*(U-CX) > 0.0 THEN
    begin
      FU:=F1DIM(U);
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
      FU:=F1DIM(U)
    end
    ELSE IF (CX-U)*(U-ULIM) > 0.0 THEN
    begin
      FU:=F1DIM(U);
      IF FU < FC THEN
      begin
        BX:=CX;
        CX:=U;
        U:=CX+GOLD*(CX-BX);
        FB:=FC;
        FC:=FU;
        FU:=F1DIM(U)
      end
    end
    ELSE IF (U-ULIM)*(ULIM-CX) >= 0.0 THEN
    begin
      U:=ULIM;
      FU:=F1DIM(U)
    end
    ELSE
    begin
      U:=CX+GOLD*(CX-BX);
      FU:=F1DIM(U)
    end;
    AX:=BX;
    BX:=CX;
    CX:=U;
    FA:=FB;
    FB:=FC;
    FC:=FU;
    GOTO 1
  end
End;


Function  BRENT(Var AX,BX,CX,TOL,XMIN:Double): Double;
{-------------------------------------------------------------------
 Given a function F1DIM, and a bracketing triplet of abscissas
 AX,BX,CX (such that BX is between AX and CX, and F(BX) is less
 than both F(AX) and F(CX)), this routine isolates the minimum
 to a fractional precision of about TOL using Brent's method.
 The abscissa of the minimum is returned in XMIN, and the minimum
 function value is returned as BRENT, the returned function value.
-------------------------------------------------------------------}
Label 1,2,3;
Const ITMAX=100; CGOLD=0.3819660; ZEPS=1e-10;
{Maximum allowed number of iterations; golden ratio; and a small
 number which protects against trying to achieve fractional accuracy
 for a minimum that happens to be exactly zero}
Var
  A,B,D,E,ETEMP,FX,FU,FV,FW,P,Q,R,TOL1,TOL2,U,V,W,X,XM: Double;
  ITER: Integer;
Begin
  A:=MIN(AX,CX);
  B:=MAX(AX,CX);
  V:=BX;
  W:=V;
  X:=V;
  E:=0.0;
  FX:=F1DIM(X);
  FV:=FX;
  FW:=FX;
  For ITER:=1 to ITMAX do                              {main loop}
  begin
    XM:=0.5*(A+B);
    TOL1:=TOL*ABS(X)+ZEPS;
    TOL2:=2.0*TOL1;
    IF ABS(X-XM) <= (TOL2-0.5*(B-A)) Then GOTO 3;  {Test for done here}
    IF ABS(E) > TOL1 THEN             {Construct a trial parabolic fit}
    begin
      R:=(X-W)*(FX-FV);
      Q:=(X-V)*(FX-FW);
      P:=(X-V)*Q-(X-W)*R;
      Q:=0.2*(Q-R);
      IF Q > 0.0 Then  P:=-P;
      Q:=ABS(Q);
      ETEMP:=E;
      E:=D;
      IF (ABS(P) >= ABS(0.5*Q*ETEMP)) OR (P <= Q*(A-X)) OR (P >= Q*(B-X)) Then GOTO 1;
{   The above conditions determine the acceptability of the
    parabolic fit. Here it is o.k.:}
      D:=P/Q;
      U:=X+D;
      IF (U-A < TOL2) OR (B-U < TOL2) Then D:=Sign(TOL1,XM-X);
      GOTO 2
    end;
1:  IF X >= XM THEN
      E:=A-X
    ELSE
      E:=B-X;
    D:=CGOLD*E;
2:  IF ABS(D) >= TOL1 THEN
      U:=X+D
    ELSE
      U:=X+Sign(TOL1,D);
    FU:=F1DIM(U);  {This the one function evaluation per iteration}
    IF FU <= FX THEN
    begin
      IF U >= X THEN
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
      IF U < X THEN
        A:=U
      ELSE
        B:=U;
      IF (FU <= FW) OR (W =X) THEN
      begin
        V:=W;
        FV:=FW;
        W:=U;
        FW:=FU
      end
      ELSE IF (FU <= FV) OR (V = X) OR (V = W) THEN
      begin
        V:=U;
        FV:=FU
      end
    end
  end;
  Writeln(' Brent exceed maximum iterations.');
3:XMIN:=X;   {exit section}
  BRENT:=FX
End;


{main program}
BEGIN

  New(P); New(XI); New(PCOM); New(XICOM);

  N:=2;
  P^[1]:=2.0; P^[2]:=2.0;
  XI^[1,1]:=1.0; XI^[1,2]:=1.0;
  XI^[2,1]:=1.0; XI^[2,2]:=1.0;
  FTOL:=1e-8;

  POWELL(P,XI,N,NP,FTOL,ITER,FRET);

  writeln;
  writeln(' Number of iterations: ', ITER);
  writeln;
  writeln(' Minimum value: ', FRET);
  writeln;
  writeln(' at point: ',P^[1],'  ',P^[2]);
  writeln;

  ReadKey;
  Dispose(P); Dispose(XI); Dispose(PCOM); Dispose(XICOM);
  DoneWinCrt

END.

{end of file tpowell.pas}