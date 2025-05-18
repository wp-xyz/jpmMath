{**********************************************************
*   Multidimensional minimization of a function FUNC(X)   *
*  where X is an NDIM-dimensional vector, by the downhill *
*  simplex method of Nelder and Mead.                     *
* ------------------------------------------------------- *
* SAMPLE RUN: Find a minimum of function F(x,y):          *
*             F=Sin(R)/R, where R = Sqrt(x*x+y*y).        *
*                                                         *
*  Number of iterations: 22                               *
*                                                         *
*  Best NDIM+1 points:                                    *
*   4.12268589437008E+0000  1.78691533207893E+0000        *
*   4.16647735983133E+0000  1.68269817531109E+0000        *
*   4.14245440904051E+0000  1.74117587320507E+0000        *
*                                                         *
*  Best NDIM+1 mimimum values:                            *
*  -2.17233626514560E-0001                                *
*  -2.17233628107148E-0001                                *
*  -2.17233627137833E-0001                                *
*                                                         *
* ------------------------------------------------------- *
* REFERENCE: "Numerical Recipes, The Art of Scientific    *
*             Computing By W.H. Press, B.P. Flannery,     *
*             S.A. Teukolsky and W.T. Vetterling,         *
*             Cambridge University Press, 1986"           *
*             [BIBLI 08].                                 *
*                                                         *
*                      TPW Release By J-P Moreau, Paris.  *
*                             (www.jpmoreau.fr)           *
**********************************************************} 
PROGRAM TEST_AMOEBA;

Uses WinCrt;

Const
      MP=21; NP=20;   {Maximum value for NDIM=20} 

Type
      MAT = Array[1..MP,1..NP] of Double;
      VEC = Array[1..MP] of Double;

Var
      P: MAT;
      Y, PT: VEC;
      I,ITER,J,NDIM: Integer;
      FTOL: Double;

{user defined function to minimize}
FUNCTION FUNC(P:VEC): Double;
Var R: double;
Begin
  R:=SQRT(P[1]*P[1]+P[2]*P[2]);
  IF ABS(R) < 1e-12 THEN
    FUNC:=1.0
  ELSE
    FUNC:=SIN(R)/R
End;


Procedure AMOEBA(Var P:MAT; Var Y:VEC; MP,NP,NDIM:Integer; FTOL:Double;
                 Var ITER:Integer);
{-------------------------------------------------------------------
! Multidimensional minimization of the function FUNC(X) where X is
! an NDIM-dimensional vector, by the downhill simplex method of
! Nelder and Mead. Input is a matrix P whose NDIM+1 rows are NDIM-
! dimensional vectors which are the vertices of the starting simplex
! (Logical dimensions of P are P(NDIM+1,NDIM); physical dimensions
! are input as P(NP,NP)). Also input is the vector Y of length NDIM
! +1, whose components must be pre-initialized to the values of FUNC
! evaluated at the NDIM+1 vertices (rows) of P; and FTOL the fractio-
! nal convergence tolerance to be achieved in the function value. On
! output, P and Y will have been reset to NDIM+1 new points all within
! FTOL of a minimum function value, and ITER gives the number of ite-
! rations taken.
!-------------------------------------------------------------------}
Label 1, 10;
Const NMAX=20; ALPHA=1.0; BETA=0.5; GAMMA=2.0; ITMAX=500;
{ Expected maximum number of dimensions, three parameters which define
  the expansions and contractions, and maximum allowed number of
  iterations. }
Var
  PR, PRR, PBAR: VEC;
  I,IHI,ILO,INHI,J,MPTS: Integer;
  RTOL,YPR,YPRR: Double;
Begin
  MPTS:=NDIM+1;
  ITER:=0;
1:ILO:=1;
  IF Y[1] > Y[2] THEN
  begin
    IHI:=1;
    INHI:=2
  end
  ELSE
  begin
    IHI:=2;
    INHI:=1
  end;
  For I:=1 to MPTS do
  begin
    IF Y[I] < Y[ILO] Then ILO:=I;
    IF Y[I] > Y[IHI] THEN
    begin
      INHI:=IHI;
      IHI:=I
    end
    ELSE IF Y[I] > Y[INHI] THEN
      IF I <> IHI Then INHI:=I
  end;
{ Compute the fractional range from highest to lowest and return if
  satisfactory. }
  RTOL:=2.0*ABS(Y[IHI]-Y[ILO])/(ABS(Y[IHI])+ABS(Y[ILO]));
  IF RTOL < FTOL Then goto 10;  {normal return}
  IF ITER = ITMAX Then
  begin
    Writeln(' Amoeba exceeding maximum iterations.');
    ReadKey;
    Goto 10  {return}
  end;
  ITER:=ITER+1;
  For J:=1 to NDIM do PBAR[J]:=0.0;
  For I:=1 to MPTS do
    IF I <> IHI THEN
      For J:=1 to NDIM do
        PBAR[J]:=PBAR[J] + P[I,J];
  For J:=1 to NDIM do
  begin
    PBAR[J]:=PBAR[J]/(1.0*NDIM);
    PR[J]:=(1.0+ALPHA)*PBAR[J] - ALPHA*P[IHI,J]
  end;
  YPR:=FUNC(PR);
  IF YPR <= Y[ILO] THEN
  begin
    For J:=1 to NDIM do
      PRR[J]:=GAMMA*PR[J] + (1.0-GAMMA)*PBAR[J];
    YPRR:=FUNC(PRR);
    IF YPRR < Y[ILO] THEN
    begin
      For J:=1 to NDIM do P[IHI,J]:=PRR[J];
      Y[IHI]:=YPRR
    end
    ELSE
    begin
      For J:=1 to NDIM do P[IHI,J]:=PR[J];
      Y[IHI]:=YPR
    end
  end
  ELSE IF YPR >= Y[INHI] THEN
  begin
    IF YPR < Y[IHI] THEN
    begin
      For J:=1 to NDIM do P[IHI,J]:=PR[J];
      Y[IHI]:=YPR
    end;
    For J:=1 to NDIM do PRR[J]:=BETA*P[IHI,J] + (1.0-BETA)*PBAR[J];
    YPRR:=FUNC(PRR);
    IF YPRR < Y[IHI] THEN
    begin
      For J:=1 to NDIM do P[IHI,J]:=PRR[J];
      Y[IHI]:=YPRR
    end
    ELSE
      For I:=1 to MPTS do
        IF I <> ILO THEN
        begin
          For J:=1 to NDIM do
          begin
            PR[J]:=0.5*(P[I,J] + P[ILO,J]);
	    P[I,J]:=PR[J]
	  end;
          Y[I]:=FUNC(PR)
	end
  end
  ELSE
  begin
    For J:=1 to NDIM do P[IHI,J]:=PR[J];
    Y[IHI]:=YPR
  end;
  GOTO 1;
10: End;

{main program}
BEGIN

  NDIM:=2;      { 2 variables }
  FTOL:=1e-8;   { User given tolerance }

  {define NDIM+1 initial vertices (one by row) }
  P[1,1]:= 1.0; P[1,2]:=2.0;
  P[2,1]:=-2.0; P[2,2]:=-3.0;
  P[3,1]:= 4.0; P[3,2]:=2.0;

  {Initialize Y to the values of FUNC evaluated 
   at the NDIM+1 vertices (rows] of P }
  For I:=1 to NDIM+1 do
  begin
    PT[1]:=P[I,1]; PT[2]:=P[I,2];
    Y[I]:=FUNC(PT)
  end;

  {call main subroutine}
  AMOEBA(P,Y,MP,NP,NDIM,FTOL,ITER);

  {print results}
  writeln;
  writeln(' Number of iterations: ', ITER);
  writeln;
  writeln(' Best NDIM+1 points:');
  For I:=1 to NDIM+1 do
  begin
    For J:=1 to NDIM do write(' ', P[I,J]);
    writeln
  end;
  writeln;
  writeln(' Best NDIM+1 mimimum values:');
  For I:=1 to NDIM+1 do writeln(' ',Y[I]);
  writeln;
  ReadKey;
  DoneWinCrt

END.
    
{end of file tamoeba.pas}