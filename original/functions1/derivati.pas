{*****************************************************************
*   Program to demonstrate Lagrange derivative interpolation     *
* -------------------------------------------------------------- *
*     Ref.: Basic Scientific Subroutines Vol. II page 314        *
*    By F.R. Ruckdeschel. BYTE/McGRAW-HILL, 1981 [BIBLI 01]      *
*                                                                *
*                         Pascal version by J-P Moreau, Paris    *
*                                  (www.jpmoreau.fr)             *
* -------------------------------------------------------------- *
* SAMPLE RUN:                                                    *
*                                                                *
*   X    2COS(2X)     YP         YP1        ERROR 1     ERROR 2  *
* -------------------------------------------------------------- *
* 0.00   2.000000   1.999961   2.006643  -0.0000394   0.0066434  *
* 0.05   1.990008   1.989970   1.996569  -0.0000385   0.0065603  *
* 0.10   1.960133   1.960096   1.966545  -0.0000373   0.0064118  *
* 0.15   1.910673   1.910637   1.916872  -0.0000357   0.0061991  *
* 0.20   1.842122   1.842088   1.848047  -0.0000337   0.0059246  *
* 0.25   1.755165   1.755134   1.760756  -0.0000314   0.0055908  *
* 0.30   1.650671   1.650642   1.655872  -0.0000288   0.0052011  *
* 0.35   1.529684   1.529659   1.534444  -0.0000259   0.0047595  *
* 0.40   1.393413   1.393391   1.397684  -0.0000227   0.0042704  *
* 0.45   1.243220   1.243201   1.246958  -0.0000193   0.0037386  *
* 0.50   1.080605   1.080589   1.083774  -0.0000157   0.0031694  *
* 0.55   0.907192   0.907180   0.909761  -0.0000120   0.0025685  *
* 0.60   0.724715   0.724707   0.726658  -0.0000081   0.0019420  *
* 0.65   0.534998   0.534993   0.536294  -0.0000042   0.0012961  *
* 0.70   0.339934   0.339934   0.340571  -0.0000002   0.0006372  *
* 0.75   0.141474   0.141478   0.141446   0.0000038  -0.0000280  *
* 0.80  -0.058399  -0.058391  -0.059092   0.0000078  -0.0006929  *
* 0.85  -0.257689  -0.257677  -0.259040   0.0000116  -0.0013510  *
* 0.90  -0.454404  -0.454389  -0.456400   0.0000154  -0.0019955  *
* 0.95  -0.646579  -0.646560  -0.649199   0.0000190  -0.0026201  *
* 1.00  -0.832294  -0.832271  -0.835512   0.0000224  -0.0032185  *
* -------------------------------------------------------------- *
*****************************************************************}
Program Test_derivative;
Uses WinCrt;

Const   NMAX = 100;
        pas = 0.05;

Type
        Tab = ARRAY[1..NMAX] of DOUBLE;

VAR
        X,Y : Tab;
        xx,yy,yy1 : DOUBLE;
        i,n,ndata : INTEGER;
        erreur : BOOLEAN;


  PROCEDURE Deriv(N,NL:INTEGER; X,Y:Tab; x1:DOUBLE; VAR Yp:DOUBLE; VAR error:BOOLEAN);
  {****************************************************
  * Lagrange derivative interpolation procedure Deriv *
  * NL is the level of the interpolation ( ex. NL=3 ) *
  * N is the total number of table values.            *
  * X[i], Y[i] are the coordinate table values, Y[i]  *
  * being the dependant variable. The X[i] may be     *
  * arbitrarily spaced. x1 is the interpolation point *
  * which is assumed to be in the interval with at    *
  * least one table value to the left, and NL to the  *
  * right. Yp is returned as the desired derivative.  *
  * error is set at TRUE if x1 is not in the interval.*
  ****************************************************}
  VAR i,j,k,ll:INTEGER;
      L:array[0..9] of DOUBLE;
      M:array[0..9,0..9] of DOUBLE;
  begin
    error:=FALSE;
    { x1 not in interval [1:N-NL] }
    if (x1<X[1]) OR (x1>X[N-NL]) then
    begin
      error:=TRUE;
      writeln(' STOP: x not between X[1] or X[N-3].')
    end;
    if Not error then
    begin
      i:=0;
      Repeat
        Inc(i)
      Until x1<X[i];
      Dec(i);
      for j:=0 to NL do
      begin
        L[j]:=0.0;
        for k:=0 to NL do M[j,k]:=1.0
      end;
      Yp:=0.0;
      for k:=0 to NL do
      begin
        for j:=0 to NL do
        begin
          if j<>k then
          begin
            for ll:=0 to NL do
            begin
              if ll<>k then
              begin
                if ll=j then
                  M[ll,k]:=M[ll,k]/(X[i+k]-X[i+j])
                else
                  M[ll,k]:=M[ll,k]*(x1-X[j+i])/(X[i+k]-X[i+j])
              end
            end
          end
        end;
        for ll:=0 to NL do
          if ll<>k then L[k]:=L[k]+M[ll,k];
        yp:=yp+L[k]*Y[i+k]
      end
    end
  end;

  {***************************************************
  *     Interpolation of order=2 ( parabola )        *
  *                            by J-P Moreau         *
  ***************************************************}
  PROCEDURE DERIV1(N:INTEGER;X,Y:Tab;X1:DOUBLE;VAR YP:DOUBLE;
                  VAR ERREUR:BOOLEAN);
  Var
      A,B : DOUBLE;
      I : INTEGER;

      {********************************************************
      CALCULE LES COEFFICIENTS A,B DE LA PARABOLE Y=A*X*X+B*X+C
      PASSANT PAR LES 3 POINTS : (X1,Y1),(X2,Y2) ET (X3,Y3)
      COEFFICIENT C NON UTILISE ICI.
      ---------------------------------------------------------
      Calculates coefficients, a, b of parabola Y=A*X+X+B*X=C
      passing through 3 points: (X1,Y1), (X2,Y2) and (X3,Y3).
      Coefficient c is not used here.
      ********************************************************}
      PROCEDURE PARABOLE(X1,Y1,X2,Y2,X3,Y3:DOUBLE; VAR A,B:DOUBLE);
      Var
        ALPHA,BETA,GAMMA,DELTA : DOUBLE;
      Begin
        ALPHA:=X1*X1-X2*X2;
        BETA:=X1-X2;
        GAMMA:=X2*X2-X3*X3;
        DELTA:=X2-X3;
        A:=(Y1-2.0*Y2+Y3)/(ALPHA-GAMMA);
        B:=(Y1-Y2-ALPHA*A)/BETA
      End;

    Begin { Deriv1 }
      ERREUR:=FALSE;
      IF (X1<X[1]) OR (X1>X[N-2]) THEN
      begin
	ERREUR:=TRUE;
	exit
      end;
      I:=0;
      While X1>=X[I] do Inc(I);
      Dec(I);

      PARABOLE(X[I],Y[I],X[I+1],Y[I+1],X[I+2],Y[I+2],A,B);

     { ESTIMATION DE LA DERIVEE EN X1 }

      YP:=2.0*A*X1+B

    End;


{main program: derivative of 2*sin(x)*cos(x) between 0 and 1}
BEGIN
  n:=4;       {level of Lagrange interpolation}
  ndata:=26;  {number of table points}
  {building X & Y Tables}
  for i:=1 to ndata do
  begin
    X[i]:=pas*(i-1);
    Y[i]:=2*cos(X[i])*sin(X[i])
  end;
  xx:=0.0;
  Writeln('     X    2COS(2X)      YP         YP1      ERROR 1    ERROR 2 ');  {heading}
  Writeln('   ------------------------------------------------------------');
  {main loop of derivation}
  Repeat
    Deriv(ndata,n,X,Y,xx,yy, erreur);
    Deriv1(ndata,X,Y,xx,yy1,erreur);
    if erreur then halt(1);
    write('   ',xx:4:2,'  ',2*cos(2*xx):9:6,'  ',yy:9:6,'  ',yy1:9:6);
    writeln('  ',(yy-2*cos(2*xx)):9:6,'  ',(yy1-2*cos(2*xx)):9:6);
    xx:=xx+pas
  Until xx>=1.0+pas;
  {ending section}
  Write('   ------------------------------------------------------------');
  ReadKey; DoneWinCrt
END.

{end of file derivati.pas}